#pragma once

// ISOSURFACE MESHES
// The alternative to raycasting the isosurfaces of a volume (volume::render_isosurfaces): triangle meshes
// extracted from the field itself (md_gto_iso), drawn with stochastic transparency a sample per pixel at a time
// and averaged over frames until the image has converged, when it is kept as it is until something changes.
//
// The result is what the raycaster gives: premultiplied linear radiance with coverage in alpha, blended over
// the colour target with (ONE, ONE_MINUS_SRC_ALPHA), under the same lights and shading model, so the two are
// interchangeable for the caller. Every surface is in the one composite, whatever the order of the draws and
// however many there are; no sorting and no extra passes for the caller.
//
// The frames are independent estimates of the composite; a Target keeps their mean. It restarts when anything
// it was given changes (the meshes, their colours, the view, the projection without its jitter, the shading,
// the size, desc.scene_hash), and after MAX_FRAMES frames it stops drawing the meshes and only composites what
// it has (unless desc.keep_drawing). A caller drawing into a view that is redrawn only on change (a panel) redraws
// until 'converged'.
//
// A frame is five passes over the meshes (isomesh.frag: the transmittance, two layers peeled exactly, the
// stochastic depth behind them, the radiance times visibility), each a program of its own, into a target of a sample
// per pixel, and one pass into the mean that also puts it over the caller's target (once the mean is kept, that
// alone), both over the pixels the surfaces' bounds cover. The mean antialiases the edges: it is of the caller's
// jittered projections, or when there are none (proj as proj_no_jitter) of a jitter of its own over the pixel. A
// Target holds 28 bytes per pixel for the frame and 24 for the history (the multisampled experiments: 28 per sample of
// the frame, the samples per pixel (4, 2 or 1) as many as fit MEMORY_BUDGET unless set_experiment_samples says).

#include <core/md_vec_math.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/isolines.h>

#include <stdint.h>
#include <stddef.h>

namespace isomesh {

// Frames averaged before the image is kept as it is. With a sample per pixel a frame is about twice as noisy as one of
// two samples was, and (on an M4 Pro) a fifth of its time: 128 rather than 64
static constexpr uint32_t MAX_FRAMES = 128;

// Bytes of multisampled storage a Target may take in the multisampled experiments: 4 samples per pixel up to about
// 1920 x 1200, 2 up to about 2560 x 1600, 1 above
static constexpr size_t MEMORY_BUDGET = (size_t)272 << 20;

void initialize();
void shutdown();

// Triangles (counter-clockwise seen from outside, which is where the normals point) in one vertex and index
// buffer on the GL side. Surfaces are ranges of the index buffer.
struct Mesh {
    uint32_t vao = 0;
    uint32_t vbo = 0;   // position xyz, normal xyz: 6 floats per vertex
    uint32_t ibo = 0;
    uint32_t vbo_value = 0;     // a scalar per vertex and its gradient (mesh_upload_values), 0 until there are any
    uint32_t vertex_count = 0;
    uint32_t index_count  = 0;
    vec3_t   bounds_min = {0, 0, 0};    // of the vertices, in the mesh's own coordinates
    vec3_t   bounds_max = {0, 0, 0};
    uint64_t version = 0;   // bumped on every upload
};

// position and normal 'stride' floats apart per vertex (3 or 4), xyz used
bool mesh_upload(Mesh* mesh, const float* position, const float* normal, size_t stride, uint32_t vertex_count, const uint32_t* index, uint32_t index_count);

// A scalar per vertex (a field evaluated at the vertices), for Coloring_VertexField: vertex_count values,
// interpolated across the triangles, and optionally its gradient along the surface (packed float3 per vertex, per
// unit of the mesh's coordinates; mesh_field_gradient), which the isolines take their width from. NULL values drop
// them; NULL gradients are zero. A new upload of the mesh drops them too.
bool mesh_upload_values(Mesh* mesh, const float* values, const float* gradients, uint32_t count);

// The gradient along the surface of a scalar given at the vertices of triangles (xyz packed float3, index three per
// triangle), per vertex: the gradient of its linear interpolation over each triangle, averaged over the triangles
// around the vertex by their area. 'out' is packed float3, vertex_count of them.
void mesh_field_gradient(float* out, const float* xyz, const float* values, uint32_t vertex_count, const uint32_t* index, uint32_t index_count);
void mesh_free(Mesh* mesh);

// What one view keeps between frames: the frame, the running mean, and the state it belongs to
struct Target {
    uint32_t width = 0, height = 0;
    uint32_t samples = 0;           // of the frame's textures, per pixel: 1 but in the multisampled experiments
    uint32_t tex_frame_color = 0;      // radiance and opacity times visibility (RGBA16F)
    uint32_t tex_frame_tod   = 0;      // log transmittance and optical depth (RG16F)
    uint32_t tex_frame_depth = 0;      // the stochastic depth
    uint32_t tex_peel_depth[2] = {};   // the depths of the two nearest layers (32F)
    uint32_t tex_peel_opacity  = 0;    // and their opacities (RG16F)
    uint32_t fbo_color = 0;         // tex_frame_color + tex_frame_depth
    uint32_t fbo_tod   = 0;         // tex_frame_tod
    uint32_t fbo_peel[2] = {};      // tex_peel_opacity + tex_peel_depth[i]
    uint32_t fbo       = 0;         // the history being written
    uint32_t tex_history_color[2] = {};
    uint32_t tex_history_t[2]     = {};
    uint32_t cur = 0;           // which history holds the mean
    uint32_t frame = 0;         // frames in the mean
    uint64_t hash = 0;          // of what the mean is of
    uint64_t triangles = 0;     // drawn in its last frame, of all the surfaces
    uint32_t tex_target = 0;    // of the frame's textures: GL_TEXTURE_2D, or GL_TEXTURE_2D_MULTISAMPLE (the experiments)
    uint32_t tex_moments = 0;   // the moment-based experiments: the moments of the absorbance over depth (RGBA32F)
    uint32_t samples_asked = 0; // the samples per pixel asked for (an implementation may give more)
    uint32_t seed = 0;          // target_set_seed
};

void target_free(Target* target);

enum Coloring {
    Coloring_Uniform,
    Coloring_ColorVolume,   // the colour times the colour volume
    Coloring_Field,         // the field volume through its colour map
    Coloring_VertexField,   // the values of the mesh's vertices (mesh_upload_values) through the colour map
};

struct Surface {
    const Mesh* mesh = nullptr;
    uint32_t index_offset = 0;
    uint32_t index_count  = 0;
    mat4_t   model = MD_MAT4_IDENT_INIT;      // mesh -> world
    vec4_t   color = {1, 1, 1, 1};      // alpha: the surface's opacity
    float    optical_density = 0.0f;    // of what it encloses, as the raycaster's (IsoRenderDesc::iso.optical_densities)

    Coloring coloring = Coloring_Uniform;
    uint32_t color_volume   = 0;        // 3D, RGBA
    uint32_t field_volume   = 0;        // 3D, R
    uint32_t field_colormap = 0;        // 2D, N x 1
    mat4_t   world_to_volume = MD_MAT4_IDENT_INIT;    // world -> [0,1]^3 of the colour or field volume
    float    field_range_beg = 0.0f;
    float    field_range_end = 1.0f;

    // Coloring_VertexField: lines where the field takes given values (gfx/isolines.h), their width in pixels from
    // the gradient the values were uploaded with and ray differentials (isomesh.frag)
    IsolineDesc isolines = {};
};

struct RenderDesc {
    struct {
        uint32_t depth = 0;     // the opaque scene (single sampled), optional
        uint32_t color = 0;     // RGBA16F, required
        uint32_t width = 0;
        uint32_t height = 0;
        bool clear_color = false;
    } render_target;

    struct {
        mat4_t view = {};
        mat4_t proj = {};               // as the opaque scene was drawn, jitter and all
        mat4_t proj_no_jitter = {};     // what the accumulation is keyed on
    } matrix;

    const Surface* surfaces = nullptr;
    size_t num_surfaces = 0;

    // As volume::IsoRenderDesc::shading
    struct {
        float roughness = 0.4f;
        float ior = 1.5f;
        postprocess_pipeline::ShadingMode mode = postprocess_pipeline::ShadingMode_PhysicallyBased;
        float ambient = 0.3f;
        float specular = 0.25f;
        postprocess_pipeline::ColorSpace colors = postprocess_pipeline::ColorSpace_Linear;
        float light = 1.0f;
        float environment = 1.0f;
    } shading;

    // Whatever else the image depends on that the desc does not show: the opaque scene behind the surfaces,
    // the contents of the colour and field volumes
    uint64_t scene_hash = 0;

    // For profiling a still view: the meshes drawn every frame, past MAX_FRAMES too, the mean then of about the last
    // MAX_FRAMES frames (each new one weighs 1 / MAX_FRAMES). It never converges.
    bool keep_drawing = false;
};

// Returns true when the colour target holds this pass's result (cleared and/or drawn into). 'converged' (optional)
// says whether the image is final: nothing changes it until the desc does.
bool render(Target* target, const RenderDesc& desc, bool* converged = nullptr);

// EXPERIMENTS
// For timing the passes and measuring their noise (Windows > Profiler, Sweep): what render() does differently, for
// every Target, until it is set back. Those marked 'timing only' do not draw the right image.
//
// The moment-based ones (Muenstermann, Kolb, Kolb and Peters, "Moment-Based Order-Independent Transparency", 2018):
// the transmittance pass also sums, per pixel, the absorbance -log(1 - opacity) of the fragments times the powers 1 to 4
// of their depth (warped to [-1, 1] over the surfaces' bounds); the accumulation takes the transmittance in front of
// each fragment from the bound those moments give, instead of from the peels and the stochastic depth. Deterministic,
// so the frames differ by their jitter alone; approximate where few, opaque layers lie close in depth.
//
// The multisampled ones: per sample into a multisampled target (the transmittance at pixel rate, MSAA coverage putting
// it into the samples). Under sample shading Apple's OpenGL gives gl_FragCoord.z at the pixel's centre, not at the
// sample the depth buffer holds, so there the passes have to write the depth they compare (gl_FragDepth), which takes
// the depth test from before the shading to after it. The spec, Mesa, NVIDIA and AMD give the sample's depth: there
// the rasterizer's will do (fragcoord_z_per_sample tells which), and the first peel, which compares with nothing, runs at
// pixel rate as plain MSAA.
enum Experiment {
    Experiment_None = 0,
    Experiment_Multisampled,            // as it was: one program for the passes after the transmittance, writing
                                        // gl_FragDepth
    Experiment_MultisampledPerPass,     // a program per pass, writing gl_FragDepth
    Experiment_MultisampledNoFragDepth, // a program per pass, the rasterizer's depth, the first peel at pixel rate:
                                        // timing only where fragcoord_z_per_sample is false (Apple's OpenGL)
    Experiment_SeparateResolve,         // the mean and the composite in passes of their own, over all of the target
    Experiment_Moments,                 // moment-based: the transmittance in front of a fragment reconstructed from four
                                        // moments of the absorbance over depth, gathered with the transmittance; no
                                        // peels, no stochastic depth: two passes over the meshes
    Experiment_MomentsPeel,             // the nearest layer peeled exactly, the moments of those behind it: three passes
    Experiment_Count
};
void set_experiment(Experiment e);
Experiment experiment();
const char* experiment_name(Experiment e);     // a short one, for tables
bool experiment_multisampled(Experiment e);

// The samples per pixel of the multisampled experiments: 0 as many as fit MEMORY_BUDGET (at most 4), or n (at most what
// the implementation has)
void set_experiment_samples(uint32_t n);
uint32_t experiment_samples();

// Two Targets of different seeds draw independent frames: the thresholds of the stochastic depth and the jitter of
// render()'s own (proj as proj_no_jitter) are offset by it (Cranley-Patterson rotations; 0 leaves them as they are).
// The mean restarts.
void target_set_seed(Target* target, uint32_t seed);

// The composite of the mean as it is (what render() puts over the caller's target): premultiplied linear radiance and
// coverage, width x height RGBA floats from the bottom row up. False when there is no mean.
bool target_read_mean(const Target* target, float* rgba);

// Whether gl_FragCoord.z under sample shading is the depth of the sample (as the rasterizer writes it into the depth
// buffer) rather than of the pixel's centre: tested once, a steep triangle into four samples, gl_FragCoord.z written
// as a colour and compared with the depth buffer sample by sample. 'max_diff' (optional) gets the largest difference
// found. False also when it could not be tested ('max_diff' then negative).
bool fragcoord_z_per_sample(float* max_diff = nullptr);

}  // namespace isomesh
