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
// Target holds 28 bytes per pixel for the frame (40 with the caller's geometry) and 24 for the history.
//
// The caller's geometry (RenderDesc::geometry) is drawn in the same passes, into the same targets and the same mean:
// md_gl's representations drawn transparent, through output snippets of md_gl's shaders (geometry_snippet) that do in
// each pass what isomesh.frag does, from the fragment md_gl hands them (its view position, normal, colour), shaded as
// the surfaces are. In four of the passes (geometry_drawn_in): its transmittance is written by its accumulation,
// which compares the stochastic depth in its shader rather than in the depth test, so a pass over the geometry fewer.
// With geometry the peels' colour target holds the opacities and the depths of both peels (o0, o1, z0, z1), so that a
// snippet's pass reads no more than md_gl's four textures of the caller's, and the stochastic depth is 32F.

#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/isolines.h>

#include <stdint.h>
#include <stddef.h>

namespace isomesh {

// The layers peeled exactly, at most (RenderDesc::peels)
static constexpr uint32_t PEELS_MAX = 3;

// Frames averaged before the image is kept as it is. With a sample per pixel a frame is about twice as noisy as one of
// two samples was, and (on an M4 Pro) a fifth of its time: 128 rather than 64
static constexpr uint32_t MAX_FRAMES = 128;

void initialize();
void shutdown();

// The passes of a frame, in the order render() draws them (isomesh.frag's PASS_*)
enum Pass {
    Pass_Transmittance = 0,     // the coverage of everything together, and the optical depth of what the surfaces enclose
    Pass_Peel0         = 1,     // the nearest layer, exactly
    Pass_Peel1         = 2,     // the one behind it
    Pass_Depth         = 3,     // the stochastic depth of the layers behind those two
    Pass_Accumulate    = 4,     // the radiance of everything times its visibility
    Pass_Count
};

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
    uint32_t tex_frame_color = 0;      // radiance and opacity times visibility (RGBA16F)
    uint32_t tex_frame_tod   = 0;      // log transmittance and optical depth (RG16F)
    uint32_t tex_frame_depth = 0;      // the stochastic depth (32F with the caller's geometry, which reads it)
    uint32_t tex_peel_depth[2] = {};   // the depth tests of the peels, a pass each in turn (32F): pass k tests against
                                       // [k % 2] and reads the one before from the other
    uint32_t tex_peel_z = 0;           // the peeled layers' depths (RGBA32F: z0, z1, z2, a layer not peeled the depth of
                                       // the one in front), and with the caller's geometry their opacities packed in a
    uint32_t tex_peel_o = 0;           // their opacities (RGBA16F: o0, o1, o2)
    uint32_t fbo_color = 0;         // tex_frame_color + tex_frame_depth
    uint32_t fbo_tod   = 0;         // tex_frame_tod
    uint32_t fbo_peel[2] = {};      // tex_peel_z + tex_peel_o + tex_peel_depth[i]
    uint32_t fbo_pack = 0;          // tex_peel_z: the opacities packed into its a (the caller's geometry)
    uint32_t fbo_geometry = 0;      // tex_frame_color + tex_frame_tod: the caller's geometry's accumulation (peel_depths)
    uint32_t fbo       = 0;         // the history being written
    bool     peel_depths = false;   // made for the caller's geometry (fbo_geometry, the stochastic depth 32F)
    uint32_t tex_history_color[2] = {};
    uint32_t tex_history_t[2]     = {};
    uint32_t cur = 0;           // which history holds the mean
    uint32_t frame = 0;         // frames in the mean
    uint64_t hash = 0;          // of what the mean is of
    uint64_t triangles = 0;     // drawn in its last frame, of all the surfaces
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

    // The layers peeled exactly in front of the stochastic ones (PEELS_MAX at most): each is a pass over all of the
    // surfaces and the geometry a frame, and takes a layer from the stochastic estimate, the noisy part. With none every
    // layer is stochastic.
    uint32_t peels = 2;

    // The stochastic thresholds of another seed: the same mean in expectation, its frames other independent estimates
    // (for measuring the noise against a reference, as the profiler's sweep does)
    uint32_t seed = 0;

    // Geometry of the caller's in the same passes: in each it is drawn in (geometry_drawn_in), after the surfaces, with
    // the pass's targets bound and its state set, draw(pass, user), drawing with md_gl's shaders made of
    // geometry_snippet(pass), and geometry_set_draw before each md_gl_draw; what the snippets read is bound at
    // MD_GL_USER_TEXTURE_UNIT + 0..3. With it the mean is of all of the target rather than of the surfaces' bounds, and
    // render() draws without surfaces too. What it draws belongs in scene_hash.
    struct {
        void (*draw)(Pass pass, void* user) = nullptr;
        void* user = nullptr;
    } geometry;
};

// Returns true when the colour target holds this pass's result (cleared and/or drawn into). 'converged' (optional)
// says whether the image is final: nothing changes it until the desc does.
bool render(Target* target, const RenderDesc& desc, bool* converged = nullptr);

// The passes the caller's geometry is drawn in: all but Pass_Transmittance, whose part its accumulation does
bool geometry_drawn_in(Pass pass);

// The output snippet of md_gl's shaders (md_gl_shaders_create) for the caller's geometry in a pass it is drawn in. A
// caller's defines may come before it: DIMMED_SATURATION, the saturation kept by atoms dimmed while others are selected
// (1 by default), as the G-buffer's snippet has it.
str_t geometry_snippet(Pass pass);

// Before each draw of the caller's geometry in a pass: the opacity of what it draws (its colours' alpha times this), its
// index of refraction (the reflection that adds to the opacity toward its silhouette, as RenderDesc::shading.ior is the
// surfaces'; 1 none), and an id that tells its thresholds from those of another draw of the same atoms (another
// representation's), the same in every pass of a frame. Binds its parameters for the snippets (MD_GL_USER_TEXTURE_UNIT
// + 1), written once per frame and id, never while queued draws may still read them.
void geometry_set_draw(float opacity, float ior, uint32_t id);

}  // namespace isomesh
