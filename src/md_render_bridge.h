#pragma once

// The main view, drawn on md_gpu (mdlib): md_render draws the representations (the glyphs among them) and the shapes
// of the dipoles and the cell (their arrows, its edges); md_compose shades what it drew (the ambient occlusion, the
// shadows, the transparent layers) with the electronic structure representations' volumes (gpu_volume.h: their raycast
// isosurfaces, volume renderings and planes; their isosurfaces drawn as meshes are md_render's triangles); md_post
// makes the image of it (the display transform, the temporal antialiasing at the size of the image shown) and presents
// it with the marks of the selection and the highlight, and the overlay over those (viamd.h: what the script shows,
// the planes' outlines, the critical points), into a texture shared with OpenGL, which draws it into the framebuffer
// bound, under the GUI. OpenGL draws nothing of the scene. See render_scene and render (main.cpp).
//
// Needs a device that draws (md_gpu's graphics) and shares textures with the GL context (the same GPU,
// GL_EXT_memory_object or IOSurface); otherwise md_render_bridge_status says why, and the main view is empty.

#include <stddef.h>
#include <stdint.h>
#include <core/md_vec_math.h>

struct ApplicationState;
struct Representation;
struct PickingHit;
struct MdRenderBridge;
struct MdRenderPanels;
struct Overlay;
struct OverlayItem;
struct IsoMeshes;
struct md_bitset_t;
struct InteractionSurfaceState;

namespace postprocess_pipeline { struct Settings; }

// Made (if it is not yet), at start. False when it cannot be had; then md_render_bridge_status says why
bool md_render_bridge_create(ApplicationState* state);
void md_render_bridge_destroy(ApplicationState* state);
// What it is doing, or why it is not
const char* md_render_bridge_status(const ApplicationState& state);

// Picking: md_render is asked about the pixel under the mouse (x, y in the pixels of the scene drawn, rows from the
// top), and its answer, a frame or two later, is the hit (hit): an atom, a bond, a segment of the cartoon, a ring, a
// gap, a dipole, a glyph, in their domains as md_gl's were, or nothing (the background, the cell's edges). An item of
// the overlay drawn last that picks is the hit at once, wherever it is over the scene. Its surface_coord
// and depth are as the GL picking's were (rows from the bottom, GL's window depth); the caller sets the rest (source,
// screen_coord, world_pos). False: no answer (yet).
void md_render_bridge_pick(ApplicationState* state, float x, float y);
bool md_render_bridge_hit(ApplicationState* state, PickingHit* hit);

// The colours of a representation's atoms as md_gl gets them, which also say which it shows (alpha 0: not shown)
void md_render_bridge_rep_colors(ApplicationState* state, const Representation& rep, const uint32_t* colors, size_t count);

// The pairs of a Bonds representation as md_gl gets them (representation_set_weak_bonds): their atoms and their bonds'
// flags (md_bond_flags_t, for the kind; null: covalent)
void md_render_bridge_rep_pairs(ApplicationState* state, const Representation& rep, const uint32_t* atom_a, const uint32_t* atom_b,
                                const uint32_t* flags, size_t count);

// From update_md_buffers: what changed (MolBit_ flags), taken up when the scene is next made
void md_render_bridge_buffers(ApplicationState* state, uint32_t dirty);

// The scene (render_scene, after update_view_param): drawn by md_render at the G-buffer's size, composed and made into
// the image kept, at the size of the image shown (gbuffer.out_width, out_height): the temporal antialiasing's history
// when it is on (settings.taa), else the frame itself. The settings' shading, ambient occlusion, display transform,
// background and TAA are used; the rest are not done.
void md_render_bridge_render(ApplicationState* state, const postprocess_pipeline::Settings& settings);

// The image kept, with the marks over it (settings.marks; null: none) and the overlay over those (viamd.h; null: none),
// into the framebuffer bound at the size of the image shown (the viewport the caller set). When the marks changed
// without the scene being made, md_render makes its masks again first. The overlay is kept to be picked
// (md_render_bridge_hit). False: there is no image (nothing was made yet).
bool md_render_bridge_present(ApplicationState* state, const postprocess_pipeline::Settings* marks, const Overlay* overlay);

// ---- Panels -------------------------------------------------------------------------------------------------------
// The component windows' views (the Orbital Grid's, the Transition Analysis', the Charge Transfer side views, the
// Density Volume's): panels of one size and one camera, each the molecule (representations of the system's atoms)
// with surfaces of its own: volumes (gpu_volume.h) raycast by md_compose, as isosurfaces or through a transfer function,
// or the sheets of meshes on md_gpu (iso_meshes.h, gpu_mesh), md_render's triangles. The panels share one md_render view
// while none has meshes (the molecule rasterised once), else each has its own. Composed, displayed and settled
// (md_render's samples of a still view, a jittered raster each frame, averaged by md_post's resolve until there are
// enough, after which nothing is drawn until something changes), presented with the marks of the selection and the
// highlight and an overlay, side by side in an image of them all: a GL texture for ImGui. The system's positions and
// marks are the main view's (taken up by md_render_bridge_render): the panels follow them a frame later.

constexpr int PANEL_TRANSFER_RES = 256;     // a transfer function's texels (PanelVolume::transfer)

// A volume of a panel (gpu_volume.h), within a box of its texture space: its isosurfaces, raycast, and its emission and
// absorption through a transfer function (direct volume rendering); either or both
struct PanelVolume {
    uint32_t           volume = 0;
    uint32_t           iso_count = 0;               // at most 8
    const float*       iso_values = nullptr;
    const md_float4_t* iso_colors = nullptr;        // sRGB, alpha: coverage
    const uint32_t*    transfer = nullptr;          // PANEL_TRANSFER_RES RGBA8: sRGB colour, the opacity of a slab 1/150 of
                                                    // the box thick; null: none
    float              transfer_range[2] = {0.0f, 1.0f};   // the values at its ends
    md_float3_t        clip_min = {0.0f, 0.0f, 0.0f}, clip_max = {1.0f, 1.0f, 1.0f};
};

// The sheets of a mesh shown in a panel: a range of its triangles each, in its colour (sRGB, alpha: coverage)
struct PanelSheet {
    uint32_t    index_offset = 0, index_count = 0;
    md_float4_t color = {};
};

struct PanelSurfaces {
    const PanelVolume* volumes = nullptr;           // at most 4
    uint32_t           volume_count = 0;
    const IsoMeshes*   meshes = nullptr;            // on md_gpu (gpu_mesh); null: none
    const PanelSheet*  sheets = nullptr;
    uint32_t           sheet_count = 0;
};

// A representation of the molecule in the panels: of the atoms of `atoms` (null: all) those not transparent in its
// colours (PanelsDesc::colors[colors]), under a transform (null: none). Of rep its type and parameters (spacefill,
// licorice, ball and stick, ribbons, cartoon), not its colours nor its filter.
struct PanelRep {
    const Representation* rep = nullptr;
    uint32_t              colors = 0;
    const md_bitset_t*    atoms = nullptr;
    const md_float4x4_t*  transform = nullptr;
};

struct PanelsDesc {
    md_float4x4_t view = {}, proj = {};         // world to view, view to clip (OpenGL's, as camera.h makes them)
    uint32_t width = 0, height = 0;             // each panel's, pixels
    uint32_t count = 0;                         // at most 16
    uint32_t columns = 1;                       // the image's: panel i at column i % columns, row i / columns
    const PanelSurfaces* panels = nullptr;      // count of them
    const uint32_t* const* colors = nullptr;    // sets of colours, RGBA8 per atom of the system, alpha 0: not shown
    uint32_t color_count = 0;
    const PanelRep* reps = nullptr;             // the molecule, in every panel
    uint32_t rep_count = 0;
    const OverlayItem* overlay = nullptr;       // over every panel (viamd.h), not picked
    uint32_t overlay_count = 0;
    const postprocess_pipeline::Settings* settings = nullptr;  // the shading, tone mapping, background; the marks
    bool  marks = true;                         // shown
};

// Brings the panels' image up to date (made at the first call, *panels null, and again at another size or count):
// what changed drawn again, another sample of what settles, nothing of what has. All of the desc is read at the call.
// False: no image (the bridge does not run, or has no system yet).
bool md_render_bridge_panels_draw(ApplicationState* state, MdRenderPanels** panels, const PanelsDesc& desc);
// The image: its GL texture (GL_TEXTURE_2D, RGBA8, rows from the bottom) and panel i's corners in it, as ImGui's
// AddImage takes them (uv0 its top left, uv1 its bottom right). False: none (yet).
bool md_render_bridge_panels_image(const MdRenderPanels* panels, uint32_t i, uint32_t* gl_texture, md_float2_t* uv0, md_float2_t* uv1);
// Picking in panel i (x, y: its pixels, rows from the top), and the answer of a frame or two later, as the main view's
// (md_render_bridge_pick, md_render_bridge_hit): an atom or a bond, or nothing (a mesh, the background)
void md_render_bridge_panels_pick(ApplicationState* state, MdRenderPanels* panels, uint32_t i, float x, float y);
bool md_render_bridge_panels_hit(ApplicationState* state, MdRenderPanels* panels, PickingHit* hit);
// The two for the mouse over panel i, an interaction surface (its points, as ImGui has them) of the size of a panel,
// width x height pixels: the hit (*hit, empty for none) as a surface's, of the source, where it is in the world
// (clip_to_world: the panels' camera's)
void md_render_bridge_panels_surface_hit(PickingHit* hit, ApplicationState* state, MdRenderPanels* panels, uint32_t i,
                                         const InteractionSurfaceState& surface, uint64_t source, const md_float4x4_t& clip_to_world,
                                         int width, int height);
void md_render_bridge_panels_destroy(ApplicationState* state, MdRenderPanels** panels);
