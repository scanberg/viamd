#pragma once

// The main view, drawn on md_gpu (mdlib): md_render draws the representations it can and the shapes of the dipoles and
// the cell (their arrows, its edges); md_compose shades what it drew (the ambient occlusion, the shadows, the
// transparent layers) with the electronic structure representations' volumes (gpu_volume.h: their raycast isosurfaces,
// volume renderings and planes; their isosurfaces drawn as meshes are md_render's triangles); md_post makes the image of
// it (the display transform, the temporal antialiasing at the size of the image shown) and presents it with the marks
// of the selection and the highlight into a texture shared with OpenGL, which draws it into the framebuffer bound,
// under the GUI. OpenGL draws nothing of the scene.
//
// What is not in the main view yet: the glyphs, what a script shows, the components' drawings (the critical points).
// See render_scene (main.cpp).
//
// Needs a device that draws (md_gpu's graphics) and shares textures with the GL context (the same GPU,
// GL_EXT_memory_object or IOSurface); otherwise md_render_bridge_status says why, and the main view is empty.

#include <stddef.h>
#include <stdint.h>

struct ApplicationState;
struct Representation;
struct PickingHit;
struct MdRenderBridge;

namespace postprocess_pipeline { struct Settings; }

// Made (if it is not yet), at start. False when it cannot be had; then md_render_bridge_status says why
bool md_render_bridge_create(ApplicationState* state);
void md_render_bridge_destroy(ApplicationState* state);
// What it is doing, or why it is not
const char* md_render_bridge_status(const ApplicationState& state);

// Picking: md_render is asked about the pixel under the mouse (x, y in the pixels of the scene drawn, rows from the
// top), and its answer, a frame or two later, is the hit (hit): an atom, a bond, a segment of the cartoon, a ring, a
// gap, a dipole, in their domains as md_gl's were, or nothing (the background, the cell's edges). Its surface_coord
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

// The image kept, with the marks over it (settings.marks; null: none), into the framebuffer bound at the size of the
// image shown (the viewport the caller set). When the marks changed without the scene being made, md_render makes its
// masks again first. False: there is no image (nothing was made yet).
bool md_render_bridge_present(ApplicationState* state, const postprocess_pipeline::Settings* marks);

// The depth of the scene drawn last as OpenGL reads it (for what is drawn over the image behind the scene: the planes'
// outlines): a texture of the G-buffer's size, GL's window depth (1 the far plane), rows from the bottom. 0: none.
uint32_t md_render_bridge_depth_gl(ApplicationState* state);
