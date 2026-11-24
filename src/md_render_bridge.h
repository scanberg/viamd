#pragma once

// The molecular representations of the main view drawn by md_render (mdlib) on md_gpu, while the rest of viamd still
// draws with OpenGL. A step on the way to md_render drawing everything: what it draws is merged into the G-buffer, so
// post-processing, picking and the marks see one scene.
//
// Each frame the scene is made (render_scene), md_render draws the representations it can (md_render_bridge_draws),
// with the camera and the jitter of md_gl's, and exports its targets into textures shared with the GL context
// (md_render_view_export, md_gpu_gl). After md_gl's opaque draws, a full screen pass writes them into the G-buffer
// where they are nearer than what md_gl drew: the colour (dimmed as md_gl dims what is not selected), the normal and
// the velocity in the G-buffer's encodings, and the depth. md_gl draws the rest, and the marks.
//
// Experimental, off unless enabled (the Debug window, or VIAMD_MD_RENDER=1 in the environment). Needs a device that
// draws (md_gpu's graphics) and shares textures with the GL context (the same GPU, GL_EXT_memory_object or IOSurface);
// otherwise md_render_bridge_status says why, and md_gl draws everything.

#include <stddef.h>
#include <stdint.h>

struct ApplicationState;
struct Representation;
struct PickingHit;
struct MdRenderBridge;

// On (made if it is not yet) or off. False when it cannot be had; then md_render_bridge_status says why
bool md_render_bridge_enable(ApplicationState* state, bool enable);
void md_render_bridge_destroy(ApplicationState* state);
// What it is doing, or why it is not
const char* md_render_bridge_status(const ApplicationState& state);

// Whether md_render draws (what it can) this frame
bool md_render_bridge_active(const ApplicationState& state);
// Whether md_render draws the representation this frame: then md_gl leaves it (draw_representations_opaque)
bool md_render_bridge_draws(const ApplicationState& state, const Representation& rep);

// Picking. Where md_render's surface is nearest, the G-buffer's picking target holds the one index of
// PickingDomain_MdRender. Each frame the GL picking reads a pixel, md_render is asked about the same one (pick, with the
// pixel as the GL readback has it: of the G-buffer, from the bottom), and a hit of that index read back is resolved to
// what md_render found there (resolve_hit): an atom, a bond, a segment of the cartoon, a ring, a gap, in their domains,
// as md_gl's would be. Both arrive a frame or two later; a hit md_render has no answer for yet is none.
void md_render_bridge_pick(ApplicationState* state, float x, float y);
void md_render_bridge_resolve_hit(ApplicationState* state, PickingHit* hit);

// The colours of a representation's atoms as md_gl gets them, which also say which it shows (alpha 0: not shown)
void md_render_bridge_rep_colors(ApplicationState* state, const Representation& rep, const uint32_t* colors, size_t count);

// The pairs of a Bonds representation as md_gl gets them (representation_set_weak_bonds): their atoms and their bonds'
// flags (md_bond_flags_t, for the kind; null: covalent)
void md_render_bridge_rep_pairs(ApplicationState* state, const Representation& rep, const uint32_t* atom_a, const uint32_t* atom_b,
                                const uint32_t* flags, size_t count);

// From update_md_buffers: what changed (MolBit_ flags), taken up when the scene is next made
void md_render_bridge_buffers(ApplicationState* state, uint32_t dirty);

// The frame: drawn by md_render and exported (render_scene, after update_view_param)
void md_render_bridge_render(ApplicationState* state);

// Into the G-buffer bound, after md_gl's opaque draws: its colour, normal, velocity and picking targets drawn, its depth
// tested and written
void md_render_bridge_merge(ApplicationState* state);

// The marks of the selection and the highlight (with its groups) of what md_render draws, into md_gl's marks pass
// (after draw_representations_marks, which leaves those representations): md_render's masks, in view where its surface
// is the scene's. When the marks changed without the scene being made, md_render makes its masks again first.
void md_render_bridge_marks(ApplicationState* state);
