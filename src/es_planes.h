#pragma once

// The planes of the electronic structure representations: the cuts of a volume rendering, and the plane a field is
// drawn on (ElectronicStructureRepresentation::planes, Plane in plane.h). Shown in the view as outlines where they
// cut their representation's box (hidden behind what of the scene is in front of them), a click on an outline (or on
// its row in the representation's panel) picking it for the gizmo, a click on nothing letting go of it (main.cpp)
// (es_planes.cpp): a ring to turn its normal as a trackball, and the normal as an arrow to slide it along. A slice's
// normal is kept on the camera's side of it: it turns over as the plane is seen edge on.
//
// And what the volume rendering of a representation is drawn with: its planes in the space of its volume, and its
// transfer function baked, the colour map with the opacity curve over it.

#include <core/md_vec_math.h>
#include <plane.h>

struct ApplicationState;
struct ElectronicStructureRepresentation;
struct Plane;
struct PlaneFitInfo;


// The box of a representation's volume in world space, the unit cube of its texture: false when it has none yet
bool es_volume_box(mat4_t* box_to_world, const ElectronicStructureRepresentation& es);

// The outlines of the planes of the representations shown, the picking of one by a click on its outline, and the
// gizmo on the picked one. Called in the main view's canvas, before its interaction surface, with the rectangle the
// scene is drawn over (screen coordinates, the whole viewport): true when the mouse is the planes' this frame (over the gizmo, using it, or
// picking an outline), in which case the view is not to turn or pick with it.
bool es_planes_gizmo(ApplicationState* state, vec2_t canvas_min, vec2_t canvas_size);

// The outlines es_planes_gizmo found this frame, drawn into the framebuffer bound (fb_width x fb_height pixels, the
// window's: the image presented, before the GUI): hidden where the scene, depth_tex (the G-buffer's, of its size), is
// in front of them. They write no depth. px_per_point: framebuffer pixels per ImGui point.
void es_planes_draw_outlines(uint32_t depth_tex, int depth_width, int depth_height, int fb_width, int fb_height, float px_per_point);
void es_planes_shutdown();

// A plane through the middle of the representation's box, facing away from the camera: as a cut it removes the half
// in front, showing the inside. Returns its index, or -1 when there is no room (ES_MAX_PLANES) or no box yet.
int es_planes_add(ApplicationState* state, ElectronicStructureRepresentation* es, uint32_t rep_id);

// A plane put through the selected atoms (plane_fit: three exactly, more by least squares, turned as little as they
// allow and kept on its side), at their positions in the frame shown, the nearest images of the first where the
// system is periodic, and how (info, optional). The plane is returned, and left as it was in 'plane' with only
// 'info' wanted (nullptr: it is set); kind None when no atoms are selected.
Plane es_plane_fit_selection(const ApplicationState* state, const Plane& plane, PlaneFitInfo* info = nullptr);

// The enabled cuts of a representation (its volume display) in the space of its volume's texture ([0,1]^3), as n . x + d >= 0 for the
// side kept: what the volume rendering clips to. Returns how many (at most ES_MAX_PLANES).
int es_volume_clip_planes(vec4_t* out, const ElectronicStructureRepresentation& es);

// A plane of a representation in the space of its volume's texture, as above
vec4_t es_volume_plane(const ElectronicStructureRepresentation& es, const Plane& plane);

// The enabled slices of a representation (its planes display) of some opacity, in the space of its volume's texture,
// and their opacities: what the slice rendering draws. Returns how many (at most ES_MAX_PLANES).
int es_volume_slices(vec4_t* out, float* out_opacity, const ElectronicStructureRepresentation& es, int* out_index = nullptr);

// The images of planes (the planes display). A plane's frame: the plane with its normal one of its two ways (the
// same frame either way it is given), u and v along it, set by the box (not the view: the image stays put as the view
// turns), and the polygon where it cuts the box (box_to_world: the unit cube to the world) with the rectangle around it
// in (u, v), Angstrom. False where the plane misses the box.
struct EsPlaneFrame {
    Plane  plane;
    vec3_t u, v;
    float  rect[4];     // s0, s1, t0, t1: along u, along v
    vec3_t poly[6];
    int    num_poly;
};
bool es_plane_frame(EsPlaneFrame* f, const Plane& plane, const mat4_t& box_to_world);

// The grid of an image over rect (in the frame) at spacing Angstrom a pixel: one layer whose voxel centres are on the
// plane, in bohr as the evaluations take it
void es_plane_image_grid(struct md_grid_t* grid, const EsPlaneFrame& f, const float rect[4], float spacing);

// Whether the gizmo is dragging this plane of this representation now
bool es_planes_dragging(uint32_t rep_id, int plane);

// The transfer function of a representation's volume rendering (ElectronicStructureRepresentation::volume.tf) baked
// into a 2D texture of 256 x 1 RGBA32F texels, made again when it changes.
uint32_t es_volume_transfer_function(ElectronicStructureRepresentation* es);
