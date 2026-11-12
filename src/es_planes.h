#pragma once

// The planes of the electronic structure representations: the cuts of a volume rendering, and the slices a field is
// drawn on (ElectronicStructureRepresentation::cuts and ::slices, Plane in plane.h). A cut is the plane through its
// representation's box; a slice, a rectangle on its plane of its own (centre, direction, size), anywhere, past the box
// too. Shown in the view: a cut's outline where it cuts the box, a slice's edge along its rectangle (drawn with the
// slices, behind what of them is in front), hidden behind what of the scene is in front of them. A click on an
// outline (or on its row in the representation's panel) picks it for the gizmo, a click on nothing lets go of it
// (main.cpp), Delete removes it (es_planes.cpp): a ring to turn its normal as a trackball, the normal as an arrow to
// slide it along, and a slice's rectangle's handles, to move it along its plane, resize it and turn it about its
// normal. A slice's normal is kept on the camera's side of it: it turns over as the plane is seen edge on.
//
// And what the volume rendering of a representation is drawn with: its planes in the space of its volume, and its
// transfer function baked, the colour map with the opacity curve over it.

#include <core/md_vec.h>
#include <plane.h>

struct ApplicationState;
struct ElectronicStructureRepresentation;
struct ElectronicStructurePlane;
struct Representation;
struct Plane;
struct PlaneFitInfo;


// The box of a representation's volume in world space, the unit cube of its texture: false when it has none yet
bool es_volume_box(md_float4x4_t* box_to_world, const ElectronicStructureRepresentation& es);

// The outlines of the planes of the representations shown, the picking of one by a click on its outline, and the
// gizmo on the picked one. Called in the main view's canvas, before its interaction surface, with the rectangle the
// scene is drawn over (screen coordinates, the whole viewport): true when the mouse is the planes' this frame (over the gizmo, using it, or
// picking an outline), in which case the view is not to turn or pick with it.
bool es_planes_gizmo(ApplicationState* state, md_float2_t canvas_min, md_float2_t canvas_size);

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
// 'info' wanted (nullptr: it is set); kind None when no atoms are selected. out_centre (optional): the atoms' centre,
// where there are any.
Plane es_plane_fit_selection(const ApplicationState* state, const Plane& plane, PlaneFitInfo* info = nullptr, md_float3_t* out_centre = nullptr);

// The enabled cuts of a representation (its volume display) in the space of its volume's texture ([0,1]^3), as n . x + d >= 0 for the
// side kept: what the volume rendering clips to. Returns how many (at most ES_MAX_PLANES).
int es_volume_clip_planes(md_float4_t* out, const ElectronicStructureRepresentation& es);

// A plane of a representation in the space of its volume's texture, as above
md_float4_t es_volume_plane(const ElectronicStructureRepresentation& es, const Plane& plane);

// The slices of a representation (its planes display), every one (at most ES_MAX_PLANES): their planes in the space of
// its volume's texture, their opacities (0 where hidden) and the matrix of each from that space to its rectangle's
// unit square (x, y in [0, 1] on it; es_slice_rect_from_world): what the slice rendering draws. Returns how many.
int es_volume_slices(md_float4_t* out_plane, float* out_opacity, md_float4x4_t* out_rect_from_model, const ElectronicStructureRepresentation& es);

// Gives each slice of a representation its rectangle where it has none yet (the rectangle around where its plane cuts
// the box, once the box is there), and keeps each on its plane as the plane is moved: the centre projected onto it,
// the first side's direction turned into it. Idempotent; called where the slices are read each frame.
void es_slices_conform(ElectronicStructureRepresentation* es);

// A slice's frame, of its plane with the normal one of its two ways (the same frame whichever way it faces: a slice's
// is turned to the camera, which changes nothing of it), and of its rectangle: the centre, u and v = n x u along it,
// and its half sizes. False when it has no rectangle yet.
struct EsSliceFrame {
    md_float3_t n, u, v, c;
    float  hu, hv;
};
bool es_slice_frame(EsSliceFrame* f, const ElectronicStructurePlane& slice);

// World to the rectangle's unit square: x along u, y along v, both [0, 1] over it, z the distance from the plane
// (Angstrom), w 1
md_float4x4_t es_slice_rect_from_world(const EsSliceFrame& f);

// The grid of an image of the rectangle, covering it exactly (its pixels' edges on the rectangle's), with pixels of
// at most 'spacing' Angstrom a side: one layer whose voxel centres are on the plane, in bohr as the evaluations take it
void es_slice_image_grid(struct md_grid_t* grid, const EsSliceFrame& f, float spacing);

// Removes plane i of the representation's display (a cut or a slice; a slice's image with it), keeping the gizmo on
// the plane it was on
void es_planes_remove(ApplicationState* state, Representation* rep, int i);

// Whether the gizmo is dragging this plane of this representation now
bool es_planes_dragging(uint32_t rep_id, int plane);

// The transfer function of a representation's volume rendering (ElectronicStructureRepresentation::volume.tf) baked
// into a 2D texture of 256 x 1 RGBA32F texels, made again when it changes.
uint32_t es_volume_transfer_function(ElectronicStructureRepresentation* es);
