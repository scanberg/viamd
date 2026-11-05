#pragma once

// The planes of the electronic structure representations: the cuts of a volume rendering, and the plane a field is
// drawn on (ElectronicStructureRepresentation::planes, Plane in plane.h). Shown in the view as outlines where they
// cut their representation's box, a click on an outline picking it for the gizmo (es_planes.cpp): a ring to turn
// its normal as a trackball, and the normal as an arrow to slide it along.
//
// And what the volume rendering of a representation is drawn with: its planes in the space of its volume, and its
// transfer function baked, the colour map with the opacity curve over it.

#include <core/md_vec_math.h>

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
int es_volume_slices(vec4_t* out, float* out_opacity, const ElectronicStructureRepresentation& es);

// The transfer function of a representation's volume rendering (ElectronicStructureRepresentation::volume.tf) baked
// into a 2D texture of 256 x 1 RGBA32F texels, made again when it changes.
uint32_t es_volume_transfer_function(ElectronicStructureRepresentation* es);
