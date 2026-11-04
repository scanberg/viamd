#pragma once

// The planes of the electronic structure representations: the cuts of a volume rendering, and the plane a field is
// drawn on (ElectronicStructureRepresentation::planes, Plane in viamd.h). Shown in the view as outlines where they
// cut their representation's box, a click on an outline picking it for the gizmo, which turns it about the axes in
// it and moves it along its normal.
//
// And what the volume rendering of a representation is drawn with: its planes in the space of its volume, and its
// transfer function, the colour map with an alpha ramp.

#include <core/md_vec_math.h>

struct ApplicationState;
struct ElectronicStructureRepresentation;
struct Plane;


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
// system is periodic. False, and the plane left, when no atoms are selected.
bool es_plane_fit_selection(const ApplicationState* state, Plane* plane);

// The enabled planes of a representation in the space of its volume's texture ([0,1]^3), as n . x + d >= 0 for the
// side kept: what the volume rendering clips to. Returns how many (at most ES_MAX_PLANES).
int es_volume_clip_planes(vec4_t* out, const ElectronicStructureRepresentation& es);

// Plane i of a representation in the space of its volume's texture, as above, whether it is enabled or not
vec4_t es_volume_plane(const ElectronicStructureRepresentation& es, int i);

// The transfer function of a representation's volume rendering, its colour map with the alpha ramp
// (ElectronicStructureRepresentation::volume), as a 2D texture of 256 x 1 RGBA32F texels, remade when they change.
uint32_t es_volume_transfer_function(ElectronicStructureRepresentation* es);
