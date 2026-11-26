#pragma once

#include <stdint.h>
#include <stddef.h>

#include <system/md_attributes.h>
#include <color_scale.h>

struct md_system_t;
struct md_system_state_t;
struct ViewParam;
struct VdwSurfaceData;

// VDW SURFACE
//
// The surface of a union of spheres, one per atom the representation shows: the van der Waals surface, or a solvation
// cavity grown from it (a GOSTSHYP or C-PCM cavity: the spheres the file gives the QM atoms, 1.5 and 1.2 times the
// van der Waals radii in VeloxChem). Exact per pixel: each sphere ray cast, the depth test keeping what is outside the
// others. It picks as the atoms.
//
// RADII: the atoms' van der Waals radii (md_system's, with the element remappings) times the scale, or a radius series
// of the QM atoms (a cavity's G/sphere_radius, vdw_surface_radii_sets), in which case only the QM atoms that have one.
//
// COLOUR: the atoms' colours (the representation's colouring), or a QUANTITY of a point set on the surface (glyphs.h:
// a cavity's tesserae, with their decay length, GOSTSHYP energy per area, surface charge, ...), in two steps:
//   SMOOTHING, on the surface: the values diffused over the points for a diffusion length l (the width of the heat
//     kernel, sqrt(2 D t)), along a graph of the points next to each other (each with those within twice its spacing),
//     the fluxes weighted by the areas: sum_j a_j v_j is kept, to rounding. The integral of a density (the GOSTSHYP
//     energy of its energy per area) is the same however smooth; nothing leaks across a gap. 0: the values as they are.
//   DISPLAY, at a pixel x: the values under the points' own Gaussians, as the solvation model smears its tesserae,
//         f(x) = sum_j a_j v_j exp(-w_j |x - r_j|^2) / sum_j a_j exp(-w_j |x - r_j|^2)
//     a_j the point's area (G/area, 1 without), w_j its exponent (G/exponent, a GOSTSHYP tessera's omega_j; without,
//     1 / (2 h^2) of its spacing h: the square root of its area, else the mean distance to the nearest point). It reaches about a point and no
//     further: continuous, conservative (within its neighbours' values), seamless across the spheres.
// Points of no area (a tessera its switching buried) take part in neither. What a point set marks as a total of its
// point rather than a density (MD_ATTRIBUTE_FLAG_EXTENSIVE: an energy, a charge) is not offered: on the surface its
// pattern would be that of the areas. Nor is what is of the discretization (MD_ATTRIBUTE_FLAG_SAMPLING: a quadrature
// weight, a switching, an exponent from the area), whose pattern would be that of the tessellation.
//
// UNITS: Angstrom, as the system's.

struct VdwSurfaceRepresentation {
    float scale = 1.0f;                                     // of the van der Waals radii
    md_attribute_id_t radii    = MD_ATTRIBUTE_INVALID;      // a radius series of the QM atoms in its place
    md_attribute_id_t quantity = MD_ATTRIBUTE_INVALID;      // a point set's quantity it is coloured by; else the atoms'
    ColorScale color_scale = {COLOR_SCALE_COLORMAP_VIRIDIS, 0.0f, 1.0f, false, true, false};
    float smoothing = 0.0f;                                 // the diffusion length over the surface, Angstrom

    VdwSurfaceData* data = nullptr;                         // owned, made by vdw_surface_update
};

// The radius series of the QM atoms the system has (attributes G/sphere_radius: a length per QM atom), in path order;
// returns how many there are, writes at most cap
size_t vdw_surface_radii_sets(md_attribute_id_t out[], size_t cap, const md_system_t& sys);

// The quantities it can be coloured by: those of the point sets of the system (glyphs_quantities) but their indices
// and what is marked a total (MD_ATTRIBUTE_FLAG_EXTENSIVE) or of the discretization (MD_ATTRIBUTE_FLAG_SAMPLING)
size_t vdw_surface_quantities(md_attribute_id_t out[], size_t cap, const md_system_t& sys);

// As a new one of a system with a cavity starts out: its radii, coloured by its principal quantity (the decay length
// of a GOSTSHYP cavity, MD_ATTRIBUTE_FLAG_PRINCIPAL). Returns whether radii was one.
bool vdw_surface_select_cavity(VdwSurfaceRepresentation* v, const md_system_t& sys, md_attribute_id_t radii);

// The colour by a quantity (or none: the atoms'), its scale set up for its values
void vdw_surface_select_quantity(VdwSurfaceRepresentation* v, const md_system_t& sys, md_attribute_id_t quantity);

// The atoms' colours, as the representation computed them for the atoms it shows (alpha 0: not shown)
void vdw_surface_set_colors(VdwSurfaceRepresentation* v, const uint32_t* colors, size_t count);

// Up to date with the system and the frame: the spheres where the atoms are, the points where they are and their
// values. Cheap when nothing changed. Returns false when there is nothing to draw.
bool vdw_surface_update(VdwSurfaceRepresentation* v, const md_system_t& sys, const md_system_state_t& state);

void vdw_surface_invalidate(VdwSurfaceRepresentation* v);
void vdw_surface_free(VdwSurfaceRepresentation* v);

// What the colour's quantity spans over its points, for the scale
ColorScaleSpan vdw_surface_span(const VdwSurfaceRepresentation& v);

// A hash of what the image depends on, for the lazy rendering
uint64_t vdw_surface_hash(const VdwSurfaceRepresentation& v);

// Into the G-buffer bound for drawing (the opaque pass), picking as the atoms (atom_base + the atom's index)
void vdw_surface_draw(const VdwSurfaceRepresentation& v, const ViewParam& view, uint32_t atom_base);
