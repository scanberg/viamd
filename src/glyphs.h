#pragma once

#include <stdint.h>
#include <stddef.h>

#include <core/md_unit.h>
#include <system/md_attributes.h>
#include <color_scale.h>

struct md_system_t;
struct md_system_state_t;
struct GBuffer;
struct ViewParam;
struct GlyphData;

// GLYPHS
//
// A glyph per point of a POINT SET of the system's attributes, sized and coloured by the point set's own quantities.
// What can be drawn is decided by the shape of the attributes alone, not by where they came from:
//   a point set is an attribute group G holding G/coordinate, a series (rank 1) of 3-vectors in a length;
//   it is ORIENTED when G/normal is a series of 3-vectors of the same length: then the glyphs can be cylinders;
//   its QUANTITIES are G's other numeric series of one component and the same length (G/atom_index among them).
// The tesserae of a GOSTSHYP cavity are one (vlx/scf/gostshyp/tessera, with what mdlib makes of them: the wall's energy
// per area, the density it sees), a C-PCM cavity's would be another, the QM atoms with their charges a third (spheres).
//
// WHERE THE POINTS ARE: a point set with G/system_index (an atom of the system per point) is at its atoms as they are
// in the state; one with G/atom_index (a QM atom per point) is in the frame of the QM atoms as the file had them
// (qm/atom/coordinate) and follows them rigidly (qm_placement.h); any other as it is.
//
// THE GLYPHS stand on their point along its normal and pick as the point (PickingDomain_Glyph: the tooltip lists the
// point's quantities, and marks it). Two shapes:
//   Sphere     centred on the point: its radius and its colour
//   Cylinder   from the point along its normal: its radius, its height and its colour
// Each from a quantity (or none: a constant), mapped by its UNIT:
//   a length as it is (in Angstrom) times the scale;
//   an area as the radius of a disc of it, sqrt(a / pi), times the scale;
//   anything else as its magnitude over the 99th percentile of the magnitudes, at the scale in Angstrom - a radius by
//   its square root, so that the cross section goes with the value, a height linearly.
// A height keeps its sign: a negative value stands inward. The colour through a colour scale over the quantity's
// spread (the 1st to 99th percentile), symmetric about 0 for one of both signs. A glyph is never thinner than a pixel
// or two, so the smallest can be seen and hovered.
//
// UNITS: positions in Angstrom as the system's.

enum class GlyphShape : int { Cylinder, Sphere, Count };

inline const char* glyph_shape_str[(int)GlyphShape::Count] = {"Cylinder", "Sphere"};

// A size from a quantity: its attribute (MD_ATTRIBUTE_INVALID: a constant), and a scale - a factor for a length or an
// area, else Angstrom
struct GlyphMapping {
    md_attribute_id_t quantity = MD_ATTRIBUTE_INVALID;
    float scale = 0.15f;
};

struct GlyphRepresentation {
    md_attribute_id_t points = MD_ATTRIBUTE_INVALID;  // the point set: its G/coordinate
    GlyphShape   shape = GlyphShape::Cylinder;
    GlyphMapping radius = {MD_ATTRIBUTE_INVALID, 0.15f};
    GlyphMapping height = {MD_ATTRIBUTE_INVALID, 0.1f};    // a cylinder's
    md_attribute_id_t color = MD_ATTRIBUTE_INVALID;        // MD_ATTRIBUTE_INVALID: uniform
    ColorScale   scale = {COLOR_SCALE_COLORMAP_VIRIDIS, 0.0f, 1.0f, false, true, false};
    bool         unlit = false;     // the colour as it is, not shaded: shaded by default, the light shows the shape

    GlyphData*   data = nullptr;    // owned, made by glyphs_update
};

// The point sets of the system (their G/coordinate), in path order; returns how many there are, writes at most cap
size_t glyphs_point_sets(md_attribute_id_t out[], size_t cap, const md_system_t& sys);

// Whether id is a point set, and whether it is oriented
bool glyphs_is_point_set(const md_system_t& sys, md_attribute_id_t points);
bool glyphs_is_oriented(const md_system_t& sys, md_attribute_id_t points);

// Its quantities, in path order; returns how many there are, writes at most cap
size_t glyphs_quantities(md_attribute_id_t out[], size_t cap, const md_system_t& sys, md_attribute_id_t points);

// The point set shown (an oriented one first), and what it is mapped by, as a new representation starts out: a disc of
// each point's area where it has one, of a constant height, coloured by the set's principal quantity, what its
// publisher says it is about (MD_ATTRIBUTE_FLAG_PRINCIPAL: a GOSTSHYP cavity's decay length, a surface charge), or
// else uniform
void glyphs_select(GlyphRepresentation* g, const md_system_t& sys, md_attribute_id_t points);

// The colour by a quantity (or none), the scale set up for its values: symmetric and diverging for both signs
void glyphs_select_color(GlyphRepresentation* g, const md_system_t& sys, md_attribute_id_t quantity);

// A scale for a size from a quantity, as the mapping is switched to it
float glyphs_default_scale(const md_system_t& sys, md_attribute_id_t quantity, bool radius);

// Whether a quantity's scale is a factor (a length, an area) rather than Angstrom
bool glyphs_scale_is_factor(const md_system_t& sys, md_attribute_id_t quantity);

// Brings the representation up to date with the system and the frame: reads the points when they are new, sizes the
// glyphs when their mapping changed, and places them where they are in the state. Cheap when nothing changed. Returns
// false when there is nothing to draw.
bool glyphs_update(GlyphRepresentation* g, const md_system_t& sys, const md_system_state_t& state);

// Forgets what was read, for a system loaded under the representation
void glyphs_invalidate(GlyphRepresentation* g);
void glyphs_free(GlyphRepresentation* g);

// What the colour's quantity spans over the points (the 1st to 99th percentile), for the scale
ColorScaleSpan glyphs_span(const GlyphRepresentation& g);

// A hash of what the image depends on, for the lazy rendering
uint64_t glyphs_hash(const GlyphRepresentation& g);

// The picking indices it needs: one per point
uint32_t glyphs_count(const GlyphRepresentation& g);

// The glyphs into the G-buffer bound for drawing (the opaque pass), picking as the points (base + the index of each)
void glyphs_draw(const GlyphRepresentation& g, const ViewParam& view, uint32_t base);

// The glyph of a point into the marks of the highlight (the marks target bound, gbuffer_marks_begin): there where it
// is, in view where nothing is in front of it
void glyphs_draw_marks(const GlyphRepresentation& g, uint32_t point, const GBuffer& gbuf, const ViewParam& view);

// The system's atom a point is on or of (G/system_index, or G/atom_index through qm/atom/system_index), or UINT32_MAX
uint32_t glyphs_point_atom(const md_system_t& sys, md_attribute_id_t points, uint32_t point);

// The value of an attribute (rank 1) at a point, in its own unit; NAN where it has none
double glyphs_value(const md_system_t& sys, md_attribute_id_t quantity, uint32_t point);

// How a quantity's values are shown: a geometric unit (a length, an area, an inverse area: no other dimension) in the
// display units (display_units.h), any other as published - its publisher chose it (Ha/bohr^2, e/bohr^3, kJ/mol), and
// the display units would rebuild a compound of an energy from mass and time. Returns whether it follows the display
// units; *factor (may be NULL) takes a value to how it is shown, buf its unit's spelling.
bool glyphs_unit_print(char* buf, size_t cap, double* factor, md_unit_t unit);

// The group a point set is: G of G/coordinate (empty when it is none). A view into the attribute table.
md_str_t glyphs_group(const md_system_t& sys, md_attribute_id_t points);
