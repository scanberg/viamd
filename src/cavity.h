#pragma once

#include <stdint.h>
#include <stddef.h>

#include <core/md_unit.h>
#include <color_scale.h>

struct md_system_t;
struct md_system_state_t;
struct GBuffer;
struct ViewParam;
struct CavityData;

// CAVITY
//
// The tesserae of the GOSTSHYP cavity of a VeloxChem run under pressure (gostshyp.h), drawn as what they are: points
// with a normal and quantities of their own, one glyph each. The fields the wall makes over R^3 - its potential, the
// potential the electrons feel, its energy density - are drawn by the electronic structure representation (a source
// of its own, and fields an isosurface is coloured by): there is no surface of the cavity here, as GOSTSHYP has none.
//
// THE GLYPHS stand on their tessera's point along its normal and pick as the tessera (PickingDomain_Tessera): hovering
// one shows its tooltip and marks it, and follows it in the GOSTSHYP window's table. Two shapes:
//   Sphere     centred on the point: its radius and its colour
//   Cylinder   from the point outward: its radius, its height and its colour
// Each of those is mapped from a quantity of the tessera (CavityAttribute), chosen per dimension:
//   a size from a length as it is (the decay length lambda_j; the area as the radius sqrt(a_j / pi), where the
//   tessera's Gaussian is at half its peak - zeta_j a_j = pi ln 2 - and the disc has the tessera's area), times the
//   scale; from another quantity as its value over the 99th percentile of the tesserae's, at the scale in Angstrom -
//   a radius by its square root, so that the cross section goes with the value, a height linearly; or constant, the
//   scale in Angstrom
//   the colour from the quantity through a colour scale over its spread (the 1st to 99th percentile)
// By default a disc of the tessera's area, of a constant height, coloured by the wall's energy per area. A glyph is
// never thinner than a pixel or two, so the smallest tesserae can be seen and hovered.
//
// UNITS: positions in Angstrom as the system's.

enum class CavityGlyphShape : int { Cylinder, Sphere, Count };

enum class CavityAttribute : int {
    Constant,
    Area,               // a_j, Angstrom^2 (switched: w_j R^2 F_j); as a radius, sqrt(a_j / pi)
    DecayLength,        // lambda_j, Angstrom
    WallEnergyPerArea,  // P lambda_j, kJ/mol/Angstrom^2
    WallEnergy,         // P a_j lambda_j, kJ/mol
    WallDensity,        // <rho|g_j> = P a_j lambda_j / p_j, e/bohr^3
    Amplitude,          // p_j, a.u. (Eh bohr^3)
    AmplitudePerArea,   // p_j / a_j, a.u./Angstrom^2
    Switching,          // F_j
    LebedevWeight,      // w_j, the point's solid angle (4 pi over a sphere)
    Count
};

inline const char* cavity_glyph_shape_str[(int)CavityGlyphShape::Count] = {"Cylinder", "Sphere"};

inline const char* cavity_attribute_str[(int)CavityAttribute::Count] = {
    "Constant",
    "Area",
    "Decay length",
    "Wall energy per area",
    "Wall energy",
    "Wall density",
    "Amplitude",
    "Amplitude per area",
    "Switching",
    "Lebedev weight",
};

// Their units in words: an energy per area is a force per length to md_unit, which it would print as kg/s^2
inline const char* cavity_attribute_unit_str[(int)CavityAttribute::Count] = {
    "",
    "Å²",
    "Å",
    "kJ/mol/Å²",
    "kJ/mol",
    "e/bohr³",
    "a.u.",
    "a.u./Å²",
    "",
    "sr",
};

// Whether a quantity is a length as a size: the area as the radius of its disc, the decay length. Its scale is then a
// factor; for the others a size in Angstrom at the 99th percentile of their values
static inline bool cavity_attribute_is_length(CavityAttribute a) {
    return a == CavityAttribute::Area || a == CavityAttribute::DecayLength;
}

struct CavityMapping {
    CavityAttribute attribute = CavityAttribute::Constant;
    float scale = 1.0f;         // a factor for a length (cavity_attribute_is_length), else Angstrom
};

struct CavityRepresentation {
    CavityGlyphShape shape = CavityGlyphShape::Cylinder;
    CavityMapping radius = {CavityAttribute::Area, 1.0f};
    CavityMapping height = {CavityAttribute::Constant, 0.1f};   // a cylinder's
    CavityAttribute color = CavityAttribute::WallEnergyPerArea;
    ColorScale scale = {COLOR_SCALE_COLORMAP_VIRIDIS, 0.0f, 1.0f, false, true, false};
    bool  unlit = false;        // the colour as it is, not shaded: shaded by default, the light shows the shape

    CavityData* data = nullptr; // owned, made by cavity_update
};

// Whether the system has a cavity to draw
bool cavity_available(const md_system_t& sys);

// Brings the representation up to date with the system and the frame: reads the tesserae when they are new, sizes the
// glyphs when their mapping changed, and places them on the atoms as they are in the state (rigidly: the tesserae are
// of the geometry of the file). Cheap when nothing changed. Returns false when there is nothing to draw.
bool cavity_update(CavityRepresentation* cav, const md_system_t& sys, const md_system_state_t& state);

// Forgets what was read, for a system loaded under the representation
void cavity_invalidate(CavityRepresentation* cav);
void cavity_free(CavityRepresentation* cav);

// What the quantity coloured by spans over the tesserae (the 1st to 99th percentile), for the scale
ColorScaleSpan cavity_span(const CavityRepresentation& cav);

// A hash of what the image of the cavity depends on, for the lazy rendering
uint64_t cavity_hash(const CavityRepresentation& cav);

// The picking indices it needs: one per tessera
uint32_t cavity_tessera_count(const CavityRepresentation& cav);

// The glyphs into the G-buffer bound for drawing (the opaque pass), picking as the tesserae (tessera_base + the index
// of each)
void cavity_draw_glyphs(const CavityRepresentation& cav, const ViewParam& view, uint32_t tessera_base);

// The glyph of a tessera into the marks of the highlight (the marks target bound, gbuffer_marks_begin): there where it
// is, in view where nothing is in front of it
void cavity_draw_marks(const CavityRepresentation& cav, uint32_t tessera, const GBuffer& gbuf, const ViewParam& view);

// The value of a quantity at a tessera, in its unit (cavity_attribute_unit_str); 1 for Constant
double cavity_attribute_value(const CavityRepresentation& cav, CavityAttribute a, uint32_t tessera);

// What a tessera is, for the tooltip
struct CavityTesseraInfo {
    uint32_t atom = 0;          // the QM atom whose sphere it is on
    uint32_t sys_atom = 0;      // and that atom in the system
    float area = 0;             // Angstrom^2
    float decay_length = 0;     // Angstrom
    float wall_energy = 0;      // kJ/mol/Angstrom^2
    float wall_density = 0;     // e/bohr^3
    float amplitude = 0;        // a.u.
    float switching = 0;        // F_j
};
bool cavity_tessera_info(CavityTesseraInfo* info, const CavityRepresentation& cav, uint32_t tessera);

// What the cavity is as a whole, for the representation's properties
struct CavitySummary {
    float pressure_gpa = 0;
    uint32_t num_tesserae = 0;
    uint32_t num_spheres = 0;
    float tessera_area = 0;     // sum a_j, Angstrom^2
    float wall_energy = 0;      // sum_j P a_j lambda_j, kJ/mol
    float effective_volume = 0; // sum_j a_j lambda_j, Angstrom^3: the wall energy over the pressure
    uint32_t num_buried = 0;    // tesserae inside the sphere of another atom: the switching keeps them where spheres overlap
    float buried_area = 0;      // their share of sum a_j
};
bool cavity_summary(CavitySummary* out, const CavityRepresentation& cav);
