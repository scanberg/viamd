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
// A tessellated cavity around the QM atoms: points (tesserae) on the union of scaled van der Waals spheres, each with
// a normal, an area and quantities of its own. Today that is the GOSTSHYP cavity of a VeloxChem run under pressure
// (vlx/scf/gostshyp/*, md_vlx.h): a normalized 3D Gaussian on each tessera, with its exponent zeta_j (zeta_j a_j is
// constant), its amplitude p_j and the decay length lambda_j = <rho|g_j> / (-d/dn <rho|g_j>). The amplitudes hold
// the pressure: p_j (-d/dn <rho|g_j>) = P a_j, so the energy of the wall at a tessera, p_j <rho|g_j>, is P a_j
// lambda_j, and its energy per area is P lambda_j.
//
// HOW IT IS DRAWN
// The spheres are ray cast into the G-buffer by a pass of their own (cavity_draw_geometry), whose fragments pick as the
// sphere they are on (PickingDomain_CavitySphere): their union. The surface is GOSTSHYP's own where it is switched:
// VeloxChem's 'swig' weights the points of each sphere by a smoothstep across a band about each neighbouring sphere's
// surface, so that the areas of the tesserae are their Lebedev weights times that switching. The surface is the union
// where the switching is at least 1/2 (on 5helicene.h5 306.7 Å^2 against the tesserae's 312.5, the union 361.5); the
// bands along the seams where it is less are drawn neutral and are of no tessera. They are not left open: F >= 1/2 is
// not a closed surface, and through the bands one would see into the cavity. A full screen pass after the opaque geometry (cavity_resolve) then
// looks at every pixel of the G-buffer that a sphere of the cavity won, intersects the view ray with that sphere again
// for the exact point, finds the tesserae near it in a uniform grid, and writes the colour and what the pixel picks
// as: the tessera nearest to the point (PickingDomain_Tessera). The surface is split exactly per pixel between the
// tesserae, with no map in between. The tessera hovered is marked as its patch by a second full screen pass into the
// marks (cavity_draw_marks).
//
// COLOURING
//   Nearest  each tessera's own value over its patch: the points nearest to it on the surface, among the tesserae
//            facing the same way as the surface there. Nearest by power, |x - c_j|^2 - w_j, with the weights fitted
//            so that each patch has the area of its tessera (fit_areas), or by distance alone (a Voronoi split, where
//            the tesserae the switching kept inside other spheres crowd the seams with slivers)
//   Kernel   the values spread over the surface by a Gaussian kernel per tessera, of the width of its own Gaussian
//            (sigma_j = 1 / sqrt(2 zeta_j), times kernel_scale), normalized over the exposed surface (the surface
//            drawn: of the union of the spheres, where the switching is at least 1/2 when switched):
//              K_j(x) = exp(-|x - c_j|^2 / 2 sigma_j^2) where n(x).n_j > 0 and |x - c_j| < 3 sigma_j, else 0
//              N_j    = integral of K_j over the exposed surface (by quadrature on the spheres)
//            A quantity that is a total of the tessera (the energy of the wall at it, its amplitude, its area) is
//            spread conservatively, f(x) = sum_j v_j K_j(x) / N_j, so that the integral of f over the surface is
//            sum_j v_j exactly. One per tessera (a decay length, a density) is averaged by area,
//            f(x) = sum_j a_j v_j K_j(x) / N_j / sum_j a_j K_j(x) / N_j.
//
// UNITS: positions in Ångström as the system's.

enum class CavityField : int {
    WallEnergy,         // energy of the wall per area, P lambda_j: kJ/mol/Å^2 (conservative)
    DecayLength,        // lambda_j: Å (by area)
    WallDensity,        // <rho|g_j> = P a_j lambda_j / p_j: e/bohr^3 (by area)
    Amplitude,          // amplitude per area p_j / a_j (conservative)
    AreaCoverage,       // the areas of the tesserae spread over the surface: 1 where they account for it (conservative)
    Count
};

enum class CavityColoring : int {
    Nearest,
    Kernel,
    Count
};

inline const char* cavity_field_str[(int)CavityField::Count] = {
    "Wall energy per area",
    "Decay length",
    "Wall density",
    "Amplitude per area",
    "Area coverage",
};

// The units of the fields as they are shown. An energy per area is a force per length to md_unit (J/m^2 = kg/s^2),
// which it would print as that: the fields carry their units in words, and no md_unit.
inline const char* cavity_field_unit_str[(int)CavityField::Count] = {
    "kJ/mol/Å²",
    "Å",
    "e/bohr³",
    "a.u./Å²",
    "",
};

inline const char* cavity_coloring_str[(int)CavityColoring::Count] = {
    "Nearest tessera",
    "Kernel",
};

struct CavityRepresentation {
    CavityField    field    = CavityField::WallEnergy;
    CavityColoring coloring = CavityColoring::Kernel;
    float kernel_scale = 2.0f;  // of the width of each tessera's Gaussian: at 1 the tesserae are two widths apart and
                                // the spread ripples between them; at 2, one
    bool  conserve = true;      // Kernel: totals spread conservatively (else averaged by area, as values per tessera)
    bool  unlit = false;        // the colour as it is, not shaded. Shaded by default: the light and the occlusion are
                                // what show the shape of the surface, and they leave its colours readable
    bool  edges = true;         // Nearest: the borders between the patches
    bool  switched = true;      // GOSTSHYP's own surface: each sphere where the switching of its points is at least 1/2,
                                // not the union of the spheres (when the file's discretization is 'swig')
    bool  fit_areas = true;     // Split the surface by a power diagram whose patches have the areas of the tesserae
                                // (relative to the surface drawn), not by the nearest tessera alone
    ColorScale scale = {COLOR_SCALE_COLORMAP_VIRIDIS, 0.0f, 1.0f, false, true, false};

    CavityData* data = nullptr; // owned, made by cavity_update
};

// Whether the system has a cavity to draw
bool cavity_available(const md_system_t& sys);

// Brings the representation up to date with the system and the frame: reads the tesserae when they are new, makes
// the normalization of the kernel when its width changed, the values of the field when it changed, and places the
// cavity on the atoms as they are in the state (rigidly: the tesserae are of the geometry of the file). Cheap when
// nothing changed. Returns false when there is nothing to draw.
bool cavity_update(CavityRepresentation* cav, const md_system_t& sys, const md_system_state_t& state);

// Forgets what was read, for a system loaded under the representation
void cavity_invalidate(CavityRepresentation* cav);
void cavity_free(CavityRepresentation* cav);

// What the field spans on the surface (the 1st to 99th percentile of it over the exposed surface), for the scale
ColorScaleSpan cavity_span(const CavityRepresentation& cav);
md_unit_t      cavity_field_unit(CavityField field);  // Å for the decay length, none for the others (cavity_field_unit_str)

// A hash of what the image of the cavity depends on, for the lazy rendering
uint64_t cavity_hash(const CavityRepresentation& cav);

// The picking indices it needs: one per sphere and one per tessera
uint32_t cavity_sphere_count(const CavityRepresentation& cav);
uint32_t cavity_tessera_count(const CavityRepresentation& cav);

// The spheres into the G-buffer bound for drawing (the opaque pass): ray cast, depth tested (their union), picking as
// their spheres (sphere_base + the QM atom of each)
void cavity_draw_geometry(const CavityRepresentation& cav, const ViewParam& view, uint32_t sphere_base);

// Whether the cavity can be drawn as GOSTSHYP's switched surface (the file's discretization is 'swig'), and is
bool cavity_has_switching(const CavityRepresentation& cav);
bool cavity_switched(const CavityRepresentation& cav);

// The switching at the point of the sphere of a QM atom in the direction of a point (world space); 1 when not switched
float cavity_sphere_switching(const CavityRepresentation& cav, uint32_t sphere, const float world[3]);

// The full screen pass over the G-buffer bound for drawing: where a sphere of this cavity was drawn (its picking
// index in [sphere_base, sphere_base + cavity_sphere_count)), writes the colour and the picking index of the tessera
// (tessera_base + its index). Leaves the depth, the normals and the velocities as they are. The picking target is
// read from a copy made first, so the G-buffer is written in place.
void cavity_resolve(const CavityRepresentation& cav, const GBuffer& gbuf, const ViewParam& view, uint32_t sphere_base, uint32_t tessera_base);

// The patch of a tessera into the marks of the highlight (the marks target bound, gbuffer_marks_begin): there where it
// is, in view where nothing is in front of it
void cavity_draw_marks(const CavityRepresentation& cav, uint32_t tessera, const GBuffer& gbuf, const ViewParam& view);

// What a tessera is, for the tooltip
struct CavityTesseraInfo {
    uint32_t atom = 0;          // the QM atom whose sphere it is on
    uint32_t sys_atom = 0;      // and that atom in the system
    float area = 0;             // Å^2
    float decay_length = 0;     // Å
    float wall_energy = 0;      // kJ/mol/Å^2
    float wall_density = 0;     // e/bohr^3
    float amplitude = 0;        // a.u.
    float field_here = 0;       // the field shown, where the cursor is
    bool  has_field_here = false;
};
bool cavity_tessera_info(CavityTesseraInfo* info, const CavityRepresentation& cav, uint32_t tessera, const float* world_pos);

// What the cavity is as a whole, for the representation's properties
struct CavitySummary {
    float pressure_gpa = 0;
    uint32_t num_tesserae = 0;
    uint32_t num_spheres = 0;
    float tessera_area = 0;     // sum a_j, Å^2
    float exposed_area = 0;     // of the surface drawn, Å^2
    bool  switched = false;     // which: GOSTSHYP's switched surface, or the union of the spheres
    float wall_energy = 0;      // sum_j P a_j lambda_j, kJ/mol
    float max_kernel_error = 0; // largest relative difference of the integral of a conservative field from its sum
    uint32_t num_buried = 0;    // tesserae inside the sphere of another atom: the switching keeps them where spheres overlap
    float buried_area = 0;      // their share of sum a_j
    float unreached_energy = 0; // the share of the wall energy on tesserae no point of the exposed surface is within reach of
};
bool cavity_summary(CavitySummary* out, const CavityRepresentation& cav);
