#pragma once

#include <stddef.h>
#include <stdint.h>
#include <vector>

#include <core/md_vec_math.h>
#include <md_gto.h>
#include <qm_placement.h>

struct md_system_t;

// GOSTSHYP: the cavity of a VeloxChem run under pressure, as mdlib publishes it (vlx/scf/gostshyp/*, md_vlx.h), and
// the fields over R^3 it makes.
//
// The cavity is a set of tesserae: points c_j on the scaled van der Waals sphere of one QM atom each, with an outward
// normal n_j, an area a_j = w_j R^2 F_j (its Lebedev weight, the sphere's radius squared, its switching), and a
// normalized Gaussian g_j(r) = (zeta_j / pi)^(3/2) exp(-zeta_j |r - c_j|^2) of amplitude p_j, zeta_j a_j = pi ln 2.
// The amplitudes hold the pressure, p_j f_j = P a_j with f_j = -<rho | 2 zeta_j n_j.(r - c_j) g_j> the decrease of
// the Gaussian's overlap with the density as it moves outward, and the decay length is lambda_j = <rho|g_j> / f_j.
//
// THE FIELDS (VeloxChem's gostshypdriver.py, checked against its source)
//   WallPotential       V1(r) = sum_j p_j g_j(r): the Gaussians as a potential. The energy of the wall is its
//                       expectation value, E = <rho|V1> = sum_j p_j <rho|g_j>. Hartree, >= 0
//   EffectivePotential  v(r) = sum_j p_j g_j(r) (1 + 2 zeta_j lambda_j n_j.(r - c_j)): what the electrons feel. The
//                       amplitudes depend on the density, and the Fock matrix is V1 - V2, V2 the derivative of that
//                       dependence: a p-type term of prefactor p_j lambda_j. To first order each Gaussian moved
//                       outward by its decay length; negative more than 1 / (2 zeta_j lambda_j) inward of a
//                       tessera. Hartree
// and rho(r) V1(r), where the energy of the wall is (integrating to E, Hartree / bohr^3), which takes the density too:
// a surface field of its own (surface_field.h, GostshypEnergyDensity), V1 times the density where the wall reaches.
//
// UNITS: the tesserae as published, in Angstrom; the fields over positions in bohr, as the grids of the volumes.

// What the wall is made of: per tessera its point and normal (Angstrom, the frame of the file: qm/atom/coordinate),
// its Gaussian's amplitude and exponent and its decay length; the QM atoms the frame is of. What else the tesserae
// are (their areas, the wall's energy over them, ...) are attributes of their own, published by mdlib beside these
// (md_vlx.h), and drawn as glyphs (glyphs.h).
struct GostshypTesserae {
    uint32_t num = 0;                   // tesserae
    uint32_t num_atoms = 0;             // QM atoms
    std::vector<double> pos, nrm;       // 3 per tessera: Angstrom (the frame of the file), unit
    std::vector<double> amp;            // Eh bohr^3
    std::vector<double> zeta;           // Angstrom^-2
    std::vector<double> lambda;         // Angstrom
    std::vector<double> qm_xyz;         // 3 per QM atom: Angstrom, where the file had it
};

// Whether the system has a GOSTSHYP cavity
bool gostshyp_available(const md_system_t& sys);

// Reads the wall; false with a log line when it is incomplete
bool gostshyp_read(GostshypTesserae* out, const md_system_t& sys);

// The tesserae placed where the QM atoms are (qm_placement.h): xyz the QM atoms' positions, `scale` of them per Angstrom
typedef QmPlacement GostshypPlacement;
static inline bool gostshyp_place(GostshypPlacement* out, const GostshypTesserae& t, const md_float3_t* xyz, size_t count, const uint32_t* idx, double scale) {
    return qm_place(out, t.qm_xyz.data(), t.num_atoms, xyz, count, idx, scale);
}

enum class GostshypField : int {
    WallPotential,
    EffectivePotential,
    Count
};

// The Gaussians, placed in bohr and binned for evaluation at many points
struct GostshypWall;
GostshypWall* gostshyp_wall_create(const GostshypTesserae& t, const GostshypPlacement& bohr_placement);
void          gostshyp_wall_destroy(GostshypWall* w);

// V1 or v at n points (bohr, packed float3) into out (Hartree). Thread safe.
void gostshyp_wall_evaluate(const GostshypWall* w, GostshypField field, const float* xyz, size_t n, double* out);

// Whether V1 is zero at the point to the cutoff of the Gaussians: nothing of the wall reaches it
bool gostshyp_wall_reaches(const GostshypWall* w, const float xyz[3]);

// A field of the wall as what it exactly is, an orbital over a basis of its own (md_gto.h): a shell on each tessera of
// one primitive, its exponent zeta_j and coefficient 1 (an s or p Cartesian function is not normalized further, so
// the AO is exp(-zeta_j |r - c_j|^2) times 1 or a component of r - c_j). V1 is the s shells with coefficients
// p_j (zeta_j / pi)^(3/2); v adds the p shells with that times 2 zeta_j lambda_j n_j. Exact, no cutoff: what the
// isosurfaces of the field are meshed from (iso_meshes.h). Tesserae of no amplitude (buried) are left out.
// The arrays the basis points into are the struct's own: not to be copied.
struct GostshypOrbital {
    md_gto_basis_t basis = {};
    std::vector<md_gto_shell_t> shells;    // the s shells of all tesserae, then the p shells (md_gto.h's order)
    std::vector<float> alpha, coeff;       // bohr^-2; 1
    std::vector<md_float3_t> centres;      // bohr, the basis' atoms: one per tessera
    std::vector<double> ao;                // the field's coefficient of each AO, Hartree
};
bool gostshyp_wall_orbital(GostshypOrbital* out, const GostshypWall* w, GostshypField field);
