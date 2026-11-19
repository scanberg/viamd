#pragma once

#include <stddef.h>
#include <stdint.h>
#include <vector>

#include <core/md_vec_math.h>

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

struct GostshypTesserae {
    uint32_t num = 0;                   // tesserae
    uint32_t num_atoms = 0;             // QM atoms, one sphere each
    std::vector<double> pos, nrm;       // 3 per tessera: Angstrom (the frame of the file), unit
    std::vector<double> area;           // Angstrom^2, switched
    std::vector<double> amp;            // Eh bohr^3
    std::vector<double> zeta;           // Angstrom^-2
    std::vector<double> lambda;         // Angstrom
    std::vector<double> switching;      // F_j
    std::vector<double> lebedev;        // w_j, the point's solid angle: 4 pi over a sphere
    std::vector<uint32_t> atom;         // the QM atom of each tessera
    std::vector<double> qm_xyz;         // 3 per QM atom: Angstrom, where the file had it
    std::vector<uint32_t> sys_idx;      // the system's atom of each QM atom
    std::vector<double> radius;         // per QM atom, Angstrom
    double pressure_gpa = 0.0;
};

// Whether the system has a GOSTSHYP cavity, and a hash of what it is read from
bool     gostshyp_available(const md_system_t& sys);
uint64_t gostshyp_source_hash(const md_system_t& sys);

// Reads the cavity; false with a log line when it is incomplete
bool gostshyp_read(GostshypTesserae* out, const md_system_t& sys);

// The wall energy per area at a tessera, P lambda_j (kJ/mol/Angstrom^2), its energy P a_j lambda_j (kJ/mol), and the
// density its Gaussian sees, <rho|g_j> = P a_j lambda_j / p_j (e/bohr^3)
double gostshyp_energy_per_area(const GostshypTesserae& t, uint32_t j);
double gostshyp_energy(const GostshypTesserae& t, uint32_t j);
double gostshyp_wall_density(const GostshypTesserae& t, uint32_t j);

// The rigid motion taking the QM atoms where the file had them onto where they are: x = R x_file + t (column major R).
// xyz holds the positions, `scale` of them per Angstrom (1: Angstrom, 1 / 0.529...: bohr); idx the entry of each QM
// atom in xyz (NULL: the QM atom's own index). The identity when they have not moved.
struct GostshypPlacement {
    double R[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    double t[3] = {0, 0, 0};
    double scale = 1.0;                 // the units of x per Angstrom: R and t act on the file's positions so scaled
};
bool gostshyp_place(GostshypPlacement* out, const GostshypTesserae& t, const md_f32x3_t* xyz, size_t count, const uint32_t* idx, double scale);

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
