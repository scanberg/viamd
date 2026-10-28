#include <surface_field.h>

#include <task_system.h>
#include <gfx/gl.h>
#include <gfx/gl_utils.h>

#include <md_system.h>
#include <md_attributes.h>
#include <md_gto.h>
#include <md_gto_int.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <atomic>
#include <math.h>
#include <string.h>

#define ANGSTROM_TO_BOHR_D 1.8897261246257702

// The spacing the electrostatic potential is sampled at, bohr. See surface_field_grid.
#define ESP_GRID_SPACING 0.4

// Screening of the QM density's gaussians: one whose potential cannot exceed this anywhere is
// dropped (hartree/e). That costs a few 1e-6 au on a molecule of tens of atoms and ~1e-4 au on C60
// (md_gto_int.h), against a spread of ~0.1 au over an isodensity surface, and drops a quarter to a
// third of the gaussians.
#define ESP_SCREENING_THRESHOLD 1.0e-6

static const str_t ESP_TOTAL_DENSITY_PATH    = STR_INIT("orbital/total/density");
static const str_t ESP_ALPHA_DENSITY_PATH    = STR_INIT("orbital/alpha/density");
static const str_t ESP_BETA_DENSITY_PATH     = STR_INIT("orbital/beta/density");
static const str_t ESP_NUCLEAR_CHARGE_PATH   = STR_INIT("qm/atom/nuclear_charge");
static const str_t ESP_MOLECULAR_CHARGE_PATH = STR_INIT("vlx/molecular_charge");

// Read on the bits: viamd is built with fast math, where NAN need not compare as itself
static inline bool is_finite_f64(double v) {
    uint64_t u;
    MEMCPY(&u, &v, sizeof(u));
    return (u & 0x7FF0000000000000ull) != 0x7FF0000000000000ull;
}

md_unit_t surface_field_unit(SurfaceFieldKind kind) {
    switch (kind) {
    case SurfaceFieldKind::CpcmChargeDensity:
        return md_unit_div(md_unit_elementary_charge(), md_unit_pow(md_unit_bohr_radius(), 3));
    case SurfaceFieldKind::GostshypPotential:
        return md_unit_hartree();
    case SurfaceFieldKind::EmbeddingPotential:
    case SurfaceFieldKind::ElectrostaticPotential:
    case SurfaceFieldKind::CpcmPotential:
    default:
        return md_unit_div(md_unit_hartree(), md_unit_elementary_charge());
    }
}

// ---------------------------------------------------------------------------
// Embedding potential: the environment's classical sites
// ---------------------------------------------------------------------------

// The charges of the system's atoms, NAN for an atom that has none
static const md_attribute_t* find_atom_charges(const md_system_t& sys) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, STR_LIT("atom/charge"));
    if (!attr || md_attribute_element_count(&attr->format) != md_system_atom_count(&sys)) {
        return nullptr;
    }
    return attr;
}

// The embedding potential is the environment's: the QM region is what the surfaces are drawn from,
// never a source of the field on them. A QM atom carries a charge too when the system came from a
// topology a QM calculation supplements (the force field's), and that charge sits inside the very
// surface it would colour, so the QM flag excludes an atom whatever it carries.
static inline bool atom_in_environment(const md_system_t& sys, size_t i) {
    return (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_QM) == 0;
}

// Whether any atom is outside the QM region: what the 'environment' keyword requires as well
static bool has_environment(const md_system_t& sys) {
    const size_t num_atoms = md_system_atom_count(&sys);
    for (size_t i = 0; i < num_atoms; ++i) {
        if (atom_in_environment(sys, i)) return true;
    }
    return false;
}

// A per atom column of 'components' values in 'unit', or nothing when the system has none
static bool extract_atom_column(double* dst, const md_system_t& sys, str_t path, uint32_t components, md_unit_t unit) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, path);
    const size_t num_atoms = md_system_atom_count(&sys);
    if (!attr || attr->format.components != components || md_attribute_value_count(&attr->format) != num_atoms) {
        return false;
    }
    const size_t n = num_atoms * components;
    return md_attribute_extract_f64(dst, n, attr, md_attribute_slice_all(), unit) == n;
}

// The dipoles a polarizable embedding induced on the environment's sites, over the system's atoms
// (md_vlx publishes them for a run that kept them): zero on a site that is not polarizable, NAN on
// the QM atoms. Part of the embedding potential - the QM region's Fock operator has the field of
// the induced dipoles in it beside that of the charges, as bare point dipoles - so they are added to
// the sites' permanent dipoles rather than offered as a field of their own: the permanent part alone
// is not a potential the QM density was ever in, and not that of a non-polarizable run either.
static const str_t EMBEDDING_INDUCED_DIPOLE_PATH = STR_INIT("atom/induced_dipole");

// Environment sites whose `components` values in `path` are finite and not all zero; 0 without the column
static size_t count_environment_sites(const md_system_t& sys, double* scratch, str_t path, uint32_t components, md_unit_t unit) {
    if (!extract_atom_column(scratch, sys, path, components, unit)) {
        return 0;
    }
    const size_t num_atoms = md_system_atom_count(&sys);
    size_t count = 0;
    for (size_t i = 0; i < num_atoms; ++i) {
        if (!atom_in_environment(sys, i)) continue;
        const double* v = scratch + (size_t)components * i;
        bool finite = true, nonzero = false;
        for (uint32_t k = 0; k < components; ++k) {
            finite  &= is_finite_f64(v[k]);
            nonzero |= v[k] != 0.0;
        }
        count += finite && nonzero;
    }
    return count;
}

SurfaceFieldEmbeddingContent surface_field_embedding_content(const md_system_t& sys) {
    SurfaceFieldEmbeddingContent c = {};
    const size_t num_atoms = md_system_atom_count(&sys);
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    double* scratch = (double*)md_temp_alloc(temp, sizeof(double) * 6 * MAX(num_atoms, (size_t)1));
    const md_unit_t e_bohr  = md_unit_elementary_charge_bohr();
    const md_unit_t e_bohr2 = md_unit_mul(md_unit_elementary_charge(), md_unit_pow(md_unit_bohr_radius(), 2));
    c.charges        = count_environment_sites(sys, scratch, STR_LIT("atom/charge"),        1, md_unit_elementary_charge());
    c.dipoles        = count_environment_sites(sys, scratch, STR_LIT("atom/dipole"),        3, e_bohr);
    c.quadrupoles    = count_environment_sites(sys, scratch, STR_LIT("atom/quadrupole"),    6, e_bohr2);
    c.induced        = count_environment_sites(sys, scratch, EMBEDDING_INDUCED_DIPOLE_PATH, 3, e_bohr);
    c.polarizable    = count_environment_sites(sys, scratch, STR_LIT("atom/polarizability"), 1, md_unit_none());
    return c;
}

// The environment's classical sites as points with their charge, and their dipole and quadrupole
// where the potential has them. A site's dipole is its permanent one and the one induced on it.
static bool build_embedding_charges(md_gto_int_charges_t* out, const SurfaceFieldDesc& desc, md_allocator_i* alloc) {
    const md_system_t& sys = *desc.sys;
    const md_attribute_t* attr = find_atom_charges(sys);
    const size_t num_atoms = md_system_atom_count(&sys);
    if (!attr || !desc.state || desc.state->num_atoms < num_atoms || !desc.state->xyz) {
        return false;
    }
    const md_system_state_t& state = *desc.state;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    defer { md_temp_end(temp); };

    double* q = (double*)md_temp_alloc(temp, sizeof(double) * num_atoms);
    if (md_attribute_extract_f64(q, num_atoms, attr, md_attribute_slice_all(), md_unit_elementary_charge()) != num_atoms) {
        return false;
    }
    double* mu = (double*)md_temp_alloc(temp, sizeof(double) * 3 * num_atoms);
    double* Q  = (double*)md_temp_alloc(temp, sizeof(double) * 6 * num_atoms);
    const md_unit_t e_bohr2 = md_unit_mul(md_unit_elementary_charge(), md_unit_pow(md_unit_bohr_radius(), 2));
    double* mu_ind = (double*)md_temp_alloc(temp, sizeof(double) * 3 * num_atoms);
    const bool has_mu  = extract_atom_column(mu, sys, STR_LIT("atom/dipole"),     3, md_unit_elementary_charge_bohr());
    const bool has_ind = extract_atom_column(mu_ind, sys, EMBEDDING_INDUCED_DIPOLE_PATH, 3, md_unit_elementary_charge_bohr());
    const bool has_Q   = extract_atom_column(Q,  sys, STR_LIT("atom/quadrupole"), 6, e_bohr2);

    float*  xyz    = (float*)md_temp_alloc(temp, sizeof(float) * 3 * num_atoms);
    double* charge = (double*)md_temp_alloc(temp, sizeof(double) * num_atoms);
    double* dipole = (has_mu || has_ind) ? (double*)md_temp_alloc(temp, sizeof(double) * 3 * num_atoms) : nullptr;
    double* quad   = has_Q  ? (double*)md_temp_alloc(temp, sizeof(double) * 6 * num_atoms) : nullptr;
    size_t  n = 0;
    for (size_t i = 0; i < num_atoms; ++i) {
        // An atom of the QM region is not a site, whatever it carries
        if (!atom_in_environment(sys, i)) continue;
        // An atom with no charge at all (NAN) is not a site; a site is one whatever it carries, as an
        // expansion point can have a dipole or a quadrupole and no charge
        if (!is_finite_f64(q[i])) continue;
        bool any = q[i] != 0.0;
        double m[9] = {};
        for (int k = 0; has_mu && k < 3; ++k) m[k]     = is_finite_f64(mu[3 * i + k]) ? mu[3 * i + k] : 0.0;
        for (int k = 0; has_ind && k < 3; ++k) m[k]   += is_finite_f64(mu_ind[3 * i + k]) ? mu_ind[3 * i + k] : 0.0;
        for (int k = 0; has_Q  && k < 6; ++k) m[3 + k] = is_finite_f64(Q[6 * i + k])  ? Q[6 * i + k]  : 0.0;
        for (int k = 0; k < 9; ++k) any |= m[k] != 0.0;
        if (!any) continue;

        xyz[3 * n + 0] = (float)(state.xyz[i].x * ANGSTROM_TO_BOHR_D);
        xyz[3 * n + 1] = (float)(state.xyz[i].y * ANGSTROM_TO_BOHR_D);
        xyz[3 * n + 2] = (float)(state.xyz[i].z * ANGSTROM_TO_BOHR_D);
        charge[n] = q[i];
        if (dipole) MEMCPY(dipole + 3 * n, m, sizeof(double) * 3);
        if (quad)   MEMCPY(quad + 6 * n, m + 3, sizeof(double) * 6);
        n += 1;
    }
    if (n == 0) {
        return false;
    }

    md_gto_int_charges_desc_t cd = {};
    cd.point_xyz        = xyz;
    cd.point_charge     = charge;
    cd.num_points       = n;
    cd.point_dipole     = dipole;
    cd.point_quadrupole = quad;
    cd.far_field_theta  = MD_GTO_INT_FAR_FIELD_THETA;   // an environment can be tens of thousands of sites
    return md_gto_int_charges_init(out, &cd, alloc);
}

// ---------------------------------------------------------------------------
// Electrostatic potential: the QM region's nuclei and electron density
// ---------------------------------------------------------------------------

// The total density; or, where the calculation has no beta channel and md_qm publishes no total, the
// alpha density, which then carries the whole occupation.
static const md_attribute_t* find_esp_density(const md_system_t& sys) {
    if (const md_attribute_t* total = md_attributes_find(&sys.attributes, ESP_TOTAL_DENSITY_PATH)) {
        return total;
    }
    if (!md_attributes_find(&sys.attributes, ESP_BETA_DENSITY_PATH)) {
        return md_attributes_find(&sys.attributes, ESP_ALPHA_DENSITY_PATH);
    }
    return nullptr;
}

static bool esp_available(const md_system_t& sys) {
    return md_attributes_find(&sys.attributes, STR_LIT("basis/shell/atom_index")) != nullptr &&
           md_attributes_find(&sys.attributes, STR_LIT("basis/primitive/exponent")) != nullptr &&
           md_attributes_find(&sys.attributes, ESP_NUCLEAR_CHARGE_PATH) != nullptr &&
           find_esp_density(sys) != nullptr;
}

// The nuclear charges are what the readers publish as qm/atom/nuclear_charge: the effective ones where
// the format states them (TREXIO), the atomic numbers where it does not. A calculation with effective
// core potentials read from a format of the second kind leaves the density short of the core
// electrons the atomic numbers count, and the potential off by a long range Q / r - on a molecule
// with one iodine, 46 e / r. That cannot be corrected here, only noticed: the net charge of nuclei
// and density against the molecular charge the calculation states, and failing that against what a
// molecule plausibly carries.
static void esp_check_net_charge(const md_gto_int_charges_t* q, const md_system_t& sys) {
    const double Q = md_gto_int_charges_moments(q, nullptr).charge;

    double stated = 0.0;
    bool has_stated = false;
    if (const md_attribute_t* attr = md_attributes_find(&sys.attributes, ESP_MOLECULAR_CHARGE_PATH)) {
        has_stated = md_attribute_extract_f64(&stated, 1, attr, md_attribute_slice_all(), md_unit_none()) == 1 && is_finite_f64(stated);
    }

    // Said once per distinct mismatch, not on every frame of a trajectory (main thread only)
    static double last_reported = NAN;
    const bool mismatch = has_stated ? fabs(Q - stated) > 0.05 : (fabs(Q) > 3.5 || fabs(Q - round(Q)) > 0.05);
    if (mismatch && is_finite_f64(last_reported) && fabs(Q - last_reported) < 1.0e-3) {
        MD_LOG_DEBUG("Electrostatic potential: net charge of nuclei and density %.4f e (reported before)", Q);
        return;
    }
    if (mismatch) last_reported = Q;

    if (has_stated && fabs(Q - stated) > 0.05) {
        MD_LOG_INFO("Electrostatic potential: the nuclei and the electron density add up to %.2f e, and the calculation states a charge of %.0f e. "
                    "With effective core potentials the file does not state, the nuclear charges count core electrons the density does not have, and the potential is off by the difference over r.",
                    Q, stated);
    } else if (!has_stated && (fabs(Q) > 3.5 || fabs(Q - round(Q)) > 0.05)) {
        MD_LOG_INFO("Electrostatic potential: the nuclei and the electron density add up to %.2f e. "
                    "If the calculation used effective core potentials the file does not state, the nuclear charges count core electrons the density does not have, and the potential is off by the difference over r.",
                    Q);
    } else {
        MD_LOG_DEBUG("Electrostatic potential: net charge of nuclei and density %.4f e", Q);
    }
}

// The QM region as a charge distribution: the basis and the density over it, and the nuclei at the
// basis atoms with their nuclear charges. Positions are the caller's (desc.basis_atom_xyz).
static bool build_esp_charges(md_gto_int_charges_t* out, const SurfaceFieldDesc& desc, md_allocator_i* alloc) {
    const md_system_t& sys = *desc.sys;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    defer { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        MD_LOG_ERROR("Electrostatic potential: the system publishes no usable basis");
        return false;
    }
    const size_t num_atoms = md_gto_basis_num_atoms(&basis);
    const size_t num_ao    = md_gto_basis_num_ao(&basis);
    if (!desc.basis_atom_xyz || desc.num_basis_atoms < num_atoms) {
        MD_LOG_ERROR("Electrostatic potential: the basis spans %zu atoms and %zu positions were given", num_atoms, desc.num_basis_atoms);
        return false;
    }

    const md_attribute_t* d_attr = find_esp_density(sys);
    if (!d_attr || d_attr->format.rank != 2 || d_attr->format.components != 1 ||
        d_attr->format.shape[0] != num_ao || d_attr->format.shape[1] != num_ao) {
        MD_LOG_ERROR("Electrostatic potential: no ground state density over the %zu atomic orbitals of the basis", num_ao);
        return false;
    }
    const size_t num_d = num_ao * num_ao;
    double* D = (double*)md_temp_alloc(temp, sizeof(double) * num_d);
    if (!D || md_attribute_extract_f64(D, num_d, d_attr, md_attribute_slice_all(), md_unit_none()) != num_d) {
        MD_LOG_ERROR("Electrostatic potential: could not read '" STR_FMT "'", STR_ARG(d_attr->path));
        return false;
    }

    // The nuclear charges are of the QM atoms, which is what the basis' atom indices are too
    const md_attribute_t* z_attr = md_attributes_find(&sys.attributes, ESP_NUCLEAR_CHARGE_PATH);
    const size_t num_z = z_attr ? md_attribute_element_count(&z_attr->format) : 0;
    if (num_z < num_atoms) {
        MD_LOG_ERROR("Electrostatic potential: the basis spans %zu atoms and '" STR_FMT "' holds %zu", num_atoms, STR_ARG(ESP_NUCLEAR_CHARGE_PATH), num_z);
        return false;
    }
    double* Z = (double*)md_temp_alloc(temp, sizeof(double) * num_z);
    if (!Z || md_attribute_extract_f64(Z, num_z, z_attr, md_attribute_slice_all(), md_unit_elementary_charge()) != num_z) {
        return false;
    }

    md_gto_int_charges_desc_t cd = {};
    cd.basis            = &basis;
    cd.atom_xyz         = desc.basis_atom_xyz[0].elem;
    cd.atom_xyz_stride  = sizeof(vec3_t);
    cd.density_matrix   = D;
    cd.density_scale    = -1.0;
    cd.point_xyz        = desc.basis_atom_xyz[0].elem;
    cd.point_xyz_stride = sizeof(vec3_t);
    cd.point_charge     = Z;
    cd.num_points       = num_atoms;
    cd.threshold        = ESP_SCREENING_THRESHOLD;
    // Beyond the reach of its gaussians, a part of the density is its multipole expansion: on a
    // 1152 atom protein ~25x less work per point near its surface, within ~5e-5 au of the exact sum
    // (md_gto_int.h). Without it every gaussian is evaluated at every point, ~12 M steps a point there.
    cd.far_field_theta  = MD_GTO_INT_FAR_FIELD_THETA;
    if (!md_gto_int_charges_init(out, &cd, alloc)) {
        MD_LOG_ERROR("Electrostatic potential: could not build the charge distribution");
        return false;
    }
    esp_check_net_charge(out, sys);
    return true;
}

// ---------------------------------------------------------------------------
// Surface gaussians: CPCM and GOSTSHYP
// ---------------------------------------------------------------------------
// md_vlx publishes the gaussians these environment models put on a surface around the molecule as
// <group>/{center,exponent,normalization,amplitude}: g_k(r) = a_k N_k exp(-w_k |r - s_k|^2), with
// N_k = (w_k/pi)^(3/2) so that g_k integrates to a_k. For CPCM a_k is a surface charge, and the
// gaussians are the solvent's polarisation smeared out: the field is either their potential (the
// reaction potential, Coulomb, through md_gto_int) or the charge density itself. For GOSTSHYP the
// sum of the gaussians IS the field, the potential that exerts the pressure on the electrons.
//
// The centres are where the calculation had the molecule. The density is drawn where the system has
// it now, so the gaussians are carried along with the QM atoms: shifted by how far the centroid of
// the QM atoms' current positions is from that of their positions in the calculation
// (qm/atom/coordinate) - a recentred or unwrapped system keeps its cavity around it.

static const str_t QM_ATOM_COORDINATE_PATH = STR_INIT("qm/atom/coordinate");

static str_t surface_gaussians_group(SurfaceFieldKind kind) {
    switch (kind) {
    case SurfaceFieldKind::CpcmPotential:
    case SurfaceFieldKind::CpcmChargeDensity: return STR_LIT("cpcm");
    case SurfaceFieldKind::GostshypPotential: return STR_LIT("gostshyp");
    default: return (str_t){0};
    }
}

// The field is the sum of the gaussians as they are, not a potential md_gto_int evaluates
static bool kind_is_gaussian_sum(SurfaceFieldKind kind) {
    return kind == SurfaceFieldKind::CpcmChargeDensity || kind == SurfaceFieldKind::GostshypPotential;
}

static const md_attribute_t* surface_gaussians_leaf(const md_system_t& sys, str_t group, const char* leaf) {
    char path[64];
    snprintf(path, sizeof(path), STR_FMT "/%s", STR_ARG(group), leaf);
    return md_attributes_find(&sys.attributes, str_from_cstr(path));
}

// How many gaussians the group has: 0 unless all four columns are there, the centre 3-vectors and
// the others scalars, all of the same length
static size_t surface_gaussians_count(const md_system_t& sys, str_t group) {
    if (str_empty(group)) return 0;
    const md_attribute_t* c = surface_gaussians_leaf(sys, group, "center");
    const md_attribute_t* w = surface_gaussians_leaf(sys, group, "exponent");
    const md_attribute_t* n = surface_gaussians_leaf(sys, group, "normalization");
    const md_attribute_t* a = surface_gaussians_leaf(sys, group, "amplitude");
    if (!c || !w || !n || !a || c->format.components != 3) return 0;
    const size_t num = md_attribute_value_count(&c->format);
    for (const md_attribute_t* x : { w, n, a }) {
        if (x->format.components != 1 || md_attribute_value_count(&x->format) != num) return 0;
    }
    return num;
}

// The unit the amplitudes are read in, which makes amplitude x normalization the field's unit
static md_unit_t surface_gaussians_amplitude_unit(SurfaceFieldKind kind) {
    const md_unit_t bohr3 = md_unit_pow(md_unit_bohr_radius(), 3);
    return kind == SurfaceFieldKind::GostshypPotential ? md_unit_mul(md_unit_hartree(), bohr3) : md_unit_elementary_charge();
}

// The shift that carries the calculation's geometry to where the QM atoms are now (bohr); zero
// without the positions of both
static void surface_gaussians_shift(double out[3], const md_system_t& sys, const vec3_t* basis_atom_xyz, size_t num_basis_atoms, md_temp_scope_t temp) {
    out[0] = out[1] = out[2] = 0.0;
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, QM_ATOM_COORDINATE_PATH);
    if (!attr || attr->format.components != 3 || !basis_atom_xyz) return;
    const size_t num = md_attribute_value_count(&attr->format);
    if (num == 0 || num != num_basis_atoms) return;
    double* ref = (double*)md_temp_alloc(temp, sizeof(double) * 3 * num);
    if (md_attribute_extract_f64(ref, 3 * num, attr, md_attribute_slice_all(), md_unit_bohr_radius()) != 3 * num) return;
    for (size_t i = 0; i < num; ++i) {
        for (int a = 0; a < 3; ++a) out[a] += ((double)basis_atom_xyz[i].elem[a] - ref[3 * i + a]) / (double)num;
    }
}

// The group's gaussians: centres (bohr, carried along), exponents (bohr^-2), normalisations
// (bohr^-3) and amplitudes (surface_gaussians_amplitude_unit), each [count]
static size_t surface_gaussians_extract(double** xyz, double** w, double** n, double** a, const SurfaceFieldDesc& desc, md_temp_scope_t temp) {
    const md_system_t& sys = *desc.sys;
    const str_t group = surface_gaussians_group(desc.kind);
    const size_t num = surface_gaussians_count(sys, group);
    if (num == 0) return 0;
    *xyz = (double*)md_temp_alloc(temp, sizeof(double) * 3 * num);
    *w   = (double*)md_temp_alloc(temp, sizeof(double) * num);
    *n   = (double*)md_temp_alloc(temp, sizeof(double) * num);
    *a   = (double*)md_temp_alloc(temp, sizeof(double) * num);
    const md_unit_t bohr = md_unit_bohr_radius();
    const md_attribute_slice_t all = md_attribute_slice_all();
    if (md_attribute_extract_f64(*xyz, 3 * num, surface_gaussians_leaf(sys, group, "center"),        all, bohr) != 3 * num ||
        md_attribute_extract_f64(*w,   num, surface_gaussians_leaf(sys, group, "exponent"),          all, md_unit_pow(bohr, -2)) != num ||
        md_attribute_extract_f64(*n,   num, surface_gaussians_leaf(sys, group, "normalization"),     all, md_unit_pow(bohr, -3)) != num ||
        md_attribute_extract_f64(*a,   num, surface_gaussians_leaf(sys, group, "amplitude"),         all, surface_gaussians_amplitude_unit(desc.kind)) != num) {
        MD_LOG_ERROR("Surface field: the " STR_FMT " gaussians could not be read in the units they are evaluated in", STR_ARG(group));
        return 0;
    }
    double shift[3];
    surface_gaussians_shift(shift, sys, desc.basis_atom_xyz, desc.num_basis_atoms, temp);
    for (size_t k = 0; k < num; ++k) {
        for (int c = 0; c < 3; ++c) (*xyz)[3 * k + c] += shift[c];
    }
    return num;
}

// CPCM's reaction potential: the surface charges as the normalised spherical gaussians they are,
// whose potential md_gto_int evaluates exactly, sum_k q_k erf(zeta_k r) / r - on the GPU as well.
// The normalisation is (w/pi)^(3/2) by construction (md_vlx says so where a file disagrees), which
// is what md_gto_int's gaussian charges are.
static bool build_cpcm_charges(md_gto_int_charges_t* out, const SurfaceFieldDesc& desc, md_allocator_i* alloc) {
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    defer { md_temp_end(temp); };
    double *xyz = nullptr, *w = nullptr, *n = nullptr, *q = nullptr;
    const size_t num = surface_gaussians_extract(&xyz, &w, &n, &q, desc, temp);
    if (num == 0) return false;
    md_gto_int_charges_desc_t cd = {};
    cd.gaussian_xyz         = xyz;
    cd.gaussian_exponent    = w;
    cd.gaussian_charge      = q;
    cd.num_gaussian_charges = num;
    cd.far_field_theta      = num > 2000 ? MD_GTO_INT_FAR_FIELD_THETA : 0.0;   // a protein's cavity is tens of thousands
    return md_gto_int_charges_init(out, &cd, alloc);
}

// A plain sum of gaussians, sum_k c_k exp(-w_k |r - s_k|^2), the field evaluated as it is
struct SurfaceGaussians {
    size_t  count = 0;
    double* xyz = nullptr;        // [count * 3] bohr
    double* exponent = nullptr;   // [count] bohr^-2
    double* coeff = nullptr;      // [count] amplitude x normalisation, in the field's unit
};

static void surface_gaussians_destroy(SurfaceGaussians* g) {
    if (!g) return;
    free(g->xyz);
    free(g->exponent);
    free(g->coeff);
    delete g;
}

static SurfaceGaussians* build_surface_gaussians(const SurfaceFieldDesc& desc) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    double *xyz = nullptr, *w = nullptr, *n = nullptr, *a = nullptr;
    const size_t num = surface_gaussians_extract(&xyz, &w, &n, &a, desc, temp);
    if (num == 0) return nullptr;
    SurfaceGaussians* g = new SurfaceGaussians();
    g->xyz      = (double*)malloc(sizeof(double) * 3 * num);
    g->exponent = (double*)malloc(sizeof(double) * num);
    g->coeff    = (double*)malloc(sizeof(double) * num);
    if (!g->xyz || !g->exponent || !g->coeff) {
        surface_gaussians_destroy(g);
        return nullptr;
    }
    size_t m = 0;
    for (size_t k = 0; k < num; ++k) {
        const double c = a[k] * n[k];
        if (!(w[k] > 0.0) || c == 0.0 || !is_finite_f64(c)) continue;
        MEMCPY(g->xyz + 3 * m, xyz + 3 * k, sizeof(double) * 3);
        g->exponent[m] = w[k];
        g->coeff[m]    = c;
        m += 1;
    }
    g->count = m;
    return g;
}

// exp(-x) below this relative to the largest term is not added: 1e-16 of it
#define SURFACE_GAUSSIAN_CUTOFF 36.8

static void surface_gaussians_evaluate(double* out, const float* xyz, size_t num_xyz, const SurfaceGaussians* g) {
    for (size_t i = 0; i < num_xyz; ++i) {
        const double x = xyz[3 * i + 0], y = xyz[3 * i + 1], z = xyz[3 * i + 2];
        double v = 0.0;
        for (size_t k = 0; k < g->count; ++k) {
            const double dx = x - g->xyz[3 * k + 0], dy = y - g->xyz[3 * k + 1], dz = z - g->xyz[3 * k + 2];
            const double e = g->exponent[k] * (dx * dx + dy * dy + dz * dz);
            if (e < SURFACE_GAUSSIAN_CUTOFF) v += g->coeff[k] * exp(-e);
        }
        out[i] = v;
    }
}

// ---------------------------------------------------------------------------
// The surface gaussians as points
// ---------------------------------------------------------------------------

static str_t surface_points_group_name(SurfacePointsGroup group) {
    return group == SurfacePointsGroup::Gostshyp ? STR_LIT("gostshyp") : STR_LIT("cpcm");
}

size_t surface_points_count(SurfacePointsGroup group, const md_system_t& sys) {
    return surface_gaussians_count(sys, surface_points_group_name(group));
}

md_unit_t surface_points_unit(SurfacePointsGroup group) {
    return surface_gaussians_amplitude_unit(group == SurfacePointsGroup::Gostshyp ? SurfaceFieldKind::GostshypPotential : SurfaceFieldKind::CpcmPotential);
}

bool surface_points_signed(SurfacePointsGroup group) {
    return group != SurfacePointsGroup::Gostshyp;
}

size_t surface_points_extract(vec3_t* out_xyz, float* out_amplitude, size_t cap, SurfacePointsGroup group, const md_system_t& sys,
                              const vec3_t* basis_atom_xyz, size_t num_basis_atoms) {
    const str_t name = surface_points_group_name(group);
    const size_t num = surface_gaussians_count(sys, name);
    if (num == 0 || num > cap) return 0;
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    double* xyz = (double*)md_temp_alloc(temp, sizeof(double) * 3 * num);
    double* a   = (double*)md_temp_alloc(temp, sizeof(double) * num);
    if (md_attribute_extract_f64(xyz, 3 * num, surface_gaussians_leaf(sys, name, "center"), md_attribute_slice_all(), md_unit_bohr_radius()) != 3 * num ||
        md_attribute_extract_f64(a, num, surface_gaussians_leaf(sys, name, "amplitude"), md_attribute_slice_all(), surface_points_unit(group)) != num) {
        return 0;
    }
    double shift[3];
    surface_gaussians_shift(shift, sys, basis_atom_xyz, num_basis_atoms, temp);
    for (size_t k = 0; k < num; ++k) {
        for (int c = 0; c < 3; ++c) out_xyz[k].elem[c] = (float)((xyz[3 * k + c] + shift[c]) / ANGSTROM_TO_BOHR_D);
        out_amplitude[k] = (float)a[k];
    }
    return num;
}

// ---------------------------------------------------------------------------
// Kinds
// ---------------------------------------------------------------------------

bool surface_field_available(SurfaceFieldKind kind, const md_system_t& sys) {
    switch (kind) {
    case SurfaceFieldKind::EmbeddingPotential:
        return find_atom_charges(sys) != nullptr && has_environment(sys);
    case SurfaceFieldKind::ElectrostaticPotential:
        return esp_available(sys);
    case SurfaceFieldKind::CpcmPotential:
    case SurfaceFieldKind::CpcmChargeDensity:
    case SurfaceFieldKind::GostshypPotential:
        return surface_gaussians_count(sys, surface_gaussians_group(kind)) > 0;
    default:
        return false;
    }
}

bool surface_field_needs_basis_atoms(SurfaceFieldKind kind) {
    // The ESP evaluates the density there; the surface gaussians are carried along with them
    return kind == SurfaceFieldKind::ElectrostaticPotential || !str_empty(surface_gaussians_group(kind));
}

bool surface_field_signed(SurfaceFieldKind kind) {
    return kind != SurfaceFieldKind::GostshypPotential;
}

SurfaceFieldKind surface_field_first_available(const md_system_t& sys) {
    for (int k = 0; k < (int)SurfaceFieldKind::Count; ++k) {
        if (surface_field_available((SurfaceFieldKind)k, sys)) return (SurfaceFieldKind)k;
    }
    return SurfaceFieldKind::Count;
}

static bool build_charges(md_gto_int_charges_t* out, const SurfaceFieldDesc& desc, md_allocator_i* alloc) {
    switch (desc.kind) {
    case SurfaceFieldKind::EmbeddingPotential:     return build_embedding_charges(out, desc, alloc);
    case SurfaceFieldKind::ElectrostaticPotential: return build_esp_charges(out, desc, alloc);
    case SurfaceFieldKind::CpcmPotential:          return build_cpcm_charges(out, desc, alloc);
    default: return false;
    }
}

md_grid_t surface_field_grid(SurfaceFieldKind kind, const md_grid_t& density_grid) {
    if (kind != SurfaceFieldKind::ElectrostaticPotential) {
        return density_grid;
    }
    // The same box (origin, orientation, extent) at about ESP_GRID_SPACING, never finer than the
    // density's: the field texture is sampled with the density's normalised coordinates
    md_grid_t grid = density_grid;
    for (int a = 0; a < 3; ++a) {
        const double extent = (double)density_grid.spacing.elem[a] * density_grid.dim[a];
        const int dim = (int)ceil(extent / ESP_GRID_SPACING - 1.0e-3);
        grid.dim[a] = CLAMP(dim, MIN(2, density_grid.dim[a]), density_grid.dim[a]);
        grid.spacing.elem[a] = (float)(extent / grid.dim[a]);
    }
    return grid;
}

// Whether two grids sample the same points, to float rounding: a grid derived again from the same box
// can differ in the last bits of its spacing, which is not a reason to evaluate again
static bool grid_same(const md_grid_t& a, const md_grid_t& b) {
    for (int i = 0; i < 3; ++i) {
        if (a.dim[i] != b.dim[i]) return false;
        if (fabsf(a.spacing.elem[i] - b.spacing.elem[i]) > 1.0e-6f * MAX(fabsf(a.spacing.elem[i]), 1.0f)) return false;
        if (fabsf(a.origin.elem[i] - b.origin.elem[i]) > 1.0e-4f) return false;
        for (int j = 0; j < 3; ++j) {
            if (fabsf(a.orientation.elem[i][j] - b.orientation.elem[i][j]) > 1.0e-6f) return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// The band and the volume
// ---------------------------------------------------------------------------

static inline bool bit_get(const uint64_t* bits, size_t i) { return (bits[i >> 6] >> (i & 63)) & 1; }
static inline void bit_set(uint64_t* bits, size_t i)       { bits[i >> 6] |= (1ull << (i & 63)); }

// ---------------------------------------------------------------------------
// What is under way, for a progress display
// ---------------------------------------------------------------------------
struct FieldProgress {
    const char* label = "";
    float       fraction = 0.0f;
    bool        preparing = false;
    bool        active = false;
    md_tick_t   t0 = 0;
};

static constexpr uint32_t MAX_PROGRESS = 64;
static FieldProgress* g_progress[MAX_PROGRESS];
static uint32_t       g_num_progress = 0;

static void progress_set(FieldProgress* p, const char* label, float fraction, bool preparing) {
    if (!p->active && g_num_progress < MAX_PROGRESS) {
        p->active = true;
        p->t0 = md_tick_now();
        g_progress[g_num_progress++] = p;
    }
    p->label = label;
    p->fraction = CLAMP(fraction, 0.0f, 1.0f);
    p->preparing = preparing;
}

static void progress_clear(FieldProgress* p) {
    if (!p->active) return;
    for (uint32_t i = 0; i < g_num_progress; ++i) {
        if (g_progress[i] == p) {
            g_progress[i] = g_progress[--g_num_progress];
            break;
        }
    }
    p->active = false;
}

size_t surface_field_progress(SurfaceFieldProgress* out, size_t cap) {
    const md_tick_t now = md_tick_now();
    size_t n = 0;
    for (uint32_t i = 0; i < g_num_progress && n < cap; ++i) {
        const FieldProgress* p = g_progress[i];
        if (md_tick_to_seconds(now - p->t0) < 0.25) continue;   // what is done within a few frames shows nothing
        out[n].label = p->label;
        out[n].fraction = p->fraction;
        out[n].preparing = p->preparing;
        n += 1;
    }
    return n;
}

// ---------------------------------------------------------------------------
// The charge distribution of a source, made on a worker thread
// ---------------------------------------------------------------------------
// For a QM region that is reading the density matrix and making md_gto_int's distribution of it:
// seconds on a protein, so it is a pool task and the frame goes on. What the build reads is copied
// when it is asked for (the positions), or does not change while pool tasks run (the system's
// attributes: a load waits for the pool, interrupt_async_tasks; the on demand density reader takes the
// HDF5 lock). One build at a time per source: a request for another while one runs is kept, and
// started when it ends, unless a newer request has replaced it.
struct FieldBuild {
    uint64_t           hash = 0;
    SurfaceFieldKind   kind = SurfaceFieldKind::EmbeddingPotential;
    const md_system_t* sys = nullptr;
    vec3_t*            basis_atom_xyz = nullptr;
    size_t             num_basis_atoms = 0;
    vec3_t*            state_xyz = nullptr;
    size_t             num_state_atoms = 0;
    // made by the task: the one or the other (kind_is_gaussian_sum)
    md_gto_int_charges_t* charges = nullptr;
    SurfaceGaussians*  gaussians = nullptr;
    bool               ran = false;     // an interrupted pool task never runs its function
    bool               ok = false;
    double             seconds = 0.0;
};

static void charges_destroy(md_gto_int_charges_t* q) {
    if (!q) return;
    md_allocator_i* heap = md_get_heap_allocator();
    md_gto_int_charges_free(q, heap);
    md_free(heap, q, sizeof(md_gto_int_charges_t));
}

static void build_destroy(FieldBuild* b) {
    if (!b) return;
    charges_destroy(b->charges);
    surface_gaussians_destroy(b->gaussians);
    free(b->basis_atom_xyz);
    free(b->state_xyz);
    delete b;
}

static FieldBuild* build_create(const SurfaceFieldDesc& desc) {
    FieldBuild* b = new FieldBuild();
    b->hash = desc.source_hash;
    b->kind = desc.kind;
    b->sys  = desc.sys;
    if (desc.basis_atom_xyz && desc.num_basis_atoms) {
        b->basis_atom_xyz = (vec3_t*)malloc(sizeof(vec3_t) * desc.num_basis_atoms);
        if (b->basis_atom_xyz) {
            MEMCPY(b->basis_atom_xyz, desc.basis_atom_xyz, sizeof(vec3_t) * desc.num_basis_atoms);
            b->num_basis_atoms = desc.num_basis_atoms;
        }
    }
    if (desc.state && desc.state->xyz && desc.state->num_atoms) {
        b->state_xyz = (vec3_t*)malloc(sizeof(vec3_t) * desc.state->num_atoms);
        if (b->state_xyz) {
            MEMCPY(b->state_xyz, desc.state->xyz, sizeof(vec3_t) * desc.state->num_atoms);
            b->num_state_atoms = desc.state->num_atoms;
        }
    }
    return b;
}

// On a worker thread
static void build_run(FieldBuild* b) {
    b->ran = true;
    const md_tick_t t0 = md_tick_now();
    md_system_state_t state = {};
    state.num_atoms = b->num_state_atoms;
    state.xyz       = b->state_xyz;
    SurfaceFieldDesc desc = {};
    desc.kind            = b->kind;
    desc.sys             = b->sys;
    desc.state           = b->state_xyz ? &state : nullptr;
    desc.basis_atom_xyz  = b->basis_atom_xyz;
    desc.num_basis_atoms = b->num_basis_atoms;
    desc.source_hash     = b->hash;
    if (kind_is_gaussian_sum(b->kind)) {
        b->gaussians = build_surface_gaussians(desc);
        b->ok = b->gaussians != nullptr;
        b->seconds = md_tick_to_seconds(md_tick_now() - t0);
        return;
    }
    md_allocator_i* heap = md_get_heap_allocator();
    md_gto_int_charges_t* q = (md_gto_int_charges_t*)md_alloc(heap, sizeof(md_gto_int_charges_t));
    if (q) {
        MEMSET(q, 0, sizeof(*q));
        if (build_charges(q, desc, heap)) {
            b->charges = q;
            b->ok = true;
        } else {
            md_free(heap, q, sizeof(md_gto_int_charges_t));
        }
    }
    b->seconds = md_tick_to_seconds(md_tick_now() - t0);
}

struct FieldSource {
    uint64_t              hash = 0;          // what charges (or gaussians) is of
    md_gto_int_charges_t* charges = nullptr;
    SurfaceGaussians*     gaussians = nullptr;   // instead, for a field that is a sum of gaussians
    void*                 gpu = nullptr;     // md_gto_int_gpu_charges_t, made when first wanted
    bool                  gpu_failed = false;
    uint64_t              failed_hash = 0;   // a source whose build failed: not tried again
    FieldBuild*           build = nullptr;   // under way
    task_system::ID       task = task_system::INVALID_ID;
    FieldBuild*           next = nullptr;    // asked for while it ran
};

static void source_start(FieldSource* fs, FieldBuild* b) {
    fs->build = b;
    // and the worker's temp arenas given back after: the charges of a protein's density read its orbitals there
    fs->task  = task_system::create_pool_task(STR_LIT("##Surface field charges"), [b]() { build_run(b); md_temp_trim(); });
    task_system::enqueue_task(fs->task);
}

// The charges, and their copy on the device. Whatever reads them on the device must have ended.
static void source_drop_charges(FieldSource* fs, md_gpu_stream_t stream) {
#if MD_ENABLE_GPU
    if (fs->gpu && stream) md_gto_int_gpu_charges_destroy(stream, (md_gto_int_gpu_charges_t)fs->gpu);
#else
    (void)stream;
#endif
    fs->gpu = nullptr;
    fs->gpu_failed = false;
    charges_destroy(fs->charges);
    fs->charges = nullptr;
    surface_gaussians_destroy(fs->gaussians);
    fs->gaussians = nullptr;
    fs->hash = 0;
}

static bool source_ready(const FieldSource& fs, uint64_t hash)  { return (fs.charges || fs.gaussians) && fs.hash == hash; }
static bool source_failed(const FieldSource& fs, uint64_t hash) { return hash != 0 && fs.failed_hash == hash; }
static bool source_building(const FieldSource& fs)              { return fs.build != nullptr || fs.next != nullptr; }

// Asks for the charges of desc's source, unless they are there, being made, or failed to be
static void source_request(FieldSource* fs, const SurfaceFieldDesc& desc) {
    const uint64_t h = desc.source_hash;
    if (source_ready(*fs, h) || source_failed(*fs, h)) return;
    if (fs->build && fs->build->hash == h) {
        build_destroy(fs->next);    // what was asked for after it is not wanted any more
        fs->next = nullptr;
        return;
    }
    if (fs->next && fs->next->hash == h) return;
    FieldBuild* b = build_create(desc);
    if (!fs->build) {
        source_start(fs, b);
    } else {
        build_destroy(fs->next);
        fs->next = b;
    }
}

// Takes a finished build: its charges replace the source's if they are of the source wanted. Whatever
// reads the old ones on the device must have ended (they go when the new ones come). True when the
// charges changed.
static bool source_poll(FieldSource* fs, uint64_t wanted, md_gpu_stream_t stream, const char* what) {
    bool changed = false;
    if (fs->build && !task_system::task_is_running(fs->task)) {
        FieldBuild* b = fs->build;
        fs->build = nullptr;
        fs->task  = task_system::INVALID_ID;
        if (!b->ran && b->hash == wanted && !fs->next) {
            source_start(fs, b);   // interrupted before it ran (interrupt_async_tasks): again
            return false;
        }
        if (b->ran && b->hash == wanted) {
            if (b->ok) {
                source_drop_charges(fs, stream);
                fs->charges   = b->charges;
                fs->gaussians = b->gaussians;
                fs->hash      = b->hash;
                b->charges    = nullptr;
                b->gaussians  = nullptr;
                changed = true;
                if (fs->charges) {
                    MD_LOG_DEBUG("Surface field '%s': %u gaussians and %u points, %llu steps per point, made in %.1f ms on a worker thread", what,
                                 fs->charges->num_gaussians, fs->charges->num_points, (unsigned long long)md_gto_int_charges_work_per_point(fs->charges),
                                 b->seconds * 1.0e3);
                } else {
                    MD_LOG_DEBUG("Surface field '%s': a sum of %zu gaussians, made in %.1f ms on a worker thread", what, fs->gaussians->count, b->seconds * 1.0e3);
                }
            } else {
                fs->failed_hash = b->hash;
                MD_LOG_ERROR("Surface field: the system has nothing to make '%s' of", what);
            }
        }
        build_destroy(b);
    }
    if (!fs->build && fs->next) {
        FieldBuild* n = fs->next;
        fs->next = nullptr;
        if (n->hash == wanted) source_start(fs, n);
        else build_destroy(n);
    }
    return changed;
}

// Waits for a build under way (a few seconds at most, and only when a field goes while it is made)
static void source_free(FieldSource* fs, md_gpu_stream_t stream) {
    if (fs->build) {
        task_system::task_wait_for(fs->task);
        build_destroy(fs->build);
        fs->build = nullptr;
        fs->task  = task_system::INVALID_ID;
    }
    build_destroy(fs->next);
    fs->next = nullptr;
    source_drop_charges(fs, stream);
    fs->failed_hash = 0;
}

#if MD_ENABLE_GPU
// The charges on the device; NULL where md_gto_int cannot evaluate them there (point multipoles)
static md_gto_int_gpu_charges_t source_gpu(FieldSource* fs, md_gpu_stream_t stream) {
    if (!fs->charges || !stream) return nullptr;
    if (fs->charges->point_dipole || fs->charges->point_quadrupole) return nullptr;
    if (!fs->gpu && !fs->gpu_failed) {
        fs->gpu = md_gto_int_gpu_charges_create(stream, fs->charges);
        fs->gpu_failed = fs->gpu == nullptr;
    }
    return (md_gto_int_gpu_charges_t)fs->gpu;
}
#endif

// What a CPU evaluation reads: the source's charges, whose potential the field is, or its gaussians,
// whose sum it is
struct FieldValues {
    const md_gto_int_charges_t* charges = nullptr;
    const SurfaceGaussians*     gaussians = nullptr;
};

// ---------------------------------------------------------------------------
// The volume's work
// ---------------------------------------------------------------------------
struct SurfaceFieldWork {
    md_gpu_stream_t stream = nullptr;
    SurfaceFieldKind kind = SurfaceFieldKind::EmbeddingPotential;
    FieldSource   source;
    FieldProgress progress;

    uint8_t* want = nullptr;        // [num_voxels] the voxels the surfaces of the last update can sample
    size_t   num_left = 0;          // of those, holding no value yet
    float*   samples = nullptr;     // [3 * num_samples] where the surfaces cross the density grid's edges,
    size_t   num_samples = 0;       // in the field grid's index space (surface_statistics)
    size_t   cap_samples = 0;

    // The evaluation under way
    enum Stage { Idle, Gpu, Readback, Cpu } stage = Idle;
    uint8_t* got = nullptr;         // [num_voxels] what it fills
    bool     gpu_failed = false;    // the device could not serve this source: the CPU does
    md_tick_t t0 = 0;
#if MD_ENABLE_GPU
    md_gto_int_gpu_job_t* job = nullptr;
    md_gpu_texture_t tex = nullptr; // the field grid, written by the job
    int      tex_dim[3] = {};
    md_gpu_mem_t  rb = {};          // its readback
    md_gpu_sync_t rb_sync = {};
#endif
    task_system::ID cpu_task = task_system::INVALID_ID;
    std::atomic<bool> cpu_cancel{false};   // looked at every few points: a cancel waits for no more than that
    uint32_t* cpu_index = nullptr;
    float*    cpu_xyz = nullptr;
    double*   cpu_out = nullptr;
    size_t    cpu_count = 0;
    FieldValues cpu_values;
};

// The field at xyz[beg, end) into out, a few points at a time until cancel is set. A pool task's
// interrupt is only looked at between the ranges it hands out, which can be a good part of the task.
static void evaluate_points_range(double* out, const float* xyz, uint32_t beg, uint32_t end, FieldValues v, const std::atomic<bool>& cancel) {
    for (uint32_t b = beg; b < end && !cancel.load(std::memory_order_relaxed); b += 8) {
        const uint32_t e = MIN(b + 8, end);
        if (v.charges) {
            md_gto_int_potential_xyz(out + b, NULL, xyz + 3 * (size_t)b, e - b, 0, v.charges);
        } else if (v.gaussians) {
            surface_gaussians_evaluate(out + b, xyz + 3 * (size_t)b, e - b, v.gaussians);
        }
    }
}

// Work of one point, in md_gto_int's unit (a few flops): for a sum of gaussians one per gaussian, ~10
static uint64_t field_work_per_point(FieldValues v) {
    if (v.charges)   return md_gto_int_charges_work_per_point(v.charges);
    if (v.gaussians) return 10 * (uint64_t)v.gaussians->count;
    return 1;
}

static void work_cpu_free(SurfaceFieldWork* w) {
    free(w->cpu_index);
    free(w->cpu_xyz);
    free(w->cpu_out);
    w->cpu_index = nullptr;
    w->cpu_xyz = nullptr;
    w->cpu_out = nullptr;
    w->cpu_count = 0;
    w->cpu_task = task_system::INVALID_ID;
}

// Stops the evaluation under way; what it has made is lost
static void work_cancel(SurfaceFieldWork* w) {
#if MD_ENABLE_GPU
    if (w->job) {
        md_gto_int_gpu_job_end(w->job);
        w->job = nullptr;
    }
    if (w->rb.gpu) {
        md_gpu_free(w->stream, w->rb.gpu);   // after the copy in flight, stream ordered
        w->rb = {};
    }
#endif
    if (w->cpu_task) {
        w->cpu_cancel = true;
        task_system::task_interrupt_and_wait_for(w->cpu_task);
    }
    work_cpu_free(w);
    w->stage = SurfaceFieldWork::Idle;
    progress_clear(&w->progress);
}

static void work_destroy(SurfaceFieldWork* w) {
    if (!w) return;
    work_cancel(w);
    source_free(&w->source, w->stream);
#if MD_ENABLE_GPU
    if (w->tex) md_gpu_texture_destroy(w->tex);
#endif
    free(w->want);
    free(w->got);
    free(w->samples);
    delete w;
}

static void field_volume_reset(SurfaceFieldVolume* vol, const md_grid_t& grid) {
    const size_t n = md_grid_num_points(&grid);
    md_allocator_i* heap = md_get_heap_allocator();
    SurfaceFieldWork* w = vol->work;
    if (n != vol->num_voxels) {
        if (vol->values) md_free(heap, vol->values, sizeof(float) * vol->num_voxels);
        if (vol->done)   md_free(heap, vol->done, sizeof(uint64_t) * ((vol->num_voxels + 63) / 64));
        vol->values = (float*)md_alloc(heap, sizeof(float) * n);
        vol->done   = (uint64_t*)md_alloc(heap, sizeof(uint64_t) * ((n + 63) / 64));
        vol->num_voxels = n;
        free(w->want);
        free(w->got);
        w->want = (uint8_t*)malloc(n);
        w->got  = (uint8_t*)malloc(n);
    }
    MEMSET(vol->values, 0, sizeof(float) * n);
    MEMSET(vol->done, 0, sizeof(uint64_t) * ((n + 63) / 64));
    MEMSET(w->want, 0, n);
    w->num_left = 0;
    w->num_samples = 0;
    vol->grid = grid;
    vol->num_evaluated = 0;
    vol->num_uploaded = 0;
    vol->complete = false;
    vol->on_gpu = false;
    vol->num_surface_samples = 0;
}

// The voxels of the density grid a surface can sample: the 8 corners of every cell that some
// isovalue passes through, grown by one voxel, because the raycaster places a hit by interpolating
// along the ray between two samples, which can put it just across the face of the cell the surface
// actually crosses.
static void band_mask(uint8_t* mask, const float* d, const int dim[3], const float* iso, size_t num_iso) {
    const size_t nx = dim[0], ny = dim[1], nz = dim[2];
    const size_t sx = 1, sy = nx, sz = nx * ny;
    MEMSET(mask, 0, nx * ny * nz);

    for (size_t k = 0; k + 1 < nz; ++k) {
        for (size_t j = 0; j + 1 < ny; ++j) {
            for (size_t i = 0; i + 1 < nx; ++i) {
                const size_t c = k * sz + j * sy + i * sx;
                const float v[8] = {
                    d[c], d[c + sx], d[c + sy], d[c + sx + sy],
                    d[c + sz], d[c + sx + sz], d[c + sy + sz], d[c + sx + sy + sz],
                };
                float lo = v[0], hi = v[0];
                for (int m = 1; m < 8; ++m) { lo = MIN(lo, v[m]); hi = MAX(hi, v[m]); }
                bool crossed = false;
                for (size_t t = 0; t < num_iso; ++t) {
                    if (lo <= iso[t] && iso[t] <= hi && lo < hi) { crossed = true; break; }
                }
                if (!crossed) continue;
                mask[c] = mask[c + sx] = mask[c + sy] = mask[c + sx + sy] = 1;
                mask[c + sz] = mask[c + sx + sz] = mask[c + sy + sz] = mask[c + sx + sy + sz] = 1;
            }
        }
    }

    // Grow by one voxel, one axis at a time: bit 1 is the band so far, bit 2 what this pass adds,
    // so that a pass does not spread further than one voxel along its own axis
    for (int a = 0; a < 3; ++a) {
        const size_t step = a == 0 ? sx : a == 1 ? sy : sz;
        for (size_t k = 0; k < nz; ++k) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t i = 0; i < nx; ++i) {
                    const size_t idx = k * sz + j * sy + i;
                    if (!(mask[idx] & 1)) continue;
                    const size_t coord = a == 0 ? i : a == 1 ? j : k;
                    const size_t len   = a == 0 ? nx : a == 1 ? ny : nz;
                    if (coord > 0)       mask[idx - step] |= 2;
                    if (coord + 1 < len) mask[idx + step] |= 2;
                }
            }
        }
        for (size_t idx = 0; idx < nx * ny * nz; ++idx) {
            mask[idx] = mask[idx] ? 1 : 0;
        }
    }
}

// Where density voxel i falls on the field grid's axis, as a texture sampler sees it: the two field
// voxels its linear interpolation reads (equal at the clamped edges). On a grid of the same size that
// is the voxel itself.
static inline void field_support(int* lo, int* hi, int i, int ddim, int fdim) {
    if (ddim == fdim) {
        *lo = *hi = i;
        return;
    }
    const float u = ((float)i + 0.5f) * ((float)fdim / (float)ddim) - 0.5f;
    const int   f = (int)floorf(u);
    *lo = CLAMP(f,     0, fdim - 1);
    *hi = CLAMP(f + 1, 0, fdim - 1);
}

// The field voxels a surface of the band can sample: the support of every density voxel of the band.
// A point on a surface lies in a density cell whose corners are all in the band, and as the field grid
// is no finer than the density's, the support of any point between two density voxel centres lies
// within the union of theirs.
static void band_to_field(uint8_t* fmask, const uint8_t* dmask, const int ddim[3], const int fdim[3], md_temp_scope_t temp) {
    MEMSET(fmask, 0, (size_t)fdim[0] * fdim[1] * fdim[2]);
    int* lo[3];
    int* hi[3];
    for (int a = 0; a < 3; ++a) {
        lo[a] = (int*)md_temp_alloc(temp, sizeof(int) * ddim[a]);
        hi[a] = (int*)md_temp_alloc(temp, sizeof(int) * ddim[a]);
        for (int i = 0; i < ddim[a]; ++i) field_support(&lo[a][i], &hi[a][i], i, ddim[a], fdim[a]);
    }
    for (int k = 0; k < ddim[2]; ++k) {
        for (int j = 0; j < ddim[1]; ++j) {
            const size_t row = ((size_t)k * ddim[1] + j) * ddim[0];
            for (int i = 0; i < ddim[0]; ++i) {
                if (!dmask[row + i]) continue;
                for (int z = lo[2][k]; z <= hi[2][k]; ++z) {
                    for (int y = lo[1][j]; y <= hi[1][j]; ++y) {
                        uint8_t* dst = fmask + ((size_t)z * fdim[1] + y) * fdim[0];
                        dst[lo[0][i]] = 1;
                        dst[hi[0][i]] = 1;
                    }
                }
            }
        }
    }
}

// The field at continuous field grid index u, interpolated as the texture sampler does (linear,
// clamped to the edge). False if a voxel it needs holds no value.
static bool field_sample(float* out, const SurfaceFieldVolume& vol, const float u[3]) {
    const int* dim = vol.grid.dim;
    int   i0[3], i1[3];
    float f[3];
    for (int a = 0; a < 3; ++a) {
        const float c = CLAMP(u[a], 0.0f, (float)(dim[a] - 1));
        i0[a] = MIN((int)c, dim[a] - 1);
        i1[a] = MIN(i0[a] + 1, dim[a] - 1);
        f[a]  = c - (float)i0[a];
    }
    float acc = 0.0f;
    for (int c = 0; c < 8; ++c) {
        const float w = ((c & 1) ? f[0] : 1.0f - f[0]) * ((c & 2) ? f[1] : 1.0f - f[1]) * ((c & 4) ? f[2] : 1.0f - f[2]);
        if (w == 0.0f) continue;
        const size_t idx = ((size_t)((c & 4) ? i1[2] : i0[2]) * dim[1] + ((c & 2) ? i1[1] : i0[1])) * dim[0] + ((c & 1) ? i1[0] : i0[0]);
        if (!bit_get(vol.done, idx)) return false;
        acc += w * vol.values[idx];
    }
    *out = acc;
    return true;
}

// Where an isosurface crosses an edge of the density grid, in the field grid's index space: the
// samples of the field ON the surface, at the density grid's resolution. Kept with the band, since
// the values they read arrive later.
static void surface_samples_collect(SurfaceFieldWork* w, const SurfaceFieldVolume& vol, const int ddim[3], const float* d, const uint8_t* dmask,
                                    const float* iso, size_t num_iso) {
    const size_t nx = ddim[0], ny = ddim[1], nz = ddim[2];
    const size_t step[3] = { 1, nx, nx * ny };
    const size_t len[3]  = { nx, ny, nz };
    const float  ratio[3] = {
        (float)vol.grid.dim[0] / (float)ddim[0],
        (float)vol.grid.dim[1] / (float)ddim[1],
        (float)vol.grid.dim[2] / (float)ddim[2],
    };
    w->num_samples = 0;
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                const size_t idx = (k * ny + j) * nx + i;
                if (!dmask[idx]) continue;
                const size_t coord[3] = { i, j, k };
                for (int a = 0; a < 3; ++a) {
                    if (coord[a] + 1 >= len[a]) continue;
                    const size_t nb = idx + step[a];
                    if (!dmask[nb]) continue;
                    const float d0 = d[idx], d1 = d[nb];
                    for (size_t t = 0; t < num_iso; ++t) {
                        if ((d0 - iso[t]) * (d1 - iso[t]) >= 0.0f) continue;
                        const float s = (iso[t] - d0) / (d1 - d0);
                        if (w->num_samples == w->cap_samples) {
                            const size_t cap = MAX(w->cap_samples * 2, (size_t)4096);
                            float* bigger = (float*)realloc(w->samples, sizeof(float) * 3 * cap);
                            if (!bigger) return;
                            w->samples = bigger;
                            w->cap_samples = cap;
                        }
                        float* u = w->samples + 3 * w->num_samples++;
                        for (int b = 0; b < 3; ++b) {
                            u[b] = ((float)coord[b] + (b == a ? s : 0.0f) + 0.5f) * ratio[b] - 0.5f;
                        }
                    }
                }
            }
        }
    }
}

// The spread of the field over the samples that have values
static void surface_statistics(SurfaceFieldVolume* vol, const SurfaceFieldWork& w) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    float* samples = (float*)md_temp_alloc(temp, sizeof(float) * MAX(w.num_samples, (size_t)1));
    size_t count = 0;
    for (size_t i = 0; samples && i < w.num_samples; ++i) {
        float v;
        if (field_sample(&v, *vol, w.samples + 3 * i)) samples[count++] = v;
    }
    vol->num_surface_samples = count;
    if (count == 0) {
        vol->surface_min = vol->surface_max = vol->surface_lo = vol->surface_hi = 0.0f;
        return;
    }
    const size_t lo_idx = (size_t)(0.01 * (double)(count - 1));
    const size_t hi_idx = (size_t)(0.99 * (double)(count - 1));
    std::nth_element(samples, samples + lo_idx, samples + count);
    vol->surface_lo = samples[lo_idx];
    std::nth_element(samples, samples + hi_idx, samples + count);
    vol->surface_hi = samples[hi_idx];
    vol->surface_min = *std::min_element(samples, samples + count);
    vol->surface_max = *std::max_element(samples, samples + count);
}

static void count_left(SurfaceFieldVolume* vol, SurfaceFieldWork* w) {
    size_t n = 0;
    for (size_t i = 0; i < vol->num_voxels; ++i) n += (w->want[i] && !bit_get(vol->done, i));
    w->num_left = n;
}

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------

// The voxels wanted that hold no value yet, on the CPU's worker threads: a pool task over copies of
// their positions
static bool start_cpu(SurfaceFieldVolume* vol, SurfaceFieldWork* w) {
    const size_t n = w->num_left;
    w->cpu_index = (uint32_t*)malloc(sizeof(uint32_t) * n);
    w->cpu_xyz   = (float*)malloc(sizeof(float) * 3 * n);
    w->cpu_out   = (double*)malloc(sizeof(double) * n);
    if (!w->cpu_index || !w->cpu_xyz || !w->cpu_out) {
        work_cpu_free(w);
        return false;
    }
    const md_grid_t& grid = vol->grid;
    const size_t nx = grid.dim[0], ny = grid.dim[1];
    size_t m = 0;
    for (size_t idx = 0; idx < vol->num_voxels && m < n; ++idx) {
        if (!w->want[idx] || bit_get(vol->done, idx)) continue;
        const size_t i = idx % nx, j = (idx / nx) % ny, k = idx / (nx * ny);
        const float l[3] = { ((float)i + 0.5f) * grid.spacing.x, ((float)j + 0.5f) * grid.spacing.y, ((float)k + 0.5f) * grid.spacing.z };
        for (int a = 0; a < 3; ++a) {
            // orientation columns are the grid axes in world space
            w->cpu_xyz[3 * m + a] = grid.origin.elem[a] + grid.orientation.elem[0][a] * l[0] + grid.orientation.elem[1][a] * l[1] + grid.orientation.elem[2][a] * l[2];
        }
        w->cpu_index[m++] = (uint32_t)idx;
    }
    w->cpu_count   = m;
    w->cpu_values = { w->source.charges, w->source.gaussians };
    // A point of a QM density's potential is ~0.1-1 ms of work, of a few point charges microseconds
    const uint64_t work = field_work_per_point(w->cpu_values);
    const uint32_t grain = (uint32_t)CLAMP(1000000 / MAX(work, (uint64_t)1), (uint64_t)8, (uint64_t)1024);
    w->cpu_cancel = false;
    w->cpu_task = task_system::create_pool_task(STR_LIT("##Surface Field"), (uint32_t)m, [w](uint32_t beg, uint32_t end, uint32_t) {
        evaluate_points_range(w->cpu_out, w->cpu_xyz, beg, end, w->cpu_values, w->cpu_cancel);
    }, grain);
    task_system::enqueue_task(w->cpu_task);
    w->stage = SurfaceFieldWork::Cpu;
    return true;
}

#if MD_ENABLE_GPU
// The voxels wanted that hold no value yet, as a job on the device into the work's own texture. With a
// far field only the 4x4x4 bricks holding such a voxel are evaluated, and every voxel of those is kept
// (they are right, and the band moves with the isovalues); without one the whole grid is. False when
// the device cannot serve it - point multipoles (which md_gto_int evaluates on the CPU only), no
// texture, no job - and the CPU takes over then.
static bool start_gpu(SurfaceFieldVolume* vol, SurfaceFieldWork* w) {
    md_gpu_stream_t stream = w->stream;
    md_gto_int_gpu_charges_t gq = source_gpu(&w->source, stream);
    if (!gq) return false;
    const md_grid_t& grid = vol->grid;
    const int nx = grid.dim[0], ny = grid.dim[1], nz = grid.dim[2];
    if (!w->tex || w->tex_dim[0] != nx || w->tex_dim[1] != ny || w->tex_dim[2] != nz) {
        if (w->tex) md_gpu_texture_destroy(w->tex);
        md_gpu_texture_desc_t td = {};
        td.type   = MD_GPU_TEX_3D;
        td.format = MD_GPU_FORMAT_R32_FLOAT;
        td.usage  = MD_GPU_TEX_STORAGE;
        td.width  = (uint32_t)nx;
        td.height = (uint32_t)ny;
        td.depth_or_layers = (uint32_t)nz;
        td.label  = "surface field";
        w->tex = md_gpu_texture_create(stream, &td);
        w->tex_dim[0] = nx; w->tex_dim[1] = ny; w->tex_dim[2] = nz;
        if (!w->tex) {
            MD_LOG_ERROR("Surface field: no %dx%dx%d texture on the device", nx, ny, nz);
            w->tex_dim[0] = w->tex_dim[1] = w->tex_dim[2] = 0;
            return false;
        }
    }

    // What to evaluate, and what will have been
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    uint8_t* mask = (uint8_t*)md_temp_alloc(temp, vol->num_voxels);
    if (!mask) return false;
    for (size_t i = 0; i < vol->num_voxels; ++i) mask[i] = (w->want[i] && !bit_get(vol->done, i)) ? 1 : 0;
    if (w->source.charges->num_nodes > 0) {
        const int bx = (nx + 3) / 4, by = (ny + 3) / 4, bz = (nz + 3) / 4;
        uint8_t* brick = (uint8_t*)md_temp_alloc(temp, (size_t)bx * by * bz);
        if (!brick) return false;
        MEMSET(brick, 0, (size_t)bx * by * bz);
        for (int z = 0; z < nz; ++z) for (int y = 0; y < ny; ++y) for (int x = 0; x < nx; ++x) {
            if (mask[((size_t)z * ny + y) * nx + x]) brick[((size_t)(z / 4) * by + y / 4) * bx + x / 4] = 1;
        }
        for (int z = 0; z < nz; ++z) for (int y = 0; y < ny; ++y) for (int x = 0; x < nx; ++x) {
            w->got[((size_t)z * ny + y) * nx + x] = brick[((size_t)(z / 4) * by + y / 4) * bx + x / 4];
        }
    } else {
        MEMSET(w->got, 1, vol->num_voxels);
    }

    md_gto_int_gpu_potential_desc_t pd = {};
    pd.charges          = gq;
    pd.out_tex          = w->tex;
    pd.grid             = &grid;
    pd.sample_offset[0] = pd.sample_offset[1] = pd.sample_offset[2] = 0.5f;   // voxel centres, as the GL texture
    pd.op               = MD_GTO_OP_SET;
    pd.voxel_mask       = mask;     // copied by the job
    w->job = md_gto_int_gpu_potential_job_begin(stream, &pd);
    if (!w->job) return false;
    w->stage = SurfaceFieldWork::Gpu;
    return true;
}

// The job is done: the texture to the host, looked for later
static bool start_readback(SurfaceFieldVolume* vol, SurfaceFieldWork* w) {
    const size_t bytes = sizeof(float) * vol->num_voxels;
    w->rb = md_gpu_malloc(w->stream, MD_GPU_MEM_HOST_READ, bytes);
    if (!w->rb.cpu) {
        MD_LOG_ERROR("Surface field: failed to allocate %zu bytes of readback", bytes);
        return false;
    }
    md_gpu_tex_region_t region = {};
    region.extent[0] = (uint32_t)vol->grid.dim[0];
    region.extent[1] = (uint32_t)vol->grid.dim[1];
    region.extent[2] = (uint32_t)vol->grid.dim[2];
    if (!md_gpu_copy_from_texture(w->stream, w->rb.gpu, w->tex, &region)) {
        md_gpu_free(w->stream, w->rb.gpu);
        w->rb = {};
        return false;
    }
    w->rb_sync = md_gpu_stream_record(w->stream);
    w->stage = SurfaceFieldWork::Readback;
    return true;
}
#endif

static void field_upload(SurfaceFieldVolume* vol) {
    const int* dim = vol->grid.dim;
    if (!vol->tex_id) {
        gl::init_texture_3D(&vol->tex_id, dim[0], dim[1], dim[2], GL_R32F);
    } else {
        int tex_dim[3] = {};
        gl::get_texture_dim(tex_dim, vol->tex_id);
        if (tex_dim[0] != dim[0] || tex_dim[1] != dim[1] || tex_dim[2] != dim[2]) {
            gl::init_texture_3D(&vol->tex_id, dim[0], dim[1], dim[2], GL_R32F);
        }
    }
    const size_t bytes = sizeof(float) * vol->num_voxels;
    if (void* dst = gl::pbo_upload_begin(bytes)) {
        MEMCPY(dst, vol->values, bytes);
        if (!gl::pbo_upload_end_texture_3D(vol->tex_id, 0, GL_R32F)) {
            gl::set_texture_3D_data(vol->tex_id, 0, vol->values, GL_R32F);
        }
    } else {
        gl::set_texture_3D_data(vol->tex_id, 0, vol->values, GL_R32F);
    }
    vol->num_uploaded = vol->num_evaluated;
}

// The values, the span over the surfaces, the texture: after new values or a new band
static void field_publish(SurfaceFieldVolume* vol, SurfaceFieldWork* w) {
    surface_statistics(vol, *w);
    if (!vol->tex_id || vol->num_uploaded != vol->num_evaluated) {
        field_upload(vol);
    }
}

// What the evaluation made, into the values
static void merge(SurfaceFieldVolume* vol, SurfaceFieldWork* w, const float* v, bool on_gpu) {
    size_t n = 0;
    for (size_t i = 0; i < vol->num_voxels; ++i) {
        if (!w->got[i] || bit_get(vol->done, i)) continue;
        vol->values[i] = v[i];
        bit_set(vol->done, i);
        n += 1;
    }
    vol->num_evaluated += n;
    vol->complete = vol->num_evaluated == vol->num_voxels;
    vol->on_gpu = on_gpu;
}

bool surface_field_poll(SurfaceFieldVolume* vol) {
    ASSERT(vol);
    SurfaceFieldWork* w = vol->work;
    if (!w) return false;
    const char* name = surface_field_kind_str[(int)w->kind];
    source_poll(&w->source, vol->source_hash, w->stream, name);

    bool changed = false;
    switch (w->stage) {
#if MD_ENABLE_GPU
    case SurfaceFieldWork::Gpu: {
        // waited for in its first frames only: a short evaluation ends within them, a long one costs
        // a frame no more than a few milliseconds
        const double age_ms = md_tick_to_milliseconds(md_tick_now() - w->t0);
        const md_gto_int_gpu_job_status_t st = md_gto_int_gpu_job_poll(w->job, age_ms < 50.0 ? 4.0 : 1.0);
        if (st == MD_GTO_INT_GPU_JOB_RUNNING) {
            progress_set(&w->progress, name, md_gto_int_gpu_job_progress(w->job), false);
            break;
        }
        md_gto_int_gpu_job_end(w->job);
        w->job = nullptr;
        if (st != MD_GTO_INT_GPU_JOB_DONE || !start_readback(vol, w)) {
            MD_LOG_ERROR("Surface field '%s': the evaluation on the GPU failed, the CPU takes over", name);
            w->gpu_failed = true;
            w->stage = SurfaceFieldWork::Idle;
        }
        break;
    }
    case SurfaceFieldWork::Readback:
        if (!md_gpu_sync_is_complete(w->rb_sync)) break;
        merge(vol, w, (const float*)w->rb.cpu, true);
        md_gpu_free(w->stream, w->rb.gpu);
        w->rb = {};
        w->stage = SurfaceFieldWork::Idle;
        changed = true;
        MD_LOG_DEBUG("Surface field '%s': %zu of %zu voxels evaluated on the GPU in %.1f ms (over frames)", name, vol->num_evaluated, vol->num_voxels,
                     md_tick_to_milliseconds(md_tick_now() - w->t0));
        break;
#endif
    case SurfaceFieldWork::Cpu: {
        if (task_system::task_is_running(w->cpu_task)) {
            progress_set(&w->progress, name, task_system::task_fraction_complete(w->cpu_task), false);
            break;
        }
        // the values for those voxels, as a field the size of the grid would be (merge reads got)
        MEMSET(w->got, 0, vol->num_voxels);
        for (size_t m = 0; m < w->cpu_count; ++m) {
            vol->values[w->cpu_index[m]] = (float)w->cpu_out[m];
        }
        size_t n = 0;
        for (size_t m = 0; m < w->cpu_count; ++m) {
            if (!bit_get(vol->done, w->cpu_index[m])) {
                bit_set(vol->done, w->cpu_index[m]);
                n += 1;
            }
        }
        vol->num_evaluated += n;
        vol->complete = vol->num_evaluated == vol->num_voxels;
        vol->on_gpu = false;
        MD_LOG_DEBUG("Surface field '%s': %zu voxels evaluated on the CPU in %.1f ms (on worker threads)", name, w->cpu_count,
                     md_tick_to_milliseconds(md_tick_now() - w->t0));
        work_cpu_free(w);
        w->stage = SurfaceFieldWork::Idle;
        changed = true;
        break;
    }
    default:
        break;
    }

    if (changed) count_left(vol, w);

    // What is left of the band, once the distribution is there
    if (w->stage == SurfaceFieldWork::Idle && w->num_left > 0 && source_ready(w->source, vol->source_hash)) {
        w->t0 = md_tick_now();
        bool started = false;
#if MD_ENABLE_GPU
        if (!w->gpu_failed && w->stream) {
            started = start_gpu(vol, w);
            if (!started) w->gpu_failed = true;
        }
#endif
        if (!started) start_cpu(vol, w);
    }

    if (changed) {
        field_publish(vol, w);
        // Everything there is to evaluate is: the distribution has no further use
        if (vol->complete && w->stage == SurfaceFieldWork::Idle) source_drop_charges(&w->source, w->stream);
    }

    if (w->stage != SurfaceFieldWork::Idle) {
        // set above
    } else if (source_building(w->source) && w->num_left > 0) {
        progress_set(&w->progress, name, 0.0f, true);
    } else {
        progress_clear(&w->progress);
    }
    return changed;
}

bool surface_field_busy(const SurfaceFieldVolume& vol) {
    const SurfaceFieldWork* w = vol.work;
    if (!w) return false;
    if (w->stage != SurfaceFieldWork::Idle || source_building(w->source)) return true;
    // values wanted that will come once the distribution does
    return w->num_left > 0 && !source_failed(w->source, vol.source_hash) && source_ready(w->source, vol.source_hash);
}

bool surface_field_update(SurfaceFieldVolume* vol, const SurfaceFieldDesc& desc) {
    ASSERT(vol);
    ASSERT(desc.sys);
    ASSERT(desc.grid);
    ASSERT(desc.density);

    if (!vol->work) vol->work = new SurfaceFieldWork();
    SurfaceFieldWork* w = vol->work;
    if (desc.gpu_stream) w->stream = desc.gpu_stream;
    w->kind = desc.kind;

    const md_grid_t& dgrid = *desc.grid;
    const md_grid_t  fgrid = surface_field_grid(desc.kind, dgrid);

    if (desc.source_hash != vol->source_hash || !vol->values || !grid_same(vol->grid, fgrid)) {
        // What is under way is of the old source; the texture keeps its values until new ones come
        work_cancel(w);
        field_volume_reset(vol, fgrid);
        vol->source_hash = desc.source_hash;
        vol->band_hash = 0;
        w->gpu_failed = false;      // that was of the old source (point multipoles, say)
    }
    if (source_failed(w->source, desc.source_hash)) {
        vol->band_hash = desc.band_hash;   // nothing to do for it: not asked again
        return false;
    }
    if (desc.band_hash == vol->band_hash) {
        return true;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    // The band on the density grid, which the surfaces are drawn from, and the voxels of the field grid
    // it needs
    const size_t num_dvox = md_grid_num_points(&dgrid);
    uint8_t* dmask = (uint8_t*)md_temp_alloc(temp, num_dvox);
    band_mask(dmask, desc.density, dgrid.dim, desc.iso_values, desc.num_iso);
    const bool same_dim = vol->grid.dim[0] == dgrid.dim[0] && vol->grid.dim[1] == dgrid.dim[1] && vol->grid.dim[2] == dgrid.dim[2];
    if (same_dim) {
        MEMCPY(w->want, dmask, vol->num_voxels);
    } else {
        band_to_field(w->want, dmask, dgrid.dim, vol->grid.dim, temp);
    }
    surface_samples_collect(w, *vol, dgrid.dim, desc.density, dmask, desc.iso_values, desc.num_iso);
    count_left(vol, w);
    vol->band_hash = desc.band_hash;

    if (w->num_left > 0) source_request(&w->source, desc);
    surface_field_poll(vol);    // starts what can be at once
    field_publish(vol, w);      // the span over the new band, of what values there are
    return true;
}

void surface_field_invalidate(SurfaceFieldVolume* vol) {
    ASSERT(vol);
    vol->source_hash = 0;
    vol->band_hash = 0;
    if (vol->work) {
        work_cancel(vol->work);
        source_free(&vol->work->source, vol->work->stream);
    }
}

ColorScaleSpan surface_field_span(const SurfaceFieldVolume& vol) {
    ColorScaleSpan span;
    if (vol.num_surface_samples == 0) return span;
    span.valid = true;
    span.min = vol.surface_min;
    span.max = vol.surface_max;
    span.lo  = vol.surface_lo;
    span.hi  = vol.surface_hi;
    return span;
}

uint32_t surface_field_colormap_texture(SurfaceFieldVolume* vol, int colormap) {
    ASSERT(vol);
    if (vol->colormap_tex && vol->colormap == colormap) {
        return vol->colormap_tex;
    }
    const int res = 256;
    uint32_t pixels[256];
    for (int i = 0; i < res; ++i) {
        ImVec4 col = ImPlot::SampleColormap((float)i / (float)(res - 1), colormap);
        col.w = 1.0f;
        pixels[i] = ImGui::ColorConvertFloat4ToU32(col);
    }
    if (!vol->colormap_tex) {
        gl::init_texture_2D(&vol->colormap_tex, res, 1, GL_RGBA8);
    }
    gl::set_texture_2D_data(vol->colormap_tex, 0, pixels, GL_RGBA8);
    vol->colormap = colormap;
    return vol->colormap_tex;
}

void surface_field_free(SurfaceFieldVolume* vol) {
    ASSERT(vol);
    md_allocator_i* heap = md_get_heap_allocator();
    work_destroy(vol->work);
    if (vol->values) md_free(heap, vol->values, sizeof(float) * vol->num_voxels);
    if (vol->done)   md_free(heap, vol->done, sizeof(uint64_t) * ((vol->num_voxels + 63) / 64));
    if (vol->tex_id)       gl::free_texture(&vol->tex_id);
    if (vol->colormap_tex) gl::free_texture(&vol->colormap_tex);
    *vol = SurfaceFieldVolume{};
}

// ---------------------------------------------------------------------------
// Field at points
// ---------------------------------------------------------------------------
struct SurfaceFieldPointsWork {
    md_gpu_stream_t  stream = nullptr;
    SurfaceFieldKind kind = SurfaceFieldKind::EmbeddingPotential;
    FieldSource      source;
    FieldProgress    progress;
    uint64_t         source_hash = 0;   // what the evaluation under way is of

    enum Stage { None, Charges, Gpu, Readback, Cpu, Done, Failed } stage = None;
    md_tick_t t0 = 0;
    float*    xyz = nullptr;            // the points on the host, packed float3, when given there
    size_t    count = 0;
    float*    values = nullptr;         // [count] once done
#if MD_ENABLE_GPU
    md_gpu_mem_t gpu_xyz = {};          // the points on the device, when given there
    uint32_t     gpu_stride = 0;
    md_gpu_mem_t gpu_out = {};
    md_gpu_mem_t rb = {};
    md_gpu_sync_t rb_sync = {};
    md_gto_int_gpu_job_t* job = nullptr;
#endif
    task_system::ID cpu_task = task_system::INVALID_ID;
    std::atomic<bool> cpu_cancel{false};
    double*   cpu_out = nullptr;
    bool      on_gpu = false;
};

// Stops the evaluation under way and drops the points
static void points_cancel(SurfaceFieldPointsWork* w) {
#if MD_ENABLE_GPU
    if (w->job) {
        md_gto_int_gpu_job_end(w->job);
        w->job = nullptr;
    }
    // stream ordered: after whatever reads or writes them
    if (w->rb.gpu)      md_gpu_free(w->stream, w->rb.gpu);
    if (w->gpu_out.gpu) md_gpu_free(w->stream, w->gpu_out.gpu);
    if (w->gpu_xyz.gpu) md_gpu_free(w->stream, w->gpu_xyz.gpu);
    w->rb = {};
    w->gpu_out = {};
    w->gpu_xyz = {};
#endif
    if (w->cpu_task) {
        w->cpu_cancel = true;
        task_system::task_interrupt_and_wait_for(w->cpu_task);
    }
    w->cpu_task = task_system::INVALID_ID;
    free(w->cpu_out);
    w->cpu_out = nullptr;
    free(w->xyz);
    w->xyz = nullptr;
    w->stage = SurfaceFieldPointsWork::None;
    progress_clear(&w->progress);
}

static bool points_start_cpu(SurfaceFieldPointsWork* w) {
    if (!w->xyz) return false;
    w->cpu_out = (double*)malloc(sizeof(double) * w->count);
    if (!w->cpu_out) return false;
    const FieldValues v = { w->source.charges, w->source.gaussians };
    const uint64_t work = field_work_per_point(v);
    const uint32_t grain = (uint32_t)CLAMP(1000000 / MAX(work, (uint64_t)1), (uint64_t)8, (uint64_t)1024);
    w->cpu_cancel = false;
    w->cpu_task = task_system::create_pool_task(STR_LIT("##Surface Field Points"), (uint32_t)w->count, [w, v](uint32_t beg, uint32_t end, uint32_t) {
        evaluate_points_range(w->cpu_out, w->xyz, beg, end, v, w->cpu_cancel);
    }, grain);
    task_system::enqueue_task(w->cpu_task);
    w->stage = SurfaceFieldPointsWork::Cpu;
    w->on_gpu = false;
    return true;
}

#if MD_ENABLE_GPU
static bool points_start_gpu(SurfaceFieldPointsWork* w) {
    md_gto_int_gpu_charges_t gq = source_gpu(&w->source, w->stream);
    if (!gq) return false;
    md_gpu_addr_t xyz = w->gpu_xyz.gpu;
    uint32_t stride = w->gpu_stride;
    if (!xyz) {
        if (!w->xyz) return false;
        // not on the device: up they go
        const size_t bytes = sizeof(float) * 3 * w->count;
        w->gpu_xyz = md_gpu_malloc(w->stream, MD_GPU_MEM_DEVICE, bytes);
        void* dst = w->gpu_xyz.gpu ? md_gpu_upload_begin(w->stream, w->gpu_xyz.gpu, bytes) : nullptr;
        if (!dst) return false;
        MEMCPY(dst, w->xyz, bytes);
        md_gpu_upload_end(w->stream);
        xyz = w->gpu_xyz.gpu;
        stride = 0;
    }
    if (!w->gpu_out.gpu) w->gpu_out = md_gpu_malloc(w->stream, MD_GPU_MEM_DEVICE, sizeof(float) * w->count);
    if (!w->gpu_out.gpu) return false;
    md_gto_int_gpu_potential_xyz_desc_t pd = {};
    pd.charges       = gq;
    pd.xyz           = xyz;
    pd.xyz_stride    = stride;
    pd.num_xyz       = (uint32_t)w->count;
    pd.out_potential = w->gpu_out.gpu;
    w->job = md_gto_int_gpu_potential_xyz_job_begin(w->stream, &pd);
    if (!w->job) return false;
    w->stage = SurfaceFieldPointsWork::Gpu;
    w->on_gpu = true;
    return true;
}
#endif

// The distribution is there: on the device where it can be, else on the CPU
static void points_start(SurfaceFieldPointsWork* w) {
#if MD_ENABLE_GPU
    if (w->stream && points_start_gpu(w)) return;
#endif
    if (!points_start_cpu(w)) {
        MD_LOG_ERROR("Surface field: the points cannot be evaluated where they are");
        w->stage = SurfaceFieldPointsWork::Failed;
    }
}

bool surface_field_points_begin(SurfaceFieldPoints* fp, const SurfaceFieldPointSet& points, const SurfaceFieldDesc& desc) {
    ASSERT(fp);
    ASSERT(desc.sys);
    if (!fp->work) fp->work = new SurfaceFieldPointsWork();
    SurfaceFieldPointsWork* w = fp->work;
    if (desc.gpu_stream) w->stream = desc.gpu_stream;
    w->kind = desc.kind;
    points_cancel(w);
    free(w->values);
    w->values = nullptr;
    w->count = points.count;
    w->source_hash = desc.source_hash;
    w->t0 = md_tick_now();
    if (points.count == 0) {
        w->stage = SurfaceFieldPointsWork::Done;
        return true;
    }
    if (source_failed(w->source, desc.source_hash)) return false;

    // The points, copied: the mesh they come from may go before the evaluation is done
    if (points.xyz) {
        w->xyz = (float*)malloc(sizeof(float) * 3 * points.count);
        if (w->xyz) MEMCPY(w->xyz, points.xyz, sizeof(float) * 3 * points.count);
    }
#if MD_ENABLE_GPU
    if (points.gpu_xyz && w->stream) {
        const uint32_t stride = points.gpu_stride ? points.gpu_stride : 12u;
        const size_t bytes = (size_t)stride * points.count;
        w->gpu_xyz = md_gpu_malloc(w->stream, MD_GPU_MEM_DEVICE, bytes);
        if (w->gpu_xyz.gpu && md_gpu_copy(w->stream, w->gpu_xyz.gpu, points.gpu_xyz, bytes)) {
            w->gpu_stride = stride;
        } else if (w->gpu_xyz.gpu) {
            md_gpu_free(w->stream, w->gpu_xyz.gpu);
            w->gpu_xyz = {};
        }
    }
    const bool have = w->xyz || w->gpu_xyz.gpu;
#else
    const bool have = w->xyz != nullptr;
#endif
    if (!have) return false;

    source_request(&w->source, desc);
    w->stage = SurfaceFieldPointsWork::Charges;
    return surface_field_points_poll(fp, nullptr, nullptr) != SurfaceFieldPointsStatus::Failed;   // starts what can be at once
}

SurfaceFieldPointsStatus surface_field_points_poll(SurfaceFieldPoints* fp, const float** values, size_t* count) {
    ASSERT(fp);
    SurfaceFieldPointsWork* w = fp->work;
    if (!w) return SurfaceFieldPointsStatus::None;
    const char* name = surface_field_kind_str[(int)w->kind];
    source_poll(&w->source, w->source_hash, w->stream, name);

    switch (w->stage) {
    case SurfaceFieldPointsWork::Charges:
        if (source_failed(w->source, w->source_hash)) {
            w->stage = SurfaceFieldPointsWork::Failed;
        } else if (source_ready(w->source, w->source_hash)) {
            points_start(w);
        } else if (!source_building(w->source)) {
            w->stage = SurfaceFieldPointsWork::Failed;   // nothing will come
        }
        break;
#if MD_ENABLE_GPU
    case SurfaceFieldPointsWork::Gpu: {
        const double age_ms = md_tick_to_milliseconds(md_tick_now() - w->t0);
        const md_gto_int_gpu_job_status_t st = md_gto_int_gpu_job_poll(w->job, age_ms < 50.0 ? 4.0 : 1.0);
        if (st == MD_GTO_INT_GPU_JOB_RUNNING) break;
        md_gto_int_gpu_job_end(w->job);
        w->job = nullptr;
        bool ok = st == MD_GTO_INT_GPU_JOB_DONE;
        if (ok) {
            w->rb = md_gpu_malloc(w->stream, MD_GPU_MEM_HOST_READ, sizeof(float) * w->count);
            ok = w->rb.cpu && md_gpu_copy(w->stream, w->rb.gpu, w->gpu_out.gpu, sizeof(float) * w->count);
        }
        if (ok) {
            w->rb_sync = md_gpu_stream_record(w->stream);
            w->stage = SurfaceFieldPointsWork::Readback;
        } else {
            MD_LOG_ERROR("Surface field '%s': the evaluation at the points failed on the GPU, the CPU takes over", name);
            if (!points_start_cpu(w)) w->stage = SurfaceFieldPointsWork::Failed;
        }
        break;
    }
    case SurfaceFieldPointsWork::Readback:
        if (!md_gpu_sync_is_complete(w->rb_sync)) break;
        w->values = (float*)malloc(sizeof(float) * MAX(w->count, (size_t)1));
        if (w->values) MEMCPY(w->values, w->rb.cpu, sizeof(float) * w->count);
        md_gpu_free(w->stream, w->rb.gpu);
        w->rb = {};
        w->stage = w->values ? SurfaceFieldPointsWork::Done : SurfaceFieldPointsWork::Failed;
        break;
#endif
    case SurfaceFieldPointsWork::Cpu:
        if (task_system::task_is_running(w->cpu_task)) break;
        w->values = (float*)malloc(sizeof(float) * MAX(w->count, (size_t)1));
        if (w->values) for (size_t i = 0; i < w->count; ++i) w->values[i] = (float)w->cpu_out[i];
        free(w->cpu_out);
        w->cpu_out = nullptr;
        w->cpu_task = task_system::INVALID_ID;
        w->stage = w->values ? SurfaceFieldPointsWork::Done : SurfaceFieldPointsWork::Failed;
        break;
    default:
        break;
    }

    // progress
    switch (w->stage) {
    case SurfaceFieldPointsWork::Charges:
        progress_set(&w->progress, name, 0.0f, true);
        break;
#if MD_ENABLE_GPU
    case SurfaceFieldPointsWork::Gpu:
        progress_set(&w->progress, name, md_gto_int_gpu_job_progress(w->job), false);
        break;
#endif
    case SurfaceFieldPointsWork::Cpu:
        progress_set(&w->progress, name, task_system::task_fraction_complete(w->cpu_task), false);
        break;
    default:
        progress_clear(&w->progress);
        break;
    }

    if (w->stage == SurfaceFieldPointsWork::Done) {
        // delivered once
        fp->on_gpu  = w->on_gpu;
        fp->seconds = md_tick_to_seconds(md_tick_now() - w->t0);
        if (values) *values = w->values;
        if (count)  *count = w->values ? w->count : 0;
        if (values || count) {
            w->stage = SurfaceFieldPointsWork::None;
#if MD_ENABLE_GPU
            // the points and the output on the device have done their part
            if (w->gpu_out.gpu) md_gpu_free(w->stream, w->gpu_out.gpu);
            if (w->gpu_xyz.gpu) md_gpu_free(w->stream, w->gpu_xyz.gpu);
            w->gpu_out = {};
            w->gpu_xyz = {};
#endif
            free(w->xyz);
            w->xyz = nullptr;
        }
        return SurfaceFieldPointsStatus::Done;
    }
    if (w->stage == SurfaceFieldPointsWork::Failed) {
        w->stage = SurfaceFieldPointsWork::None;
        progress_clear(&w->progress);
        return SurfaceFieldPointsStatus::Failed;
    }
    return w->stage == SurfaceFieldPointsWork::None ? SurfaceFieldPointsStatus::None : SurfaceFieldPointsStatus::Running;
}

void surface_field_points_free(SurfaceFieldPoints* fp) {
    ASSERT(fp);
    if (SurfaceFieldPointsWork* w = fp->work) {
        points_cancel(w);
        source_free(&w->source, w->stream);
        free(w->values);
        delete w;
    }
    *fp = SurfaceFieldPoints{};
}

ColorScaleSpan surface_field_span_of(const float* values, size_t count) {
    ColorScaleSpan span;
    if (!values || count == 0) return span;
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    float* v = (float*)md_temp_alloc(temp, sizeof(float) * count);
    MEMCPY(v, values, sizeof(float) * count);
    const size_t lo_idx = (size_t)(0.01 * (double)(count - 1));
    const size_t hi_idx = (size_t)(0.99 * (double)(count - 1));
    span.valid = true;
    span.min = *std::min_element(v, v + count);
    span.max = *std::max_element(v, v + count);
    std::nth_element(v, v + lo_idx, v + count);
    span.lo = v[lo_idx];
    std::nth_element(v, v + hi_idx, v + count);
    span.hi = v[hi_idx];
    return span;
}
