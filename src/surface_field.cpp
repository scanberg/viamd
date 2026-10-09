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
    case SurfaceFieldKind::EmbeddingPotential:
    case SurfaceFieldKind::ElectrostaticPotential:
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

// The environment's classical sites as points with their charge, and their dipole and quadrupole
// where the potential has them
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
    const bool has_mu = extract_atom_column(mu, sys, STR_LIT("atom/dipole"),     3, md_unit_elementary_charge_bohr());
    const bool has_Q  = extract_atom_column(Q,  sys, STR_LIT("atom/quadrupole"), 6, e_bohr2);

    float*  xyz    = (float*)md_temp_alloc(temp, sizeof(float) * 3 * num_atoms);
    double* charge = (double*)md_temp_alloc(temp, sizeof(double) * num_atoms);
    double* dipole = has_mu ? (double*)md_temp_alloc(temp, sizeof(double) * 3 * num_atoms) : nullptr;
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
    if (!md_gto_int_charges_init(out, &cd, alloc)) {
        MD_LOG_ERROR("Electrostatic potential: could not build the charge distribution");
        return false;
    }
    esp_check_net_charge(out, sys);
    return true;
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
    default:
        return false;
    }
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

static void field_charges_free(SurfaceFieldVolume* vol) {
    if (vol->charges) {
        md_allocator_i* heap = md_get_heap_allocator();
        md_gto_int_charges_free(vol->charges, heap);
        md_free(heap, vol->charges, sizeof(md_gto_int_charges_t));
        vol->charges = nullptr;
    }
}

static void field_volume_reset(SurfaceFieldVolume* vol, const md_grid_t& grid) {
    const size_t n = md_grid_num_points(&grid);
    md_allocator_i* heap = md_get_heap_allocator();
    if (n != vol->num_voxels) {
        if (vol->values) md_free(heap, vol->values, sizeof(float) * vol->num_voxels);
        if (vol->done)   md_free(heap, vol->done, sizeof(uint64_t) * ((vol->num_voxels + 63) / 64));
        vol->values = (float*)md_alloc(heap, sizeof(float) * n);
        vol->done   = (uint64_t*)md_alloc(heap, sizeof(uint64_t) * ((n + 63) / 64));
        vol->num_voxels = n;
    }
    MEMSET(vol->values, 0, sizeof(float) * n);
    MEMSET(vol->done, 0, sizeof(uint64_t) * ((n + 63) / 64));
    vol->grid = grid;
    vol->num_evaluated = 0;
    vol->num_uploaded = 0;
    vol->complete = false;
    vol->on_gpu = false;
    vol->num_surface_samples = 0;
    field_charges_free(vol);
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

// The field interpolated to where an isosurface crosses an edge of the density grid: samples of the
// field ON the surface, at the density grid's resolution
static void surface_statistics(SurfaceFieldVolume* vol, const int ddim[3], const float* d, const uint8_t* dmask, const float* iso, size_t num_iso) {
    const size_t nx = ddim[0], ny = ddim[1], nz = ddim[2];
    const size_t step[3] = { 1, nx, nx * ny };
    const size_t len[3]  = { nx, ny, nz };
    const float  ratio[3] = {
        (float)vol->grid.dim[0] / (float)ddim[0],
        (float)vol->grid.dim[1] / (float)ddim[1],
        (float)vol->grid.dim[2] / (float)ddim[2],
    };

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t cap = 1024, count = 0;
    float* samples = (float*)md_temp_alloc(temp, sizeof(float) * cap);

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
                        float u[3];
                        for (int b = 0; b < 3; ++b) {
                            u[b] = ((float)coord[b] + (b == a ? s : 0.0f) + 0.5f) * ratio[b] - 0.5f;
                        }
                        float v;
                        if (!field_sample(&v, *vol, u)) continue;
                        if (count == cap) {
                            float* bigger = (float*)md_temp_alloc(temp, sizeof(float) * cap * 2);
                            MEMCPY(bigger, samples, sizeof(float) * cap);
                            samples = bigger;
                            cap *= 2;
                        }
                        samples[count++] = v;
                    }
                }
            }
        }
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

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------

// The voxels of fmask that hold no value yet, on the CPU's worker threads
static void evaluate_band_cpu(SurfaceFieldVolume* vol, const uint8_t* fmask, size_t num_new) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    uint32_t* index = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_new);
    float*    xyz   = (float*)md_temp_alloc(temp, sizeof(float) * 3 * num_new);
    double*   out   = (double*)md_temp_alloc(temp, sizeof(double) * num_new);

    const md_grid_t& grid = vol->grid;
    const size_t nx = grid.dim[0], ny = grid.dim[1];
    size_t n = 0;
    for (size_t idx = 0; idx < vol->num_voxels; ++idx) {
        if (!fmask[idx] || bit_get(vol->done, idx)) continue;
        const size_t i = idx % nx, j = (idx / nx) % ny, k = idx / (nx * ny);
        const float l[3] = { ((float)i + 0.5f) * grid.spacing.x, ((float)j + 0.5f) * grid.spacing.y, ((float)k + 0.5f) * grid.spacing.z };
        for (int a = 0; a < 3; ++a) {
            // orientation columns are the grid axes in world space
            xyz[3 * n + a] = grid.origin.elem[a] + grid.orientation.elem[0][a] * l[0] + grid.orientation.elem[1][a] * l[1] + grid.orientation.elem[2][a] * l[2];
        }
        index[n++] = (uint32_t)idx;
    }
    ASSERT(n == num_new);

    struct Payload {
        const md_gto_int_charges_t* charges;
        const float* xyz;
        double* out;
    } payload = { vol->charges, xyz, out };

    // A point of a QM density's potential is ~0.1-1 ms of work, of a few point charges microseconds
    const uint64_t work = md_gto_int_charges_work_per_point(vol->charges);
    const uint32_t grain = (uint32_t)CLAMP(1000000 / MAX(work, (uint64_t)1), (uint64_t)8, (uint64_t)1024);
    task_system::ID task = task_system::create_pool_task(STR_LIT("## Surface Field"), (uint32_t)num_new, [p = &payload](uint32_t beg, uint32_t end, uint32_t) {
        md_gto_int_potential_xyz(p->out + beg, NULL, p->xyz + 3 * (size_t)beg, end - beg, 0, p->charges);
    }, grain);
    task_system::enqueue_task(task);
    task_system::task_wait_for(task);

    for (size_t m = 0; m < num_new; ++m) {
        vol->values[index[m]] = (float)out[m];
        bit_set(vol->done, index[m]);
    }
    vol->num_evaluated += num_new;
}

#if MD_ENABLE_GPU
// The whole field grid on the device, read back. False when the device cannot serve it - no device,
// a grid larger than the scratch, point multipoles (which md_gto_int evaluates on the CPU only) - or
// when it fails, and the CPU takes over then.
static bool evaluate_all_gpu(SurfaceFieldVolume* vol, const SurfaceFieldDesc& desc) {
    md_gpu_stream_t stream = desc.gpu_stream;
    if (!stream || !desc.gpu_scratch || !vol->charges) return false;
    if (vol->charges->point_dipole || vol->charges->point_quadrupole) return false;
    const md_grid_t& grid = vol->grid;
    for (int a = 0; a < 3; ++a) {
        if (grid.dim[a] <= 0 || (uint32_t)grid.dim[a] > desc.gpu_scratch_dim[a]) return false;
    }

    md_gto_int_gpu_charges_t gq = md_gto_int_gpu_charges_create(stream, vol->charges);
    if (!gq) return false;

    md_gto_int_gpu_potential_desc_t pd = {};
    pd.charges          = gq;
    pd.out_tex          = desc.gpu_scratch;
    pd.grid             = &grid;
    pd.sample_offset[0] = pd.sample_offset[1] = pd.sample_offset[2] = 0.5f;   // voxel centres, as the GL texture
    pd.op               = MD_GTO_OP_SET;

    bool ok = false;
    if (md_gto_int_gpu_potential_launch(stream, &pd)) {
        const size_t bytes = sizeof(float) * vol->num_voxels;
        md_gpu_mem_t rb = md_gpu_malloc(stream, MD_GPU_MEM_HOST_READ, bytes);
        if (rb.cpu) {
            md_gpu_tex_region_t region = {};
            region.extent[0] = (uint32_t)grid.dim[0];
            region.extent[1] = (uint32_t)grid.dim[1];
            region.extent[2] = (uint32_t)grid.dim[2];
            if (md_gpu_copy_from_texture(stream, rb.gpu, desc.gpu_scratch, &region)) {
                md_gpu_stream_sync(stream);
                MEMCPY(vol->values, rb.cpu, bytes);
                ok = true;
            }
            md_gpu_free(stream, rb.gpu);
        } else {
            MD_LOG_ERROR("Surface field: failed to allocate %zu bytes of readback", bytes);
        }
    }
    md_gto_int_gpu_charges_destroy(stream, gq);
    if (!ok) return false;

    MEMSET(vol->done, 0xFF, sizeof(uint64_t) * ((vol->num_voxels + 63) / 64));
    vol->num_evaluated = vol->num_voxels;
    vol->complete = true;
    vol->on_gpu = true;
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

bool surface_field_update(SurfaceFieldVolume* vol, const SurfaceFieldDesc& desc) {
    ASSERT(vol);
    ASSERT(desc.sys);
    ASSERT(desc.grid);
    ASSERT(desc.density);

    const md_grid_t& dgrid = *desc.grid;
    const md_grid_t  fgrid = surface_field_grid(desc.kind, dgrid);

    if (desc.source_hash != vol->source_hash || !vol->values || !grid_same(vol->grid, fgrid)) {
        field_volume_reset(vol, fgrid);
        vol->source_hash = desc.source_hash;
        vol->band_hash = 0;
    }
    if (desc.band_hash == vol->band_hash && vol->tex_id) {
        return true;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    // The band on the density grid, which the surfaces are drawn from
    const size_t num_dvox = md_grid_num_points(&dgrid);
    uint8_t* dmask = (uint8_t*)md_temp_alloc(temp, num_dvox);
    band_mask(dmask, desc.density, dgrid.dim, desc.iso_values, desc.num_iso);

    if (!vol->complete) {
        // ...and the voxels of the field grid it needs
        const bool same_dim = vol->grid.dim[0] == dgrid.dim[0] && vol->grid.dim[1] == dgrid.dim[1] && vol->grid.dim[2] == dgrid.dim[2];
        uint8_t* fmask = dmask;
        if (!same_dim) {
            fmask = (uint8_t*)md_temp_alloc(temp, vol->num_voxels);
            band_to_field(fmask, dmask, dgrid.dim, vol->grid.dim, temp);
        }

        size_t num_new = 0;
        for (size_t i = 0; i < vol->num_voxels; ++i) {
            num_new += (fmask[i] && !bit_get(vol->done, i));
        }

        if (num_new > 0) {
            if (!vol->charges) {
                const md_tick_t t0 = md_tick_now();
                md_allocator_i* heap = md_get_heap_allocator();
                vol->charges = (md_gto_int_charges_t*)md_alloc(heap, sizeof(md_gto_int_charges_t));
                MEMSET(vol->charges, 0, sizeof(md_gto_int_charges_t));
                if (!build_charges(vol->charges, desc, heap)) {
                    md_free(heap, vol->charges, sizeof(md_gto_int_charges_t));
                    vol->charges = nullptr;
                    MD_LOG_ERROR("Surface field: the system has nothing to make '%s' of", surface_field_kind_str[(int)desc.kind]);
                    field_volume_reset(vol, fgrid);
                    return false;
                }
                MD_LOG_DEBUG("Surface field '%s': %u gaussians and %u points, %llu steps per voxel, built in %.1f ms", surface_field_kind_str[(int)desc.kind],
                             vol->charges->num_gaussians, vol->charges->num_points, (unsigned long long)md_gto_int_charges_work_per_point(vol->charges),
                             md_tick_to_milliseconds(md_tick_now() - t0));
            }

            const md_tick_t t0 = md_tick_now();
            bool on_gpu = false;
#if MD_ENABLE_GPU
            on_gpu = evaluate_all_gpu(vol, desc);
#endif
            if (!on_gpu) {
                evaluate_band_cpu(vol, fmask, num_new);
            }
            MD_LOG_DEBUG("Surface field '%s': %zu voxels of %dx%dx%d on the %s, %.1f ms", surface_field_kind_str[(int)desc.kind],
                         on_gpu ? vol->num_voxels : num_new, vol->grid.dim[0], vol->grid.dim[1], vol->grid.dim[2], on_gpu ? "GPU" : "CPU",
                         md_tick_to_milliseconds(md_tick_now() - t0));

            // Everything there is to evaluate is: the distribution has no further use
            if (vol->complete) {
                field_charges_free(vol);
            }
        }
    }

    surface_statistics(vol, dgrid.dim, desc.density, dmask, desc.iso_values, desc.num_iso);

    if (!vol->tex_id || vol->num_uploaded != vol->num_evaluated) {
        field_upload(vol);
    }

    vol->band_hash = desc.band_hash;
    return true;
}

void surface_field_invalidate(SurfaceFieldVolume* vol) {
    ASSERT(vol);
    vol->source_hash = 0;
    vol->band_hash = 0;
    field_charges_free(vol);
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
    field_charges_free(vol);
    if (vol->values) md_free(heap, vol->values, sizeof(float) * vol->num_voxels);
    if (vol->done)   md_free(heap, vol->done, sizeof(uint64_t) * ((vol->num_voxels + 63) / 64));
    if (vol->tex_id)       gl::free_texture(&vol->tex_id);
    if (vol->colormap_tex) gl::free_texture(&vol->colormap_tex);
    *vol = SurfaceFieldVolume{};
}
