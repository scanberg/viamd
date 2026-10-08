#include <surface_field.h>

#include <task_system.h>
#include <gfx/gl.h>
#include <gfx/gl_utils.h>

#include <md_system.h>
#include <md_attributes.h>
#include <md_gto_int.h>
#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>

#include <imgui.h>
#include <implot.h>

#include <algorithm>
#include <math.h>
#include <string.h>

#define ANGSTROM_TO_BOHR_D 1.8897261246257702

// Read on the bits: viamd is built with fast math, where NAN need not compare as itself
static inline bool is_finite_f64(double v) {
    uint64_t u;
    MEMCPY(&u, &v, sizeof(u));
    return (u & 0x7FF0000000000000ull) != 0x7FF0000000000000ull;
}

md_unit_t surface_field_unit(SurfaceFieldKind kind) {
    switch (kind) {
    case SurfaceFieldKind::EmbeddingPotential:
    default:
        return md_unit_div(md_unit_hartree(), md_unit_elementary_charge());
    }
}

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

bool surface_field_available(SurfaceFieldKind kind, const md_system_t& sys) {
    switch (kind) {
    case SurfaceFieldKind::EmbeddingPotential:
        return find_atom_charges(sys) != nullptr && has_environment(sys);
    default:
        return false;
    }
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

// The field as a charge distribution md_gto_int evaluates: the environment's classical sites as points
// with their charge, and their dipole and quadrupole where the potential has them. The QM density and nuclei,
// for the full electrostatic potential, are the same structure with a basis.
static bool build_charges(md_gto_int_charges_t* out, SurfaceFieldKind kind, const md_system_t& sys, const md_system_state_t& state, md_allocator_i* alloc) {
    ASSERT(kind == SurfaceFieldKind::EmbeddingPotential);
    (void)kind;

    const md_attribute_t* attr = find_atom_charges(sys);
    const size_t num_atoms = md_system_atom_count(&sys);
    if (!attr || state.num_atoms < num_atoms || !state.xyz) {
        return false;
    }

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

    md_gto_int_charges_desc_t desc = {};
    desc.point_xyz        = xyz;
    desc.point_charge     = charge;
    desc.num_points       = n;
    desc.point_dipole     = dipole;
    desc.point_quadrupole = quad;
    return md_gto_int_charges_init(out, &desc, alloc);
}

static inline bool bit_get(const uint64_t* bits, size_t i) { return (bits[i >> 6] >> (i & 63)) & 1; }
static inline void bit_set(uint64_t* bits, size_t i)       { bits[i >> 6] |= (1ull << (i & 63)); }

static void field_volume_reset(SurfaceFieldVolume* vol, const int dim[3]) {
    const size_t n = (size_t)dim[0] * (size_t)dim[1] * (size_t)dim[2];
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
    MEMCPY(vol->dim, dim, sizeof(vol->dim));
    vol->num_evaluated = 0;
}

// The voxels a surface can sample: the 8 corners of every cell that some isovalue passes through,
// grown by one voxel, because the raycaster places a hit by interpolating along the ray between two
// samples, which can put it just across the face of the cell the surface actually crosses.
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

// The field interpolated to where an isosurface crosses an edge of the grid: samples of the field
// ON the surface, at the grid's resolution
static void surface_statistics(SurfaceFieldVolume* vol, const float* d, const float* iso, size_t num_iso) {
    const size_t nx = vol->dim[0], ny = vol->dim[1], nz = vol->dim[2];
    const size_t step[3] = { 1, nx, nx * ny };
    const size_t len[3]  = { nx, ny, nz };

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t cap = 1024, count = 0;
    float* samples = (float*)md_temp_alloc(temp, sizeof(float) * cap);

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                const size_t idx = (k * ny + j) * nx + i;
                if (!bit_get(vol->done, idx)) continue;
                const size_t coord[3] = { i, j, k };
                for (int a = 0; a < 3; ++a) {
                    if (coord[a] + 1 >= len[a]) continue;
                    const size_t nb = idx + step[a];
                    if (!bit_get(vol->done, nb)) continue;
                    const float d0 = d[idx], d1 = d[nb];
                    for (size_t t = 0; t < num_iso; ++t) {
                        if ((d0 - iso[t]) * (d1 - iso[t]) >= 0.0f) continue;
                        const float s = (iso[t] - d0) / (d1 - d0);
                        if (count == cap) {
                            float* bigger = (float*)md_temp_alloc(temp, sizeof(float) * cap * 2);
                            MEMCPY(bigger, samples, sizeof(float) * cap);
                            samples = bigger;
                            cap *= 2;
                        }
                        samples[count++] = vol->values[idx] + s * (vol->values[nb] - vol->values[idx]);
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

bool surface_field_update(SurfaceFieldVolume* vol, SurfaceFieldKind kind, const md_system_t& sys, const md_system_state_t& state,
                          const md_grid_t& grid, const float* density, const float* iso_values, size_t num_iso,
                          uint64_t source_hash, uint64_t band_hash) {
    ASSERT(vol);
    ASSERT(density);

    const bool same_grid = vol->dim[0] == grid.dim[0] && vol->dim[1] == grid.dim[1] && vol->dim[2] == grid.dim[2];
    if (source_hash != vol->source_hash || !same_grid || !vol->values) {
        field_volume_reset(vol, grid.dim);
        vol->source_hash = source_hash;
        vol->band_hash = 0;
    }
    if (band_hash == vol->band_hash && vol->tex_id) {
        return true;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    // The voxels of the band that do not hold a value yet
    uint8_t* mask = (uint8_t*)md_temp_alloc(temp, vol->num_voxels);
    band_mask(mask, density, vol->dim, iso_values, num_iso);

    size_t num_new = 0;
    for (size_t i = 0; i < vol->num_voxels; ++i) {
        num_new += (mask[i] && !bit_get(vol->done, i));
    }

    if (num_new > 0) {
        md_gto_int_charges_t charges = {};
        if (!build_charges(&charges, kind, sys, state, temp_alloc)) {
            MD_LOG_ERROR("Surface field: the system has nothing to make '%s' of", surface_field_kind_str[(int)kind]);
            field_volume_reset(vol, grid.dim);
            return false;
        }

        uint32_t* index = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_new);
        float*    xyz   = (float*)md_temp_alloc(temp, sizeof(float) * 3 * num_new);
        double*   out   = (double*)md_temp_alloc(temp, sizeof(double) * num_new);

        const size_t nx = vol->dim[0], ny = vol->dim[1];
        size_t n = 0;
        for (size_t idx = 0; idx < vol->num_voxels; ++idx) {
            if (!mask[idx] || bit_get(vol->done, idx)) continue;
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
        } payload = { &charges, xyz, out };

        task_system::ID task = task_system::create_pool_task(STR_LIT("## Surface Field"), (uint32_t)num_new, [p = &payload](uint32_t beg, uint32_t end, uint32_t) {
            md_gto_int_potential_xyz(p->out + beg, NULL, p->xyz + 3 * (size_t)beg, end - beg, 0, p->charges);
        }, 256);
        task_system::enqueue_task(task);
        task_system::task_wait_for(task);

        for (size_t m = 0; m < num_new; ++m) {
            vol->values[index[m]] = (float)out[m];
            bit_set(vol->done, index[m]);
        }
        vol->num_evaluated += num_new;
    }

    surface_statistics(vol, density, iso_values, num_iso);

    if (!vol->tex_id) {
        gl::init_texture_3D(&vol->tex_id, vol->dim[0], vol->dim[1], vol->dim[2], GL_R32F);
    } else {
        int tex_dim[3] = {};
        gl::get_texture_dim(tex_dim, vol->tex_id);
        if (tex_dim[0] != vol->dim[0] || tex_dim[1] != vol->dim[1] || tex_dim[2] != vol->dim[2]) {
            gl::init_texture_3D(&vol->tex_id, vol->dim[0], vol->dim[1], vol->dim[2], GL_R32F);
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

    vol->band_hash = band_hash;
    return true;
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
    if (vol->values) md_free(heap, vol->values, sizeof(float) * vol->num_voxels);
    if (vol->done)   md_free(heap, vol->done, sizeof(uint64_t) * ((vol->num_voxels + 63) / 64));
    if (vol->tex_id)       gl::free_texture(&vol->tex_id);
    if (vol->colormap_tex) gl::free_texture(&vol->colormap_tex);
    *vol = SurfaceFieldVolume{};
}
