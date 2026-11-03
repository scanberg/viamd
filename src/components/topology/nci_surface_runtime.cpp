#include "nci_surface_model.h"

#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <gfx/gl_utils.h>
#include <surface_field.h>

#include <cfloat>
#include <cmath>

namespace topology {

namespace {

constexpr double resolution_spacing[] = {0.3, 0.2, 0.13};
constexpr double grid_margin = 2.5;
constexpr size_t max_samples = 32u << 20;

void update_scatter_outline(NciSurface& n) {
    if (!nci_surface::scatter_density_outline(&n.outline, n.scatter, n.rho_max))
        MD_LOG_ERROR("NCI surface '%s': invalid scatter samples for density outline", n.name);
}

bool upload_surface_mesh(NciSurface& n) {
    const auto& m = n.source_mesh;
    const uint32_t nv = (uint32_t)m.value.size(), ni = (uint32_t)m.index.size();
    if (nv == 0 || ni == 0) {
        isomesh::mesh_free(&n.mesh);
        n.triangles = 0;
        n.span = {};
        return true;
    }
    if (!isomesh::mesh_upload(&n.mesh, m.position.data(), m.normal.data(), 3, nv, m.index.data(), ni) ||
        !isomesh::mesh_upload_values(&n.mesh, m.value.data(), nullptr, nv)) {
        MD_LOG_ERROR("NCI surface '%s': could not upload mesh", n.name);
        isomesh::mesh_free(&n.mesh);
        n.triangles = 0;
        n.span = {};
        return false;
    }
    n.triangles = ni / 3;
    n.span = surface_field_span_of(m.value.data(), m.value.size());
    color_scale_update_range(&n.scale, n.span);
    return true;
}

void poll_surface_job(NciSurface& n, bool context) {
    NciSurfaceJob* job = n.job;
    if (!job || !job->done.load(std::memory_order_acquire)) return;
    n.job = nullptr;
    if (job->key == n.want_key) {
        n.have_key = job->key;
        n.seconds = md_tick_to_seconds(md_tick_now() - job->start);
        n.cancelled = !job->ok && job->cancel.load(std::memory_order_relaxed) != 0;
        n.failed = !job->ok && !n.cancelled;
        if (n.failed) MD_LOG_ERROR("NCI surface '%s': extraction failed", n.name);
        n.scatter = job->ok ? std::move(job->scatter) : nci_surface::Scatter{};
        n.reference_failed = job->ok && !job->reference_ok;
        if (n.reference_failed) {
            MD_LOG_ERROR("NCI surface '%s': category-specific reference extraction failed", n.name);
            n.reference_scatters = {};
        } else {
            n.reference_scatters = job->ok ? std::move(job->reference_scatters)
                                           : std::array<nci_surface::Scatter, nci_surface::family_count>{};
        }
        update_scatter_outline(n);
        n.source_mesh = job->ok && !context ? std::move(job->mesh) : nci_surface::Mesh{};
        if (context) {
            n.span = surface_field_span_of(n.scatter.signed_rho.data(), n.scatter.signed_rho.size());
            color_scale_update_range(&n.scale, n.span);
        } else if (!upload_surface_mesh(n)) {
            n.failed = true;
        }
        if (job->ok) MD_LOG_DEBUG("NCI surface: %zu triangles in %.2f s", n.triangles, n.seconds);
    }
    nci_surface_job_destroy(job);
}

}  // namespace

double nci_resolution_spacing(int resolution) {
    return resolution_spacing[CLAMP(resolution, 0, (int)NciSurface::Res_Count - 1)];
}

bool upload_nci_surface_mesh(NciSurface& n) {
    return upload_surface_mesh(n);
}

void update_nci_surface(NciSurface& n, const NciSurfaceRuntimeInput& input, bool context) {
    poll_surface_job(n, context);
    if (!n.shown || (!context && !nci_scale_member(n))) {
        nci_surface_job_destroy(n.job);
        n.job = nullptr;
        n.want_key = 0;
        if (!n.saved && !n.family_preview && !n.full_system_family) n.anchors.clear();
        return;
    }
    if (n.editing) return;

    std::vector<double> near = context || n.saved || n.family_preview || n.full_system_family ? n.anchors : std::vector<double>{};
    if (!n.saved && !n.family_preview && !n.full_system_family && !context) {
        n.families = 0;
        n.environments.clear();
    }

    const double pad = 6.0 - grid_margin;
    double lo[3], hi[3];
    for (int a = 0; a < 3; ++a) {
        lo[a] = input.bounds_min[a] + pad;
        hi[a] = input.bounds_max[a] - pad;
    }
    if (n.scope != NciSurface::Scope_Everywhere) {
        if (near.empty()) {
            nci_surface_job_destroy(n.job);
            n.job = nullptr;
            n.want_key = 0;
            if (!n.saved && !n.family_preview && !n.full_system_family) n.anchors.clear();
            return;
        }
        double nlo[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, nhi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
        for (size_t p = 0; p < near.size(); p += 3) {
            for (int a = 0; a < 3; ++a) {
                nlo[a] = MIN(nlo[a], near[p + a]);
                nhi[a] = MAX(nhi[a], near[p + a]);
            }
        }
        if (n.radius_cap) {
            for (int a = 0; a < 3; ++a) {
                lo[a] = MAX(lo[a], nlo[a] - n.near_radius - 1.0);
                hi[a] = MIN(hi[a], nhi[a] + n.near_radius + 1.0);
            }
        }
    }
    if (!(hi[0] > lo[0] && hi[1] > lo[1] && hi[2] > lo[2])) {
        n.want_key = 0;
        return;
    }

    double h = nci_resolution_spacing(n.resolution);
    int dim[3];
    for (int pass = 0; pass < 2; ++pass) {
        size_t total = 1;
        for (int a = 0; a < 3; ++a) {
            dim[a] = MAX(8, (int)ALIGN_TO((int)std::ceil((hi[a] - lo[a]) / h), 8));
            total *= (size_t)dim[a];
        }
        if (total <= max_samples) break;
        h *= std::cbrt((double)total / (double)max_samples) * 1.01;
    }
    nci_surface::Grid grid = {};
    for (int a = 0; a < 3; ++a) {
        grid.dim[a] = dim[a];
        grid.spacing[a] = (float)((hi[a] - lo[a]) / dim[a]);
        grid.origin[a] = (float)lo[a] + 0.5 * grid.spacing[a];
    }

    uint64_t rho_key = md_hash64(&grid, sizeof(grid), n.generation);
    if (rho_key != n.rho_key) {
        nci_surface_job_destroy(n.job);
        n.job = nullptr;
        n.rho_key = rho_key;
        n.rho.reset();
        n.rho_failed = false;
        const md_tick_t t0 = md_tick_now();
        const size_t count = (size_t)dim[0] * dim[1] * dim[2];
        std::vector<float> density;
        const bool sampled = input.density_provider.evaluate &&
            input.density_provider.evaluate(input.density_provider.user_data, grid, &density);
        if (sampled && density.size() == count) {
            n.rho = std::make_shared<std::vector<float>>(std::move(density));
            n.grid = grid;
        } else {
            if (sampled) MD_LOG_ERROR("NCI surface: density provider returned %zu samples; expected %zu", density.size(), count);
            else MD_LOG_ERROR("NCI surface: the density could not be evaluated");
            n.rho_failed = true;
        }
        n.rho_seconds = md_tick_to_seconds(md_tick_now() - t0);
        if (n.rho)
            MD_LOG_DEBUG("NCI surface: density sampled on %d x %d x %d (%.3f Bohr) in %.2f s",
                         dim[0], dim[1], dim[2], h, n.rho_seconds);
    }
    if (!n.rho) {
        n.want_key = 0;
        return;
    }

    nci_surface::Params params = nci_surface::reference_params(near, n.competitors);
    params.s_iso = n.s_iso;
    params.scatter_s = n.scatter_s;
    params.rho_max = n.rho_max;
    params.near_radius = n.near_radius;
    // No tail (and no recompute) unless the plot range reaches above the scatter's s
    params.display_s = n.display_s > n.scatter_s ? MIN(n.display_s, nci_surface::scatter_s_max) : 0.0;
    params.surface = !context;
    uint64_t key = context ? rho_key : md_hash64(&params.s_iso, sizeof(params.s_iso), rho_key);
    key = md_hash64(&params.scatter_s, sizeof(params.scatter_s), key);
    key = md_hash64(&params.display_s, sizeof(params.display_s), key);
    key = md_hash64(&params.rho_max, sizeof(params.rho_max), key);
    key = md_hash64(&n.scope, sizeof(n.scope), key);
    key = md_hash64(&params.connected, sizeof(params.connected), key);
    key = md_hash64(&params.radius_cap, sizeof(params.radius_cap), key);
    key = md_hash64(&params.surface, sizeof(params.surface), key);
    if (!near.empty()) {
        key = md_hash64(&params.near_radius, sizeof(params.near_radius), key);
        key = md_hash64(near.data(), near.size() * sizeof(double), key);
        if (!n.competitors.empty()) key = md_hash64(n.competitors.data(), n.competitors.size() * sizeof(double), key);
    }
    std::array<std::vector<double>, nci_surface::family_count> reference_near;
    for (size_t family = 0; family < reference_near.size(); ++family) {
        const nci_surface::Family type = (nci_surface::Family)family;
        if (!nci_surface::has_nci(type)) continue;
        auto& anchors = reference_near[family];
        for (const NciSurface::CriticalPointMembership& cp : n.critical_points_by_family[family])
            anchors.insert(anchors.end(), cp.position.begin(), cp.position.end());
    }
    key |= 1;
    n.want_key = key;
    n.anchors = near;
    if (key == n.have_key || (n.job && n.job->key == key)) return;
    nci_surface_job_destroy(n.job);
    n.reference_failed = false;
    std::vector<double> competitors = near.empty() ? std::vector<double>{} : n.competitors;
    n.job = nci_surface_job_start(n.rho, n.grid, params, std::move(near), std::move(competitors), reference_near, key);
}

}  // namespace topology
