#include "nci_surface_job.h"

#include <core/md_common.h>

#include <cmath>

namespace topology {

NciSurfaceJob* nci_surface_job_start(std::shared_ptr<const std::vector<float>> rho,
                                     const nci_surface::Grid& grid,
                                     const nci_surface::Params& params,
                                     std::vector<double> near,
                                     const std::array<std::vector<double>, nci_surface::family_count>& reference_near,
                                     uint64_t key) {
    NciSurfaceJob* job = new NciSurfaceJob();
    job->rho = std::move(rho);
    job->grid = grid;
    job->near = std::move(near);
    job->params = params;
    job->params.near_xyz = job->near;
    job->key = key;
    job->start = md_tick_now();
    const auto category_anchors = reference_near;
    job->task = task_system::create_pool_task(STR_LIT("##NCI surface"), [job, category_anchors]() {
        job->ok = nci_surface::extract(&job->mesh, job->rho->data(), job->grid, job->params,
                                       &job->cancel, &job->progress, &job->scatter);
        if (job->ok) {
            for (size_t family = 0; family < category_anchors.size(); ++family) {
                if (category_anchors[family].empty()) continue;
                const nci_surface::Family type = (nci_surface::Family)family;
                const nci_surface::Family settings_family =
                    nci_surface::reference_settings_family(nci_surface::reference_category(type));
                const nci_surface::Defaults defaults =
                    nci_surface::family_defaults(nci_surface::family_bit(settings_family));
                const bool matches_surface_scope = category_anchors[family] == job->near &&
                    std::abs(job->params.s_iso - defaults.s_iso) < 1e-6 &&
                    std::abs(job->params.rho_max - defaults.rho_max) < 1e-6 &&
                    std::abs(job->params.near_radius - defaults.radius) < 1e-6 &&
                    job->params.connected && job->params.radius_cap && job->params.same_sign;
                if (matches_surface_scope) {
                    job->reference_scatters[family] = job->scatter;
                    continue;
                }
                nci_surface::Params params;
                params.s_iso = defaults.s_iso;
                params.rho_max = defaults.rho_max;
                params.near_radius = defaults.radius;
                params.connected = true;
                params.same_sign = true;
                params.radius_cap = true;
                params.surface = false;
                params.near_xyz = category_anchors[family];
                nci_surface::Mesh mesh;
                if (!nci_surface::extract(&mesh, job->rho->data(), job->grid, params,
                                          &job->cancel, &job->progress, &job->reference_scatters[family])) {
                    job->reference_ok = false;
                    break;
                }
            }
        }
        job->done.store(true, std::memory_order_release);
    });
    task_system::enqueue_task(job->task);
    return job;
}

void nci_surface_job_destroy(NciSurfaceJob* job) {
    if (!job) return;
    job->cancel = 1;
    task_system::task_wait_for(job->task);
    delete job;
}

}  // namespace topology
