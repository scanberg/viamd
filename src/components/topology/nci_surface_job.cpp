#include "nci_surface_job.h"

#include <core/md_common.h>

#include <algorithm>
#include <cmath>

namespace topology {

NciSurfaceJob* nci_surface_job_start(std::shared_ptr<const std::vector<float>> rho,
                                     const nci_surface::Grid& grid,
                                     const nci_surface::Params& params,
                                     std::vector<double> near,
                                     std::vector<double> competitors,
                                     const std::array<std::vector<double>, nci_surface::family_count>& reference_near,
                                     uint64_t key) {
    NciSurfaceJob* job = new NciSurfaceJob();
    job->rho = std::move(rho);
    job->grid = grid;
    job->near = std::move(near);
    job->competitors = std::move(competitors);
    job->params = params;
    job->params.near_xyz = job->near;
    job->params.other_xyz = job->competitors;
    job->params.by_anchor = true;
    job->key = key;
    job->start = md_tick_now();
    const auto category_anchors = reference_near;
    job->task = task_system::create_pool_task(STR_LIT("##NCI surface"), [job, category_anchors]() {
        job->ok = nci_surface::extract(&job->mesh, job->rho->data(), job->grid, job->params,
                                       &job->cancel, &job->progress, &job->scatter);
        if (!job->ok) {
            job->done.store(true, std::memory_order_release);
            return;
        }
        // Slice the selected extraction region by category, retaining nearest-CP ownership.
        std::vector<double> anchors;
        std::vector<uint8_t> anchor_family;
        for (size_t family = 0; family < category_anchors.size(); ++family) {
            anchors.insert(anchors.end(), category_anchors[family].begin(), category_anchors[family].end());
            anchor_family.insert(anchor_family.end(), category_anchors[family].size() / 3, (uint8_t)family);
        }
        if (anchors.empty()) {
            job->done.store(true, std::memory_order_release);
            return;
        }
        nci_surface::Params reference = job->params;
        reference.near_xyz = anchors;
        reference.display_s = 0.0;
        reference.surface = false;
        const auto same = [](double a, double b) { return std::abs(a - b) < 1e-6; };
        const bool same_scope = same(job->params.scatter_s > 0.0 ? job->params.scatter_s : job->params.s_iso, reference.scatter_s) &&
            same(job->params.rho_max, reference.rho_max) && same(job->params.near_radius, reference.near_radius) &&
            job->params.connected && job->params.radius_cap;
        nci_surface::Scatter own;
        const nci_surface::Scatter* source = &own;
        std::vector<uint8_t> family_of;
        if (same_scope && anchors.size() == job->near.size()) {
            // The surface's own scatter, if its anchors are the same set; map its anchor order to the types
            family_of.assign(job->near.size() / 3, UINT8_MAX);
            for (size_t i = 0; i < family_of.size(); ++i)
                for (size_t j = 0; j < anchor_family.size(); ++j)
                    if (std::equal(&job->near[3 * i], &job->near[3 * i] + 3, &anchors[3 * j])) { family_of[i] = anchor_family[j]; break; }
            if (std::find(family_of.begin(), family_of.end(), UINT8_MAX) == family_of.end()) source = &job->scatter;
        }
        if (source == &own) {
            family_of = anchor_family;
            nci_surface::Mesh mesh;
            if (!nci_surface::extract(&mesh, job->rho->data(), job->grid, reference, &job->cancel, &job->progress, &own)) {
                job->reference_ok = false;
                job->done.store(true, std::memory_order_release);
                return;
            }
        }
        for (size_t family = 0; family < category_anchors.size(); ++family) {
            if (category_anchors[family].empty()) continue;
            std::vector<uint8_t> keep(family_of.size());
            for (size_t i = 0; i < keep.size(); ++i) keep[i] = family_of[i] == family;
            job->reference_scatters[family] = nci_surface::scatter_subset(*source, keep);
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
