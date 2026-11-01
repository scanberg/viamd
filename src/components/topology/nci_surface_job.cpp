#include "nci_surface_job.h"

#include <core/md_common.h>

namespace topology {

NciSurfaceJob* nci_surface_job_start(std::shared_ptr<const std::vector<float>> rho,
                                     const nci_surface::Grid& grid,
                                     const nci_surface::Params& params,
                                     std::vector<double> near, uint64_t key) {
    NciSurfaceJob* job = new NciSurfaceJob();
    job->rho = std::move(rho);
    job->grid = grid;
    job->near = std::move(near);
    job->params = params;
    job->params.near_xyz = job->near;
    job->key = key;
    job->start = md_tick_now();
    job->task = task_system::create_pool_task(STR_LIT("##NCI surface"), [job]() {
        job->ok = nci_surface::extract(&job->mesh, job->rho->data(), job->grid, job->params,
                                       &job->cancel, &job->progress, &job->scatter);
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
