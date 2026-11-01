#pragma once

#include "nci_surface_core.h"

#include <core/md_os.h>
#include <task_system.h>

#include <atomic>
#include <memory>
#include <vector>

namespace topology {

struct NciSurfaceJob {
    std::shared_ptr<const std::vector<float>> rho;
    nci_surface::Grid grid = {};
    nci_surface::Params params = {};
    std::vector<double> near;
    nci_surface::Mesh mesh;
    nci_surface::Scatter scatter;

    uint64_t key = 0;
    md_tick_t start = 0;
    bool ok = false;
    std::atomic<int32_t> cancel = 0;
    std::atomic<float> progress = 0.0f;
    std::atomic<bool> done = false;
    task_system::ID task = task_system::INVALID_ID;
};

NciSurfaceJob* nci_surface_job_start(std::shared_ptr<const std::vector<float>> rho,
                                     const nci_surface::Grid& grid,
                                     const nci_surface::Params& params,
                                     std::vector<double> near, uint64_t key);
void nci_surface_job_destroy(NciSurfaceJob* job);

}  // namespace topology
