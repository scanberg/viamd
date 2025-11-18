#pragma once

#ifdef MD_TREXIO

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

struct md_trexio_t;
struct md_allocator_i;

namespace trexio_orbital_grid {

// Grid statistics after evaluation
struct GridStats {
    float min_value;
    float max_value;
    float mean_value;
    size_t num_points;
};

// Bounding box for grid evaluation
struct BoundingBox {
    float min_x, max_x;
    float min_y, max_y;
    float min_z, max_z;
};

// Grid resolution
struct GridResolution {
    int nx, ny, nz;
};

// Progress callback function type
// Returns false to request cancellation
using ProgressCallback = bool (*)(float progress, void* user_data);

// Log callback function type
using LogCallback = void (*)(const char* message, void* user_data);

// Compute orbital grid for a specific molecular orbital
// - trexio_data: TREXIO data structure with MO and basis information
// - mo_index: Index of the molecular orbital to evaluate (0-based)
// - bbox: Bounding box for grid evaluation
// - resolution: Grid resolution (number of points in each dimension)
// - out_grid: Output grid values (must be pre-allocated: nx * ny * nz floats)
// - out_stats: Output statistics (min/max/mean)
// - progress_callback: Optional progress callback (can be nullptr)
// - log_callback: Optional logging callback (can be nullptr)
// - user_data: User data passed to callbacks
// Returns true on success, false on error or cancellation
bool compute_orbital_grid(
    md_trexio_t* trexio_data,
    int mo_index,
    const BoundingBox& bbox,
    const GridResolution& resolution,
    float* out_grid,
    GridStats* out_stats,
    ProgressCallback progress_callback = nullptr,
    LogCallback log_callback = nullptr,
    void* user_data = nullptr
);

// Helper function to save grid to binary file
// File format: 
//   - 4 bytes: magic number (0x4754524F = "GTRO")
//   - 4 bytes: version (1)
//   - 12 bytes: resolution (nx, ny, nz as int32)
//   - 24 bytes: bounding box (min_x, max_x, min_y, max_y, min_z, max_z as float)
//   - grid_data: nx * ny * nz floats
bool save_grid_to_file(
    const char* filename,
    const float* grid_data,
    const GridResolution& resolution,
    const BoundingBox& bbox
);

} // namespace trexio_orbital_grid

#endif // MD_TREXIO
