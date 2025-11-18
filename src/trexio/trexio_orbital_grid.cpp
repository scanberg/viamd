// TREXIO Orbital Grid Evaluator
// Evaluates molecular orbitals on a 3D grid using GTO basis functions

#ifdef MD_TREXIO

#include "trexio_orbital_grid.h"

#include <md_trexio.h>
#include <md_gto.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_grid.h>

#include <cmath>
#include <cstring>
#include <cstdarg>
#include <algorithm>

namespace trexio_orbital_grid {

// Helper function to log messages
static void log_message(LogCallback log_cb, void* user_data, const char* fmt, ...) {
    if (!log_cb) return;
    
    char buffer[512];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buffer, sizeof(buffer), fmt, args);
    va_end(args);
    
    log_cb(buffer, user_data);
}

bool compute_orbital_grid(
    md_trexio_t* trexio_data,
    int mo_index,
    const BoundingBox& bbox,
    const GridResolution& resolution,
    float* out_grid,
    GridStats* out_stats,
    ProgressCallback progress_callback,
    LogCallback log_callback,
    void* user_data
) {
    if (!trexio_data || !out_grid || !out_stats) {
        MD_LOG_ERROR("Invalid parameters to compute_orbital_grid");
        return false;
    }
    
    // Validate MO index
    size_t num_mos = md_trexio_mo_num(trexio_data);
    if (mo_index < 0 || mo_index >= (int)num_mos) {
        MD_LOG_ERROR("Invalid MO index %d (must be 0-%zu)", mo_index, num_mos - 1);
        log_message(log_callback, user_data, "ERROR: Invalid MO index %d", mo_index);
        return false;
    }
    
    log_message(log_callback, user_data, "Starting orbital grid computation for MO %d", mo_index);
    log_message(log_callback, user_data, "Grid: %dx%dx%d = %d points", 
                resolution.nx, resolution.ny, resolution.nz, 
                resolution.nx * resolution.ny * resolution.nz);
    log_message(log_callback, user_data, "Bounds: [%.2f,%.2f] x [%.2f,%.2f] x [%.2f,%.2f]",
                bbox.min_x, bbox.max_x, bbox.min_y, bbox.max_y, bbox.min_z, bbox.max_z);
    
    // Get allocator
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Extract GTOs for this molecular orbital
    size_t num_gtos = md_trexio_mo_gto_count(trexio_data);
    if (num_gtos == 0) {
        MD_LOG_ERROR("No GTOs available for MO evaluation");
        log_message(log_callback, user_data, "ERROR: No GTOs available");
        return false;
    }
    
    log_message(log_callback, user_data, "Allocating space for %zu GTOs", num_gtos);
    
    md_gto_t* gtos = (md_gto_t*)md_alloc(alloc, num_gtos * sizeof(md_gto_t));
    if (!gtos) {
        MD_LOG_ERROR("Failed to allocate memory for GTOs");
        log_message(log_callback, user_data, "ERROR: Memory allocation failed");
        return false;
    }
    
    // Extract GTOs with a reasonable cutoff (1e-6 is typical for quantum chemistry)
    double cutoff_value = 1.0e-6;
    size_t extracted_gtos = md_trexio_mo_gto_extract(gtos, trexio_data, mo_index, cutoff_value);
    
    if (extracted_gtos == 0) {
        MD_LOG_ERROR("Failed to extract GTOs for MO %d", mo_index);
        log_message(log_callback, user_data, "ERROR: Failed to extract GTOs");
        md_free(alloc, gtos, num_gtos * sizeof(md_gto_t));
        return false;
    }
    
    log_message(log_callback, user_data, "Extracted %zu GTOs (cutoff: %.1e)", extracted_gtos, cutoff_value);
    
    // Report progress: GTO extraction complete (10%)
    if (progress_callback && !progress_callback(0.1f, user_data)) {
        log_message(log_callback, user_data, "Computation cancelled by user");
        md_free(alloc, gtos, num_gtos * sizeof(md_gto_t));
        return false;
    }
    
    // Setup grid structure
    md_grid_t grid;
    grid.dim[0] = resolution.nx;
    grid.dim[1] = resolution.ny;
    grid.dim[2] = resolution.nz;
    
    // Identity orientation (axis-aligned grid)
    grid.orientation[0][0] = 1.0f; grid.orientation[0][1] = 0.0f; grid.orientation[0][2] = 0.0f;
    grid.orientation[1][0] = 0.0f; grid.orientation[1][1] = 1.0f; grid.orientation[1][2] = 0.0f;
    grid.orientation[2][0] = 0.0f; grid.orientation[2][1] = 0.0f; grid.orientation[2][2] = 1.0f;
    
    // Calculate grid origin and spacing
    grid.origin[0] = bbox.min_x;
    grid.origin[1] = bbox.min_y;
    grid.origin[2] = bbox.min_z;
    
    float dx = (bbox.max_x - bbox.min_x) / (float)(resolution.nx - 1);
    float dy = (bbox.max_y - bbox.min_y) / (float)(resolution.ny - 1);
    float dz = (bbox.max_z - bbox.min_z) / (float)(resolution.nz - 1);
    
    grid.spacing[0] = dx;
    grid.spacing[1] = dy;
    grid.spacing[2] = dz;
    
    log_message(log_callback, user_data, "Grid spacing: (%.4f, %.4f, %.4f) Angstrom", dx, dy, dz);
    log_message(log_callback, user_data, "Evaluating orbital on grid...");
    
    // Report progress: Starting evaluation (20%)
    if (progress_callback && !progress_callback(0.2f, user_data)) {
        log_message(log_callback, user_data, "Computation cancelled by user");
        md_free(alloc, gtos, num_gtos * sizeof(md_gto_t));
        return false;
    }
    
    // Evaluate GTOs over the grid
    // Use PSI mode for molecular orbitals (not PSI_SQUARED)
    md_gto_grid_evaluate(out_grid, &grid, gtos, extracted_gtos, MD_GTO_EVAL_MODE_PSI);
    
    // Report progress: Evaluation complete (80%)
    if (progress_callback && !progress_callback(0.8f, user_data)) {
        log_message(log_callback, user_data, "Computation cancelled by user");
        md_free(alloc, gtos, num_gtos * sizeof(md_gto_t));
        return false;
    }
    
    log_message(log_callback, user_data, "Computing grid statistics...");
    
    // Calculate statistics
    size_t total_points = resolution.nx * resolution.ny * resolution.nz;
    float min_val = out_grid[0];
    float max_val = out_grid[0];
    double sum = 0.0;
    
    for (size_t i = 0; i < total_points; ++i) {
        float val = out_grid[i];
        if (val < min_val) min_val = val;
        if (val > max_val) max_val = val;
        sum += val;
    }
    
    out_stats->min_value = min_val;
    out_stats->max_value = max_val;
    out_stats->mean_value = (float)(sum / total_points);
    out_stats->num_points = total_points;
    
    log_message(log_callback, user_data, "Statistics: min=%.6e, max=%.6e, mean=%.6e",
                min_val, max_val, out_stats->mean_value);
    
    // Clean up
    md_free(alloc, gtos, num_gtos * sizeof(md_gto_t));
    
    // Report progress: Complete (100%)
    if (progress_callback) {
        progress_callback(1.0f, user_data);
    }
    
    log_message(log_callback, user_data, "Orbital grid computation completed successfully");
    
    return true;
}

bool save_grid_to_file(
    const char* filename,
    const float* grid_data,
    const GridResolution& resolution,
    const BoundingBox& bbox
) {
    if (!filename || !grid_data) {
        MD_LOG_ERROR("Invalid parameters to save_grid_to_file");
        return false;
    }
    
    FILE* file = fopen(filename, "wb");
    if (!file) {
        MD_LOG_ERROR("Failed to open file for writing: %s", filename);
        return false;
    }
    
    // Write header
    uint32_t magic = 0x4754524F; // "GTRO"
    uint32_t version = 1;
    int32_t dims[3] = { resolution.nx, resolution.ny, resolution.nz };
    float bounds[6] = { bbox.min_x, bbox.max_x, bbox.min_y, bbox.max_y, bbox.min_z, bbox.max_z };
    
    fwrite(&magic, sizeof(uint32_t), 1, file);
    fwrite(&version, sizeof(uint32_t), 1, file);
    fwrite(dims, sizeof(int32_t), 3, file);
    fwrite(bounds, sizeof(float), 6, file);
    
    // Write grid data
    size_t total_points = resolution.nx * resolution.ny * resolution.nz;
    fwrite(grid_data, sizeof(float), total_points, file);
    
    fclose(file);
    
    MD_LOG_INFO("Saved grid to file: %s (%d x %d x %d = %zu points)",
                filename, resolution.nx, resolution.ny, resolution.nz, total_points);
    
    return true;
}

} // namespace trexio_orbital_grid

#endif // MD_TREXIO
