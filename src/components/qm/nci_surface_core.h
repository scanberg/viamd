#pragma once

// NCI surfaces (Johnson et al., JACS 2010): the isosurface of the reduced density gradient
// s = |grad rho| / (2 (3 pi^2)^(1/3) rho^(4/3)) where the density is low (rho <= rho_max), coloured by sign(lambda2) rho,
// lambda2 the second eigenvalue of the density's Hessian: negative attractive, near zero van der Waals, positive repulsive.
//
// From a density sampled on a regular, axis-aligned grid: the gradient and the Hessian by central finite differences,
// the surface by Surface Nets on s. Samples denser than rho_max, and the grid's outer layer, count as s = infinity (out),
// as NCIPLOT's density cutoff does, so the low gradient regions of covalent bonds and of the atoms are not shown.

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace nci_surface {

// Sample (i, j, k) is at origin + (i, j, k) * spacing, Bohr; rho[i + dim[0] * (j + dim[1] * k)]
struct Grid {
    int    dim[3]     = {0, 0, 0};
    double origin[3]  = {0, 0, 0};
    double spacing[3] = {0, 0, 0};
};

struct Params {
    double s_iso   = 0.5;
    double rho_max = 0.05;   // au
    // When given, only the surface within radius (Bohr) of any of these points (Bohr), e.g. the BCPs of interest
    std::span<const double> near_xyz = {};   // 3 per point
    double near_radius = 3.0;
};

struct Mesh {
    std::vector<float>    position;   // xyz per vertex, Bohr
    std::vector<float>    normal;     // xyz per vertex, unit, towards higher s (out of the low gradient region)
    std::vector<float>    value;      // per vertex: sign(lambda2) rho, au
    std::vector<uint32_t> index;      // three per triangle, counter-clockwise seen from where the normals point
};

constexpr size_t scatter_max_points = 50000;
constexpr float scatter_s_max = 2.0f;

struct Scatter {
    std::vector<float> signed_rho;
    std::vector<float> s;
    size_t eligible_count = 0;
};

// Returns false when cancelled (cancel set non-zero) or the input is invalid; progress (optional) goes from 0 to 1.
// Optional scatter: a deterministic reservoir of grid samples with s <= 2, respecting rho_max and near_xyz.
// It samples the volume, not the surface; changing s_iso does not change the scatter.
bool extract(Mesh* out, const float* rho, const Grid& grid, const Params& params,
             const std::atomic<int32_t>* cancel = nullptr, std::atomic<float>* progress = nullptr,
             Scatter* scatter = nullptr);

// The reduced density gradient and sign(lambda2) rho of sample (i, j, k), not on the grid's outer layer
void sample(double* s, double* signed_rho, const float* rho, const Grid& grid, int i, int j, int k);

}  // namespace nci_surface
