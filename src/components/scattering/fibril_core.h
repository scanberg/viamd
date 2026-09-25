#pragma once

// Continuous fibril density from a slice based coarse grained bead model.
//
// A fibril is described by a sequence of slices (e.g. 7 beads: one center bead and 6 surrounding beads). From the
// beads of every slice we extract
//   - a centerline point (electron weighted centroid of the slice),
//   - the local fibril axis (tangent of a Catmull-Rom spline through the centroids),
//   - the in-plane orientation of the cross-section (2D Procrustes fit of the slice beads to the template anchors,
//     in a rotation minimizing frame transported along the fibril).
// A cross-section template (a set of 2D Gaussians in the slice frame) is then swept continuously along the
// centerline, with sample spacing <= the smallest template width, such that the resulting density is continuous
// along the fibril. The electrons of each slice are conserved: half of them are spread over each adjacent segment,
// the outer halves of the end slices over end caps of half a slice spacing.
//
// Lengths are in Ångström.

#include <stddef.h>
#include <stdint.h>
#include <vector>

struct FibrilGauss {
    float u, v;         // offset in the slice frame
    float weight;       // fraction of the slice electrons (normalized to sum 1)
    float sigma;        // Gaussian width
};

struct FibrilTemplate {
    // Bead layout in the slice frame, used to determine the orientation of each slice (index = bead index within slice).
    // Centered on the electron weighted centroid.
    std::vector<float> anchor_u;
    std::vector<float> anchor_v;
    // Cross-section density
    std::vector<FibrilGauss> gauss;
};

struct FibrilInput {
    size_t num_slices = 0;
    int    stride = 0;                      // beads per slice
    const float* x = nullptr;               // num_slices * stride bead positions (may be wrapped by the periodic box)
    const float* y = nullptr;
    const float* z = nullptr;
    const float* e = nullptr;               // electrons per bead (num_slices * stride)
    const uint32_t* chain_offset = nullptr; // num_chains + 1, chains are slice ranges [offset[c], offset[c+1])
    size_t num_chains = 0;
    double box[3] = {0, 0, 0};              // periodic lengths, 0 = not periodic
    int    ref_bead = 1;                    // bead used for the initial orientation estimate
};

struct FibrilStats {
    size_t num_points = 0;
    double electrons_in = 0.0;
    double electrons_out = 0.0;
    double mean_spacing = 0.0;              // mean distance between consecutive slice centroids
    double mean_radius = 0.0;               // RMS bead distance from the slice centroid, in the slice plane
    double mean_axial_offset = 0.0;         // RMS bead offset along the fibril axis (slice thickness indicator)
};

struct FibrilOutput {
    std::vector<float> x, y, z, w, s;
    FibrilStats stats;
};

// Mean bead layout of the slices in their own frame (orientation from ref_bead only). Returns false on failure.
bool fibril_mean_layout(std::vector<float>& out_u, std::vector<float>& out_v, std::vector<float>& out_frac, const FibrilInput& in);

// Template presets
// Beads: one Gaussian per anchor, weights = electron fractions
FibrilTemplate fibril_template_beads(const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& frac, float sigma);
// Uniform disk of radius R represented by Gaussians on a hexagonal grid with the given spacing (sigma = 0.5 * spacing)
FibrilTemplate fibril_template_disk(const std::vector<float>& anchor_u, const std::vector<float>& anchor_v, float radius, float spacing);
// Text file:
//   anchor <u> <v>                  one line per bead in slice order (optional, the mean layout is used otherwise)
//   gauss <u> <v> <weight> <sigma>  cross-section components
// '#' starts a comment. Lengths in Å.
bool fibril_template_load(FibrilTemplate* out, const char* path, char* err, size_t err_cap);

// Makes the handedness of the template anchors consistent with a reference layout (e.g. the mean MD layout): if the
// mirrored anchors (v -> -v) match the reference better after an optimal in-plane rotation, the whole template is
// mirrored. Guards against an opposite fibril direction between the reference slice and the MD. Returns true if mirrored.
// out_rms receives the RMS anchor mismatch (Å) of the chosen handedness.
bool fibril_template_match_handedness(FibrilTemplate& tmpl, const std::vector<float>& ref_u, const std::vector<float>& ref_v, double* out_rms);

// Sweep the template along all chains. ds <= 0 -> smallest template sigma.
bool fibril_sweep(FibrilOutput* out, const FibrilInput& in, const FibrilTemplate& tmpl, double ds, char* err, size_t err_cap);
