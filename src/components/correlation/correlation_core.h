#pragma once

// The statistics of the Correlation window (correlation.cpp): free of ImGui and of the application, so that they can be
// tested on their own.

#include <stddef.h>
#include <stdint.h>

namespace correlation {

// What two series are: values along a line, or angles around a circle
enum CorrKind : int {
    CorrKind_Linear = 0,    // Pearson's
    CorrKind_CircularX,     // x an angle, y linear: Mardia's circular-linear (0..1, without a sign)
    CorrKind_CircularY,     // y an angle, x linear: the same
    CorrKind_Circular,      // Both angles: Jammalamadaka and SenGupta's (-1..1)
};

static inline CorrKind corr_kind(bool x_circular, bool y_circular) {
    return x_circular ? (y_circular ? CorrKind_Circular : CorrKind_CircularX) : (y_circular ? CorrKind_CircularY : CorrKind_Linear);
}

// The ranks of v[0..n[, from 1, ties given the mean of theirs
void corr_rank(double* out, const double* v, size_t n);

// Pearson's of a against b. False (r untouched) for fewer than 3 samples, or when either does not vary.
bool corr_pearson(double* r, const double* a, const double* b, size_t n);

// Spearman's: Pearson's of the ranks
bool corr_spearman(double* rho, const double* a, const double* b, size_t n);

// Jammalamadaka and SenGupta's (1988) of angles a against angles b, in radians: about their circular means
bool corr_circular(double* r, const double* a, const double* b, size_t n);

// Mardia's (1976) of angles a, in radians, against linear values x: how well x is predicted by the cosine and sine of a
bool corr_circular_linear(double* r, const double* a, const double* x, size_t n);

// The one of kind (angles in radians)
bool corr_of(double* r, CorrKind kind, const double* x, const double* y, size_t n);

// The signed parts of an angle against a linear value: Pearson's of x against the cosine and against the sine of a (in
// radians). Mardia's is the multiple correlation of x on both. False when either is undefined.
bool corr_circular_linear_parts(double* r_cos, double* r_sin, const double* a, const double* x, size_t n);

// ## Pairing

// What series on the two axes of a layer make: none, a mismatch (several on both, not as many), or pairs
enum CorrShape : int {
    CorrShape_Empty = 0,
    CorrShape_Mismatch,
    CorrShape_Valid,
};

static inline CorrShape corr_shape(int num_x, int num_y) {
    if (num_x <= 0 || num_y <= 0) return CorrShape_Empty;
    if (num_x > 1 && num_y > 1 && num_x != num_y) return CorrShape_Mismatch;
    return CorrShape_Valid;
}

// One against one, one against several (1 x n, n x 1), or as many against as many, pairwise
static inline int corr_num_pairs(int num_x, int num_y) {
    return corr_shape(num_x, num_y) == CorrShape_Valid ? (num_x > num_y ? num_x : num_y) : 0;
}

// The series of an axis pair p is of: the only one, or the p'th
static inline int corr_pair_slot(int num, int p) { return num == 1 ? 0 : p; }

// The members two series of dx and dy values per sample pair up into: the other's, against a single value; the fewer,
// of two populations (truncated set when they differ)
static inline int corr_members(int dx, int dy, bool* truncated) {
    if (truncated) *truncated = dx > 1 && dy > 1 && dx != dy;
    return dx == 1 ? dy : dy == 1 ? dx : (dx < dy ? dx : dy);
}

// The member of a series of dim values per sample member m of a pair takes: its only one, or the m'th
static inline int corr_member_slot(int dim, int m) { return dim == 1 ? 0 : m; }

// Member k of values sampled at n places, v[i * stride + k] * scale, at fractional index fi: linear between the samples
// around it (clamped to the ends). Of periodic values (period > 0), the short way around the circle between them,
// possibly outside the period they are given in.
double corr_interp(const float* v, size_t stride, int n, double scale, double fi, int k, double period);

// ## Confidence

// The lag-1 autocorrelation of sequences interleaved in v: v[i] the next value of sequence group[i] (of num_groups),
// each about its own mean, pooled. 0 when undefined.
double corr_lag1(const double* v, const uint32_t* group, size_t n, size_t num_groups);

// The number of independent samples n samples of two series of lag-1 autocorrelations ax and ay are worth (Bretherton et
// al. 1999): n (1 - ax ay) / (1 + ax ay), within [2, n]
double corr_effective_n(double n, double ax, double ay);

// The confidence interval of a correlation r of n_eff independent samples, by Fisher's z, at z standard errors (1.96: 95
// %). False when n_eff is too few (3 or less).
bool corr_fisher_interval(double* lo, double* hi, double r, double n_eff, double z);

// Of the samples of a sequence at frames[0..n[ (increasing): those i, j with frames[j] = frames[i] + lag, into ia and
// ib (room for n). Returns how many.
size_t corr_lag_match(uint32_t* ia, uint32_t* ib, const uint32_t* frames, size_t n, int64_t lag);

}  // namespace correlation
