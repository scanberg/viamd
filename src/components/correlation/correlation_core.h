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

// Of the samples of a sequence at frames[0..n[ (increasing): those i, j with frames[j] = frames[i] + lag, into ia and
// ib (room for n). Returns how many.
size_t corr_lag_match(uint32_t* ia, uint32_t* ib, const uint32_t* frames, size_t n, int64_t lag);

}  // namespace correlation
