#pragma once

// The geometry of periodic values laid out flat (plot_periodic.h): free of ImGui, so that it can be tested on its own.

// A step between two samples longer than this, of the period, is not drawn: so far around the circle, which way the values
// went is anybody's guess (a 180 degree turn is as near one way as the other), and the samples do not resolve it.
#define PERIODIC_MAX_STEP 0.4

// x moved by whole periods of length len into [beg, beg + len[
double periodic_wrap(double x, double beg, double len);

// The place of value nearest to ref, a whole number of periods of length len away from it
static inline double periodic_near(double value, double ref, double len) {
    return periodic_wrap(value, ref - 0.5 * len, len);
}

// A periodic series as a polyline in the window [base, base + len[, for a plot to draw: sample i at (x[i], y[i * stride + k]
// * scale), moved into the window by whole periods, and each step from one sample to the next taken the short way around
// the circle. A step which crosses an edge of the window is cut there: the line runs out at one edge and in again at the
// other, at the x where it crosses, interpolated, so that a sample is seen to be joined to its neighbour across the cut. A
// step longer than max_step, and a sample which is not a number (a frame not evaluated), leave a gap. The polyline is
// written to out_x and out_y, broken where there is no line by points which are not numbers (ImPlot leaves a gap there);
// cap is at least 4 n. Returns the number of points written.
int periodic_polyline(double* out_x, double* out_y, int cap, const float* x, const float* y, int n, int stride, int k, double scale,
                      double base, double len, double max_step);

// The band of a periodic series about its centre, [lo[i], hi[i]] at x[i], hi[i] - lo[i] at most a period: as runs of samples
// along which the centre is unwrapped (each step the short way around the circle) and stays within half a period of the
// window [base, base + len[. Drawn as it is and a period below and above, the runs cover the window wherever the band goes
// around. A run breaks where the centre steps further than max_step or a sample is not a number, and where the centre
// strays from the window, after which the next run takes over a whole period back, from the sample before. out_x, out_lo
// and out_hi have room for 2 n, run_beg and run_count for n. Returns the number of runs.
int periodic_band_runs(double* out_x, double* out_lo, double* out_hi, int* run_beg, int* run_count,
                       const double* x, const double* lo, const double* hi, int n, double base, double len, double max_step);

// The bin of num_bins around a period, density over them, where it is least populated: after smoothing around the circle over
// a tenth of the period, for the least of a noisy histogram is wherever the noise dipped. -1 for no bins.
int periodic_least_populated(const float* density, int num_bins);
