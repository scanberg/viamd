#include <plot_periodic_core.h>

#include <core/md_common.h>

#include <float.h>
#include <math.h>
#include <stddef.h>

double periodic_wrap(double x, double beg, double len) {
    if (!(len > 0.0)) return x;
    const double t = x - beg;
    double w = beg + (t - len * floor(t / len));
    // Rounding may land a value just short of the end on it: the same value, at the beginning
    if (w >= beg + len) w = beg;
    return w;
}

int periodic_polyline(double* out_x, double* out_y, int cap, const float* x, const float* y, int n, int stride, int k, double scale,
                      double base, double len, double max_step) {
    ASSERT(out_x && out_y && x && y);
    const double top = base + len;
    const double nan = NAN;
    int m = 0;
    auto emit = [&](double px, double py) {
        if (m < cap) {
            out_x[m] = px;
            out_y[m] = py;
            m += 1;
        }
    };

    bool have = false;      // a line is under way
    double px = 0, py = 0;  // its last point, in the window
    double praw = 0;        // and the value it was
    for (int i = 0; i < n; ++i) {
        const double raw = (double)y[(size_t)i * stride + k] * scale;
        const double xi = x[i];
        if (!(raw == raw)) {
            if (have) emit(nan, nan);
            have = false;
            continue;
        }
        const double wrapped = periodic_wrap(raw, base, len);
        if (!have) {
            emit(xi, wrapped);
            have = true;
            px = xi; py = wrapped; praw = raw;
            continue;
        }

        // The step, the short way around
        const double d = periodic_wrap(raw - praw, -0.5 * len, len);
        if (fabs(d) > max_step) {
            emit(nan, nan);
            emit(xi, wrapped);
        } else {
            const double to = py + d;
            if (to >= top || to < base) {
                // Out at one edge and in at the other, where it crosses
                const double edge  = to >= top ? top : base;
                const double other = to >= top ? base : top;
                const double t = d != 0.0 ? (edge - py) / d : 0.0;
                const double xc = px + t * (xi - px);
                emit(xc, edge);
                emit(nan, nan);
                emit(xc, other);
            }
            emit(xi, wrapped);
        }
        px = xi; py = wrapped; praw = raw;
    }
    return m;
}

int periodic_band_runs(double* out_x, double* out_lo, double* out_hi, int* run_beg, int* run_count,
                       const double* x, const double* lo, const double* hi, int n, double base, double len, double max_step) {
    int num_runs = 0;
    int m = 0;
    bool have = false;
    double pc = 0;          // the centre of the last sample, unwrapped along the run
    double praw = 0;        // and as it was
    int prev = -1;          // the last sample, for a run which takes over from another
    auto begin_run = [&]() {
        run_beg[num_runs] = m;
        run_count[num_runs] = 0;
        num_runs += 1;
    };
    auto emit = [&](int i, double c) {
        const double h = 0.5 * (hi[i] - lo[i]);
        out_x[m] = x[i];
        out_lo[m] = c - h;
        out_hi[m] = c + h;
        m += 1;
        run_count[num_runs - 1] += 1;
    };

    for (int i = 0; i < n; ++i) {
        const double c_raw = 0.5 * (lo[i] + hi[i]);
        if (!(c_raw == c_raw) || !(hi[i] == hi[i])) {
            have = false;
            continue;
        }
        if (!have) {
            begin_run();
            pc = periodic_wrap(c_raw, base, len);
            emit(i, pc);
            have = true;
            praw = c_raw;
            prev = i;
            continue;
        }
        const double d = periodic_wrap(c_raw - praw, -0.5 * len, len);
        if (fabs(d) > max_step) {
            begin_run();
            pc = periodic_wrap(c_raw, base, len);
            emit(i, pc);
        } else {
            double c = pc + d;
            if (c < base - 0.5 * len || c >= base + 1.5 * len) {
                // Strayed from the window: taken over a period nearer it, from the sample before, so that the band
                // goes on unbroken
                const double shift = c < base ? len : -len;
                begin_run();
                emit(prev, pc + shift);
                c += shift;
            }
            emit(i, c);
            pc = c;
        }
        praw = c_raw;
        prev = i;
    }
    return num_runs;
}

int periodic_least_populated(const float* density, int num_bins) {
    if (!density || num_bins <= 0) return -1;
    const int r = MAX(1, num_bins / 20);
    int best = -1;
    double best_sum = DBL_MAX;
    for (int b = 0; b < num_bins; ++b) {
        double sum = 0.0;
        for (int i = -r; i <= r; ++i) {
            int s = (b + i) % num_bins;
            if (s < 0) s += num_bins;
            sum += density[s];
        }
        if (sum < best_sum) {
            best_sum = sum;
            best = b;
        }
    }
    return best;
}
