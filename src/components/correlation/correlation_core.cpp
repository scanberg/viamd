#include "correlation_core.h"

#include <algorithm>
#include <math.h>
#include <vector>

namespace correlation {

void corr_rank(double* out, const double* v, size_t n) {
    std::vector<uint32_t> idx(n);
    for (size_t i = 0; i < n; ++i) idx[i] = (uint32_t)i;
    std::sort(idx.begin(), idx.end(), [v](uint32_t a, uint32_t b) { return v[a] < v[b]; });
    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && v[idx[j]] == v[idx[i]]) ++j;
        const double r = 0.5 * (double)(i + j + 1);   // The mean of the ranks i + 1 .. j
        for (size_t k = i; k < j; ++k) out[idx[k]] = r;
        i = j;
    }
}

bool corr_pearson(double* r, const double* a, const double* b, size_t n) {
    if (n < 3) return false;
    double ma = 0, mb = 0;
    for (size_t i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
    ma /= (double)n;
    mb /= (double)n;
    double saa = 0, sbb = 0, sab = 0;
    for (size_t i = 0; i < n; ++i) {
        const double da = a[i] - ma, db = b[i] - mb;
        saa += da * da;
        sbb += db * db;
        sab += da * db;
    }
    if (!(saa > 0.0) || !(sbb > 0.0)) return false;
    const double v = sab / sqrt(saa * sbb);
    *r = v < -1.0 ? -1.0 : v > 1.0 ? 1.0 : v;
    return true;
}

bool corr_spearman(double* rho, const double* a, const double* b, size_t n) {
    if (n < 3) return false;
    std::vector<double> ra(n), rb(n);
    corr_rank(ra.data(), a, n);
    corr_rank(rb.data(), b, n);
    return corr_pearson(rho, ra.data(), rb.data(), n);
}

bool corr_circular(double* r, const double* a, const double* b, size_t n) {
    if (n < 3) return false;
    double sa = 0, ca = 0, sb = 0, cb = 0;
    for (size_t i = 0; i < n; ++i) {
        sa += sin(a[i]); ca += cos(a[i]);
        sb += sin(b[i]); cb += cos(b[i]);
    }
    const double ma = atan2(sa, ca), mb = atan2(sb, cb);
    double num = 0, daa = 0, dbb = 0;
    for (size_t i = 0; i < n; ++i) {
        const double da = sin(a[i] - ma);
        const double db = sin(b[i] - mb);
        num += da * db;
        daa += da * da;
        dbb += db * db;
    }
    if (!(daa > 0.0) || !(dbb > 0.0)) return false;
    const double v = num / sqrt(daa * dbb);
    *r = v < -1.0 ? -1.0 : v > 1.0 ? 1.0 : v;
    return true;
}

bool corr_circular_linear(double* r, const double* a, const double* x, size_t n) {
    if (n < 3) return false;
    const double inv_n = 1.0 / (double)n;
    double mx = 0, mc = 0, ms = 0;
    for (size_t i = 0; i < n; ++i) {
        mx += x[i]; mc += cos(a[i]); ms += sin(a[i]);
    }
    mx *= inv_n; mc *= inv_n; ms *= inv_n;
    double sxx = 0, scc = 0, sss = 0, sxc = 0, sxs = 0, scs = 0;
    for (size_t i = 0; i < n; ++i) {
        const double dx = x[i] - mx, dc = cos(a[i]) - mc, ds = sin(a[i]) - ms;
        sxx += dx * dx; scc += dc * dc; sss += ds * ds;
        sxc += dx * dc; sxs += dx * ds; scs += dc * ds;
    }
    if (!(sxx > 0.0) || !(scc > 0.0) || !(sss > 0.0)) return false;
    const double rxc = sxc / sqrt(sxx * scc);
    const double rxs = sxs / sqrt(sxx * sss);
    const double rcs = scs / sqrt(scc * sss);
    const double den = 1.0 - rcs * rcs;
    // The cosine and the sine move together: the angles are all but on a line through the circle
    if (!(den > 1e-12)) return false;
    const double r2 = (rxc * rxc + rxs * rxs - 2.0 * rxc * rxs * rcs) / den;
    *r = sqrt(r2 < 0.0 ? 0.0 : r2 > 1.0 ? 1.0 : r2);
    return true;
}

bool corr_of(double* r, CorrKind kind, const double* x, const double* y, size_t n) {
    switch (kind) {
    case CorrKind_CircularX: return corr_circular_linear(r, x, y, n);
    case CorrKind_CircularY: return corr_circular_linear(r, y, x, n);
    case CorrKind_Circular:  return corr_circular(r, x, y, n);
    default:                 return corr_pearson(r, x, y, n);
    }
}

size_t corr_lag_match(uint32_t* ia, uint32_t* ib, const uint32_t* frames, size_t n, int64_t lag) {
    size_t count = 0;
    size_t j = 0;
    for (size_t i = 0; i < n; ++i) {
        const int64_t want = (int64_t)frames[i] + lag;
        if (want < 0) continue;
        while (j < n && (int64_t)frames[j] < want) ++j;
        if (j == n) break;
        if ((int64_t)frames[j] == want) {
            ia[count] = (uint32_t)i;
            ib[count] = (uint32_t)j;
            ++count;
        }
    }
    return count;
}

}  // namespace correlation
