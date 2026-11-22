#include "heatmaps_core.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace heatmaps {

TileSpan tile_span(uint32_t rows, uint32_t cols, uint32_t level, uint32_t tx, uint32_t ty) {
    const uint64_t span = (uint64_t)TILE << level;
    TileSpan t;
    t.i0 = (uint32_t)std::min<uint64_t>(ty * span, rows);
    t.i1 = (uint32_t)std::min<uint64_t>((ty + 1) * span, rows);
    t.j0 = (uint32_t)std::min<uint64_t>(tx * span, cols);
    t.j1 = (uint32_t)std::min<uint64_t>((tx + 1) * span, cols);
    t.nr = (uint32_t)(((uint64_t)(t.i1 - t.i0) + (1ull << level) - 1) >> level);
    t.nc = (uint32_t)(((uint64_t)(t.j1 - t.j0) + (1ull << level) - 1) >> level);
    return t;
}

uint32_t tile_top(uint32_t rows, uint32_t cols) {
    uint32_t top = 0;
    while (((uint64_t)TILE << top) < std::max(rows, cols)) ++top;
    return top;
}

uint32_t tile_level(double per_px, uint32_t rows, uint32_t cols) {
    if (!(per_px > 1.0)) return 0;
    return std::min((uint32_t)std::ceil(std::log2(per_px)), tile_top(rows, cols));
}

float sparse_value(const Sparse& s, uint32_t i, uint32_t j) {
    if (!s.num_frames) return NAN;
    const uint32_t e = i * s.cols + j;
    const uint32_t* k = std::lower_bound(s.entries, s.entries + s.num, e);
    return k != s.entries + s.num && *k == e ? (float)s.counts[k - s.entries] / (float)s.num_frames : 0.0f;
}

void tile_cells(float* cell, uint32_t* known, const Sparse& s, uint32_t level, uint32_t tx, uint32_t ty) {
    const TileSpan t = tile_span(s.rows, s.cols, level, tx, ty);
    std::fill(cell, cell + TILE * TILE, 0.0f);
    std::fill(known, known + TILE * TILE, 0u);
    for (uint32_t r = 0; r < t.nr; ++r) {
        for (uint32_t c = 0; c < t.nc; ++c) known[r * TILE + c] = 1;
    }
    if (!s.num_frames) return;
    // A row's entries in the tile's columns, a binary search on from the last row's: the most frames of them, then as
    // a fraction
    const uint32_t* end = s.entries + s.num;
    const uint32_t* k = s.entries;
    for (uint32_t i = t.i0; i < t.i1; ++i) {
        const uint32_t row = i * s.cols;    // Below 2^32: rows × cols is
        k = std::lower_bound(k, end, row + t.j0);
        float* line = cell + (size_t)((i - t.i0) >> level) * TILE;
        for (; k != end && *k < row + t.j1; ++k) {
            float& v = line[(*k - row - t.j0) >> level];
            v = std::max(v, (float)s.counts[k - s.entries]);
        }
    }
    for (uint32_t c = 0; c < TILE * TILE; ++c) cell[c] /= (float)s.num_frames;
}

void tile_cells(float* cell, uint32_t* known, const float* values, uint32_t rows, uint32_t cols, uint32_t level, uint32_t tx, uint32_t ty) {
    const TileSpan t = tile_span(rows, cols, level, tx, ty);
    std::vector<double> sum(TILE * TILE, 0.0);
    std::fill(known, known + TILE * TILE, 0u);
    for (uint32_t i = t.i0; i < t.i1; ++i) {
        const float* src = values + (size_t)i * cols;
        const size_t line = (size_t)((i - t.i0) >> level) * TILE;
        for (uint32_t j = t.j0; j < t.j1; ++j) {
            if (std::isnan(src[j])) continue;
            const size_t c = line + ((j - t.j0) >> level);
            sum[c] += src[j];
            known[c] += 1;
        }
    }
    for (uint32_t c = 0; c < TILE * TILE; ++c) cell[c] = known[c] ? (float)(sum[c] / known[c]) : 0.0f;
}

}  // namespace heatmaps
