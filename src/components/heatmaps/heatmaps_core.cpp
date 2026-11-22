#include "heatmaps_core.h"

#include <core/md_pair_frames.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <script_engine.h>
#include <task_system.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
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

// The slot of entry (i, j), or false when it is never true
static bool slot_of(const Sparse& s, uint32_t i, uint32_t j, uint32_t* slot) {
    if (i >= s.rows || j >= s.cols) return false;
    const uint32_t e = i * s.cols + j;
    const uint32_t* k = std::lower_bound(s.entries, s.entries + s.num, e);
    if (k == s.entries + s.num || *k != e) return false;
    *slot = s.slots[k - s.entries];
    return true;
}

float sparse_value(const Sparse& s, uint32_t i, uint32_t j) {
    if (!s.num_frames) return NAN;
    uint32_t slot, count = 0;
    if (!slot_of(s, i, j, &slot)) return 0.0f;
    md_pair_frames_slot_counts(&count, s.pf, &slot, 1, s.frames);
    return (float)count / (float)s.num_frames;
}

md_bitset_t sparse_frames(const Sparse& s, uint32_t i, uint32_t j, md_allocator_i* alloc) {
    uint32_t slot;
    return slot_of(s, i, j, &slot) ? md_pair_frames_slot_frames(s.pf, slot, alloc) : md_bitset_empty_set();
}

// The rows of the matrix cell rows [r0, r1) of a tile hold
static void band(uint32_t* b0, uint32_t* b1, const TileSpan& t, uint32_t level, uint32_t r0, uint32_t r1) {
    *b0 = (uint32_t)std::min<uint64_t>(t.i0 + ((uint64_t)r0 << level), t.i1);
    *b1 = (uint32_t)std::min<uint64_t>(t.i0 + ((uint64_t)r1 << level), t.i1);
}

static void clear_rows(float* cell, uint32_t* known, uint32_t r0, uint32_t r1) {
    std::fill(cell + (size_t)r0 * TILE, cell + (size_t)r1 * TILE, 0.0f);
    std::fill(known + (size_t)r0 * TILE, known + (size_t)r1 * TILE, 0u);
}

void tile_cells(float* cell, uint32_t* known, const Sparse& s, uint32_t level, uint32_t tx, uint32_t ty, uint32_t r0, uint32_t r1) {
    r1 = std::min(r1, TILE);
    if (r0 >= r1) return;
    const TileSpan t = tile_span(s.rows, s.cols, level, tx, ty);
    clear_rows(cell, known, r0, r1);
    for (uint32_t r = r0; r < std::min(r1, t.nr); ++r) {
        for (uint32_t c = 0; c < t.nc; ++c) known[r * TILE + c] = 1;
    }
    if (!s.num_frames) return;
    // The slots in the band, by a binary search on from the last row's, and the cell of each; then their counts, the
    // most of a cell's as its fraction
    uint32_t b0, b1;
    band(&b0, &b1, t, level, r0, r1);
    std::vector<uint32_t> slots, cells;
    const uint32_t* end = s.entries + s.num;
    const uint32_t* k = s.entries;
    for (uint32_t i = b0; i < b1; ++i) {
        const uint32_t row = i * s.cols;    // Below 2^32: rows × cols is
        k = std::lower_bound(k, end, row + t.j0);
        const uint32_t line = ((i - t.i0) >> level) * TILE;
        for (; k != end && *k < row + t.j1; ++k) {
            slots.push_back(s.slots[k - s.entries]);
            cells.push_back(line + ((*k - row - t.j0) >> level));
        }
    }
    std::vector<uint32_t> counts(slots.size());
    md_pair_frames_slot_counts(counts.data(), s.pf, slots.data(), slots.size(), s.frames);
    for (size_t q = 0; q < slots.size(); ++q) cell[cells[q]] = std::max(cell[cells[q]], (float)counts[q]);
    for (size_t c = (size_t)r0 * TILE; c < (size_t)r1 * TILE; ++c) cell[c] /= (float)s.num_frames;
}

void tile_cells(float* cell, uint32_t* known, const float* values, uint32_t rows, uint32_t cols, uint32_t level, uint32_t tx, uint32_t ty,
                uint32_t r0, uint32_t r1) {
    r1 = std::min(r1, TILE);
    if (r0 >= r1) return;
    const TileSpan t = tile_span(rows, cols, level, tx, ty);
    clear_rows(cell, known, r0, r1);
    uint32_t b0, b1;
    band(&b0, &b1, t, level, r0, r1);
    std::vector<double> sum((size_t)(r1 - r0) * TILE, 0.0);
    for (uint32_t i = b0; i < b1; ++i) {
        const float* src = values + (size_t)i * cols;
        const size_t line = (size_t)((i - t.i0) >> level) * TILE;
        for (uint32_t j = t.j0; j < t.j1; ++j) {
            if (std::isnan(src[j])) continue;
            const size_t c = line + ((j - t.j0) >> level);
            sum[c - (size_t)r0 * TILE] += src[j];
            known[c] += 1;
        }
    }
    for (size_t c = (size_t)r0 * TILE; c < (size_t)r1 * TILE; ++c) cell[c] = known[c] ? (float)(sum[c - (size_t)r0 * TILE] / known[c]) : 0.0f;
}

Index::~Index() {
    md_array_free(entries, md_get_heap_allocator());
    md_array_free(slots, md_get_heap_allocator());
}

View::View() { arena = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(16)); }
View::~View() { md_arena_allocator_destroy(arena); }

Sparse View::sparse() const {
    return { pf, index ? index->entries : nullptr, index ? index->slots : nullptr, index ? md_array_size(index->entries) : 0,
             frames, num_frames, rows, cols };
}

float View::value(uint32_t i, uint32_t j) const {
    if (pairs) return sparse_value(sparse(), i, j);
    return values.empty() || i >= rows || j >= cols ? NAN : values[(size_t)i * cols + j];
}

static float attr_value(const md_attribute_t* a, size_t i) {
    switch (a->format.type) {
    case MD_ATTRIBUTE_TYPE_F32: return ((const float*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_F64: return (float)((const double*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_U8:  return (float)((const uint8_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_I8:  return (float)((const int8_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_U16: return (float)((const uint16_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_I16: return (float)((const int16_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_U32: return (float)((const uint32_t*)a->data)[i];
    case MD_ATTRIBUTE_TYPE_I32: return (float)((const int32_t*)a->data)[i];
    default: return NAN;
    }
}

void make_view(View& v, std::unique_ptr<Index>* made, const Index* cur, const ViewSource& src) {
    if (src.pf) {
        const md_pair_frames_t* pf = src.pf;
        v.pairs = true;
        v.pf = pf;
        v.rows = md_pair_frames_rows(pf);
        v.cols = md_pair_frames_cols(pf);
        v.lo = 0.0f;
        v.hi = 1.0f;
        v.truth = src.frame;
        v.fraction = !v.truth;
        const md_bitset_t added = md_pair_frames_added(pf, v.arena);
        v.frames = md_bitset_and(added, md_bitset_range(src.beg, src.end), v.arena);
        v.num_frames = (size_t)md_bitset_count(v.frames);
        // After the frames: it holds every entry true in them
        const size_t count = md_pair_frames_count(pf);
        if (cur && cur->pf == pf && cur->count == count) {
            v.index = cur;
        } else {
            made->reset(new Index);
            (*made)->pf = pf;
            (*made)->count = md_pair_frames_index(&(*made)->entries, &(*made)->slots, pf, md_get_heap_allocator());
            v.index = made->get();
        }
        const size_t num = md_array_size(v.index->entries);
        if (v.truth) {
            if (!v.num_frames) {
                snprintf(v.status, sizeof(v.status), "Frame %u is not evaluated yet", src.beg);
            } else {
                // How many are true at it
                std::vector<uint32_t> c(num);
                md_pair_frames_slot_counts(c.data(), pf, v.index->slots, num, v.frames);
                size_t on = 0;
                for (uint32_t x : c) on += x != 0;
                snprintf(v.status, sizeof(v.status), "Frame %u: %zu entries true", src.beg, on);
            }
        } else if (!v.num_frames) {
            snprintf(v.status, sizeof(v.status), "None of its %u frames is evaluated yet", src.end - src.beg);
        } else {
            snprintf(v.status, sizeof(v.status), "The fraction of %zu frames%s each entry is true in, of %zu entries ever true",
                     v.num_frames, v.num_frames < src.end - src.beg ? " evaluated so far" : "", num);
        }
        return;
    }

    const md_attribute_t* attr = &src.attr;
    const bool per_frame = attr->format.rank == 3;
    v.rows = attr->format.shape[per_frame ? 1 : 0];
    v.cols = attr->format.shape[per_frame ? 2 : 1];
    v.unit = attr->unit;
    const size_t n = (size_t)v.rows * v.cols;
    v.values.assign(n, NAN);
    if (!per_frame) {
        for (size_t i = 0; i < n; ++i) v.values[i] = attr_value(attr, i);
    } else {
        const uint32_t nf = attr->format.shape[0];
        const uint32_t beg = src.beg, end = std::min(src.end, nf);
        // The frames it holds: a frame not computed yet is NaN, or 0 for integers, and is left out
        const md_bitset_t held = script::eval_output_frames(src.eval, md_str_from_cstr(src.name), v.arena);
        const uint64_t num = beg < end ? md_bitset_count_range(held, beg, end) : 0;
        if (num == 0) {
            v.values.clear();
            snprintf(v.status, sizeof(v.status), src.frame ? "Frame %u is not evaluated yet" : "None of its frames is evaluated yet", beg);
            return;
        }
        std::vector<uint32_t> count(n, 0);
        std::vector<double> sum(n, 0.0);
        md_bitset_iter_t it = md_bitset_iter_range(held, beg, end);
        while (md_bitset_iter_next(&it)) {
            if (task_system::cancelled()) return;
            const size_t base = (size_t)it.idx * n;
            for (size_t i = 0; i < n; ++i) {
                const float x = attr_value(attr, base + i);
                if (std::isnan(x)) continue;
                sum[i] += x;
                count[i] += 1;
            }
        }
        for (size_t i = 0; i < n; ++i) v.values[i] = count[i] ? (float)(sum[i] / count[i]) : NAN;
        v.num_frames = (size_t)num;
        if (src.frame) snprintf(v.status, sizeof(v.status), "Frame %u", beg);
        else snprintf(v.status, sizeof(v.status), "The mean over %llu frames", (unsigned long long)num);
    }
    // Its scale, in the display unit
    float lo = INFINITY, hi = -INFINITY;
    for (float& x : v.values) {
        if (std::isnan(x)) continue;
        x = (float)(x * src.factor);
        lo = std::min(lo, x);
        hi = std::max(hi, x);
    }
    if (!(hi >= lo)) { lo = 0.0f; hi = 1.0f; }
    if (hi == lo) hi = lo + 1.0f;
    v.lo = lo;
    v.hi = hi;
}

bool make_tile(uint32_t* rgba, const View& v, uint32_t level, uint32_t tx, uint32_t ty, const uint32_t lut[256], uint32_t unknown) {
    std::vector<float> cell(TILE * TILE);
    std::vector<uint32_t> known(TILE * TILE);
    const Sparse s = v.sparse();
    task_system::parallel_for(TILE, 16, [&](uint32_t r0, uint32_t r1) {
        if (v.pairs) tile_cells(cell.data(), known.data(), s, level, tx, ty, r0, r1);
        else tile_cells(cell.data(), known.data(), v.values.data(), v.rows, v.cols, level, tx, ty, r0, r1);
    });
    if (task_system::cancelled()) return false;
    const TileSpan span = tile_span(v.rows, v.cols, level, tx, ty);
    const float range = v.hi - v.lo;
    for (uint32_t r = 0; r < TILE; ++r) {
        for (uint32_t c = 0; c < TILE; ++c) {
            const size_t idx = (size_t)r * TILE + c;
            if (r >= span.nr || c >= span.nc) { rgba[idx] = 0; continue; }
            if (!known[idx]) { rgba[idx] = unknown; continue; }
            const float u = std::clamp((cell[idx] - v.lo) / range, 0.0f, 1.0f);
            rgba[idx] = lut[(int)(u * 255.0f + 0.5f)];
        }
    }
    return true;
}

}  // namespace heatmaps
