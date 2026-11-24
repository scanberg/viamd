#include "utest.h"

#include "../heatmaps_core.h"

#include <core/md_pair_frames.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <task_system.h>

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>

/* Tiles of a heatmap against the whole matrix: every tile of every level, the cells as the entries they hold make them,
   whole and a band of cell rows at a time. */

using namespace heatmaps;

static inline uint64_t min64(uint64_t a, uint64_t b) { return a < b ? a : b; }

// A matrix of truths over 10 frames, where (i, j) true in some of them: near the diagonal, and a row and a column
static uint32_t truth_count(uint32_t i, uint32_t j) {
    const uint32_t d = i > j ? i - j : j - i;
    if (d < 4) return 10 - d;
    if (i == 77 || j == 301) return 1 + (i + j) % 9;
    uint64_t h = (uint64_t)i * 2654435761u ^ (uint64_t)j * 40503u;
    h ^= h >> 13; h *= 0x5bd1e995; h ^= h >> 15;
    return h % 97 == 0 ? (uint32_t)(h % 10) + 1 : 0;
}

UTEST(heatmaps_core, levels) {
    EXPECT_EQ(0u, tile_top(100, 128));
    EXPECT_EQ(1u, tile_top(129, 3));
    EXPECT_EQ(9u, tile_top(65536, 1000));
    // A cell a pixel or a little more
    EXPECT_EQ(0u, tile_level(0.3, 5000, 5000));
    EXPECT_EQ(0u, tile_level(1.0, 5000, 5000));
    EXPECT_EQ(1u, tile_level(1.01, 5000, 5000));
    EXPECT_EQ(1u, tile_level(2.0, 5000, 5000));
    EXPECT_EQ(3u, tile_level(7.5, 5000, 5000));
    EXPECT_EQ(tile_top(5000, 5000), tile_level(1e6, 5000, 5000));
    // Past the matrix, empty; at its end, the cells it has
    const TileSpan t = tile_span(300, 1000, 1, 3, 1);
    EXPECT_EQ(256u, t.i0); EXPECT_EQ(300u, t.i1); EXPECT_EQ(768u, t.j0); EXPECT_EQ(1000u, t.j1);
    EXPECT_EQ(22u, t.nr); EXPECT_EQ(116u, t.nc);
    EXPECT_EQ(0u, tile_span(300, 1000, 1, 3, 2).nr);
}

// Of frames [f0, f1): how many (i, j) is true in, true in frames [0, truth_count)
static uint32_t truth_in(uint32_t i, uint32_t j, uint32_t f0, uint32_t f1) {
    const uint32_t c = truth_count(i, j);
    return c > f0 ? (c < f1 ? c : f1) - f0 : 0;
}

UTEST(heatmaps_core, sparse_against_dense) {
    const uint32_t rows = 450, cols = 700, nf = 10;
    md_allocator_i* heap = md_get_heap_allocator();
    md_allocator_i* arena = md_arena_allocator_create(heap, MD_MEGABYTES(1));
    md_pair_frames_t* pf = md_pair_frames_create(rows, cols, nf, heap);
    ASSERT_TRUE(pf != NULL);
    // Frames added from the last, so that slots are not in the order of entries
    for (uint32_t f = nf; f-- > 0;) {
        std::vector<uint32_t> on;
        for (uint32_t i = 0; i < rows; ++i) {
            for (uint32_t j = 0; j < cols; ++j) if (truth_count(i, j) > f) on.push_back(i * cols + j);
        }
        ASSERT_TRUE(md_pair_frames_add(pf, f, on.data(), on.size()));
    }
    uint32_t* entries = NULL;
    uint32_t* slots = NULL;
    const size_t num = md_pair_frames_index(&entries, &slots, pf, heap);

    // Over every frame, some of them, and one
    const uint32_t ranges[3][2] = { {0, nf}, {2, 7}, {4, 5} };
    std::vector<float> cell(TILE * TILE), part(TILE * TILE);
    std::vector<uint32_t> known(TILE * TILE), known_part(TILE * TILE);
    for (const auto& fr : ranges) {
        const Sparse s = { pf, entries, slots, num, md_bitset_range(fr[0], fr[1]), fr[1] - fr[0], rows, cols };
        const float n = (float)(fr[1] - fr[0]);
        for (uint32_t i = 0; i < rows; i += 7) {
            for (uint32_t j = 0; j < cols; ++j) EXPECT_EQ((float)truth_in(i, j, fr[0], fr[1]) / n, sparse_value(s, i, j));
        }
        EXPECT_EQ((uint64_t)truth_in(5, 6, 0, nf), md_bitset_count(sparse_frames(s, 5, 6, arena)));
        for (uint32_t level = 0; level <= tile_top(rows, cols); ++level) {
            const uint64_t span = (uint64_t)TILE << level;
            for (uint32_t ty = 0; ty * span < rows; ++ty) {
                for (uint32_t tx = 0; tx * span < cols; ++tx) {
                    tile_cells(cell.data(), known.data(), s, level, tx, ty);
                    // In bands of cell rows, as threads make it: the same
                    for (uint32_t r = 0; r < TILE; r += 24) tile_cells(part.data(), known_part.data(), s, level, tx, ty, r, r + 24);
                    EXPECT_TRUE(cell == part);
                    EXPECT_TRUE(known == known_part);
                    const TileSpan t = tile_span(rows, cols, level, tx, ty);
                    for (uint32_t r = 0; r < TILE; ++r) {
                        for (uint32_t c = 0; c < TILE; ++c) {
                            // The highest of the entries the cell holds
                            uint32_t top = 0;
                            bool in = false;
                            for (uint64_t i = t.i0 + ((uint64_t)r << level); i < min64(t.i1, t.i0 + ((uint64_t)(r + 1) << level)); ++i) {
                                for (uint64_t j = t.j0 + ((uint64_t)c << level); j < min64(t.j1, t.j0 + ((uint64_t)(c + 1) << level)); ++j) {
                                    in = true;
                                    const uint32_t v = truth_in((uint32_t)i, (uint32_t)j, fr[0], fr[1]);
                                    top = top > v ? top : v;
                                }
                            }
                            EXPECT_EQ(in ? 1u : 0u, known[r * TILE + c]);
                            EXPECT_EQ((float)top / n, cell[r * TILE + c]);
                        }
                    }
                }
            }
        }
    }
    const Sparse none = { pf, entries, slots, num, md_bitset_empty_set(), 0, rows, cols };
    EXPECT_TRUE(isnan(sparse_value(none, 3, 3)));
    md_array_free(entries, heap);
    md_array_free(slots, heap);
    md_pair_frames_free(pf);
    md_arena_allocator_destroy(arena);
}

// task_system owns a global pool with no re-entrant initialize: brought up once, as test_task_system.cpp does
static void ensure_task_system() {
    static bool initialized = false;
    if (!initialized) {
        if (!task_system::executor()) {
            task_system::initialize(4);
            atexit([]() { task_system::shutdown(); });
        }
        initialized = true;
    }
}

// A view of a matrix of truths, and its tiles made over the threads, as the cells made alone color them; its index
// taken over while nothing has come, and made again when an entry has
UTEST(heatmaps_core, views_and_tiles) {
    ensure_task_system();
    const uint32_t rows = 300, cols = 520, nf = 10;
    md_allocator_i* heap = md_get_heap_allocator();
    md_pair_frames_t* pf = md_pair_frames_create(rows, cols, nf, heap);
    ASSERT_TRUE(pf != NULL);
    for (uint32_t f = 0; f < nf - 1; ++f) {     // The last frame added later
        std::vector<uint32_t> on;
        for (uint32_t i = 0; i < rows; ++i) {
            for (uint32_t j = 0; j < cols; ++j) if (truth_count(i, j) > f) on.push_back(i * cols + j);
        }
        ASSERT_TRUE(md_pair_frames_add(pf, f, on.data(), on.size()));
    }
    ViewSource src;
    src.pf = pf;
    src.beg = 2;
    src.end = nf;
    View v;
    std::unique_ptr<Index> made;
    make_view(v, &made, nullptr, src);
    ASSERT_TRUE(made != nullptr);
    EXPECT_TRUE(v.index == made.get());
    EXPECT_TRUE(v.pairs && v.fraction && !v.truth);
    EXPECT_EQ(7u, v.num_frames);        // 2 to 8: the last is not added
    EXPECT_EQ(rows, v.rows);
    EXPECT_EQ(cols, v.cols);
    EXPECT_EQ(md_pair_frames_count(pf), md_array_size(v.index->entries));

    // Tiles over the threads, one job's worth of them at once, at two levels
    uint32_t lut[256];
    for (uint32_t i = 0; i < 256; ++i) lut[i] = 0xFF000000u | i;
    const uint32_t unknown = 0xFF5A5A5Au;
    for (uint32_t level = 0; level <= 2; level += 2) {
        const uint64_t span = (uint64_t)TILE << level;
        const uint32_t ntx = (uint32_t)((cols + span - 1) / span), nty = (uint32_t)((rows + span - 1) / span);
        std::vector<std::vector<uint32_t>> rgba(ntx * nty, std::vector<uint32_t>(TILE * TILE));
        task_system::parallel_for(ntx * nty, 1, [&](uint32_t beg, uint32_t end) {
            for (uint32_t k = beg; k < end; ++k) make_tile(rgba[k].data(), v, level, k % ntx, k / ntx, lut, unknown);
        });
        std::vector<float> cell(TILE * TILE);
        std::vector<uint32_t> known(TILE * TILE);
        for (uint32_t k = 0; k < ntx * nty; ++k) {
            tile_cells(cell.data(), known.data(), v.sparse(), level, k % ntx, k / ntx);
            const TileSpan t = tile_span(rows, cols, level, k % ntx, k / ntx);
            for (uint32_t r = 0; r < TILE; ++r) {
                for (uint32_t c = 0; c < TILE; ++c) {
                    const size_t q = (size_t)r * TILE + c;
                    const uint32_t want = r >= t.nr || c >= t.nc ? 0 : lut[(int)(cell[q] * 255.0f + 0.5f)];
                    EXPECT_EQ(want, rgba[k][q]);
                }
            }
        }
    }
    // Hover: the exact fraction
    EXPECT_EQ((float)truth_in(10, 12, 2, nf - 1) / 7.0f, v.value(10, 12));

    // Nothing come since: the index taken over
    View same;
    std::unique_ptr<Index> none;
    make_view(same, &none, made.get(), src);
    EXPECT_TRUE(none == nullptr);
    EXPECT_TRUE(same.index == made.get());
    // The last frame, with an entry never true before: an index of its own, with it
    std::vector<uint32_t> on;
    for (uint32_t i = 0; i < rows; ++i) {
        for (uint32_t j = 0; j < cols; ++j) if (truth_count(i, j) > nf - 1) on.push_back(i * cols + j);
    }
    uint32_t fresh = 400;
    while (truth_count(5, fresh) != 0) ++fresh;
    on.push_back(5 * cols + fresh);
    ASSERT_TRUE(md_pair_frames_add(pf, nf - 1, on.data(), on.size()));
    View later;
    std::unique_ptr<Index> again;
    make_view(later, &again, made.get(), src);
    ASSERT_TRUE(again != nullptr);
    EXPECT_EQ(md_array_size(made->entries) + 1, md_array_size(again->entries));
    EXPECT_EQ(8u, later.num_frames);
    EXPECT_EQ(1.0f / 8.0f, later.value(5, fresh));

    // At a frame: truths
    ViewSource at = src;
    at.frame = true;
    at.beg = 4;
    at.end = 5;
    View fv;
    std::unique_ptr<Index> fi;
    make_view(fv, &fi, again.get(), at);
    EXPECT_TRUE(fv.truth);
    EXPECT_EQ(1u, fv.num_frames);
    EXPECT_EQ(truth_count(0, 0) > 4 ? 1.0f : 0.0f, fv.value(0, 0));

    md_pair_frames_free(pf);
}

// A static matrix of numbers: its values in the display unit, and its scale
UTEST(heatmaps_core, dense_view) {
    const uint32_t rows = 3, cols = 4;
    float data[rows * cols];
    for (uint32_t k = 0; k < rows * cols; ++k) data[k] = k == 5 ? NAN : (float)k;
    ViewSource src;
    src.attr.format.type = MD_ATTRIBUTE_TYPE_F32;
    src.attr.format.components = 1;
    src.attr.format.rank = 2;
    src.attr.format.shape[0] = rows;
    src.attr.format.shape[1] = cols;
    src.attr.data = data;
    src.factor = 10.0;
    View v;
    make_view(v, nullptr, nullptr, src);
    EXPECT_FALSE(v.pairs);
    EXPECT_EQ(rows, v.rows);
    EXPECT_EQ(cols, v.cols);
    EXPECT_EQ(0.0f, v.lo);
    EXPECT_EQ(110.0f, v.hi);
    EXPECT_EQ(70.0f, v.value(1, 3));
    EXPECT_TRUE(isnan(v.value(1, 1)));
}

UTEST(heatmaps_core, dense_means) {
    const uint32_t rows = 300, cols = 260;
    std::vector<float> values((size_t)rows * cols);
    for (uint32_t i = 0; i < rows; ++i) {
        for (uint32_t j = 0; j < cols; ++j) values[(size_t)i * cols + j] = (i * 7 + j * 3) % 11 == 0 ? NAN : (float)(i * 0.5 - j * 0.25);
    }
    std::vector<float> cell(TILE * TILE);
    std::vector<uint32_t> known(TILE * TILE);
    for (uint32_t level = 0; level <= tile_top(rows, cols); ++level) {
        const uint64_t span = (uint64_t)TILE << level;
        for (uint32_t ty = 0; ty * span < rows; ++ty) {
            for (uint32_t tx = 0; tx * span < cols; ++tx) {
                tile_cells(cell.data(), known.data(), values.data(), rows, cols, level, tx, ty);
                std::vector<float> part(TILE * TILE);
                std::vector<uint32_t> known_part(TILE * TILE);
                for (uint32_t r = 0; r < TILE; r += 40) tile_cells(part.data(), known_part.data(), values.data(), rows, cols, level, tx, ty, r, r + 40);
                EXPECT_TRUE(cell == part);
                EXPECT_TRUE(known == known_part);
                const TileSpan t = tile_span(rows, cols, level, tx, ty);
                for (uint32_t r = 0; r < TILE; ++r) {
                    for (uint32_t c = 0; c < TILE; ++c) {
                        double sum = 0;
                        uint32_t n = 0;
                        for (uint64_t i = t.i0 + ((uint64_t)r << level); i < min64(t.i1, t.i0 + ((uint64_t)(r + 1) << level)); ++i) {
                            for (uint64_t j = t.j0 + ((uint64_t)c << level); j < min64(t.j1, t.j0 + ((uint64_t)(c + 1) << level)); ++j) {
                                const float v = values[i * cols + j];
                                if (isnan(v)) continue;
                                sum += v;
                                n += 1;
                            }
                        }
                        EXPECT_EQ(n, known[r * TILE + c]);
                        if (n) EXPECT_NEAR((float)(sum / n), cell[r * TILE + c], 1e-4f);
                    }
                }
            }
        }
    }
}
