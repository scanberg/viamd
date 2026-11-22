#include "utest.h"

#include "../heatmaps_core.h"

#include <math.h>
#include <stdint.h>
#include <vector>

/* Tiles of a heatmap against the whole matrix: every tile of every level, the cells as the entries they hold make them. */

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

UTEST(heatmaps_core, sparse_against_dense) {
    const uint32_t rows = 450, cols = 700;
    std::vector<uint32_t> entries, counts;
    for (uint32_t i = 0; i < rows; ++i) {
        for (uint32_t j = 0; j < cols; ++j) {
            const uint32_t c = truth_count(i, j);
            if (c) { entries.push_back(i * cols + j); counts.push_back(c); }
        }
    }
    const Sparse s = { entries.data(), counts.data(), entries.size(), 10, rows, cols };
    for (uint32_t i = 0; i < rows; i += 7) {
        for (uint32_t j = 0; j < cols; ++j) EXPECT_EQ((float)truth_count(i, j) / 10.0f, sparse_value(s, i, j));
    }
    const Sparse none = { entries.data(), counts.data(), entries.size(), 0, rows, cols };
    EXPECT_TRUE(isnan(sparse_value(none, 3, 3)));

    std::vector<float> cell(TILE * TILE);
    std::vector<uint32_t> known(TILE * TILE);
    for (uint32_t level = 0; level <= tile_top(rows, cols); ++level) {
        const uint64_t span = (uint64_t)TILE << level;
        for (uint32_t ty = 0; ty * span < rows; ++ty) {
            for (uint32_t tx = 0; tx * span < cols; ++tx) {
                tile_cells(cell.data(), known.data(), s, level, tx, ty);
                const TileSpan t = tile_span(rows, cols, level, tx, ty);
                for (uint32_t r = 0; r < TILE; ++r) {
                    for (uint32_t c = 0; c < TILE; ++c) {
                        // The highest of the entries the cell holds
                        uint32_t top = 0;
                        bool in = false;
                        for (uint64_t i = t.i0 + ((uint64_t)r << level); i < min64(t.i1, t.i0 + ((uint64_t)(r + 1) << level)); ++i) {
                            for (uint64_t j = t.j0 + ((uint64_t)c << level); j < min64(t.j1, t.j0 + ((uint64_t)(c + 1) << level)); ++j) {
                                in = true;
                                top = top > truth_count((uint32_t)i, (uint32_t)j) ? top : truth_count((uint32_t)i, (uint32_t)j);
                            }
                        }
                        EXPECT_EQ(in ? 1u : 0u, known[r * TILE + c]);
                        EXPECT_EQ((float)top / 10.0f, cell[r * TILE + c]);
                    }
                }
            }
        }
    }
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
