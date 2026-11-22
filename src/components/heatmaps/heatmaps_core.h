#pragma once

// Tiles of a heatmap: a matrix drawn a part at a time, at the resolution the screen shows it at, from a matrix too large
// to hold whole as cells.
//
// A tile is TILE × TILE cells, a cell 2^level entries a side: tile (level, tx, ty) is the rows from ty TILE 2^level
// and the columns from tx TILE 2^level. A matrix of truths over frames comes as what is true (md_pair_frames_counts):
// its cells are the highest fraction of the entries they hold, so one contact is not lost from far away. Any other comes
// whole, and its cells are the mean of the entries known.
//
// Free of ImGui, OpenGL and ApplicationState: tested in tests/.

#include <stddef.h>
#include <stdint.h>

namespace heatmaps {

constexpr uint32_t TILE = 128;          // A tile's cells a side

// Of tile (level, tx, ty), in a matrix of rows × cols: its rows [i0, i1) and columns [j0, j1), and the cells they
// make, nr × nc, the last of a row or a column fewer entries where the matrix ends. Empty past it.
struct TileSpan {
    uint32_t i0 = 0, i1 = 0, j0 = 0, j1 = 0;
    uint32_t nr = 0, nc = 0;
};
TileSpan tile_span(uint32_t rows, uint32_t cols, uint32_t level, uint32_t tx, uint32_t ty);

// The level of one tile over the whole matrix; the level a view of per_px entries a pixel is drawn at: a cell a pixel
// or a little more, up to two, so what is drawn is never more than the screen shows
uint32_t tile_top(uint32_t rows, uint32_t cols);
uint32_t tile_level(double per_px, uint32_t rows, uint32_t cols);

// A matrix of truths over num_frames frames, as what is true: entries i * cols + j, ascending, entry k true in counts[k]
struct Sparse {
    const uint32_t* entries = nullptr;
    const uint32_t* counts = nullptr;
    size_t   num = 0;
    size_t   num_frames = 0;
    uint32_t rows = 0, cols = 0;
};

// The fraction of the frames entry (i, j) is true in: 0 when in none, NaN when there are no frames
float sparse_value(const Sparse& s, uint32_t i, uint32_t j);

// The cells of tile (level, tx, ty), row by row, into cell[TILE * TILE], and how many entries known each holds into
// known[TILE * TILE]; past the matrix 0 and 0. Of truths, every entry is known and a cell is the highest fraction; of
// values[rows * cols], NaN where not known, it is their mean.
void tile_cells(float* cell, uint32_t* known, const Sparse& s, uint32_t level, uint32_t tx, uint32_t ty);
void tile_cells(float* cell, uint32_t* known, const float* values, uint32_t rows, uint32_t cols, uint32_t level, uint32_t tx, uint32_t ty);

}  // namespace heatmaps
