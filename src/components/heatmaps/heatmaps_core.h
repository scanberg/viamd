#pragma once

// Tiles of a heatmap: a matrix drawn a part at a time, at the resolution the screen shows it at, from a matrix too large
// to go through whole for what is in view.
//
// A tile is TILE × TILE cells, a cell 2^level entries a side: tile (level, tx, ty) is the rows from ty TILE 2^level
// and the columns from tx TILE 2^level. A matrix of truths over frames is read as md_pair_frames keeps it, through an
// index of its slots by entry (md_pair_frames_index), and counted over the frames shown: a tile costs what is true in
// it. Its cells are the highest fraction of the entries they hold, so one contact is not lost from far away. Any other
// matrix comes whole, and its cells are the mean of the entries known. A tile's cells can be made a band of cell rows
// at a time, on as many threads.
//
// What a heatmap shows is a View, made on a job's thread (make_view) and then read, unchanged, by the main thread and by
// the jobs making its tiles (make_tile).
//
// Free of ImGui, OpenGL and ApplicationState: tested in tests/.

#include <core/md_bitset.h>
#include <core/md_unit.h>
#include <md_attributes.h>

#include <stddef.h>
#include <stdint.h>
#include <memory>
#include <vector>

struct md_pair_frames_t;
struct md_allocator_i;

namespace script { struct Eval; }

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

// A matrix of truths as md_pair_frames keeps it, over some of its frames: the entries ever true, i * cols + j,
// ascending, with their slots (md_pair_frames_index), and the frames counted, num_frames of them, all added
struct Sparse {
    const md_pair_frames_t* pf = nullptr;
    const uint32_t* entries = nullptr;
    const uint32_t* slots = nullptr;
    size_t      num = 0;
    md_bitset_t frames = {};
    size_t      num_frames = 0;
    uint32_t    rows = 0, cols = 0;
};

// The fraction of the frames entry (i, j) is true in: 0 when in none, NaN when no frames are counted. The frames of
// all it is true in, words from alloc, empty when none.
float       sparse_value(const Sparse& s, uint32_t i, uint32_t j);
md_bitset_t sparse_frames(const Sparse& s, uint32_t i, uint32_t j, md_allocator_i* alloc);

// The cell rows [r0, r1) of tile (level, tx, ty), into those rows of cell[TILE * TILE], and how many entries known
// each holds into known[TILE * TILE]; past the matrix 0 and 0. Of truths, every entry is known and a cell is the
// highest fraction; of values[rows * cols], NaN where not known, it is their mean. Other rows are left as they are.
void tile_cells(float* cell, uint32_t* known, const Sparse& s, uint32_t level, uint32_t tx, uint32_t ty, uint32_t r0 = 0, uint32_t r1 = TILE);
void tile_cells(float* cell, uint32_t* known, const float* values, uint32_t rows, uint32_t cols, uint32_t level, uint32_t tx, uint32_t ty,
                uint32_t r0 = 0, uint32_t r1 = TILE);

// Where the entries of a matrix of truths are kept, by entry (md_pair_frames_index), as of count of them
struct Index {
    const md_pair_frames_t* pf = nullptr;
    uint32_t* entries = nullptr;    // md_arrays from the heap
    uint32_t* slots = nullptr;
    size_t    count = 0;
    Index() = default;
    Index(const Index&) = delete;
    Index& operator=(const Index&) = delete;
    ~Index();
};

// What a heatmap shows
struct View {
    uint64_t gen = 0;               // Of the tiles made of it
    uint32_t rows = 0, cols = 0;
    bool     pairs = false;         // A matrix of truths per frame
    bool     truth = false;         // Truths, of the frame shown
    bool     fraction = false;      // Fractions of frames, of truths
    size_t   num_frames = 0;        // Of a fraction or a mean: over how many frames
    md_unit_t unit = {};
    float    lo = 0.0f, hi = 1.0f;  // Its scale, in the display unit
    char     status[192] = "";
    // A matrix of truths: the frames counted, all added, words in arena, and the index
    const md_pair_frames_t* pf = nullptr;
    md_bitset_t  frames = {};
    const Index* index = nullptr;
    // Any other: the matrix, row by row, in the display unit, NaN where nothing is known
    std::vector<float> values;
    md_allocator_i* arena = nullptr;

    View();
    ~View();
    View(const View&) = delete;
    View& operator=(const View&) = delete;

    bool   shown() const { return pairs ? num_frames > 0 : !values.empty(); }
    Sparse sparse() const;
    float  value(uint32_t i, uint32_t j) const;     // Entry (i, j) as it is shown
};

// What a view is made from, taken where it does not change while it is made: a matrix of truths, or a copy of the
// description of an attribute (its frames held from eval, when it is per frame)
struct ViewSource {
    const md_pair_frames_t* pf = nullptr;
    const script::Eval* eval = nullptr;
    md_attribute_t attr = {};
    char     name[64] = "";
    bool     frame = false;         // Over the frame shown, beg
    uint32_t beg = 0, end = 0;      // The frames it is over
    double   factor = 1.0;          // To the display unit
};

// The view of src into v. A matrix of truths: the frames added among those it is over, taken once so that every tile
// counts the same, and an index, cur when no entry has come since it was made, else made again into *made.
void make_view(View& v, std::unique_ptr<Index>* made, const Index* cur, const ViewSource& src);

// Tile (level, tx, ty) of v into rgba[TILE * TILE], RGBA8, its first row at the top: a cell's value on the scale
// [lo, hi] colored by lut, unknown where nothing is known, 0 past the matrix. Its cell rows in bands over the threads
// (task_system). False when the job it is made in was cancelled.
bool make_tile(uint32_t* rgba, const View& v, uint32_t level, uint32_t tx, uint32_t ty, const uint32_t lut[256], uint32_t unknown);

}  // namespace heatmaps
