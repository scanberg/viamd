#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/math_utils.h>

namespace spatialhash {

struct Entry {
    vec3 position{};
    int index = 0;
};

struct Cell {
    int offset = 0;
    int count = 0;
};

struct Frame {
    vec3 min_box{};
    vec3 max_box{};
    vec3 cell_ext{};
    ivec3 cell_count{};

    DynamicArray<Cell> cells{};
    DynamicArray<Entry> entries{};
};

inline int compute_cell_idx(const Frame& frame, ivec3 cell_coord) {
    ASSERT(cell_coord.x < frame.cell_count.x);
    ASSERT(cell_coord.y < frame.cell_count.y);
    ASSERT(cell_coord.z < frame.cell_count.z);

    return cell_coord.z * frame.cell_count.x * frame.cell_count.y + cell_coord.y * frame.cell_count.x + cell_coord.x;
}

inline ivec3 compute_cell_coord(const Frame& frame, vec3 coord) {
    return math::clamp(ivec3((coord - frame.min_box) / frame.cell_ext), ivec3(0), frame.cell_count - 1);
}

inline int compute_cell_idx(const Frame& frame, vec3 coord) { return compute_cell_idx(frame, compute_cell_coord(frame, coord)); }

inline Cell get_cell(const Frame& frame, ivec3 cell_coord) {
    int idx = compute_cell_idx(frame, cell_coord);
    ASSERT(idx < frame.cells.count);
    return frame.cells[idx];
}

inline Array<const Entry> get_cell_entries(const Frame& frame, ivec3 cell_coord) {
    Cell cell = get_cell(frame, cell_coord);
    return {frame.entries.beg() + cell.offset, cell.count};
}

inline DynamicArray<int> query_indices(const Frame& frame, vec3 coord, float radius) {
    const float r2 = radius * radius;
    ivec3 min_cc = compute_cell_coord(frame, coord - radius);
    ivec3 max_cc = compute_cell_coord(frame, coord + radius);
    ivec3 cc;
    DynamicArray<int> res;
    for (cc.z = min_cc.z; cc.z <= max_cc.z; cc.z++) {
        for (cc.y = min_cc.y; cc.y <= max_cc.y; cc.y++) {
            for (cc.x = min_cc.x; cc.x <= max_cc.x; cc.x++) {
                for (const auto& e : get_cell_entries(frame, cc)) {
                    if (math::distance2(coord, e.position) < r2) {
                        res.push_back(e.index);
                    }
                }
            }
        }
    }
    return res;
}

template <typename Callback>
void for_each_within(const Frame& frame, vec3 coord, float radius, Callback cb) {
    const float r2 = radius * radius;
    const ivec3 min_cc = compute_cell_coord(frame, coord - radius);
    const ivec3 max_cc = compute_cell_coord(frame, coord + radius);
    ivec3 cc;
    for (cc.z = min_cc.z; cc.z <= max_cc.z; cc.z++) {
        for (cc.y = min_cc.y; cc.y <= max_cc.y; cc.y++) {
            for (cc.x = min_cc.x; cc.x <= max_cc.x; cc.x++) {
                for (const auto& e : get_cell_entries(frame, cc)) {
                    if (math::distance2(coord, e.position) < r2) {
                        cb(e.index, e.position);
                    }
                }
            }
        }
    }
}

Frame compute_frame(Array<const vec3> positions, vec3 cell_ext);
void compute_frame(Frame* frame, Array<const vec3> positions, vec3 cell_ext);

Frame compute_frame(Array<const vec3> positions, vec3 cell_ext, vec3 min_box, vec3 max_box);
void compute_frame(Frame* frame, Array<const vec3> positions, vec3 cell_ext, vec3 min_box, vec3 max_box);

}  // namespace spatialhash
