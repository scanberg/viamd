#pragma once

#include <core/types.h>
#include <core/array.h>

namespace spatialhash {

struct Cell {
	int offset;
	int count;
};

struct Frame {
	vec3 min_box;
	vec3 max_box;
	vec3 cell_ext;
	ivec3 cell_count;

	DynamicArray<vec3> positions;
	DynamicArray<Cell> cells;
};

inline int compute_cell_idx(const Frame& frame, ivec3 cell_coord) {
	ASSERT(cell_coord.x < frame.cell_count.x);
	ASSERT(cell_coord.y < frame.cell_count.y);
	ASSERT(cell_coord.z < frame.cell_count.z);

	return cell_coord.z * frame.cell_count.x * frame.cell_count.y + cell_coord.y * frame.cell_count.x + cell_coord.x;
}

inline ivec3 compute_cell_coord(const Frame& frame, vec3 coord) {
	return ivec3((coord - frame.min_box) / frame.cell_ext);
}

inline int compute_cell_idx(const Frame& frame, vec3 coord) {
	return compute_cell_idx(frame, compute_cell_coord(frame, coord));
}

inline Cell get_cell(const Frame& frame, ivec3 cell_coord) {
	int idx = compute_cell_idx(frame, cell_coord);
	ASSERT(idx < frame.cells.count);
	return frame.cells[idx];
}

/*
inline Array<const int> get_cell_indices(const Frame& frame, ivec3 cell_coord) {
	Cell cell = get_cell(frame, cell_coord);
	return { frame.indices.beg() + cell.offset, cell.count };
}

inline DynamicArray<int> query_indices(const Frame& frame, vec3 coord, float radius) {
	DynamicArray<int> res;
	ivec3 min_cc = compute_cell_coord(frame, coord - radius);
	ivec3 max_cc = compute_cell_coord(frame, coord + radius);
	ivec3 cc;
	for (cc.z = min_cc.z; cc.z <= max_cc.z; cc.z++) {
		for (cc.y = min_cc.y; cc.y <= max_cc.y; cc.y++) {
			Cell beg_cell = get_cell(frame, ivec3(cc.x, cc.y, min_cc.x));
			Cell end_cell = get_cell(frame, ivec3(cc.x, cc.y, max_cc.x));
			Array<const int> arr(frame.indices.beg() + beg_cell.index, frame.indices.beg() + end_cell.index + end_cell.count);
			res.append(get_cell_indices(frame, cc));
		}
	}

	return res;
}
*/

Frame compute_frame(Array<vec3> positions, vec3 cell_ext);

Frame compute_frame(Array<vec3> positions, vec3 cell_ext, vec3 min_box, vec3 max_box);

}