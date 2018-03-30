#pragma once

#include <core/types.h>
#include <core/array.h>

namespace spatialhash {

struct Cell {
	int beg_idx;
	int count;
};

struct Frame {
	vec3 min_box;
	vec3 max_box;
	vec3 cell_ext;
	ivec3 cell_count;

	DynamicArray<int> indices;
	DynamicArray<Cell> cells;
};

inline int compute_cell_idx(const Frame& frame, ivec3 cell_coord) {
	ASSERT(cell_coord.x < frame.cell_count.x);
	ASSERT(cell_coord.y < frame.cell_count.y);
	ASSERT(cell_coord.z < frame.cell_count.z);

	return cell_coord.z * frame.cell_count.x * frame.cell_count.y + cell_coord.y * frame.cell_count.x + cell_coord.x;
}

inline int compute_cell_idx(const Frame& frame, vec3 coord) {
	ivec3 cell_coord = ivec3((coord - frame.min_box) / frame.cell_ext);
	return compute_cell_idx(frame, cell_coord);
}

inline Cell get_cell(const Frame& frame, ivec3 cell_coord) {
	int idx = compute_cell_idx(frame, cell_coord);
	ASSERT(idx < frame.cells.count);
	return frame.cells[idx];
}

inline Array<int> get_cell_indices(Frame* frame, ivec3 cell_coord) {
	auto cell = get_cell(*frame, cell_coord);
	return { frame->indices.beg() + cell.beg_idx, frame->indices.beg() + cell.beg_idx + cell.count };
}

Frame compute_frame(Array<vec3> positions, vec3 cell_ext);

Frame compute_frame(Array<vec3> positions, vec3 cell_ext, vec3 min_box, vec3 max_box);

}