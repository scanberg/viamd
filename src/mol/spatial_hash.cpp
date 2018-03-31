#include "spatial_hash.h"
#include <core/common.h>
#include <core/math_utils.h>

namespace spatialhash {

Frame spatialhash::compute_frame(Array<vec3> positions, vec3 cell_ext) {
	if (positions.count == 0) return {};
	vec3 min_box(FLT_MAX);
	vec3 max_box(-FLT_MAX);
	for (const auto& p : positions) {
		min_box = math::min(min_box, p);
		max_box = math::max(max_box, p);
	}

	return compute_frame(positions, cell_ext, min_box, max_box);
}

Frame spatialhash::compute_frame(Array<vec3> positions, vec3 cell_ext, vec3 min_box, vec3 max_box) {
	Frame frame;
	frame.min_box = min_box;
	frame.max_box = max_box;
	frame.cell_count = ivec3((max_box - min_box) / cell_ext);
	frame.cell_ext = (max_box - min_box) / (vec3)frame.cell_count;
	frame.cells.resize(frame.cell_count.x * frame.cell_count.y * frame.cell_count.z);
	frame.positions = positions;

	if (positions.count == 0) return frame;

	int* l_idx = (int*)TMP_MALLOC(positions.count * sizeof(int));
	int* g_idx = (int*)TMP_MALLOC(positions.count * sizeof(int));

	for (int i = 0; i < positions.count; i++) {
		int cell_idx = compute_cell_idx(frame, positions[i]);
		l_idx[i] = frame.cells[cell_idx].count++;
		g_idx[i] = cell_idx;
	}

	for (int i = 1; i < frame.cells.count; i++) {
		frame.cells[i].offset = frame.cells[i-1].offset + frame.cells[i-1].count;
	}

	for (int i = 0; i < frame.positions.count; i++) {
		int dst = frame.cells[g_idx[i]].offset + l_idx[i];
		//frame.positions[dst] = i;
	}

	TMP_FREE(g_idx);
	TMP_FREE(l_idx);

	return frame;
}

}  // namespace spatialhash
