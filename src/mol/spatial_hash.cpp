#include "spatial_hash.h"
#include <core/common.h>
#include <core/math_utils.h>

#include <atomic>

namespace spatialhash {

void compute_frame(Frame* frame, Array<const vec3> positions, vec3 cell_ext) {
    ASSERT(frame);
    if (positions.count == 0) {
        *frame = {};
    };
    vec3 min_box(FLT_MAX);
    vec3 max_box(-FLT_MAX);
    for (const auto& p : positions) {
        min_box = math::min(min_box, p);
        max_box = math::max(max_box, p);
    }
    min_box -= 1.f;
    max_box += 1.f;
    compute_frame(frame, positions, cell_ext, min_box, max_box);
}

Frame compute_frame(Array<const vec3> positions, vec3 cell_ext) {
    Frame frame;
    compute_frame(&frame, positions, cell_ext);
    return frame;
}

Frame compute_frame(Array<const vec3> positions, vec3 cell_ext, vec3 min_box, vec3 max_box) {
    Frame frame;
    compute_frame(&frame, positions, cell_ext, min_box, max_box);
    return frame;
}

void compute_frame(Frame* frame, Array<const vec3> positions, vec3 cell_ext, vec3 min_box, vec3 max_box) {
    ASSERT(frame);
    if (positions.count == 0) return;

    frame->min_box = min_box;
    frame->max_box = max_box;
    frame->cell_count = math::max(ivec3(1), ivec3((max_box - min_box) / cell_ext));
    frame->cell_ext = (max_box - min_box) / (vec3)frame->cell_count;
    frame->cells.resize(frame->cell_count.x * frame->cell_count.y * frame->cell_count.z);
    frame->entries.resize(positions.count);
    memset(frame->cells.data, 0, frame->cells.count * sizeof(Cell));

    unsigned int num_points = positions.count;
    unsigned int num_cells = frame->cells.count;

    uint32* l_idx = (uint32*)TMP_MALLOC(num_points * sizeof(uint32));
    uint32* g_idx = (uint32*)TMP_MALLOC(num_points * sizeof(uint32));
    defer {
        TMP_FREE(l_idx);
        TMP_FREE(g_idx);
    };

    for (int i = 0; i < positions.count; i++) {
        int cell_idx = compute_cell_idx(*frame, positions[i]);
        l_idx[i] = frame->cells[cell_idx].count++;
        g_idx[i] = cell_idx;
    }

    for (int i = 1; i < frame->cells.count; i++) {
        frame->cells[i].offset = frame->cells[i - 1].offset + frame->cells[i - 1].count;
    }

    for (int i = 0; i < frame->entries.count; i++) {
        int dst = frame->cells[g_idx[i]].offset + l_idx[i];
        frame->entries[dst].position = positions[i];
        frame->entries[dst].index = i;
    }
}

}  // namespace spatialhash
