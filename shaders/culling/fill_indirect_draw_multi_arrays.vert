#version 430

layout(std430, binding=0) buffer visiblity_buffer {
  int visible[];
};

layout(std430, binding=1) buffer cmd_buffer {
    uint cmd_data[];
};

layout(binding = 2) uniform atomic_uint cmd_count;

in uint in_draw_offset;
in uint in_draw_count;

/*
    Layout of cmd_buffer:
    uint  count;
    uint  instanceCount;
    uint  first;
    uint  baseInstance;
*/

bool is_visible(int idx) {
    return visible[idx] == 1;
}

void main() {
    int idx = gl_VertexID;
    if (is_visible(idx)) {
        uint slot = atomicCounterIncrement(cmd_count);
        cmd_data[slot * 4 + 0] = in_draw_count;
        cmd_data[slot * 4 + 1] = 1;
        cmd_data[slot * 4 + 2] = in_draw_offset;
        cmd_data[slot * 4 + 3] = 0;
    }
}