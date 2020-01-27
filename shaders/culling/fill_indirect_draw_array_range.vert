#version 430

layout(std430, binding=0) buffer visiblity_buffer {
  int visible[];
};

layout(std430, binding=1) buffer cmd_buffer {
    uint cmd_data;
}

layout(binding = 3) uniform atomic_uint cmd_count;

in uint in_draw_offset;
in uint in_draw_count;

/*
    Layout of cmd_buffer:
    uint  count;
    uint  instanceCount;
    uint  first;
    uint  baseInstance;
*/

void main() {
    int idx = gl_VertexID;
    if (visible[i] == 1) {
        uint slot = atomicCounterIncrement(cmd_count);
        cmd_data[slot * 4U + 0U] = in_draw_count;
        cmd_data[slot * 4U + 1U] = 1U;
        cmd_data[slot * 4U + 2U] = in_draw_offset;
        cmd_data[slot * 4U + 3U] = 0U;
    }
}