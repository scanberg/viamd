#version 430

// @TODO: This is bogus, this should fill an element index buffer instead.

layout(std430, binding=0) buffer visiblity_buffer {
  int visible[];
};

layout(std430, binding=1) buffer cmd_buffer {
    uint cmd_data[];
};

layout(binding = 2) uniform atomic_uint cmd_count;

/*
    Layout of cmd_buffer:
    uint  count;
    uint  instanceCount;
    uint  first;
    uint  baseInstance;
*/

void main() {
    int idx = gl_VertexID;
    if (visible[idx] == 1) {
        uint slot = atomicCounterIncrement(cmd_count);
        cmd_data[slot * 4 + 0] = 1;
        cmd_data[slot * 4 + 1] = 1;
        cmd_data[slot * 4 + 2] = uint(idx);
        cmd_data[slot * 4 + 3] = 0;
    }
}