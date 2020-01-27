#version 430

layout(early_fragment_tests) in;
layout(std430,binding=0) buffer visiblity_buffer {
  int visible[];
};

flat in uint id;
out vec4 out_frag;

void main() {
    visible[id] = 1;
    out_frag = unpackUnorm4x8(id);
}