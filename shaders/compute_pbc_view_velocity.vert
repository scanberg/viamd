#version 150

uniform mat4 u_view_mat;
uniform vec3 u_box_ext;

layout (location = 0) in vec3 in_position;
layout (location = 0) in vec3 in_old_position;

out vec3 out_view_velocity;

vec3 de_periodize(vec3 pos, vec3 ref_pos, vec3 ext) {
    vec3 delta = pos - ref_pos;
    vec3 signed_mask = sign(delta) * step(ext * 0.5, abs(delta));
    return pos - full_ext * signed_mask;
}

void main() {
    vec3 cur = in_position;
    vec3 old = de_periodize(in_old_position, in_position, u_box_ext);
    out_view_velocity = cur - old;
}