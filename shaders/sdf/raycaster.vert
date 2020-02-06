#version 150 core

mat4 u_model_view_proj_mat;
vec4 u_model_eye_pos;

in vec3 in_pos;

out vec3 model_pos;
out vec3 model_eye;

void main() {
    model_pos = in_pos;
    model_eye = u_model_eye_pos.xyz;
    gl_Position = u_model_view_proj_mat * vec4(model_pos, 1);
}