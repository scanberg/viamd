#version 410 core

// Isosurface meshes (isomesh_utils.cpp): positions in the mesh's own units (bohr), placed by u_model.

layout(location = 0) in vec3 a_pos;
layout(location = 1) in vec3 a_nrm;
layout(location = 2) in float a_value;  // a field at the vertex (isomesh::mesh_upload_values), 0 without
layout(location = 3) in vec3  a_grad;   // its gradient along the surface, per mesh unit, 0 without

uniform mat4 u_model;       // mesh -> world (a uniform scale and a translation)
uniform mat4 u_view;
uniform mat4 u_proj;

out vec3 v_world;
out vec3 v_view_pos;
out vec3 v_view_nrm;
out float v_value;
out vec3  v_grad;       // a_grad in view space, per world unit

invariant gl_Position;

void main() {
    vec4 w  = u_model * vec4(a_pos, 1.0);
    vec4 vp = u_view * w;
    v_world    = w.xyz;
    v_view_pos = vp.xyz;
    v_view_nrm = mat3(u_view) * (mat3(u_model) * a_nrm);
    v_value    = a_value;
    // A gradient goes with the inverse transpose of u_model, which for its uniform scale s is u_model / s^2
    float s2   = dot(u_model[0].xyz, u_model[0].xyz);
    v_grad     = mat3(u_view) * (mat3(u_model) * a_grad) / max(s2, 1.0e-30);
    gl_Position = u_proj * vp;
}
