#version 330 core

/*
    This shader can be optimized alot furhter.
    Especially for Intel integrated GPUs, where we want to minimize the workload in the geometry shader
    Move matrix multiplications and transformations to vertex shader
    Output 'Ready' corner points for the box that is the body for the ribbon-segment
    Output 'Ready' normals for the edges of the boxes
    In geometry shader: Pick the right indices and construct geometry
*/

uniform mat4 u_view_proj_mat;
uniform vec2 u_scale = vec2(1.0, 0.1);

layout(lines) in;
layout(triangle_strip, max_vertices = 24) out;

in Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    int discard_vert;
} in_vert[];

void main() {
    if (in_vert[0].discard_vert != 0 || in_vert[1].discard_vert != 0) {
        EndPrimitive();
        return;
    }

    vec4 p[2];
    vec4 x[2];
    vec4 y[2];
    vec4 z[2];
    
    p[0] = in_vert[0].control_point;
    p[1] = in_vert[1].control_point;
    x[0] = in_vert[0].support_vector;
    x[1] = in_vert[1].support_vector * sign(dot(in_vert[0].support_vector, in_vert[1].support_vector));
    z[0] = in_vert[0].support_tangent;
    z[1] = in_vert[1].support_tangent;
    y[0] = vec4(cross(z[0].xyz, x[0].xyz), 0); // To maintain right-handedness
    y[1] = vec4(cross(z[1].xyz, x[1].xyz), 0);

#if 0
    // This is to possible fix tesselation direction so it follows the curvature of the segment
    float flip_sign = sign(dot(x[1].xyz, y[0].xyz));
    x[0].xyz *= flip_sign;
    x[1].xyz *= flip_sign;
    y[0].xyz *= flip_sign;
    y[1].xyz *= flip_sign;
    z[0].xyz *= flip_sign;
    z[1].xyz *= flip_sign;
#endif

    x[0] *= u_scale.x;
    x[1] *= u_scale.x;
    y[0] *= u_scale.y;
    y[1] *= u_scale.y;
    float sz = 1.f;

    mat4 M[2];
    M[0] = mat4(x[0], y[0], z[0], p[0]);
    M[1] = mat4(x[1], y[1], z[1], p[1]);

    vec4 v0[4];
    v0[0] = M[0] * vec4(vec2(-1,-1), 0, 1);
    v0[1] = M[0] * vec4(vec2( 1,-1), 0, 1);
    v0[2] = M[0] * vec4(vec2(-1, 1), 0, 1);
    v0[3] = M[0] * vec4(vec2( 1, 1), 0, 1);

    vec4 curr_clip0[4];
    curr_clip0[0] = u_view_proj_mat * v0[0];
    curr_clip0[1] = u_view_proj_mat * v0[1];
    curr_clip0[2] = u_view_proj_mat * v0[2];
    curr_clip0[3] = u_view_proj_mat * v0[3];

    vec4 v1[4];
    v1[0] = M[1] * vec4(vec2(-1,-1), 0, 1);
    v1[1] = M[1] * vec4(vec2( 1,-1), 0, 1);
    v1[2] = M[1] * vec4(vec2(-1, 1), 0, 1);
    v1[3] = M[1] * vec4(vec2( 1, 1), 0, 1);

    vec4 curr_clip1[4];
    curr_clip1[0] = u_view_proj_mat * v1[0];
    curr_clip1[1] = u_view_proj_mat * v1[1];
    curr_clip1[2] = u_view_proj_mat * v1[2];
    curr_clip1[3] = u_view_proj_mat * v1[3];
 
    // BOTTOM
    gl_Position = curr_clip0[0]; EmitVertex();
    gl_Position = curr_clip0[1]; EmitVertex();
    gl_Position = curr_clip1[0]; EmitVertex();
    gl_Position = curr_clip1[1]; EmitVertex();
    EndPrimitive();

    // TOP
    gl_Position = curr_clip0[3]; EmitVertex();
    gl_Position = curr_clip0[2]; EmitVertex();
    gl_Position = curr_clip1[3]; EmitVertex();
    gl_Position = curr_clip1[2]; EmitVertex();
    EndPrimitive();

    // LEFT
    gl_Position = curr_clip0[2]; EmitVertex();
    gl_Position = curr_clip0[0]; EmitVertex();
    gl_Position = curr_clip1[2]; EmitVertex();
    gl_Position = curr_clip1[0]; EmitVertex();
    EndPrimitive();

    // RIGHT
    gl_Position = curr_clip0[1]; EmitVertex();
    gl_Position = curr_clip0[3]; EmitVertex();
    gl_Position = curr_clip1[1]; EmitVertex();
    gl_Position = curr_clip1[3]; EmitVertex();
    EndPrimitive();

    // FRONT
    gl_Position = curr_clip1[0]; EmitVertex();
    gl_Position = curr_clip1[1]; EmitVertex();
    gl_Position = curr_clip1[2]; EmitVertex();
    gl_Position = curr_clip1[3]; EmitVertex();
    EndPrimitive();

    // BACK
    gl_Position = curr_clip0[1]; EmitVertex();
    gl_Position = curr_clip0[0]; EmitVertex();
    gl_Position = curr_clip0[3]; EmitVertex();
    gl_Position = curr_clip0[2]; EmitVertex();
    EndPrimitive();
}