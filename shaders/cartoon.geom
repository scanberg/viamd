#version 330 core

uniform mat4 u_normal_mat;
uniform mat4 u_view_proj_mat;
uniform vec2 u_scale = vec2(1.0, 1.0);

layout(lines) in;
layout(triangle_strip, max_vertices = 26) out;

in Vertex {
    vec4 control_point;
    vec4 support_vector;
    vec4 support_tangent;
    vec4 color;
    vec4 picking_color;
    vec4 weights;
} in_vert[];

out Fragment {
    smooth vec4 color;
    smooth vec4 view_normal;
    flat vec4 picking_color;
} out_frag;

void emit_vertex(in vec4 pos, in vec4 normal, in int idx) {
    out_frag.color = in_vert[idx].color;
    out_frag.picking_color = in_vert[idx].picking_color;
    out_frag.view_normal = normal;
    gl_Position = pos;
    EmitVertex();
}

vec2 elipse_normal(in vec2 p, in float a, in float b) {
    return normalize(vec2(p.x * a, p.y * b));
}

void main() {
    if (in_vert[0].color.a == 0.0f || in_vert[1].color.a == 0.0f) return;

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

    const float TWO_PI = 2.0 * 3.14159265;

    mat4 M[2];
    M[0] = mat4(x[0], y[0], z[0], p[0]);
    M[1] = mat4(x[1], y[1], z[1], p[1]);

    mat4 N[2];
    N[0] = u_normal_mat * M[0];
    N[1] = u_normal_mat * M[1];

    mat4 V[2];
    V[0] = u_view_proj_mat * M[0];
    V[1] = u_view_proj_mat * M[1];

    vec4 w0 = in_vert[0].weights;
    vec4 w1 = in_vert[1].weights;

    const vec2 coil_scale = vec2(0.2, 0.2);
    const vec2 sheet_scale = vec2(1.5, 0.1);
    const vec2 helix_scale = vec2(1.2, 0.2);

    vec2 s0 = u_scale * (w0.x * coil_scale + w0.y * sheet_scale + w0.z * helix_scale);
    vec2 s1 = u_scale * (w1.x * coil_scale + w1.y * sheet_scale + w1.z * helix_scale);

    vec4 v0[12];
    vec4 v1[12];
    
    vec4 n0[12];
    vec4 n1[12];

    for (int i = 0; i < 12; i++) {
        float t = float(i) / 12.0 * TWO_PI;
        vec2 x = vec2(cos(t), sin(t));
        vec2 v = x * s0;
        vec2 n = x / s0;
        v0[i] = V[0] * vec4(x * s0,0,1);
        n0[i] = N[0] * vec4(x / s0,0,0);
        v1[i] = V[1] * vec4(x * s1,0,1);
        n1[i] = N[1] * vec4(x / s1,0,0);
    }
    
    for (int i = 0; i < 12; i++) {
        emit_vertex(v1[i], n1[i], 1);
        emit_vertex(v0[i], n0[i], 0);
    }
    emit_vertex(v1[0], n1[0], 1);
    emit_vertex(v0[0], n0[0], 0);
}