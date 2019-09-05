#version 150 core

uniform mat4 u_proj_mat;
uniform float u_radius = 1.0;

layout (lines) in;
layout (triangle_strip, max_vertices = 24) out;

in Vertex {
    flat int discard_vert;
} in_vert[];

out Fragment {
    smooth vec4 view_coord;
    flat vec4 capsule_center_radius;
    flat vec4 capsule_axis_length;
} out_frag;

vec4 view_vertices[8];
vec4 proj_vertices[8];

void emit_vertex(int i){
    out_frag.view_coord = view_vertices[i];
    gl_Position = proj_vertices[i];
    EmitVertex();
}

void emit(int a, int b, int c, int d)
{
    emit_vertex(a);
    emit_vertex(b);
    emit_vertex(c);
    emit_vertex(d);
    EndPrimitive(); 
}

vec3 get_ortho_vec(vec3 v, vec3 A, vec3 B){
    if(abs(1-dot(v,A))>0.001){
        return normalize(cross(v,A));
    }else{
        return normalize(cross(v,B));
    }
}

void main()
{
    if (in_vert[0].discard_vert != 0 || in_vert[1].discard_vert != 0) {
        EndPrimitive();
        return;
    }

    // Compute orientation vectors for the two connecting faces:
    vec3 p0 = gl_in[0].gl_Position.xyz;
    vec3 p1 = gl_in[1].gl_Position.xyz;
    float r = u_radius;
    float l = distance(p0, p1);
    vec3 a = (p1 - p0) / l;
    vec3 c = (p0 + p1) * 0.5;

    out_frag.capsule_center_radius = vec4(c, r);
    out_frag.capsule_axis_length = vec4(a, l);

    // Extend end points to properly fit the sphere caps
    p0 -= a * r;
    p1 += a * r;

    vec3 B = vec3(0,0,1);
    vec3 A = vec3(1,0,0);
    vec3 o = get_ortho_vec(a,A,B);

    // Declare scratch variables for basis vectors:
    vec3 i,j,k;

    // Compute vertices of prismoid:
    j = a; i = o; k = normalize(cross(i, j)); i = normalize(cross(k, j)); ; i *= r; k *= r;
    view_vertices[0] = vec4(p0 + i + k, 1);
    view_vertices[1] = vec4(p0 + i - k, 1);
    view_vertices[2] = vec4(p0 - i - k, 1);
    view_vertices[3] = vec4(p0 - i + k, 1);
    view_vertices[4] = vec4(p1 + i + k, 1);
    view_vertices[5] = vec4(p1 + i - k, 1);
    view_vertices[6] = vec4(p1 - i - k, 1);
    view_vertices[7] = vec4(p1 - i + k, 1);

    proj_vertices[0] = u_proj_mat * view_vertices[0];
    proj_vertices[1] = u_proj_mat * view_vertices[1];
    proj_vertices[2] = u_proj_mat * view_vertices[2];
    proj_vertices[3] = u_proj_mat * view_vertices[3];
    proj_vertices[4] = u_proj_mat * view_vertices[4];
    proj_vertices[5] = u_proj_mat * view_vertices[5];
    proj_vertices[6] = u_proj_mat * view_vertices[6];
    proj_vertices[7] = u_proj_mat * view_vertices[7];

    // Emit the six faces of the prismoid:
    emit(0,1,3,2); emit(5,4,6,7);
    emit(4,5,0,1); emit(3,2,7,6);
    emit(0,3,4,7); emit(2,1,6,5);
}