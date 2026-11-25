#version 410 core

// The spheres of a VdW Surface (vdw_surface.h), one instance each: a quad facing the camera, as large as the
// sphere's silhouette, ray cast in vdw_surface.frag

uniform samplerBuffer u_sphere;         // xyz: the centre (world, Angstrom), w: the radius
uniform mat4 u_world_to_view;
uniform mat4 u_view_to_clip;
uniform int  u_ortho;

out vec3 v_view;
flat out int  v_sphere;
flat out vec3 v_c;
flat out float v_r;

void main() {
    vec4 S = texelFetch(u_sphere, gl_InstanceID);
    vec3 c = (u_world_to_view * vec4(S.xyz, 1.0)).xyz;
    float R = S.w;
    float d = length(c);
    v_sphere = gl_InstanceID; v_c = c; v_r = R;
    if (u_ortho == 0 && d <= R * 1.0001) {
        // the camera inside it: nothing of it to draw from here
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
        v_view = vec3(0.0);
        v_sphere = -1;
        return;
    }
    vec2 corner = vec2(float(gl_VertexID & 1) * 2.0 - 1.0, float((gl_VertexID >> 1) & 1) * 2.0 - 1.0);
    vec3 p;
    if (u_ortho != 0) {
        p = c + vec3(corner * R * 1.01, R);
    } else {
        vec3 w = c / d;
        vec3 u = normalize(abs(w.y) < 0.99 ? cross(w, vec3(0, 1, 0)) : cross(w, vec3(1, 0, 0)));
        vec3 v = cross(u, w);
        float ext = R * d / sqrt(d * d - R * R) * 1.01;
        p = c + (corner.x * u + corner.y * v) * ext;
    }
    v_view = p;
    gl_Position = u_view_to_clip * vec4(p, 1.0);
}
