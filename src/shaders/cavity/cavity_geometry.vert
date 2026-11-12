#version 410 core

// The spheres of a cavity (cavity.h), one instance each: a quad facing the camera through the centre of the sphere, as
// large as the cone of the rays that touch it cuts that plane (R d / sqrt(d^2 - R^2), the silhouette exactly)

uniform samplerBuffer u_spheres;        // xyz: centre (the frame of the file), w: radius, 0 for none
uniform mat4  u_cavity_to_view;
uniform mat4  u_view_to_clip;
uniform int   u_ortho;

out vec3 v_view;
flat out int  v_sphere;
flat out vec4 v_center;                 // view space centre, radius

void main() {
    int  i = gl_InstanceID;
    vec4 s = texelFetch(u_spheres, i);
    vec3 c = (u_cavity_to_view * vec4(s.xyz, 1.0)).xyz;
    float R = s.w;
    float d = length(c);
    v_sphere = i;
    v_center = vec4(c, R);
    // None, or the camera within it: nothing to draw
    if (R <= 0.0 || (u_ortho == 0 && d <= R * 1.0001)) {
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
        v_view = vec3(0);
        return;
    }
    vec2 corner = vec2(float(gl_VertexID & 1) * 2.0 - 1.0, float((gl_VertexID >> 1) & 1) * 2.0 - 1.0);
    vec3 p;
    if (u_ortho != 0) {
        p = c + vec3(corner * R * 1.01, 0.0);
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
