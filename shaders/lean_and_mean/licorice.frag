#version 150 core
#extension GL_ARB_conservative_depth : enable
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_proj_mat;
uniform vec4 u_color;

in Fragment {
    smooth vec4 view_coord;
    flat vec4 capsule_center_radius;
    flat vec4 capsule_axis_length;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
layout(location = 0) out vec4 out_color;

// Source from Ingo Quilez (https://www.shadertoy.com/view/Xt3SzX)
float intersect_capsule(in vec3 ro, in vec3 rd, in vec3 cc, in vec3 ca, float cr,
                      float ch)  // cc center, ca orientation axis, cr radius, ch height
{
    vec3 oc = ro - cc;
    ch *= 0.5;

    float card = dot(ca, rd);
    float caoc = dot(ca, oc);

    float a = 1.0 - card * card;
    float b = dot(oc, rd) - caoc * card;
    float c = dot(oc, oc) - caoc * caoc - cr * cr;
    float h = b * b - a * c;
    if (h < 0.0) return -1.0;
    float t = (-b - sqrt(h)) / a;

    float y = caoc + t * card;

    // body
    if (abs(y) < ch) {
        return t;
    }

    // caps
    float sy = sign(y);
    oc = ro - (cc + sy * ca * ch);
    b = dot(rd, oc);
    c = dot(oc, oc) - cr * cr;
    h = b * b - c;
    if (h > 0.0) {
        t = -b - sqrt(h);
        return t;
    }

    return -1.0;
}

void main() {
    vec3 ro = vec3(0);
    vec3 rd = normalize(in_frag.view_coord.xyz);
    vec3 cc = in_frag.capsule_center_radius.xyz;
    float cr = in_frag.capsule_center_radius.w;
    vec3 ca = in_frag.capsule_axis_length.xyz;
    float ch = in_frag.capsule_axis_length.w;

    float t = intersect_capsule(ro, rd, cc, ca, cr, ch);
    if (t < 0.0) {
        discard;
        return;
    }

    vec4 view_coord = vec4(rd * t, 1);
    vec4 curr_clip_coord = u_proj_mat * view_coord;

    gl_FragDepth = (curr_clip_coord.z / curr_clip_coord.w) * 0.5 + 0.5;
    out_color = u_color;
}