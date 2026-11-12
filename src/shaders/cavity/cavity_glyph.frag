#version 410 core

// The glyph of a tessera (cavity_glyph.vert), ray cast: a cylinder capped at both ends, or a capsule. Into the G-buffer
// as the opaque geometry writes it (the colour, the normal, the velocity of the camera's motion and the picking index of
// the tessera), or, for the highlight, into the marks as cavity_marks.frag marks a patch.

in vec3 v_view;
flat in int   v_tess;
flat in vec3  v_p0;
flat in vec3  v_p1;
flat in float v_r;

uniform samplerBuffer u_tess_b;
uniform samplerBuffer u_glyph;
uniform sampler2D     u_colormap;
uniform sampler2D     u_tex_depth;      // marks: the G-buffer's depth (jittered, at its render scale)

uniform mat4  u_view_to_clip;
uniform mat4  u_curr_view_to_prev_clip;
uniform mat4  u_clip_to_view;           // marks: of the G-buffer's depth
uniform vec4  u_jitter_uv;
uniform int   u_ortho;
uniform int   u_capsule;
uniform int   u_color_mode;             // 0: the field over u_range, 1: the weight, 2: uniform
uniform vec2  u_range;
uniform int   u_unlit;
uniform uint  u_tessera_base;
uniform int   u_marks;
uniform float u_proj_11;
uniform vec2  u_inv_size;

layout(location = 0) out vec4 out_color;    // marks: the marks
layout(location = 1) out vec4 out_normal;   // marks: the group
layout(location = 2) out vec4 out_velocity;
layout(location = 3) out vec4 out_picking;

vec2 encode_normal(vec3 n) {
    float p = sqrt(n.z * 8.0 + 8.0);
    return n.xy / p + 0.5;
}

vec4 encode_index(uint index) {
    return vec4(
        (index & 0x000000FFU) >> 0U,
        (index & 0x0000FF00U) >> 8U,
        (index & 0x00FF0000U) >> 16U,
        (index & 0xFF000000U) >> 24U) / 255.0;
}

// The nearest hit of the ray with the capped cylinder from pa to pb of radius ra: t and the normal, t < 0 for none
vec4 cylinder(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra) {
    vec3  ba = pb - pa;
    vec3  oc = ro - pa;
    float baba = dot(ba, ba);
    float bard = dot(ba, rd);
    float baoc = dot(ba, oc);
    float k2 = baba - bard * bard;
    float k1 = baba * dot(oc, rd) - baoc * bard;
    float k0 = baba * dot(oc, oc) - baoc * baoc - ra * ra * baba;
    float h = k1 * k1 - k2 * k0;
    if (h < 0.0) return vec4(-1.0);
    h = sqrt(h);
    float t = (-k1 - h) / k2;
    float y = baoc + t * bard;
    if (y > 0.0 && y < baba) return vec4(t, (oc + t * rd - ba * y / baba) / ra);
    t = (((y < 0.0) ? 0.0 : baba) - baoc) / bard;
    if (abs(k1 + k2 * t) < h) return vec4(t, ba * sign(y) / sqrt(baba));
    return vec4(-1.0);
}

// And with the capsule
vec4 capsule(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float ra) {
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;
    float baba = dot(ba, ba);
    float bard = dot(ba, rd);
    float baoa = dot(ba, oa);
    float a = baba - bard * bard;
    float b = baba * dot(rd, oa) - baoa * bard;
    float c = baba * dot(oa, oa) - baoa * baoa - ra * ra * baba;
    float h = b * b - a * c;
    float t = -1.0;
    if (h >= 0.0) {
        t = (-b - sqrt(h)) / a;
        float y = baoa + t * bard;
        if (y <= 0.0 || y >= baba) {
            vec3 oc = (y <= 0.0) ? oa : ro - pb;
            b = dot(rd, oc);
            c = dot(oc, oc) - ra * ra;
            h = b * b - c;
            t = (h > 0.0) ? -b - sqrt(h) : -1.0;
        }
    }
    if (t < 0.0) return vec4(-1.0);
    vec3  pa_x = ro + rd * t - pa;
    float s = clamp(dot(pa_x, ba) / baba, 0.0, 1.0);
    return vec4(t, (pa_x - s * ba) / ra);
}

void main() {
    if (v_tess < 0) discard;
    vec3 ro = (u_ortho != 0) ? vec3(v_view.xy, 0.0) : vec3(0.0);
    vec3 rd = (u_ortho != 0) ? vec3(0.0, 0.0, -1.0) : normalize(v_view);
    vec4 hit = (u_capsule != 0) ? capsule(ro, rd, v_p0, v_p1, v_r) : cylinder(ro, rd, v_p0, v_p1, v_r);
    if (hit.x < 0.0) discard;
    vec3 x = ro + rd * hit.x;
    vec3 n = normalize(hit.yzw);

    vec4 clip = u_view_to_clip * vec4(x, 1.0);
    gl_FragDepth = (clip.z / clip.w) * 0.5 + 0.5;

    if (u_marks != 0) {
        // In view: no closer than the scene within a tolerance, as cavity_marks.frag decides it
        vec2  uv = gl_FragCoord.xy * u_inv_size;
        ivec2 size  = textureSize(u_tex_depth, 0);
        ivec2 texel = clamp(ivec2(uv * vec2(size)), ivec2(0), size - 1);
        float depth = texelFetch(u_tex_depth, texel, 0).r;
        vec4  scene = u_clip_to_view * vec4(uv * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
        float scene_z = -scene.z / scene.w;
        float z = -x.z;
        float pixel  = (u_ortho != 0) ? 2.0 / (u_proj_11 * float(size.y)) : 2.0 * z / (u_proj_11 * float(size.y));
        float facing = (u_ortho != 0) ? abs(n.z) : abs(dot(n, normalize(x)));
        float slope  = sqrt(max(1.0 - facing * facing, 0.0)) / max(facing, 0.125);
        bool  in_view = z <= scene_z + 2.0 * pixel * slope + 2e-3 * z;
        out_color  = vec4(0.0, 0.0, in_view ? 1.0 : 0.0, 1.0);
        out_normal = vec4((1.0 + (in_view ? 8.0 : 0.0)) / 255.0, 0.0, 0.0, 0.0);
        return;
    }

    vec3 color = vec3(0.75);
    if (u_color_mode != 2) {
        float param;
        if (u_color_mode == 0) {
            float ext = u_range.y - u_range.x;
            param = abs(ext) > 0.0 ? clamp((texelFetch(u_tess_b, v_tess).w - u_range.x) / ext, 0.0, 1.0) : 0.5;
        } else {
            param = clamp(texelFetch(u_glyph, v_tess).z, 0.0, 1.0);
        }
        color = texture(u_colormap, vec2(param * (255.0 / 256.0) + 0.5 / 256.0, 0.5)).rgb;
    }

    vec4 prev = u_curr_view_to_prev_clip * vec4(x, 1.0);
    vec2 ss_vel = (clip.xy / clip.w - prev.xy / prev.w) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);

    out_color    = vec4(color, u_unlit != 0 ? 1.0 / 255.0 : 1.0);
    out_normal   = vec4(encode_normal(n), 0.0, 0.0);
    out_velocity = vec4(ss_vel, 0.0, 0.0);
    out_picking  = encode_index(u_tessera_base + uint(v_tess));
}
