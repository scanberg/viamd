#version 410 core

// Temporal accumulation of the half-res AO, between the AO pass and the spatial blur.
//
// The blur can only remove the interleaved sampling pattern where it finds a full period of it on the same surface.
// In crevasses, on small atoms and along silhouettes the bilateral weights cut it short and the raw per-pixel estimate
// shows through. There is no spatial neighbourhood left to borrow samples from there, so this pass borrows from earlier
// frames instead: ssao.frag permutes the pattern every frame and the history averages each pixel over all of it.
//
// Unlike the TAA resolve, history is not clamped to the colour of the current neighbourhood (AO noise is much coarser
// than a 3x3 window). It is accepted or rejected on geometry: a history texel is used when it lies on the tangent plane
// of the current surface point, both expressed in the previous frame's view space.
//
// The history is a jitter average, like the TAA history: texels are reprojected with the jitter-free velocity, so with a
// still camera every texel reads exactly its own history while the camera jitter moves the sample point under it.
#ifndef AO_PERSPECTIVE
#define AO_PERSPECTIVE 1
#endif
#define N_STILL     32.0    // history length (frames) for still content: two periods of the 16 frame pattern cycle
#define N_MOVING    4.0     // history length for content that moves MOVING_PX or more per frame. Screen-space AO
#define MOVING_PX   1.0     // changes with the view, so a long history lags behind a rotation (full-res px per frame)
#define TOL_PX      1.5     // tangent-plane distance, in full-res pixel footprints, below which history is fully trusted
#define MIN_NZ      0.25    // |N.z| floor for the plane test, so a tap far behind a silhouette is not on its plane
#define VIS_PACK    0.5     // history.x = n + VIS_PACK * visibility (frame count and visibility share one float)

uniform sampler2D u_tex_ao;         // RG32F (visibility, linear depth) of this frame, half res
uniform sampler2D u_tex_history;    // RG32F (n + VIS_PACK * visibility, linear depth) of the previous frame, half res
uniform sampler2D u_tex_normal;     // full-res G-buffer normal
uniform sampler2D u_tex_velocity;   // full-res screen-space motion in uv, jitter excluded (prev uv = uv - vel)

uniform vec4  u_proj_info;
uniform vec4  u_proj_info_prev;
uniform mat4  u_curr_to_prev;       // view space (z into the screen) of this frame -> view space of the previous frame
uniform vec2  u_jitter_delta;       // (jitter this frame - jitter previous frame), full-res px
uniform vec2  u_full_res;
uniform float u_px_scale;
uniform float u_z_max;
uniform int   u_history_valid;

layout(location = 0) out vec2 out_history;
layout(location = 1) out vec2 out_ao;

vec3 uv_to_view(vec2 uv, float z, vec4 proj_info) {
#if AO_PERSPECTIVE
    return vec3((uv * proj_info.xy + proj_info.zw) * z, z);
#else
    return vec3((uv * proj_info.xy + proj_info.zw), z);
#endif
}

vec3 decode_normal(vec2 enc) {
    vec2 fenc = enc * 4.0 - 2.0;
    float f = dot(fenc, fenc);
    float g = sqrt(1.0 - f / 4.0);
    return vec3(fenc * g, 1.0 - f / 2.0);
}

// Full-res pixel picked for half-res texel h by the rotated-grid downsample (see depth_downsample.frag)
ivec2 rep_px(ivec2 h) { return h * 2 + ivec2(h.y & 1, h.x & 1); }

void main() {
    ivec2 hp = ivec2(gl_FragCoord.xy);
    vec2  c  = texelFetch(u_tex_ao, hp, 0).rg;
    if (c.y >= u_z_max) {
        out_history = vec2(VIS_PACK, c.y);
        out_ao = c;
        return;
    }

    ivec2 fp = rep_px(hp);
    vec2  px = vec2(fp) + 0.5;
    vec3  P  = uv_to_view(px / u_full_res, c.y, u_proj_info);
    vec3  N  = decode_normal(texelFetch(u_tex_normal, fp, 0).xy) * vec3(1, 1, -1);
    vec2  vel_px = texelFetch(u_tex_velocity, fp, 0).xy * u_full_res;

#if AO_PERSPECTIVE
    float px_world = u_px_scale * c.y;
#else
    float px_world = u_px_scale;
#endif

    // This surface point in the previous frame. Its depth comes from the camera motion; its screen position from the
    // velocity, which also carries the motion of the atoms. Only atom motion along the view axis is unaccounted for.
    float zp = (u_curr_to_prev * vec4(P, 1.0)).z;
    vec3  Pp = uv_to_view((px - vel_px + u_jitter_delta) / u_full_res, zp, u_proj_info_prev);
    vec3  Np = mat3(u_curr_to_prev) * N;
    Np = normalize(vec3(Np.xy, (Np.z < 0.0 ? -1.0 : 1.0) * max(abs(Np.z), MIN_NZ)));
    float inv_tol = 1.0 / (TOL_PX * px_world);

    // Where the history of this texel is, in half-res texels, with the jitter left out. History is interpolated on
    // the regular half-res grid: still content then reads exactly its own texel (the rotated-grid picks of diagonal
    // neighbours are only sqrt(2) px apart, a filter on true positions would blend them in every frame).
    vec2 hq = vec2(hp) - 0.5 * vel_px;

    float n_prev = 0.0;
    float v_prev = 0.0;
    ivec2 half_size = max(ivec2(u_full_res) / 2, ivec2(1));
    if (u_history_valid != 0 && all(greaterThan(hq, vec2(-0.5))) && all(lessThan(hq, vec2(half_size) - 0.5))) {
        ivec2 h0 = ivec2(floor(hq));
        vec2  f  = hq - vec2(h0);
        float w_sum  = 0.0, v_sum  = 0.0, n_sum  = 0.0;   // bilinear (2x2 at h0) x geometric
        float wf_sum = 0.0, vf_sum = 0.0, nf_sum = 0.0;   // geometric only (3x3 around h0), fallback
        for (int y = -1; y <= 1; ++y) {
            for (int x = -1; x <= 1; ++x) {
                ivec2 h = h0 + ivec2(x, y);
                if (any(lessThan(h, ivec2(0))) || any(greaterThanEqual(h, half_size))) continue;
                vec2  wb2 = mix(1.0 - f, f, vec2(x, y));
                float wb  = (x >= 0 && y >= 0) ? wb2.x * wb2.y : 0.0;

                vec2 hist = texelFetch(u_tex_history, h, 0).rg;
                if (hist.y >= u_z_max) continue;    // background last frame: disoccluded
                vec3  H  = uv_to_view((vec2(rep_px(h)) + 0.5) / u_full_res, hist.y, u_proj_info_prev);
                float wg = clamp(2.0 - abs(dot(H - Pp, Np)) * inv_tol, 0.0, 1.0);
                if (wg <= 0.0) continue;
                float n = floor(hist.x);
                float v = (hist.x - n) * (1.0 / VIS_PACK);
                w_sum  += wb * wg;  v_sum  += wb * wg * v;  n_sum  += wb * wg * n;
                wf_sum += wg;       vf_sum += wg * v;       nf_sum += wg * n;
            }
        }
        if (w_sum > 1e-3) {
            v_prev = v_sum / w_sum;
            n_prev = n_sum / w_sum;
        } else if (wf_sum > 1e-3) {
            // The texel's own history is on another surface. Typical in crevasses, where the camera jitter moves the
            // picked pixel from one atom to the next: continue from the neighbours on this surface, at reduced trust.
            v_prev = vf_sum / wf_sum;
            n_prev = 0.5 * nf_sum / wf_sum;
        }
    }

    // n must stay integral: it shares history.x with the visibility (VIS_PACK)
    float n_max = floor(mix(N_STILL, N_MOVING, clamp(length(vel_px) / MOVING_PX, 0.0, 1.0)) + 0.5);
    float n = min(floor(n_prev + 1e-3) + 1.0, n_max);
    float vis = mix(v_prev, c.x, 1.0 / n);

    out_history = vec2(n + VIS_PACK * min(vis, 1.0), c.y);
    out_ao = vec2(vis, c.y);
}
