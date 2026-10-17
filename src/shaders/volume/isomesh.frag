#version 410 core

// Isosurface meshes with stochastic transparency (Enderton et al., "Stochastic Transparency", 2010) behind two
// exactly peeled layers, into a multisampled target; the passes of one frame selected by u_pass
// (isomesh_utils.cpp):
//
//   PASS_TRANSMITTANCE  every fragment, no depth test: sum log(1 - opacity) (the exact coverage of all of them
//                       together) and the optical depth of the regions they enclose (below), additively
//   PASS_PEEL0          the nearest fragment of every sample (depth test), and its opacity
//   PASS_PEEL1          the nearest one behind that (per sample), and its opacity
//   PASS_DEPTH          the stochastic depth of the fragments behind those two: each covers a sample with
//                       probability = its opacity, and the depth test keeps the nearest covering one
//   PASS_ACCUMULATE     every fragment adds its radiance times its visibility in each sample: 1 for the
//                       nearest, 1 - o0 for the second, (1 - o0)(1 - o1) where it is not behind the stochastic
//                       depth (test LEQUAL) for the rest
//
// The two nearest layers and the one behind them are exact; only the visibility of the layers behind those
// is estimated, from the layers in between. Without the peeling the nearest layer alone is exact, and where
// surfaces lie behind nearer ones (the back of a lobe, a lobe behind another) the estimate of their
// visibility depends on two or more random thresholds: the speckle that takes long to converge. Simulated per
// pixel, four layers of opacity ~0.8: the error after 64 frames goes from ~2/255 to ~0.04/255. PASS_PEEL1,
// PASS_DEPTH and PASS_ACCUMULATE run per sample (sample shading), comparing their depth with the peeled ones;
// every pass writes the depth it compares (gl_FragDepth), so a fragment always finds its own.
//
// Every fragment is shaded like a hit of the isosurface raycaster (isosurface.frag, the model of the deferred
// compose pass): a surface reflects L = a * diffuse + specular and lets 1 - opacity through, opacity =
// 1 - (1 - a)(1 - Fe) (the basic model: L = a * colour, opacity = a). The sum of L times visibility is in
// expectation the front-to-back composite, whatever the number of surfaces and the order of the draws. Unlike
// keeping only the nearest covering fragment (hashed transparency), only visibilities are estimated, not which
// surface is seen; isomesh_accumulate.frag scales the sum to the exact coverage (Enderton's alpha correction)
// and averages the frames, which converges to the composite.
//
// The coverage thresholds (PASS_DEPTH): stratified over the samples of a pixel (u + i/n for sample i), u from a
// blue noise threshold map (a void and cluster ordering, isomesh_utils.cpp) tiled over the screen, read at an
// offset per frame (the R2 sequence, so a pixel goes through the values of the map evenly) and per layer. A layer is a
// triangle - the finest there is, so two sheets that overlap in a pixel never share their thresholds (the
// draws told apart by u_surface_id, as gl_PrimitiveID starts from 0 in each) - and its offset is drawn again
// every few frames, so no two layers stay at one relative offset for long. Thresholds correlated between the
// layers of a pixel would not average to the composite; offsets of the map leave them close to independent,
// and the noise of every frame blue. An offset of the threshold VALUES per layer instead (a constant over the
// triangle) makes the error of a frame a constant over each triangle: triangles that flicker as it converges.
//
// The optical density of what the surfaces enclose (the raycaster's absorption between them): an entry (a
// front face: triangles are counter-clockwise seen from outside) adds -tau * d, an exit +tau * d, d the distance
// along the ray, so whatever the order the sum is the integral of tau over the ray inside the surfaces. A face
// behind the opaque scene counts at its depth, so the ray ends there; a ray that starts inside (the camera in a
// lobe) has an exit without an entry and counts from the eye. It dims what lies behind all the surfaces and the
// surfaces seen from inside, the back of a lobe behind its front: exact for one enclosed region along a ray;
// with several, the exits of the nearer ones are dimmed by the farther ones as well, and the entries of the
// farther ones are not dimmed by the nearer ones.
//
// Isolines (gfx/isolines.h): on a surface coloured by a field, the lines where the field takes given values are
// painted into its colour before it is shaded, their width in pixels from the screen space derivatives of the
// field. Every pass paints them alike, so the opacity they give a surface is the same in all of them.

in vec3 v_world;
in vec3 v_view_pos;
in vec3 v_view_nrm;
in float v_value;

uniform vec4  u_color;              // the surface's colour, alpha = coverage
uniform int   u_coloring;           // 0: u_color, 1: times the colour volume, 2: the field through its colour map,
                                    // 3: the value of the vertices (v_value) through the colour map
uniform mat4  u_world_to_volume;    // world -> the [0,1]^3 texture space of the colour / field volume
uniform sampler3D u_tex_color_volume;
uniform sampler3D u_tex_field;
uniform sampler2D u_tex_field_colormap;
uniform float u_field_beg;
uniform float u_field_inv_ext;

// Isolines (gfx/isolines.h), where u_coloring is a field (2, 3): lines where the field takes these values, the
// persistent ones and the transient one, u_isoline_width pixels wide
#define ISOLINES_MAX 8
uniform int   u_isoline_count;
uniform float u_isoline_values[ISOLINES_MAX];
uniform int   u_isoline_hover;
uniform float u_isoline_hover_value;
uniform float u_isoline_width;

uniform sampler2D u_tex_depth;      // the opaque scene, single sampled
uniform int   u_use_depth;

uniform int   u_pass;               // PASS_*
uniform uint  u_frame;              // index of the frame being accumulated
uniform int   u_samples;            // samples per pixel of the target
uniform uint  u_surface_id;         // tells the draws apart (gl_PrimitiveID starts from 0 in each)

uniform sampler2DMS u_tex_tod;      // PASS_ACCUMULATE: the transmittance and optical depth of this frame
uniform sampler2DMS u_tex_peel_depth0;  // the depth of the nearest layer (PASS_PEEL0)
uniform sampler2DMS u_tex_peel_depth1;  // and of the second (PASS_PEEL1)
uniform sampler2DMS u_tex_peel_opacity; // their opacities (r, g)
uniform float u_tau;                // optical density of the surface times the extinction per world unit
uniform float u_od_ref;             // distance subtracted from d, so the sums stay small (half floats)
uniform mat4  u_inv_proj;
uniform sampler2D u_tex_blue_noise; // BLUE_NOISE_SIZE^2 thresholds
uniform vec2  u_inv_res;
uniform int   u_ortho;

// The lights and the model of the compose pass (as in isosurface.frag)
uniform vec3  u_light_dir;          // toward the key light, view space
uniform vec3  u_dir_radiance;
uniform vec3  u_sky_radiance;
uniform vec3  u_ground_radiance;
uniform float u_roughness;
uniform float u_F0;
uniform float u_basic;
uniform float u_basic_ambient;
uniform float u_basic_specular;
uniform float u_linear_colors;
uniform float u_light_strength;

layout(location = 0) out vec4 out_color;    // PASS_ACCUMULATE: radiance, opacity; PASS_PEEL*: the opacity
layout(location = 1) out vec4 out_tod;      // PASS_TRANSMITTANCE: log(1 - opacity), optical depth

const int PASS_TRANSMITTANCE = 0;
const int PASS_PEEL0         = 1;
const int PASS_PEEL1         = 2;
const int PASS_DEPTH         = 3;
const int PASS_ACCUMULATE    = 4;

// A fragment at most this far behind a peeled depth is that layer: the same rasterized depth, within rounding
const float PEEL_EPS = 2.0e-7;

const int  BLUE_NOISE_SIZE = 64;      // a power of two, as isomesh_utils.cpp makes it
const uint LAYER_FRAMES    = 8u;      // frames a layer keeps its offset of the map

// Fully opaque surfaces would make the log infinite: a little light gets through all of them
const float MAX_OPACITY = 0.999;

const float PI     = 3.1415926535;
const float INV_PI = 1.0 / PI;

float FresnelSchlick(float cos_theta, float F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cos_theta, 0.0, 1.0), 5.0);
}

vec2 EnvBRDFApprox(float NdotV, float roughness) {
    const vec4 c0 = vec4(-1.0, -0.0275, -0.572, 0.022);
    const vec4 c1 = vec4( 1.0,  0.0425,  1.04, -0.04);
    vec4  r = roughness * c0 + c1;
    float a004 = min(r.x * r.x, exp2(-9.28 * NdotV)) * r.x + r.y;
    return vec2(-1.04, 1.04) * a004 + r.zw;
}

vec3 environment_irradiance(vec3 N) {
    return mix(u_ground_radiance, u_sky_radiance, 0.5 + 0.5 * N.y);
}

vec3 environment_radiance(vec3 R, float roughness) {
    float w = 0.1 + roughness;
    return mix(u_ground_radiance, u_sky_radiance, smoothstep(-w, w, R.y));
}

vec3 srgb_to_linear(vec3 c) {
    c = clamp(c, 0.0, 1.0);
    return mix(c / 12.92, pow((c + 0.055) / 1.055, vec3(2.4)), step(vec3(0.04045), c));
}

float DistributionGGX(float NdotH, float roughness) {
    float a  = roughness * roughness;
    float a2 = a * a;
    float d  = NdotH * NdotH * (a2 - 1.0) + 1.0;
    return a2 / (PI * d * d);
}

float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = roughness + 1.0;
    float k = (r * r) / 8.0;
    return NdotV / (NdotV * (1.0 - k) + k);
}

const float BASIC_SHININESS = 48.0;

vec3 shade_basic(vec3 color, vec3 N, vec3 V) {
    float w = 0.5 + 0.5 * dot(N, u_light_dir);
    vec3  c = clamp(color, 0.0, 1.0) * (u_basic_ambient + u_light_strength * (1.0 - u_basic_ambient) * w * w);
    c /= max(1.0, max(c.r, max(c.g, c.b)));
    vec3  H = normalize(u_light_dir + V);
    float s = u_light_strength * u_basic_specular * pow(max(dot(N, H), 0.0), BASIC_SHININESS);
    c = mix(c, vec3(1.0), clamp(s, 0.0, 1.0));
    return (u_linear_colors > 0.5) ? c : srgb_to_linear(c);
}

// What the surface sends toward the eye (premultiplied: L) and how much of what lies behind it it hides
void shade_surface(vec3 color, float alpha, vec3 N, vec3 V, out vec3 L, out float opacity) {
    float a = clamp(alpha, 0.0, 1.0);
    if (u_basic > 0.5) {
        L = a * shade_basic(color, N, V);
        opacity = a;
        return;
    }
    float roughness = clamp(u_roughness, 0.04, 1.0);
    float F0 = clamp(u_F0, 0.0, 1.0);
    vec3  albedo = (u_linear_colors > 0.5) ? clamp(color, 0.0, 1.0) : srgb_to_linear(color);

    float NdotV = clamp(dot(N, V), 0.0, 1.0);
    float NdotL = clamp(dot(N, u_light_dir), 0.0, 1.0);

    vec3 diffuse  = vec3(0.0);
    vec3 specular = vec3(0.0);
    if (NdotL > 0.0) {
        vec3  H     = normalize(u_light_dir + V);
        float NdotH = clamp(dot(N, H), 0.0, 1.0);
        float HdotV = clamp(dot(H, V), 0.0, 1.0);
        float D = DistributionGGX(NdotH, roughness);
        float G = GeometrySchlickGGX(NdotV, roughness) * GeometrySchlickGGX(NdotL, roughness);
        float F = FresnelSchlick(HdotV, F0);
        vec3  E = u_dir_radiance * u_light_strength;
        specular += vec3(D * G * F / (4.0 * NdotV * NdotL + 0.0001)) * E * NdotL;
        diffuse  += (1.0 - F) * albedo * INV_PI * E * NdotL;
    }
    vec2  ab = EnvBRDFApprox(NdotV, roughness);
    float Fe = F0 * ab.x + ab.y;
    diffuse  += (1.0 - Fe) * albedo * environment_irradiance(N);
    specular += Fe * environment_radiance(reflect(-V, N), roughness);

    L = a * diffuse + specular;
    opacity = 1.0 - (1.0 - a) * (1.0 - Fe);
}

// The field the surface is coloured by, at this fragment: the vertices' values interpolated, or the field volume
float field_value() {
    if (u_coloring == 3) {
        return v_value;
    }
    if (u_coloring == 2) {
        vec3 p = (u_world_to_volume * vec4(v_world, 1.0)).xyz;
        return texture(u_tex_field, p).r;
    }
    return 0.0;
}

// f: field_value()
vec4 surface_color(float f) {
    if (u_coloring == 1) {
        vec3 p = (u_world_to_volume * vec4(v_world, 1.0)).xyz;
        return u_color * texture(u_tex_color_volume, p);
    }
    if (u_coloring == 2 || u_coloring == 3) {
        float t = clamp((f - u_field_beg) * u_field_inv_ext, 0.0, 1.0);
        return u_color * vec4(texture(u_tex_field_colormap, vec2(t, 0.5)).rgb, 1.0);
    }
    return u_color;
}

// The coverage of a line w pixels wide at distance d (pixels) from its middle, box filtered over a pixel
float isoline_coverage(float d, float w) {
    return clamp(0.5 * w + 0.5 - d, 0.0, 1.0);
}

// The isolines painted onto the colour of the surface: f the field here, g its change per pixel along the
// surface, so |f - v| / g is the distance in pixels to the line of value v. The persistent lines are dark with a
// light rim and the transient one light with a dark rim, which reads on any colour of a map and tells the two
// apart. A line is opaque paint, whatever the coverage of the surface under it.
vec4 apply_isolines(vec4 c, float f, float g) {
    const vec3  DARK  = vec3(0.02);
    const vec3  LIGHT = vec3(0.95);
    const float RIM   = 2.0;            // pixels, both sides together
    float inv_g = 1.0 / max(g, 1.0e-30);
    float w = u_isoline_width;

    float core = 0.0;
    float rim  = 0.0;
    for (int i = 0; i < u_isoline_count; ++i) {
        float d = abs(f - u_isoline_values[i]) * inv_g;
        core = max(core, isoline_coverage(d, w));
        rim  = max(rim,  isoline_coverage(d, w + RIM));
    }
    c.rgb = mix(c.rgb, LIGHT, rim);
    c.rgb = mix(c.rgb, DARK, core);
    float cover = rim;

    if (u_isoline_hover != 0) {
        float d  = abs(f - u_isoline_hover_value) * inv_g;
        float hc = isoline_coverage(d, w);
        float hr = isoline_coverage(d, w + RIM);
        c.rgb = mix(c.rgb, DARK, hr);
        c.rgb = mix(c.rgb, LIGHT, hc);
        cover = max(cover, hr);
    }
    c.a = mix(c.a, 1.0, cover);
    return c;
}

uint pcg(uint v) {
    uint s = v * 747796405u + 2891336453u;
    uint w = ((s >> ((s >> 28u) + 4u)) ^ s) * 277803737u;
    return (w >> 22u) ^ w;
}

// Distance along the ray to view position p
float ray_distance(vec3 p) {
    return (u_ortho != 0) ? -p.z : length(p);
}

void main() {
    // The depth the depth tests see is the one the passes compare with the peeled depths: written, not left to
    // the rasterizer. Under sample shading gl_FragCoord.z is not everywhere the depth at the sample - Apple's
    // OpenGL (on Metal) gives the depth at the pixel's centre, while the depth buffer holds the sample's - and a
    // fragment then misses its own peeled depth by its slope times the offset of the sample, far beyond
    // PEEL_EPS: the layers of whole triangles misclassified, bands of them where the slope turns. Where
    // gl_FragCoord.z is the sample's (Mesa, NVIDIA, AMD) this is the depth the rasterizer writes, bit for bit.
    gl_FragDepth = gl_FragCoord.z;

    // The field and its change per pixel, taken here in uniform control flow: derivatives after a discard are
    // undefined. Across a triangle the vertex values are linear, so the lines are straight within it and meet at
    // its edges.
    float f = field_value();
    float g = length(vec2(dFdx(f), dFdy(f)));

    ivec2 px = ivec2(gl_FragCoord.xy);
    float z_opaque = (u_use_depth != 0) ? texelFetch(u_tex_depth, px, 0).r : 1.0;
    bool  occluded = gl_FragCoord.z >= z_opaque;
    if (occluded && u_pass != PASS_TRANSMITTANCE) discard;

    vec3 V = normalize(-v_view_pos);
    vec3 N = normalize(v_view_nrm);
    if (dot(N, V) < 0.0) N = -N;

    vec4  c = surface_color(f);
    if (u_isoline_count > 0 || u_isoline_hover != 0) {
        c = apply_isolines(c, f, g);
    }
    vec3  L;
    float opacity;
    shade_surface(c.rgb, c.a, N, V, L, opacity);
    opacity = clamp(opacity, 0.0, MAX_OPACITY);

    if (u_pass == PASS_TRANSMITTANCE) {
        // Behind the opaque scene: hides nothing, and the ray through the enclosed regions ends at the scene
        vec3 p = v_view_pos;
        if (occluded) {
            vec4 v = u_inv_proj * vec4(gl_FragCoord.xy * u_inv_res * 2.0 - 1.0, z_opaque * 2.0 - 1.0, 1.0);
            p = v.xyz / v.w;
        }
        float od = (gl_FrontFacing ? -u_tau : u_tau) * (ray_distance(p) - u_od_ref);
        out_tod = vec4(occluded ? 0.0 : log(1.0 - opacity), od, 0.0, 0.0);
        return;
    }

    if (opacity <= 0.0) discard;

    if (u_pass == PASS_PEEL0) {
        out_color = vec4(opacity);
        return;
    }

    // Per sample from here on
    float z  = gl_FragCoord.z;
    float z0 = texelFetch(u_tex_peel_depth0, px, gl_SampleID).r;
    if (u_pass == PASS_PEEL1) {
        if (z <= z0 + PEEL_EPS) discard;
        out_color = vec4(opacity);
        return;
    }
    float z1 = texelFetch(u_tex_peel_depth1, px, gl_SampleID).r;

    if (u_pass == PASS_DEPTH) {
        if (z <= z1 + PEEL_EPS) discard;
        // stratified thresholds over the samples, from the blue noise at an offset per layer and frame
        uint  layer = pcg(pcg(uint(gl_PrimitiveID)) ^ (u_surface_id * 0x9E3779B9u));
        uint  shift = pcg(layer ^ ((u_frame / LAYER_FRAMES + 1u) * 0x68E31DA4u));
        vec2  r2    = fract(vec2(0.7548776662, 0.5698402910) * float(u_frame) + 0.5);
        ivec2 q     = (px + ivec2(int(shift & 0xFFFFu), int(shift >> 16u)) + ivec2(r2 * float(BLUE_NOISE_SIZE))) & (BLUE_NOISE_SIZE - 1);
        float u0    = texelFetch(u_tex_blue_noise, q, 0).r;
        if (fract(u0 + float(gl_SampleID) / float(u_samples)) >= opacity) discard;
        out_color = vec4(0.0);
        return;
    }

    // PASS_ACCUMULATE: the transmittance in front of this fragment in this sample
    vec2  o01 = texelFetch(u_tex_peel_opacity, px, gl_SampleID).rg;
    float vis = (z <= z0 + PEEL_EPS) ? 1.0 : (z <= z1 + PEEL_EPS) ? (1.0 - o01.r) : (1.0 - o01.r) * (1.0 - o01.g);
    // Seen from inside (the back faces): behind the absorption of what the surfaces enclose
    if (!gl_FrontFacing) {
        L *= exp(-max(texelFetch(u_tex_tod, px, gl_SampleID).g, 0.0));
    }
    out_color = vec4(L, opacity) * vis;
}
