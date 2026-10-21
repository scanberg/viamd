#version 410 core

// Isosurface meshes with stochastic transparency (Enderton et al., "Stochastic Transparency", 2010) behind two
// exactly peeled layers, a sample per pixel; the passes of a frame are programs of their own (ISOMESH_PASS)
// (isomesh_utils.cpp):
//
//   PASS_TRANSMITTANCE  every fragment, no depth test: sum log(1 - opacity) (the exact coverage of all of them
//                       together) and the optical depth of the regions they enclose (below), additively
//                       (ISOMESH_TRANSMITTANCE: the opacity alone, not the colour)
//   PASS_PEEL0          the nearest fragment of every pixel (depth test), and its opacity and depth (r, b)
//   PASS_PEEL1          the nearest one behind that, and its opacity and depth (g, a): the peels' colour target holds
//                       both when md_gl's geometry is drawn in the same passes (isomesh_utils.h): its accumulation reads it
//   PASS_DEPTH          the stochastic depth of the fragments behind those two: each covers the pixel with
//                       probability = its opacity, and the depth test keeps the nearest covering one
//   PASS_ACCUMULATE     every fragment adds its radiance times its visibility: 1 for the nearest, 1 - o0 for the
//                       second, (1 - o0)(1 - o1) where it is not behind the stochastic depth (test LEQUAL) for
//                       the rest
//
// The two nearest layers and the one behind them are exact; only the visibility of the layers behind those
// is estimated, from the layers in between. Without the peeling the nearest layer alone is exact, and where
// surfaces lie behind nearer ones (the back of a lobe, a lobe behind another) the estimate of their
// visibility depends on two or more random thresholds: the speckle that takes long to converge. Simulated per
// pixel (the thresholds below), four layers of opacity ~0.8: the error after 128 frames goes from ~3/255 to
// ~0.1/255.
//
// A sample per pixel, in plain textures: the depth the passes compare with the peeled ones is the rasterizer's,
// which with one sample is the pixel's, so nothing writes gl_FragDepth and the depth tests run before the shading.
// Each pass is a program of its own (ISOMESH_PASS), without the others' discards and registers. On an M4 Pro the
// four passes cost a tenth of what they did per sample into a multisampled target of two, in one program writing
// the depth it compared (under sample shading Apple's OpenGL gives gl_FragCoord.z at the pixel's centre, not at
// the sample the depth buffer holds). The edges are antialiased by the mean of the frames (isomesh_utils.h).
//
// Every fragment is shaded like a hit of the isosurface raycaster (isosurface.frag, the model of the deferred
// compose pass): a surface reflects L = a * diffuse + specular and lets 1 - opacity through, opacity =
// 1 - (1 - a)(1 - Fe) (the basic model: L = a * colour, opacity = a). The sum of L times visibility is in
// expectation the front-to-back composite, whatever the number of surfaces and the order of the draws. Unlike
// keeping only the nearest covering fragment (hashed transparency), only visibilities are estimated, not which
// surface is seen; isomesh_accumulate.frag scales the sum to the exact coverage (Enderton's alpha correction)
// and averages the frames, which converges to the composite.
//
// The coverage thresholds (PASS_DEPTH): from a blue noise threshold map (a void and cluster ordering,
// isomesh_utils.cpp) tiled over the screen, read at an offset per frame (the R2 sequence, so a pixel goes through
// the values of the map evenly) and per layer. A layer is a
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
// Isolines (gfx/isolines.h): on a surface coloured by the values of its vertices (u_coloring 3), the lines where
// the field takes given values are painted into its colour before it is shaded, so every pass sees the same
// opacity. Their width in pixels comes from ray differentials, as in the raycaster (isosurface.frag): the rays
// of the next pixel along x and y meet the tangent plane, and the field changes along that by its gradient
// along the surface, which the vertices carry (v_grad). No screen space derivatives: the passes depend on discard,
// and on Metal (Apple's OpenGL) a shader that takes derivatives does not stop at a discard, the invocation carries on
// for the derivatives' sake.

in vec3 v_world;
in vec3 v_view_pos;
in vec3 v_view_nrm;
in float v_value;
in vec3  v_grad;                    // the gradient of the field along the surface, view space, per world unit

uniform vec4  u_color;              // the surface's colour, alpha = coverage
uniform int   u_coloring;           // 0: u_color, 1: times the colour volume, 2: the field through its colour map,
                                    // 3: the value of the vertices (v_value) through the colour map
uniform mat4  u_world_to_volume;    // world -> the [0,1]^3 texture space of the colour / field volume
uniform sampler3D u_tex_color_volume;
uniform sampler3D u_tex_field;
uniform sampler2D u_tex_field_colormap;
uniform float u_field_beg;
uniform float u_field_inv_ext;

// Isolines (gfx/isolines.h), u_coloring 3 only: the persistent lines and the transient one, u_isoline_width
// pixels wide
#define ISOLINES_MAX 8
uniform int   u_isoline_count;
uniform float u_isoline_values[ISOLINES_MAX];
uniform int   u_isoline_hover;
uniform float u_isoline_hover_value;
uniform float u_isoline_width;

uniform sampler2D u_tex_depth;      // the opaque scene, single sampled
uniform int   u_use_depth;

uniform uint  u_frame;              // index of the frame being accumulated
uniform uint  u_surface_id;         // tells the draws apart (gl_PrimitiveID starts from 0 in each)

// A program per pass: ISOMESH_TRANSMITTANCE, or ISOMESH_PASS the number of one of the others (PASS_*)
#if !defined(ISOMESH_TRANSMITTANCE) && !defined(ISOMESH_PASS)
#error "isomesh.frag is built for a pass: ISOMESH_TRANSMITTANCE or ISOMESH_PASS"
#endif

uniform sampler2D   u_tex_tod;      // PASS_ACCUMULATE: the transmittance and optical depth of this frame
uniform sampler2D   u_tex_peel_depth0;  // PASS_PEEL1: the depth of the layer peeled before
uniform sampler2D   u_tex_peel;         // PASS_DEPTH, PASS_ACCUMULATE: the peeled layers' depths (z0, z1, z2)
uniform sampler2D   u_tex_peel_opacity; // PASS_ACCUMULATE: their opacities (o0, o1, o2)
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

#ifdef ISOMESH_TRANSMITTANCE
layout(location = 0) out vec4 out_tod;      // log(1 - opacity), optical depth
#else
layout(location = 0) out vec4 out_color;    // PASS_ACCUMULATE: radiance, opacity; PASS_PEEL*: the depth (z_k and after)
#if ISOMESH_PASS == 1 || ISOMESH_PASS == 2
layout(location = 1) out vec4 out_opacity;  // PASS_PEEL*: the opacity (o_k)
#endif
#endif

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

vec4 surface_color() {
    if (u_coloring == 1) {
        vec3 p = (u_world_to_volume * vec4(v_world, 1.0)).xyz;
        return u_color * texture(u_tex_color_volume, p);
    }
    if (u_coloring == 3) {
        float t = clamp((v_value - u_field_beg) * u_field_inv_ext, 0.0, 1.0);
        return u_color * vec4(texture(u_tex_field_colormap, vec2(t, 0.5)).rgb, 1.0);
    }
    if (u_coloring == 2) {
        vec3 p = (u_world_to_volume * vec4(v_world, 1.0)).xyz;
        float f = texture(u_tex_field, p).r;
        float t = clamp((f - u_field_beg) * u_field_inv_ext, 0.0, 1.0);
        return u_color * vec4(texture(u_tex_field_colormap, vec2(t, 0.5)).rgb, 1.0);
    }
    return u_color;
}

// The view space ray through normalized device coordinates ndc: its origin on the near plane, its direction
void isoline_view_ray(vec2 ndc, out vec3 o, out vec3 d) {
    vec4 a = u_inv_proj * vec4(ndc, -1.0, 1.0);
    vec4 b = u_inv_proj * vec4(ndc,  0.0, 1.0);
    o = a.xyz / a.w;
    d = b.xyz / b.w - o;
}

// Where the ray o + s d meets the plane through P with normal N. At a grazing angle the point runs off far: the
// field changes a lot from one pixel to the next there, which is what it does on the screen.
vec3 isoline_plane_hit(vec3 o, vec3 d, vec3 P, vec3 N) {
    float dn = dot(N, d);
    float lo = 1.0e-6 * length(d);
    dn = (dn < 0.0) ? min(dn, -lo) : max(dn, lo);
    return o + d * (dot(N, P - o) / dn);
}

// The change of the field per pixel along the surface at P (view space) with normal N, G the field's gradient:
// how far the point of the surface moves from this pixel's ray to the next pixel's, along x and along y, times G
float isoline_change_per_pixel(vec3 P, vec3 N, vec3 G) {
    vec2  ndc  = gl_FragCoord.xy * u_inv_res * 2.0 - 1.0;
    vec2  step = 2.0 * u_inv_res;
    vec3  o, d;
    isoline_view_ray(ndc, o, d);
    vec3  H  = isoline_plane_hit(o, d, P, N);
    isoline_view_ray(ndc + vec2(step.x, 0.0), o, d);
    vec3  Hx = isoline_plane_hit(o, d, P, N);
    isoline_view_ray(ndc + vec2(0.0, step.y), o, d);
    vec3  Hy = isoline_plane_hit(o, d, P, N);
    return length(vec2(dot(G, Hx - H), dot(G, Hy - H)));
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

#ifdef ISOMESH_TRANSMITTANCE
// The opacity alone, for the transmittance: the alpha of surface_color with apply_isolines' lines on it, and the
// Fresnel term of shade_surface

// The cover apply_isolines gives the lines (their rims), without their colours
float isoline_cover(float f, float g) {
    const float RIM = 2.0;              // as apply_isolines'
    float inv_g = 1.0 / max(g, 1.0e-30);
    float w = u_isoline_width + RIM;
    float cover = 0.0;
    for (int i = 0; i < u_isoline_count; ++i) {
        cover = max(cover, isoline_coverage(abs(f - u_isoline_values[i]) * inv_g, w));
    }
    if (u_isoline_hover != 0) {
        cover = max(cover, isoline_coverage(abs(f - u_isoline_hover_value) * inv_g, w));
    }
    return cover;
}

float surface_alpha(vec3 N) {
    float a = u_color.a;
    if (u_coloring == 1) {
        vec3 p = (u_world_to_volume * vec4(v_world, 1.0)).xyz;
        a *= texture(u_tex_color_volume, p).a;
    }
    if (u_coloring == 3 && (u_isoline_count > 0 || u_isoline_hover != 0)) {
        a = mix(a, 1.0, isoline_cover(v_value, isoline_change_per_pixel(v_view_pos, N, v_grad)));
    }
    return a;
}

float surface_opacity(float alpha, vec3 N, vec3 V) {
    float a = clamp(alpha, 0.0, 1.0);
    if (u_basic > 0.5) return a;
    float roughness = clamp(u_roughness, 0.04, 1.0);
    float F0 = clamp(u_F0, 0.0, 1.0);
    vec2  ab = EnvBRDFApprox(clamp(dot(N, V), 0.0, 1.0), roughness);
    float Fe = F0 * ab.x + ab.y;
    return 1.0 - (1.0 - a) * (1.0 - Fe);
}
#endif

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
#ifdef ISOMESH_TRANSMITTANCE
    ivec2 px = ivec2(gl_FragCoord.xy);
    float z_opaque = (u_use_depth != 0) ? texelFetch(u_tex_depth, px, 0).r : 1.0;
    bool  occluded = gl_FragCoord.z >= z_opaque;

    // The opacity of every fragment, outside the branch below: the colour volume is sampled with derivatives
    vec3 V = normalize(-v_view_pos);
    vec3 N = normalize(v_view_nrm);
    if (dot(N, V) < 0.0) N = -N;
    float opacity = clamp(surface_opacity(surface_alpha(N), N, V), 0.0, MAX_OPACITY);

    // Behind the opaque scene: hides nothing, and the ray through the enclosed regions ends at the scene
    vec3 p = v_view_pos;
    if (occluded) {
        vec4 v = u_inv_proj * vec4(gl_FragCoord.xy * u_inv_res * 2.0 - 1.0, z_opaque * 2.0 - 1.0, 1.0);
        p = v.xyz / v.w;
    }
    float od = (gl_FrontFacing ? -u_tau : u_tau) * (ray_distance(p) - u_od_ref);
    out_tod = vec4(occluded ? 0.0 : log(1.0 - opacity), od, 0.0, 0.0);
#else
    ivec2 px = ivec2(gl_FragCoord.xy);
    float z_opaque = (u_use_depth != 0) ? texelFetch(u_tex_depth, px, 0).r : 1.0;
    if (gl_FragCoord.z >= z_opaque) discard;       // behind the opaque scene

    vec3 V = normalize(-v_view_pos);
    vec3 N = normalize(v_view_nrm);
    if (dot(N, V) < 0.0) N = -N;

    // The passes before PASS_ACCUMULATE use the opacity alone: the compiler leaves the rest out of their programs
    vec4  c = surface_color();
    if (u_coloring == 3 && (u_isoline_count > 0 || u_isoline_hover != 0)) {
        c = apply_isolines(c, v_value, isoline_change_per_pixel(v_view_pos, N, v_grad));
    }
    vec3  L;
    float opacity;
    shade_surface(c.rgb, c.a, N, V, L, opacity);
    opacity = clamp(opacity, 0.0, MAX_OPACITY);
    if (opacity <= 0.0) discard;

    float z = gl_FragCoord.z;
#if ISOMESH_PASS == 1 || ISOMESH_PASS == 2
    // A peel: the nearest layer (behind the one peeled before), its depth and opacity in the channels of the pass
    if (ISOMESH_PASS == PASS_PEEL1 && z <= texelFetch(u_tex_peel_depth0, px, 0).r + PEEL_EPS) discard;
    out_color   = vec4(z);
    out_opacity = vec4(opacity);
    return;
#else
    // The peeled layers of the pixel from here on (a layer not peeled has the depth of the one in front)
    vec3 zp = texelFetch(u_tex_peel, px, 0).rgb;

    if (ISOMESH_PASS == PASS_DEPTH) {
        if (z <= zp.z + PEEL_EPS) discard;      // behind the last peel
        // the threshold, from the blue noise at an offset per layer and frame
        uint  layer = pcg(pcg(uint(gl_PrimitiveID)) ^ (u_surface_id * 0x9E3779B9u));
        uint  shift = pcg(layer ^ ((u_frame / LAYER_FRAMES + 1u) * 0x68E31DA4u));
        vec2  r2    = fract(vec2(0.7548776662, 0.5698402910) * float(u_frame) + 0.5);
        ivec2 q     = (px + ivec2(int(shift & 0xFFFFu), int(shift >> 16u)) + ivec2(r2 * float(BLUE_NOISE_SIZE))) & (BLUE_NOISE_SIZE - 1);
        float u0    = texelFetch(u_tex_blue_noise, q, 0).r;
        if (u0 >= opacity) discard;
        out_color = vec4(0.0);
        return;
    }

    // PASS_ACCUMULATE: the transmittance in front of this fragment
    vec3  o   = texelFetch(u_tex_peel_opacity, px, 0).rgb;
    float vis = (z <= zp.x + PEEL_EPS) ? 1.0 :
                (z <= zp.y + PEEL_EPS) ? (1.0 - o.x) :
                (z <= zp.z + PEEL_EPS) ? (1.0 - o.x) * (1.0 - o.y) : (1.0 - o.x) * (1.0 - o.y) * (1.0 - o.z);
    // Seen from inside (the back faces): behind the absorption of what the surfaces enclose
    if (!gl_FrontFacing) {
        L *= exp(-max(texelFetch(u_tex_tod, px, 0).g, 0.0));
    }
    out_color = vec4(L, opacity) * vis;
#endif
#endif
}
