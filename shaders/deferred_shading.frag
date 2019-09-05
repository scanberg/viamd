#version 150 core

uniform sampler2D u_texture_depth;
uniform sampler2D u_texture_color;
uniform sampler2D u_texture_normal;

uniform mat4 u_inv_proj_mat;
uniform float u_time;

in vec2 tc;
out vec4 out_frag;

// TODO: Use linear depth instead and use uniform vec4 for unpacking to view coords.

vec4 depth_to_view_coord(vec2 tex_coord, float depth) {
    vec4 clip_coord = vec4(vec3(tex_coord, depth) * 2.0 - 1.0, 1.0);
    vec4 view_coord = u_inv_proj_mat * clip_coord;
    return view_coord / view_coord.w;
}

float fresnel(float H_dot_V) {
    const float n1 = 1.0;
    const float n2 = 1.5;
    const float R0 = pow((n1-n2)/(n1+n2), 2);

    return R0 + (1.0 - R0)*pow(1.0 - H_dot_V, 5);
}

// https://aras-p.info/texts/CompactNormalStorage.html
vec3 decode_normal(vec2 enc) {
    vec2 fenc = enc*4-2;
    float f = dot(fenc,fenc);
    float g = sqrt(1-f/4.0);
    vec3 n;
    n.xy = fenc*g;
    n.z = 1-f/2.0;
    return n;
}

vec4 rand4(vec2 n) {
    return fract(sin(dot(n.xy, vec2(12.9898, 78.233)))* vec4(43758.5453, 28001.8384, 50849.4141, 12996.89));
}

vec4 srand4(vec2 n) {
    return rand4(n) * 2.0 - 1.0;
}

const vec3 env_radiance = vec3(5.0);
const vec3 dir_radiance = vec3(10.0);
const vec3 L = normalize(vec3(1,1,1));
const float spec_exp = 100.0;

vec3 lambert(in vec3 radiance) {
    const float ONE_OVER_PI = 1.0 / 3.1415926535;
    return radiance * ONE_OVER_PI;
}

vec3 shade(vec3 color, vec3 V, vec3 N) {
    vec3 H = normalize(L + V);
    float H_dot_V = max(0.0, dot(H, V));
    float N_dot_H = max(0.0, dot(N, H));
    float N_dot_L = max(0.0, dot(N, L));
    float fr = fresnel(H_dot_V);

    vec3 diffuse = color.rgb * lambert(env_radiance + N_dot_L * dir_radiance);
    vec3 specular = fr * (env_radiance + dir_radiance) * pow(N_dot_H, spec_exp);

    return diffuse + specular;
}

void main() {
    float depth = texelFetch(u_texture_depth, ivec2(gl_FragCoord.xy), 0).x;
    if (depth == 1.0) {
        discard;
    }
    vec4 color = texelFetch(u_texture_color, ivec2(gl_FragCoord.xy), 0);
    vec3 normal = decode_normal(texelFetch(u_texture_normal, ivec2(gl_FragCoord.xy), 0).xy);
    vec4 view_coord = depth_to_view_coord(tc, depth);

    // Add noise to reduce banding
    vec4 noise4 = srand4(tc + u_time + 0.6959174) / 10.0;
    color += color*noise4;

    vec3 N = normal;
    vec3 V = -normalize(view_coord.xyz);
    vec3 result = shade(color.rgb, V, N);
    //result = N;

    out_frag = vec4(result, color.a);
}