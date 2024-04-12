#version 150 core

layout (std140) uniform UniformData {
    mat4  u_inv_proj_mat;
    vec3  u_bg_color;
    float u_time;
    vec3  u_env_radiance;
    float u_roughness;
    vec3  u_dir_radiance;
    float u_F0;
    vec3  u_light_dir;
};

uniform sampler2D u_texture_depth;
uniform sampler2D u_texture_color;
uniform sampler2D u_texture_normal;

in vec2 tc;
out vec4 out_frag;

// TODO: Use linear depth instead and use uniform vec4 for unpacking to view coords.
vec4 depth_to_view_coord(vec2 tex_coord, float depth) {
    vec4 clip_coord = vec4(vec3(tex_coord, depth) * 2.0 - 1.0, 1.0);
    vec4 view_coord = u_inv_proj_mat * clip_coord;
    return view_coord / view_coord.w;
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

const float PI = 3.1415926535;
const float ONE_OVER_PI = 1.0 / 3.1415926535;

// Cook-Torrance model From here:
// https://learnopengl.com/PBR/Theory

float FresnelSchlick(float cosTheta, float F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

float FresnelSchlickRoughness(float cosTheta, float F0, float roughness) {
    return F0 + (max(1.0 - roughness, F0) - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}   

float DistributionGGX(float NdotH, float roughness) {
    float a      = roughness*roughness;
    float a2     = a*a;
    float NdotH2 = NdotH*NdotH;
    float num   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
    return num / denom;
}

float GeometrySchlickGGX(float NdotV, float roughness) {
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;
    float num   = NdotV;
    float denom = NdotV * (1.0 - k) + k;
    return num / denom;
}

float GeometrySmith(float NdotV, float NdotL, float roughness) {
    float ggx2  = GeometrySchlickGGX(NdotV, roughness);
    float ggx1  = GeometrySchlickGGX(NdotL, roughness);
    return ggx1 * ggx2;
}

vec3 shade(vec3 color, vec3 V, vec3 N) {
    vec3  L = u_light_dir; 
    vec3  H = normalize(L + V);
    float H_dot_V = clamp(dot(H, V), 0.0, 1.0);
    float N_dot_H = clamp(dot(N, H), 0.0, 1.0);
    float N_dot_L = clamp(dot(N, L), 0.0, 1.0);
    float N_dot_V = clamp(dot(N, V), 0.0, 1.0);

    float F0 = u_F0;
    float roughness = u_roughness;
    vec3  albedo = color;

    vec3 Lo = vec3(0);

    {
        // Add contribution from directional light
        float NDF = DistributionGGX(N_dot_H, roughness);
        float G   = GeometrySmith(N_dot_V, N_dot_L, roughness);
        float F   = FresnelSchlick(H_dot_V, F0);
        float kS = F;
        float kD = 1.0 - kS;
        float numerator   = NDF * G * F;
        float denominator = 4.0 * N_dot_V * N_dot_L + 0.0001;
        vec3 specular     = vec3(numerator / denominator);
        Lo += (kD * albedo * ONE_OVER_PI + specular) * u_dir_radiance * N_dot_L;
    }

    {
        // Add contribution from environment
        float F = FresnelSchlickRoughness(N_dot_V, F0, roughness);
        float kS = F;
        float kD = 1.0 - kS;
        Lo += (kD * albedo * ONE_OVER_PI + vec3(kS)) * u_env_radiance; // Multiply with ao here
    }

    return Lo;
}

void main() {
    vec3 result = u_bg_color.rgb;
    
    vec4 color = texelFetch(u_texture_color, ivec2(gl_FragCoord.xy), 0);
    if (color.a > 0) {
        float depth = texelFetch(u_texture_depth, ivec2(gl_FragCoord.xy), 0).x;
        vec3 normal = decode_normal(texelFetch(u_texture_normal, ivec2(gl_FragCoord.xy), 0).xy);
        vec4 view_coord = depth_to_view_coord(tc, depth);

        // Add noise to reduce banding
        // Signed random to not affect overall luminance
        color.rgb = clamp(color.rgb + color.rgb * srand4(tc + u_time).rgb * 0.15, 0.0, 1.0);

        vec3 N = normal;
        vec3 V = -normalize(view_coord.xyz);
        color.rgb = shade(color.rgb, V, N);
        result = color.rgb;
    }

    out_frag = vec4(result, 1);
}