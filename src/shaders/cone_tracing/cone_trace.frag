#version 330 core

in vec2 uv;
out vec4 frag_color;

// Textures
uniform sampler2D u_depth_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_color_alpha_texture;
uniform sampler2D u_f0_smoothness_texture;

// Voxel stuff
uniform sampler3D u_voxel_texture;
uniform vec3 u_voxel_grid_world_size;
uniform vec3 u_voxel_grid_world_min;
uniform ivec3 u_voxel_dimensions;
uniform float u_voxel_extent;

// Scaling factors
uniform float u_indirect_diffuse_scale = 1.f;
uniform float u_indirect_specular_scale = 1.f;
uniform float u_ambient_occlusion_scale = 1.f;
uniform float u_cone_angle = 0.08;

// View parameters
uniform mat4 u_inv_view_mat;
uniform mat4 u_inv_view_proj_mat;
uniform vec3 u_world_space_camera;

const float MAX_DIST = 1000.0;
const float ALPHA_THRESH = 0.95;

// 6 60 degree cone
const int NUM_CONES = 6;
const vec3 cone_directions[6] = vec3[]
(                            vec3(0, 1, 0),
                            vec3(0, 0.5, 0.866025),
                            vec3(0.823639, 0.5, 0.267617),
                            vec3(0.509037, 0.5, -0.700629),
                            vec3(-0.509037, 0.5, -0.700629),
                            vec3(-0.823639, 0.5, 0.267617)
                            );
const float cone_weights[6] = float[](0.25, 0.15, 0.15, 0.15, 0.15, 0.15);

// // 5 90 degree cones
// const int NUM_CONES = 5;
// vec3 coneDirections[5] = vec3[]
// (                            vec3(0, 1, 0),
//                             vec3(0, 0.707, 0.707),
//                             vec3(0, 0.707, -0.707),
//                             vec3(0.707, 0.707, 0),
//                             vec3(-0.707, 0.707, 0)
//                             );
// float coneWeights[5] = float[](0.28, 0.18, 0.18, 0.18, 0.18);

vec4 sample_voxels(vec3 world_position, float lod) {
    vec3 tc = (world_position - u_voxel_grid_world_min) / (u_voxel_grid_world_size);
    //vec3 d = tc - vec3(0.5);
    //if (dot(d,d) > 1) return vec4(0,0,0,1);
    //tc = tc + vec3(0.5) / vec3(u_voxel_dimensions);
    return textureLod(u_voxel_texture, tc, lod);
}

// Third argument to say how long between steps?
vec4 cone_trace(vec3 world_position, vec3 world_normal, vec3 direction, float tan_half_angle, out float occlusion) {
    
    // lod level 0 mipmap is full size, level 1 is half that size and so on
    float lod = 0.0;
    vec4 rgba = vec4(0);
    occlusion = 0.0;

    float dist = u_voxel_extent; // Start one voxel away to avoid self occlusion
    vec3 start_pos = world_position + world_normal * u_voxel_extent; // Plus move away slightly in the normal direction to avoid
                                                                     // self occlusion in flat surfaces

    while(dist < MAX_DIST && rgba.a < ALPHA_THRESH) {
        // smallest sample diameter possible is the voxel size
        float diameter = max(u_voxel_extent, 2.0 * tan_half_angle * dist);
        float lod_level = log2(diameter / u_voxel_extent);
        vec4 voxel_color = sample_voxels(start_pos + dist * direction, lod_level);
        //if (voxel_color.a < 0) break;

        // front-to-back compositing
        float a = (1.0 - rgba.a);
        rgba += a * voxel_color;
        occlusion += (a * voxel_color.a) / (1.0 + 0.03 * diameter);

        //dist += diameter * 0.5; // smoother
        dist += diameter; // faster but misses more voxels
    }

    return rgba;
}

vec4 indirect_light(in vec3 world_position, in vec3 world_normal, in mat3 tangent_to_world, out float occlusion_out) {
    vec4 color = vec4(0);
    occlusion_out = 0.0;

    for(int i = 0; i < NUM_CONES; i++) {
        float occlusion = 0.0;
        // 60 degree cones -> tan(30) = 0.577
        // 90 degree cones -> tan(45) = 1.0
        color += cone_weights[i] * cone_trace(world_position, world_normal, tangent_to_world * cone_directions[i], 0.577, occlusion);
        occlusion_out += cone_weights[i] * occlusion;
    }

    occlusion_out = 1.0 - occlusion_out;

    return color;
}

mat3 compute_ON_basis(in vec3 v1) {
    vec3 v0;
    vec3 v2;
    float d0 = dot(vec3(1,0,0), v1);
    float d1 = dot(vec3(0,0,1), v1);
    if (d0 < d1) {
        v0 = normalize(vec3(1,0,0) - v1 * d0);
    } else {
        v0 = normalize(vec3(0,0,1) - v1 * d1);
    }
    v2 = cross(v0, v1);
    return mat3(v0, v1, v2);
}

vec4 depth_to_world_coord(vec2 tex_coord, float depth) {
    vec4 clip_coord = vec4(vec3(tex_coord, depth) * 2.0 - 1.0, 1.0);
    vec4 world_coord = u_inv_view_proj_mat * clip_coord;
    return world_coord / world_coord.w;
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

vec3 fresnel(vec3 f0, float H_dot_V) {
    //const float n1 = 1.0;
    //const float n2 = 1.5;
    //const float R0 = pow((n1-n2)/(n1+n2), 2);
    return f0 + (1.0 - f0) * pow(1.0 - H_dot_V, 5.0);
}

vec3 shade(vec3 albedo, float alpha, vec3 f0, float smoothness, vec3 P, vec3 V, vec3 N) {
    const float PI_QUARTER = 3.14159265 * 0.25;
    const vec3 env_radiance = vec3(0.5);
    const vec3 dir_radiance = vec3(0.5);
    const vec3 L = normalize(vec3(1));
    const float spec_exp = 10.0;

    mat3 tangent_to_world = compute_ON_basis(N);

    float N_dot_V = max(0.0, -dot(N, V));
    vec3 R = -V + 2.0 * dot(N, V) * N;
    vec3 H = normalize(L + V);
    float H_dot_V = max(0.0, dot(H, V));
    float N_dot_H = max(0.0, dot(N, H));
    float N_dot_L = max(0.0, dot(N, L));
    vec3 fresnel_direct = fresnel(f0, H_dot_V);
    vec3 fresnel_indirect = fresnel(f0, N_dot_V);
    float tan_half_angle = tan(mix(PI_QUARTER, 0.0, smoothness));

    float diffuse_occlusion;
    vec3 direct_diffuse = env_radiance + N_dot_L * dir_radiance;
    vec3 indirect_diffuse = indirect_light(P, N, tangent_to_world, diffuse_occlusion).rgb;

    float transmissive_occlusion;
    vec3 transmissive = vec3(0);
    //alpha = 0.1;
    //if (alpha < 1.0) {
    //    transmissive = cone_trace(P, N, -V, tan_half_angle, transmissive_occlusion).rgb;
    //    transmissive *= 1.0 - alpha * albedo;
    //    direct_diffuse *= alpha;
    //    indirect_diffuse *= alpha;
    //}

    float specular_occlusion;
    vec3 direct_specular = dir_radiance * pow(N_dot_H, spec_exp);
    vec3 indirect_specular = cone_trace(P, N, R, tan_half_angle, specular_occlusion).rgb;

    vec3 result = vec3(0);
    result += albedo * (direct_diffuse + u_indirect_diffuse_scale * indirect_diffuse) * pow(diffuse_occlusion, u_ambient_occlusion_scale * 0.5);
    result += direct_specular * fresnel_direct + u_indirect_specular_scale * indirect_specular * fresnel_indirect;
    result += transmissive;

    return result;
}

void main() {
    float depth = texture(u_depth_texture, uv).r;
    if (depth == 1.0) discard;

    vec2 encoded_normal = texture(u_normal_texture, uv).rg;
    vec4 albedo_alpha = texture(u_color_alpha_texture, uv);
    vec4 f0_smoothness = texture(u_f0_smoothness_texture, uv);
    vec3 albedo = albedo_alpha.rgb;
    float alpha = albedo_alpha.a;
    vec3 f0 = f0_smoothness.rgb;
    float smoothness = f0_smoothness.a;

    vec3 world_position = depth_to_world_coord(uv, depth).xyz;
    vec3 world_normal = mat3(u_inv_view_mat) * decode_normal(encoded_normal);
    vec3 world_eye = u_world_space_camera;

    vec3 P = world_position;
    vec3 V = normalize(world_eye - world_position);
    vec3 N = world_normal;

    frag_color = vec4(shade(albedo, alpha, f0, smoothness, P, V, N), 1);
}