// Modified version of
// https://github.com/Cigg/Voxel-Cone-Tracing

#version 330 core

in vec2 uv;
out vec4 frag_color;

// Textures
uniform sampler2D u_depth_texture;
uniform sampler2D u_normal_texture;

// Voxel stuff
uniform sampler3D u_voxel_texture;
uniform vec3 u_voxel_grid_world_size;
uniform vec3 u_voxel_grid_world_min;
uniform ivec3 u_voxel_dimensions;
uniform float u_voxel_extent;

// Scaling factors
uniform float u_directional_occlusion_scale = 1.f;
uniform float u_cone_angle = 0.08;

// View parameters
uniform mat4 u_inv_view_mat;
uniform mat4 u_inv_view_proj_mat;
uniform vec3 u_world_space_camera;

const float MAX_DIST = 1000.0;

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

float sample_volume(vec3 world_position, float lod) {
    vec3 tc = (world_position - u_voxel_grid_world_min) / (u_voxel_grid_world_size);
    return textureLod(u_voxel_texture, tc, lod).r;
}

float cone_trace(in vec3 world_position, in vec3 world_normal, in vec3 direction, in float tan_half_angle) {
    float lod = 0.0;
    float occlusion = 0.0;

    float dist = u_voxel_extent * 0.5;
    vec3 start_pos = world_position + world_normal * dist;

    while(dist < MAX_DIST) {
        // smallest sample diameter possible is the voxel size
        float diameter = max(u_voxel_extent, 2.0 * tan_half_angle * dist);
        float lod_level = log2(diameter / u_voxel_extent);
        float occ = sample_volume(start_pos + dist * direction, lod_level);
        //if (voxel_color.a < 0) break;

        occlusion += (occ) / (1.0 + 0.03 * diameter);

        dist += diameter; // faster but misses more voxels
    }

    return occlusion;
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

float compute_directional_occlusion(vec3 P, vec3 V, vec3 N) {
    // 60 degree cones -> tan(30) = 0.577
    // 90 degree cones -> tan(45) = 1.0
    float tan_half_angle = 1.0;
    return cone_trace(P, N, N, tan_half_angle);
}

void main() {
    float depth = texture(u_depth_texture, uv).r;
    if (depth == 1.0) discard;
    vec2 encoded_normal = texture(u_normal_texture, uv).rg;

    vec3 world_position = depth_to_world_coord(uv, depth).xyz;
    vec3 world_normal = mat3(u_inv_view_mat) * decode_normal(encoded_normal);
    vec3 world_eye = u_world_space_camera;

    vec3 P = world_position;
    vec3 V = normalize(world_eye - world_position);
    vec3 N = world_normal;

    float res = compute_directional_occlusion(P, V, N);
    res = res * u_directional_occlusion_scale;

    frag_color = vec4(vec3(res), 1);
}