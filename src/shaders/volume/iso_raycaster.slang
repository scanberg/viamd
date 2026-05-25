// iso_raycaster.slang
// Scalar-volume isosurface renderer for the md_gpu path.

#define WG_X 8
#define WG_Y 8
#define MAX_ISO_COUNT 8

#define FLAG_USE_ABS_FIELD 1u
#define FLAG_TEMPORAL_JITTER 2u

static const float PI = 3.1415926535f;
static const float INV_PI = 1.0f / PI;
static const float EPSILON = 1.0e-6f;
static const float ERT_TRANSMITTANCE = 0.005f;

struct IsoParams {
    float4x4 inv_view_proj;
    float4x4 world_to_texture;
    float4x4 texture_to_world;

    float4 target;                 // width, height, inv_width, inv_height
    uint4  volume_dim_iso_count;   // dim.xyz, iso_count

    float4 clip_min_sampling;      // clip_min.xyz, sampling_rate
    float4 clip_max_extinction;    // clip_max.xyz, interior_extinction

    float4 env_roughness;          // env_radiance.rgb, roughness
    float4 dir_ior;                // dir_radiance.rgb, material IOR
    float4 light_exposure;         // light_dir_world.xyz, exposure
    float4 gamma_flags;            // gamma, flags, unused, unused

    float4 iso_values[2];          // values 0..3, values 4..7
    float4 iso_colors[MAX_ISO_COUNT];
};

struct RayState {
    float3 radiance;
    float transmittance;
    uint inside_mask;
};

struct IsoPushConstants {
    uint reserved;
};

[[vk::push_constant]]
ConstantBuffer<IsoPushConstants> push_constants;

[[vk::binding(0)]] StructuredBuffer<IsoParams> params_buf;

[[vk::binding(16)]] RWTexture2D<float4> out_color;
[[vk::binding(33)]] Texture3D<float> scalar_volume;
[[vk::binding(34)]] Texture2D<float> depth_buffer;
[[vk::binding(48)]] SamplerState scalar_sampler;

float3 transform_point(float4x4 mat, float3 pos) {
    float4 p = mul(mat, float4(pos, 1.0f));
    float w = (abs(p.w) > EPSILON) ? p.w : ((p.w < 0.0f) ? -EPSILON : EPSILON);
    return p.xyz / w;
}

float3 transform_vector(float4x4 mat, float3 vec) {
    return mul(mat, float4(vec, 0.0f)).xyz;
}

float get_iso_value(IsoParams params, uint idx) {
    return (idx < 4u) ? params.iso_values[0][idx] : params.iso_values[1][idx - 4u];
}

float3 get_clip_min(IsoParams params) {
    return params.clip_min_sampling.xyz;
}

float3 get_clip_max(IsoParams params) {
    return params.clip_max_extinction.xyz;
}

bool use_abs_field(IsoParams params) {
    return (uint(params.gamma_flags.y) & FLAG_USE_ABS_FIELD) != 0u;
}

float sample_density(IsoParams params, float3 tex_pos) {
    uint3 dim = max(params.volume_dim_iso_count.xyz, uint3(1, 1, 1));
    float3 sample_pos = (saturate(tex_pos) * float3(dim - uint3(1, 1, 1)) + 0.5f) / float3(dim);
    return scalar_volume.SampleLevel(scalar_sampler, sample_pos, 0.0f);
}

float sample_field(IsoParams params, float3 tex_pos) {
    float density = sample_density(params, tex_pos);
    return use_abs_field(params) ? abs(density) : density;
}

float3 sample_gradient_world(IsoParams params, float3 tex_pos) {
    uint3 dim = max(params.volume_dim_iso_count.xyz, uint3(2, 2, 2));
    float3 h = 1.0f / float3(dim);

    float dx = sample_field(params, tex_pos + float3(h.x, 0.0f, 0.0f)) - sample_field(params, tex_pos - float3(h.x, 0.0f, 0.0f));
    float dy = sample_field(params, tex_pos + float3(0.0f, h.y, 0.0f)) - sample_field(params, tex_pos - float3(0.0f, h.y, 0.0f));
    float dz = sample_field(params, tex_pos + float3(0.0f, 0.0f, h.z)) - sample_field(params, tex_pos - float3(0.0f, 0.0f, h.z));

    float3 axis_x = transform_vector(params.texture_to_world, float3(h.x, 0.0f, 0.0f));
    float3 axis_y = transform_vector(params.texture_to_world, float3(0.0f, h.y, 0.0f));
    float3 axis_z = transform_vector(params.texture_to_world, float3(0.0f, 0.0f, h.z));

    float len_x = length(axis_x);
    float len_y = length(axis_y);
    float len_z = length(axis_z);

    float3 grad = float3(0.0f);
    if (len_x > EPSILON) grad += (axis_x / len_x) * (dx / (2.0f * len_x));
    if (len_y > EPSILON) grad += (axis_y / len_y) * (dy / (2.0f * len_y));
    if (len_z > EPSILON) grad += (axis_z / len_z) * (dz / (2.0f * len_z));
    return grad;
}

bool is_inside_iso(IsoParams params, float density, uint iso_idx) {
    float value = get_iso_value(params, iso_idx);
    if (use_abs_field(params)) {
        return abs(density) >= abs(value);
    }
    return (value >= 0.0f) ? (density >= value) : (density <= value);
}

float iso_field_delta(IsoParams params, float density, uint iso_idx) {
    float value = get_iso_value(params, iso_idx);
    return use_abs_field(params) ? (abs(density) - abs(value)) : (density - value);
}

bool ray_vs_aabb(float3 ori, float3 dir, float3 box_min, float3 box_max, out float t_entry, out float t_exit) {
    float3 safe_dir = dir;
    if (abs(safe_dir.x) < EPSILON) safe_dir.x = (safe_dir.x < 0.0f) ? -EPSILON : EPSILON;
    if (abs(safe_dir.y) < EPSILON) safe_dir.y = (safe_dir.y < 0.0f) ? -EPSILON : EPSILON;
    if (abs(safe_dir.z) < EPSILON) safe_dir.z = (safe_dir.z < 0.0f) ? -EPSILON : EPSILON;
    float3 inv_dir = 1.0f / safe_dir;

    float3 t0 = (box_min - ori) * inv_dir;
    float3 t1 = (box_max - ori) * inv_dir;
    float3 tmin = min(t0, t1);
    float3 tmax = max(t0, t1);

    t_entry = max(max(tmin.x, tmin.y), tmin.z);
    t_exit = min(min(tmax.x, tmax.y), tmax.z);
    return t_exit >= max(t_entry, 0.0f);
}

float segment_length_world(IsoParams params, float3 p0_tex, float3 p1_tex) {
    float3 p0 = transform_point(params.texture_to_world, p0_tex);
    float3 p1 = transform_point(params.texture_to_world, p1_tex);
    return length(p1 - p0);
}

float3 aces_input_transform(float3 color) {
    return float3(
        dot(color, float3(0.59719f, 0.35458f, 0.04823f)),
        dot(color, float3(0.07600f, 0.90834f, 0.01566f)),
        dot(color, float3(0.02840f, 0.13383f, 0.83777f)));
}

float3 aces_output_transform(float3 color) {
    return float3(
        dot(color, float3( 1.60475f, -0.53108f, -0.07367f)),
        dot(color, float3(-0.10208f,  1.10813f, -0.00605f)),
        dot(color, float3(-0.00327f, -0.07276f,  1.07602f)));
}

float3 rrt_and_odt_fit(float3 value) {
    float3 a = value * (value + 0.0245786f) - 0.000090537f;
    float3 b = value * (0.983729f * value + 0.4329510f) + 0.238081f;
    return a / b;
}

float3 aces_fitted(float3 color) {
    color = aces_input_transform(color);
    color = rrt_and_odt_fit(color);
    color = aces_output_transform(color);
    return saturate(color);
}

float3 tonemap(IsoParams params, float3 hdr) {
    const float exposure_bias = 0.125f;
    const float white_point = 24.0f;

    float exposure = max(params.light_exposure.w, 0.0f);
    float gamma = max(params.gamma_flags.x, 0.01f);

    float3 exposed = max(hdr, float3(0.0f)) * (exposure_bias * exposure);
    float3 white_scale = max(aces_fitted(float3(white_point * exposure_bias * exposure)), float3(EPSILON));
    float3 color = saturate(aces_fitted(exposed) / white_scale);

    return pow(color, float3(1.0f / gamma));
}

float fresnel_schlick_dielectric(float cos_theta_i, float eta_i, float eta_t) {
    cos_theta_i = saturate(cos_theta_i);
    eta_i = max(eta_i, EPSILON);
    eta_t = max(eta_t, EPSILON);

    float r0 = (eta_i - eta_t) / (eta_i + eta_t);
    float f0 = r0 * r0;
    float one_minus_cos = 1.0f - cos_theta_i;
    float one_minus_cos2 = one_minus_cos * one_minus_cos;
    return saturate(f0 + (1.0f - f0) * one_minus_cos2 * one_minus_cos2 * one_minus_cos);
}

float distribution_ggx(float ndot_h, float roughness) {
    float a = roughness * roughness;
    float a2 = a * a;
    float ndot_h2 = ndot_h * ndot_h;
    float denom = ndot_h2 * (a2 - 1.0f) + 1.0f;
    return a2 / max(PI * denom * denom, EPSILON);
}

float geometry_schlick_ggx(float ndot_v, float roughness) {
    float r = roughness + 1.0f;
    float k = (r * r) * 0.125f;
    return ndot_v / max(ndot_v * (1.0f - k) + k, EPSILON);
}

float geometry_smith(float ndot_v, float ndot_l, float roughness) {
    return geometry_schlick_ggx(ndot_v, roughness) * geometry_schlick_ggx(ndot_l, roughness);
}

float3 shade_surface(IsoParams params, float3 pos_tex, float3 ray_dir_world, float4 base_color, bool entering, out float view_fresnel) {
    float3 V = normalize(-ray_dir_world);
    float3 L = normalize(params.light_exposure.xyz);

    float3 N = sample_gradient_world(params, pos_tex);
    float n_len = length(N);
    N = (n_len > EPSILON) ? (N / n_len) : V;
    if (dot(N, V) < 0.0f) {
        N = -N;
    }

    float roughness = clamp(params.env_roughness.w, 0.04f, 1.0f);
    float material_ior = max(params.dir_ior.w, 1.0f);
    float eta_i = entering ? 1.0f : material_ior;
    float eta_t = entering ? material_ior : 1.0f;
    float3 albedo = max(base_color.rgb, float3(0.0f));
    float base_alpha = saturate(base_color.a);

    float NdotL = saturate(dot(N, L));
    float NdotV = saturate(dot(N, V));
    float3 hdr = float3(0.0f);

    if (NdotL > 0.0f) {
        float3 H = normalize(L + V);
        float NdotH = saturate(dot(N, H));
        float HdotV = saturate(dot(H, V));

        float F = fresnel_schlick_dielectric(HdotV, eta_i, eta_t);
        float D = distribution_ggx(NdotH, roughness);
        float G = geometry_smith(NdotV, NdotL, roughness);

        float specular = (D * G * F) / max(4.0f * NdotV * NdotL, 1.0e-4f);
        float diffuse_weight = base_alpha * (1.0f - F);
        hdr += (diffuse_weight * albedo * INV_PI + float3(specular)) * params.dir_ior.rgb * NdotL;
    }

    float env_fresnel = fresnel_schlick_dielectric(NdotV, eta_i, eta_t);
    float env_diffuse_weight = base_alpha * (1.0f - env_fresnel);
    hdr += (env_diffuse_weight * albedo + float3(env_fresnel)) * params.env_roughness.rgb;

    view_fresnel = env_fresnel;
    return hdr;
}

void initialize_inside_mask(inout RayState state, IsoParams params, float density) {
    state.inside_mask = 0u;
    uint iso_count = min(params.volume_dim_iso_count.w, uint(MAX_ISO_COUNT));
    for (uint i = 0u; i < iso_count; ++i) {
        if (is_inside_iso(params, density, i)) {
            state.inside_mask |= (1u << i);
        }
    }
}

void accumulate_absorption(inout RayState state, IsoParams params, float3 p0_tex, float3 p1_tex) {
    if (state.inside_mask == 0u || params.clip_max_extinction.w <= 0.0f) {
        return;
    }

    float len_world = segment_length_world(params, p0_tex, p1_tex);
    if (len_world <= EPSILON) {
        return;
    }

    float transmittance = exp(-params.clip_max_extinction.w * len_world);
    state.transmittance *= saturate(transmittance);
}

void accumulate_surface(inout RayState state, IsoParams params, float3 pos_tex, float3 ray_dir_world, float4 base_color, bool entering) {
    if (state.transmittance <= ERT_TRANSMITTANCE) {
        return;
    }

    float view_fresnel = 1.0f;
    float3 hdr = shade_surface(params, pos_tex, ray_dir_world, base_color, entering, view_fresnel);
    float base_alpha = saturate(base_color.a);
    float transmission = saturate((1.0f - base_alpha) * (1.0f - view_fresnel));
    float opacity = 1.0f - transmission;

    state.radiance += state.transmittance * hdr;
    state.transmittance *= (1.0f - opacity);
}

int find_iso_hits(IsoParams params, float d0, float d1, out float hit_frac[MAX_ISO_COUNT], out uint hit_idx[MAX_ISO_COUNT]) {
    int hit_count = 0;
    uint iso_count = min(params.volume_dim_iso_count.w, uint(MAX_ISO_COUNT));

    for (uint i = 0u; i < iso_count; ++i) {
        float f0 = iso_field_delta(params, d0, i);
        float f1 = iso_field_delta(params, d1, i);

        if ((f0 < 0.0f && f1 >= 0.0f) || (f0 >= 0.0f && f1 < 0.0f)) {
            float denom = f0 - f1;
            if (abs(denom) <= EPSILON) {
                continue;
            }

            float frac = saturate(f0 / denom);
            if (frac <= EPSILON || frac >= 1.0f - EPSILON) {
                continue;
            }

            hit_frac[hit_count] = frac;
            hit_idx[hit_count] = i;
            ++hit_count;
        }
    }

    for (int i = 1; i < hit_count; ++i) {
        float frac = hit_frac[i];
        uint index = hit_idx[i];
        int j = i - 1;

        while (j >= 0 && hit_frac[j] > frac) {
            hit_frac[j + 1] = hit_frac[j];
            hit_idx[j + 1] = hit_idx[j];
            --j;
        }

        hit_frac[j + 1] = frac;
        hit_idx[j + 1] = index;
    }

    return hit_count;
}

void process_segment(inout RayState state, IsoParams params, float3 p0_tex, float3 p1_tex, float d0, float d1, float3 ray_dir_world) {
    if (distance(p0_tex, p1_tex) <= EPSILON || state.transmittance <= ERT_TRANSMITTANCE) {
        return;
    }

    float hit_frac[MAX_ISO_COUNT];
    uint hit_idx[MAX_ISO_COUNT];
    int hit_count = find_iso_hits(params, d0, d1, hit_frac, hit_idx);

    float last_frac = 0.0f;
    for (int h = 0; h < hit_count; ++h) {
        float frac = hit_frac[h];
        uint iso_idx = hit_idx[h];

        float3 seg_start = lerp(p0_tex, p1_tex, last_frac);
        float3 seg_end = lerp(p0_tex, p1_tex, frac);
        accumulate_absorption(state, params, seg_start, seg_end);

        float4 iso_color = params.iso_colors[iso_idx];
        bool entering = (state.inside_mask & (1u << iso_idx)) == 0u;
        accumulate_surface(state, params, seg_end, ray_dir_world, iso_color, entering);

        state.inside_mask ^= (1u << iso_idx);
        last_frac = frac;
    }

    accumulate_absorption(state, params, lerp(p0_tex, p1_tex, last_frac), p1_tex);
}

float rand_2d(float2 p) {
    return frac(sin(dot(p, float2(12.9898f, 78.233f))) * 43758.5453f);
}

[shader("compute")]
[numthreads(WG_X, WG_Y, 1)]
void main(uint3 pixel : SV_DispatchThreadID) {
    IsoParams params = params_buf[0];
    uint2 target_dim = uint2(uint(params.target.x), uint(params.target.y));
    if (any(pixel.xy >= target_dim)) {
        return;
    }

    if (push_constants.reserved != 0u) {
        out_color[pixel.xy] = float4(0.0f);
        return;
    }

    if (params.volume_dim_iso_count.w == 0u) {
        out_color[pixel.xy] = float4(0.0f);
        return;
    }

    float2 uv = (float2(pixel.xy) + 0.5f) * params.target.zw;
    float2 ndc = uv * 2.0f - 1.0f;

    float3 near_world = transform_point(params.inv_view_proj, float3(ndc, -1.0f));
    float3 far_world = transform_point(params.inv_view_proj, float3(ndc, 1.0f));
    float3 near_tex = transform_point(params.world_to_texture, near_world);
    float3 far_tex = transform_point(params.world_to_texture, far_world);

    float3 ray_tex = far_tex - near_tex;
    float ray_len_tex = length(ray_tex);
    if (ray_len_tex <= EPSILON) {
        out_color[pixel.xy] = float4(0.0f);
        return;
    }

    float3 ray_dir_tex = ray_tex / ray_len_tex;
    float t_entry = 0.0f;
    float t_exit = 0.0f;
    if (!ray_vs_aabb(near_tex, ray_dir_tex, get_clip_min(params), get_clip_max(params), t_entry, t_exit)) {
        out_color[pixel.xy] = float4(0.0f);
        return;
    }

    t_entry = max(t_entry, 0.0f);
    t_exit = min(t_exit, ray_len_tex);

    float depth = depth_buffer.Load(int3(pixel.xy, 0));
    if (depth < 1.0f) {
        float3 stop_world = transform_point(params.inv_view_proj, float3(ndc, depth * 2.0f - 1.0f));
        float3 stop_tex = transform_point(params.world_to_texture, stop_world);
        t_exit = min(t_exit, dot(stop_tex - near_tex, ray_dir_tex));
    }

    float ray_len = t_exit - t_entry;
    if (ray_len <= EPSILON) {
        out_color[pixel.xy] = float4(0.0f);
        return;
    }

    float3 entry_tex = near_tex + ray_dir_tex * t_entry;
    float3 exit_tex = near_tex + ray_dir_tex * t_exit;
    float3 ray_dir_world = normalize(transform_point(params.texture_to_world, exit_tex) - transform_point(params.texture_to_world, entry_tex));

    float3 voxel_ray = ray_dir_tex * ray_len * float3(params.volume_dim_iso_count.xyz);
    float num_samples = max(1.0f, params.clip_min_sampling.w * length(voxel_ray));
    float step_tex = max(ray_len / ceil(num_samples), 1.0e-5f);

    RayState state;
    state.radiance = float3(0.0f);
    state.transmittance = 1.0f;
    state.inside_mask = 0u;

    float3 prev_pos = entry_tex;
    float prev_density = sample_density(params, prev_pos);
    initialize_inside_mask(state, params, prev_density);

    bool temporal = (uint(params.gamma_flags.y) & FLAG_TEMPORAL_JITTER) != 0u;
    float t = temporal ? (rand_2d(float2(pixel.xy)) * step_tex) : step_tex;

    while (t < ray_len && state.transmittance > ERT_TRANSMITTANCE) {
        float3 curr_pos = entry_tex + ray_dir_tex * t;
        float curr_density = sample_density(params, curr_pos);

        process_segment(state, params, prev_pos, curr_pos, prev_density, curr_density, ray_dir_world);

        prev_pos = curr_pos;
        prev_density = curr_density;
        t += step_tex;
    }

    if (state.transmittance > ERT_TRANSMITTANCE) {
        float exit_density = sample_density(params, exit_tex);
        process_segment(state, params, prev_pos, exit_tex, prev_density, exit_density, ray_dir_world);
    }

    float opacity = saturate(1.0f - state.transmittance);
    out_color[pixel.xy] = float4(tonemap(params, state.radiance), opacity);
}