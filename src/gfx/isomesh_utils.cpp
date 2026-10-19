#include "isomesh_utils.h"

#include <gfx/gl.h>
#include <gfx/gl_utils.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_hash.h>
#include <core/md_log.h>

#include <math.h>
#include <float.h>

#include <shaders.inl>

#include <gfx/profiler.h>

namespace isomesh {

static constexpr str_t v_shader_src_fs_triangle = STR_LIT(R"(
#version 410 core
void main() {
    uint idx = uint(gl_VertexID) % 3U;
    gl_Position = vec4(float(idx & 1U) * 4.0 - 1.0, float((idx >> 1U) & 1U) * 4.0 - 1.0, 0.0, 1.0);
}
)");

// The mean of the frames over the target (isomesh_accumulate.frag): the radiance scaled to the exact coverage
// of the surfaces (Enderton's alpha correction), and the coverage of the surfaces and what they enclose
static constexpr str_t f_shader_src_composite = STR_LIT(R"(
#version 410 core
uniform sampler2D u_tex_color;
uniform sampler2D u_tex_t;
layout(location = 0) out vec4 out_color;
void main() {
    ivec2 px = ivec2(gl_FragCoord.xy);
    vec4 c = texelFetch(u_tex_color, px, 0);
    vec2 t = texelFetch(u_tex_t, px, 0).rg;
    float coverage = 1.0 - t.r;
    vec3 L = (c.a > 1.0e-4) ? c.rgb * (coverage / c.a) : c.rgb;
    out_color = vec4(L, 1.0 - t.g);
}
)");

// Extinction per world unit (Ångström) of unit optical density, as the raycaster's (ISO_OPTICAL_SCALE)
static constexpr float OPTICAL_SCALE = 150.0f;
static constexpr int   MAX_SAMPLES   = 4;
static constexpr int   BLUE_NOISE_SIZE = 64;    // BLUE_NOISE_SIZE in isomesh.frag

enum Pass {     // isomesh.frag
    Pass_Transmittance = 0,
    Pass_Peel0         = 1,
    Pass_Peel1         = 2,
    Pass_Depth         = 3,
    Pass_Accumulate    = 4,
};

// isomesh.frag, built twice: PASS_TRANSMITTANCE at pixel rate (ISOMESH_TRANSMITTANCE) and the other passes per
// sample. A uniform one of them has not got is at -1, which glUniform ignores, so both are set up alike.
struct MeshProgram {
    GLuint program = 0;
    GLint model, view, proj;
    GLint color, coloring, world_to_volume, tex_color_volume, tex_field, tex_field_colormap, field_beg, field_inv_ext;
    GLint tex_depth, use_depth, pass, frame, samples, surface_id;
    GLint tex_tod, tau, od_ref, inv_proj, inv_res, ortho, tex_blue_noise;
    GLint tex_peel_depth0, tex_peel_depth1, tex_peel_opacity;
    GLint light_dir, dir_radiance, sky_radiance, ground_radiance, roughness, F0;
    GLint basic, basic_ambient, basic_specular, linear_colors, light_strength;
    GLint isoline_count, isoline_values, isoline_hover, isoline_hover_value, isoline_width;
};

static struct {
    MeshProgram mesh;       // PASS_PEEL0 .. PASS_ACCUMULATE, per sample
    MeshProgram mesh_tod;   // PASS_TRANSMITTANCE, at pixel rate
    MeshProgram mesh_nfd;   // as mesh, without writing gl_FragDepth (Experiment_*NoFragDepth)
    struct {
        GLuint program = 0;
        GLint tex_color, tex_tod, samples, tex_history_color, tex_history_t, weight;
    } accumulate;
    struct {
        GLuint program = 0;
        GLint tex_color, tex_t;
    } composite;
    GLuint tex_blue_noise = 0;
    GLuint vao = 0;     // empty, for the fullscreen triangle
    GLuint fbo = 0;     // onto the caller's colour target
    bool   ok  = false;
} gl;

// A blue noise threshold map, n x n (a power of two), toroidal: the rank of each pixel in a void and cluster
// ordering (Ulichney 1993), as (rank + 0.5) / n^2. Each threshold level is a pattern of evenly spread pixels, so any
// threshold keeps its pixels apart, and so does the map shifted by any offset.
static void blue_noise_generate(float* out, int n) {
    const int N = n * n;
    const float sigma = 1.9f;
    // The toroidal Gaussian, by offset
    float* kernel = (float*)md_alloc(md_get_heap_allocator(), sizeof(float) * N);
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            const int dx = MIN(x, n - x), dy = MIN(y, n - y);
            kernel[y * n + x] = expf(-(float)(dx * dx + dy * dy) / (2.0f * sigma * sigma));
        }
    }
    uint8_t* on     = (uint8_t*)md_alloc(md_get_heap_allocator(), N);
    float*   energy = (float*)md_alloc(md_get_heap_allocator(), sizeof(float) * N);
    int*     rank   = (int*)md_alloc(md_get_heap_allocator(), sizeof(int) * N);
    MEMSET(on, 0, N);
    MEMSET(energy, 0, sizeof(float) * N);

    auto splat = [&](int i, float sign) {
        const int ix = i % n, iy = i / n;
        for (int y = 0; y < n; ++y) {
            const int ky = ((y - iy) & (n - 1)) * n;
            for (int x = 0; x < n; ++x) {
                energy[y * n + x] += sign * kernel[ky + ((x - ix) & (n - 1))];
            }
        }
    };
    auto tightest_cluster = [&]() {
        int best = -1;
        for (int i = 0; i < N; ++i) if (on[i] && (best < 0 || energy[i] > energy[best])) best = i;
        return best;
    };
    auto largest_void = [&]() {
        int best = -1;
        for (int i = 0; i < N; ++i) if (!on[i] && (best < 0 || energy[i] < energy[best])) best = i;
        return best;
    };

    // The initial pattern: a tenth of the pixels at random, relaxed until moving the tightest cluster to the
    // largest void changes nothing
    uint32_t rng = 0x9E3779B9u;
    const int num_initial = N / 10;
    for (int k = 0; k < num_initial;) {
        rng = rng * 1664525u + 1013904223u;
        const int i = (int)((rng >> 8) % (uint32_t)N);
        if (on[i]) continue;
        on[i] = 1;
        splat(i, 1.0f);
        ++k;
    }
    for (int iter = 0; iter < 4 * N; ++iter) {
        const int c = tightest_cluster();
        on[c] = 0;
        splat(c, -1.0f);
        const int v = largest_void();
        on[v] = 1;
        splat(v, 1.0f);
        if (v == c) break;
    }
    // Ranks below the initial: take the tightest clusters away, on a copy of the state
    uint8_t* on0     = (uint8_t*)md_alloc(md_get_heap_allocator(), N);
    float*   energy0 = (float*)md_alloc(md_get_heap_allocator(), sizeof(float) * N);
    MEMCPY(on0, on, N);
    MEMCPY(energy0, energy, sizeof(float) * N);
    for (int r = num_initial - 1; r >= 0; --r) {
        const int c = tightest_cluster();
        on[c] = 0;
        splat(c, -1.0f);
        rank[c] = r;
    }
    MEMCPY(on, on0, N);
    MEMCPY(energy, energy0, sizeof(float) * N);
    // and above: fill the largest voids
    for (int r = num_initial; r < N; ++r) {
        const int v = largest_void();
        on[v] = 1;
        splat(v, 1.0f);
        rank[v] = r;
    }
    for (int i = 0; i < N; ++i) out[i] = ((float)rank[i] + 0.5f) / (float)N;

    md_free(md_get_heap_allocator(), kernel, sizeof(float) * N);
    md_free(md_get_heap_allocator(), on, N);
    md_free(md_get_heap_allocator(), energy, sizeof(float) * N);
    md_free(md_get_heap_allocator(), rank, sizeof(int) * N);
    md_free(md_get_heap_allocator(), on0, N);
    md_free(md_get_heap_allocator(), energy0, sizeof(float) * N);
}

static GLuint link(str_t v_src, str_t f_src, const char* name, str_t f_defines = {}) {
    GLuint v = gl::compile_shader_from_source(v_src, GL_VERTEX_SHADER);
    GLuint f = gl::compile_shader_from_source(f_src, GL_FRAGMENT_SHADER, f_defines);
    GLuint p = 0;
    if (v && f) {
        p = glCreateProgram();
        const GLuint shaders[] = {v, f};
        if (!gl::attach_link_detach(p, shaders, 2)) {
            glDeleteProgram(p);
            p = 0;
        }
    }
    if (v) glDeleteShader(v);
    if (f) glDeleteShader(f);
    if (!p) MD_LOG_ERROR("Isosurface meshes: the %s program failed to build", name);
    return p;
}

void initialize() {
    const str_t vert = {(const char*)isomesh_vert, isomesh_vert_size};
    const str_t frag = {(const char*)isomesh_frag, isomesh_frag_size};
    const str_t acc  = {(const char*)isomesh_accumulate_frag, isomesh_accumulate_frag_size};

    shutdown();

    gl.mesh.program       = link(vert, frag, "mesh");
    gl.mesh_tod.program   = link(vert, frag, "transmittance", STR_LIT("#define ISOMESH_TRANSMITTANCE"));
    gl.mesh_nfd.program   = link(vert, frag, "mesh without gl_FragDepth", STR_LIT("#define ISOMESH_NO_FRAGDEPTH"));
    gl.accumulate.program = link(v_shader_src_fs_triangle, acc, "accumulate");
    gl.composite.program  = link(v_shader_src_fs_triangle, f_shader_src_composite, "composite");
    gl.ok = gl.mesh.program && gl.mesh_tod.program && gl.mesh_nfd.program && gl.accumulate.program && gl.composite.program;
    if (!gl.ok) return;

#define LOC(prog, name) glGetUniformLocation(prog, name)
    MeshProgram* const mesh_programs[] = {&gl.mesh, &gl.mesh_tod, &gl.mesh_nfd};
    for (MeshProgram* mp : mesh_programs) {
        auto& m = *mp;
        const GLuint p = m.program;
        m.model = LOC(p, "u_model"); m.view = LOC(p, "u_view"); m.proj = LOC(p, "u_proj");
        m.color = LOC(p, "u_color"); m.coloring = LOC(p, "u_coloring"); m.world_to_volume = LOC(p, "u_world_to_volume");
        m.tex_color_volume = LOC(p, "u_tex_color_volume"); m.tex_field = LOC(p, "u_tex_field"); m.tex_field_colormap = LOC(p, "u_tex_field_colormap");
        m.field_beg = LOC(p, "u_field_beg"); m.field_inv_ext = LOC(p, "u_field_inv_ext");
        m.tex_depth = LOC(p, "u_tex_depth"); m.use_depth = LOC(p, "u_use_depth"); m.pass = LOC(p, "u_pass"); m.frame = LOC(p, "u_frame");
        m.samples = LOC(p, "u_samples"); m.surface_id = LOC(p, "u_surface_id");
        m.tex_tod = LOC(p, "u_tex_tod"); m.tau = LOC(p, "u_tau"); m.od_ref = LOC(p, "u_od_ref");
        m.inv_proj = LOC(p, "u_inv_proj"); m.inv_res = LOC(p, "u_inv_res"); m.ortho = LOC(p, "u_ortho");
        m.tex_blue_noise = LOC(p, "u_tex_blue_noise");
        m.tex_peel_depth0 = LOC(p, "u_tex_peel_depth0"); m.tex_peel_depth1 = LOC(p, "u_tex_peel_depth1");
        m.tex_peel_opacity = LOC(p, "u_tex_peel_opacity");
        m.light_dir = LOC(p, "u_light_dir"); m.dir_radiance = LOC(p, "u_dir_radiance"); m.sky_radiance = LOC(p, "u_sky_radiance");
        m.ground_radiance = LOC(p, "u_ground_radiance"); m.roughness = LOC(p, "u_roughness"); m.F0 = LOC(p, "u_F0");
        m.basic = LOC(p, "u_basic"); m.basic_ambient = LOC(p, "u_basic_ambient"); m.basic_specular = LOC(p, "u_basic_specular");
        m.linear_colors = LOC(p, "u_linear_colors"); m.light_strength = LOC(p, "u_light_strength");
        m.isoline_count = LOC(p, "u_isoline_count"); m.isoline_values = LOC(p, "u_isoline_values");
        m.isoline_hover = LOC(p, "u_isoline_hover"); m.isoline_hover_value = LOC(p, "u_isoline_hover_value");
        m.isoline_width = LOC(p, "u_isoline_width");
    }
    {
        auto& a = gl.accumulate;
        const GLuint p = a.program;
        a.tex_color = LOC(p, "u_tex_color"); a.tex_tod = LOC(p, "u_tex_tod"); a.samples = LOC(p, "u_samples");
        a.tex_history_color = LOC(p, "u_tex_history_color"); a.tex_history_t = LOC(p, "u_tex_history_t"); a.weight = LOC(p, "u_weight");
    }
    gl.composite.tex_color = LOC(gl.composite.program, "u_tex_color");
    gl.composite.tex_t     = LOC(gl.composite.program, "u_tex_t");
#undef LOC

    glGenVertexArrays(1, &gl.vao);
    glGenFramebuffers(1, &gl.fbo);

    {
        const int n = BLUE_NOISE_SIZE;
        float* bn = (float*)md_alloc(md_get_heap_allocator(), sizeof(float) * n * n);
        blue_noise_generate(bn, n);
        gl::init_texture_2D(&gl.tex_blue_noise, n, n, GL_R32F);
        gl::set_texture_2D_data(gl.tex_blue_noise, 0, bn, GL_R32F);
        glBindTexture(GL_TEXTURE_2D, gl.tex_blue_noise);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glBindTexture(GL_TEXTURE_2D, 0);
        md_free(md_get_heap_allocator(), bn, sizeof(float) * n * n);
    }
}

void shutdown() {
    GLuint* programs[] = {&gl.mesh.program, &gl.mesh_tod.program, &gl.mesh_nfd.program, &gl.accumulate.program, &gl.composite.program};
    for (GLuint* p : programs) {
        if (*p) glDeleteProgram(*p);
        *p = 0;
    }
    if (gl.vao) glDeleteVertexArrays(1, &gl.vao);
    if (gl.fbo) glDeleteFramebuffers(1, &gl.fbo);
    if (gl.tex_blue_noise) glDeleteTextures(1, &gl.tex_blue_noise);
    gl.vao = 0;
    gl.fbo = 0;
    gl.tex_blue_noise = 0;
    gl.ok = false;
}

bool mesh_upload(Mesh* mesh, const float* position, const float* normal, size_t stride, uint32_t vertex_count, const uint32_t* index, uint32_t index_count) {
    ASSERT(mesh);
    if (vertex_count > 0 && (!position || !normal || stride < 3)) return false;
    if (index_count > 0 && !index) return false;

    if (!mesh->vao) glGenVertexArrays(1, &mesh->vao);
    if (!mesh->vbo) glGenBuffers(1, &mesh->vbo);
    if (!mesh->ibo) glGenBuffers(1, &mesh->ibo);

    float* interleaved = (float*)md_alloc(md_get_heap_allocator(), sizeof(float) * 6 * MAX(vertex_count, 1u));
    vec3_t lo = {FLT_MAX, FLT_MAX, FLT_MAX}, hi = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (uint32_t i = 0; i < vertex_count; ++i) {
        for (int k = 0; k < 3; ++k) {
            const float x = position[stride * i + k];
            interleaved[6 * i + k]     = x;
            interleaved[6 * i + 3 + k] = normal[stride * i + k];
            lo.elem[k] = MIN(lo.elem[k], x);
            hi.elem[k] = MAX(hi.elem[k], x);
        }
    }
    mesh->bounds_min = vertex_count ? lo : vec3_t{0, 0, 0};
    mesh->bounds_max = vertex_count ? hi : vec3_t{0, 0, 0};

    glBindVertexArray(mesh->vao);
    glBindBuffer(GL_ARRAY_BUFFER, mesh->vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * vertex_count, interleaved, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (const void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 6, (const void*)(sizeof(float) * 3));
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * index_count, index, GL_STATIC_DRAW);
    // The values were of the vertices before
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    md_free(md_get_heap_allocator(), interleaved, sizeof(float) * 6 * MAX(vertex_count, 1u));

    mesh->vertex_count = vertex_count;
    mesh->index_count  = index_count;
    mesh->version += 1;
    return true;
}

bool mesh_upload_values(Mesh* mesh, const float* values, const float* gradients, uint32_t count) {
    ASSERT(mesh);
    if (!mesh->vao) return false;
    if (values && count != mesh->vertex_count) return false;
    glBindVertexArray(mesh->vao);
    if (!values) {
        glDisableVertexAttribArray(2);
        glDisableVertexAttribArray(3);
        glVertexAttrib1f(2, 0.0f);
        glVertexAttrib3f(3, 0.0f, 0.0f, 0.0f);
    } else {
        // The value and its gradient, 4 floats per vertex
        float* data = (float*)md_alloc(md_get_heap_allocator(), sizeof(float) * 4 * MAX(count, 1u));
        for (uint32_t i = 0; i < count; ++i) {
            data[4 * i + 0] = values[i];
            data[4 * i + 1] = gradients ? gradients[3 * i + 0] : 0.0f;
            data[4 * i + 2] = gradients ? gradients[3 * i + 1] : 0.0f;
            data[4 * i + 3] = gradients ? gradients[3 * i + 2] : 0.0f;
        }
        if (!mesh->vbo_value) glGenBuffers(1, &mesh->vbo_value);
        glBindBuffer(GL_ARRAY_BUFFER, mesh->vbo_value);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4 * count, data, GL_STATIC_DRAW);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float) * 4, (const void*)0);
        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 4, (const void*)(sizeof(float)));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        md_free(md_get_heap_allocator(), data, sizeof(float) * 4 * MAX(count, 1u));
    }
    glBindVertexArray(0);
    mesh->version += 1;
    return true;
}

void mesh_field_gradient(float* out, const float* xyz, const float* values, uint32_t vertex_count, const uint32_t* index, uint32_t index_count) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(float) * 3 * vertex_count);
    if (!xyz || !values || !index || vertex_count == 0) return;
    float* weight = (float*)md_alloc(md_get_heap_allocator(), sizeof(float) * vertex_count);
    MEMSET(weight, 0, sizeof(float) * vertex_count);

    for (uint32_t t = 0; t + 2 < index_count; t += 3) {
        const uint32_t i0 = index[t], i1 = index[t + 1], i2 = index[t + 2];
        if (i0 >= vertex_count || i1 >= vertex_count || i2 >= vertex_count) continue;
        const vec3_t p0 = {xyz[3 * i0], xyz[3 * i0 + 1], xyz[3 * i0 + 2]};
        const vec3_t e1 = vec3_t{xyz[3 * i1], xyz[3 * i1 + 1], xyz[3 * i1 + 2]} - p0;
        const vec3_t e2 = vec3_t{xyz[3 * i2], xyz[3 * i2 + 1], xyz[3 * i2 + 2]} - p0;
        const vec3_t n  = vec3_cross(e1, e2);
        const float  a2 = vec3_dot(n, n);     // (twice the area)^2
        if (!(a2 > 0.0f)) continue;
        // The gradient g of the linear interpolation in the triangle's plane: g . e1 = f1 - f0, g . e2 = f2 - f0
        const float  df1 = values[i1] - values[i0];
        const float  df2 = values[i2] - values[i0];
        const vec3_t g   = (vec3_cross(e2, n) * df1 + vec3_cross(n, e1) * df2) * (1.0f / a2);
        const float  w   = sqrtf(a2);
        const uint32_t tri[3] = {i0, i1, i2};
        for (uint32_t v : tri) {
            out[3 * v + 0] += g.x * w;
            out[3 * v + 1] += g.y * w;
            out[3 * v + 2] += g.z * w;
            weight[v] += w;
        }
    }
    for (uint32_t v = 0; v < vertex_count; ++v) {
        if (weight[v] > 0.0f) {
            const float inv = 1.0f / weight[v];
            out[3 * v + 0] *= inv;
            out[3 * v + 1] *= inv;
            out[3 * v + 2] *= inv;
        }
    }
    md_free(md_get_heap_allocator(), weight, sizeof(float) * vertex_count);
}

void mesh_free(Mesh* mesh) {
    ASSERT(mesh);
    if (mesh->vao) glDeleteVertexArrays(1, &mesh->vao);
    if (mesh->vbo) glDeleteBuffers(1, &mesh->vbo);
    if (mesh->ibo) glDeleteBuffers(1, &mesh->ibo);
    if (mesh->vbo_value) glDeleteBuffers(1, &mesh->vbo_value);
    const uint64_t version = mesh->version;
    *mesh = {};
    mesh->version = version + 1;
}

void target_free(Target* t) {
    ASSERT(t);
    GLuint* fbos[] = {&t->fbo_color, &t->fbo_tod, &t->fbo, &t->fbo_peel[0], &t->fbo_peel[1]};
    for (GLuint* f : fbos) {
        if (*f) glDeleteFramebuffers(1, f);
    }
    GLuint* textures[] = {&t->tex_ms_color, &t->tex_ms_tod, &t->tex_ms_depth, &t->tex_ms_peel_depth[0], &t->tex_ms_peel_depth[1],
                          &t->tex_ms_peel_opacity, &t->tex_history_color[0], &t->tex_history_color[1],
                          &t->tex_history_t[0], &t->tex_history_t[1]};
    for (GLuint* tex : textures) {
        if (*tex) glDeleteTextures(1, tex);
    }
    *t = {};
}

static GLuint texture_ms(GLsizei samples, GLenum format, uint32_t width, uint32_t height) {
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, tex);
    glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, samples, format, (GLsizei)width, (GLsizei)height, GL_TRUE);
    glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
    return tex;
}

static bool fbo_complete(GLuint fbo, const char* name) {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    const GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE) {
        MD_LOG_ERROR("Isosurface meshes: the %s framebuffer is incomplete (0x%04X)", name, (unsigned int)status);
        return false;
    }
    return true;
}

static Experiment g_experiment = Experiment_None;

void set_experiment(Experiment e) { g_experiment = (e >= 0 && e < Experiment_Count) ? e : Experiment_None; }
Experiment experiment() { return g_experiment; }

const char* experiment_name(Experiment e) {
    switch (e) {
        case Experiment_None:                 return "baseline";
        case Experiment_NoFragDepth:          return "no-fragdepth";
        case Experiment_OneSample:            return "1-sample";
        case Experiment_OneSampleNoFragDepth: return "1-sample-no-fragdepth";
        case Experiment_GeometryOnly:         return "geometry-only";
        default:                              return "?";
    }
}

static bool target_init(Target* t, uint32_t width, uint32_t height) {
    const int cap = (g_experiment == Experiment_OneSample || g_experiment == Experiment_OneSampleNoFragDepth) ? 1 : MAX_SAMPLES;
    if (t->fbo && t->width == width && t->height == height && t->sample_cap == (uint32_t)cap) return true;
    target_free(t);

    GLint max_color = 1, max_depth = 1;
    glGetIntegerv(GL_MAX_COLOR_TEXTURE_SAMPLES, &max_color);
    glGetIntegerv(GL_MAX_DEPTH_TEXTURE_SAMPLES, &max_depth);
    // As many samples as the budget allows: colour 8, transmittance 4, three depths 12, peeled opacities 4 bytes
    const size_t bytes_per_sample = 28;
    GLsizei samples = CLAMP(MIN(max_color, max_depth), 1, cap);
    while (samples > 1 && (size_t)width * height * samples * bytes_per_sample > MEMORY_BUDGET) samples /= 2;

    t->tex_ms_color = texture_ms(samples, GL_RGBA16F, width, height);
    t->tex_ms_tod   = texture_ms(samples, GL_RG16F, width, height);
    t->tex_ms_depth = texture_ms(samples, GL_DEPTH_COMPONENT24, width, height);
    t->tex_ms_peel_depth[0] = texture_ms(samples, GL_DEPTH_COMPONENT32F, width, height);
    t->tex_ms_peel_depth[1] = texture_ms(samples, GL_DEPTH_COMPONENT32F, width, height);
    t->tex_ms_peel_opacity  = texture_ms(samples, GL_RG16F, width, height);
    GLint actual = samples;
    glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_color);
    glGetTexLevelParameteriv(GL_TEXTURE_2D_MULTISAMPLE, 0, GL_TEXTURE_SAMPLES, &actual);
    glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);

    glGenFramebuffers(1, &t->fbo_color);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, t->fbo_color);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_color, 0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,  GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_depth, 0);
    glGenFramebuffers(1, &t->fbo_tod);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, t->fbo_tod);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_tod, 0);
    for (int i = 0; i < 2; ++i) {
        glGenFramebuffers(1, &t->fbo_peel[i]);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, t->fbo_peel[i]);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_peel_opacity, 0);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,  GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_peel_depth[i], 0);
    }

    for (int i = 0; i < 2; ++i) {
        gl::init_texture_2D(&t->tex_history_color[i], (int)width, (int)height, GL_RGBA16F);
        gl::init_texture_2D(&t->tex_history_t[i],     (int)width, (int)height, GL_RG16F);
    }
    glGenFramebuffers(1, &t->fbo);

    t->width   = width;
    t->height  = height;
    t->samples = (uint32_t)MAX(actual, 1);
    t->sample_cap = (uint32_t)cap;
    t->frame   = 0;
    t->hash    = 0;
    t->cur     = 0;

    if (!fbo_complete(t->fbo_color, "multisampled colour") || !fbo_complete(t->fbo_tod, "multisampled transmittance") ||
        !fbo_complete(t->fbo_peel[0], "first peel") || !fbo_complete(t->fbo_peel[1], "second peel")) {
        target_free(t);
        return false;
    }
    return true;
}

// The draw target and the state the passes change, put back afterwards
struct SavedState {
    GLint fbo = 0;
    GLint read_fbo = 0;
    GLint viewport[4] = {};
    GLint draw_buffer[8] = {};
    GLint draw_buffer_count = 0;
    GLboolean depth_test, blend, cull_face, scissor_test, depth_mask;
    GLboolean color_mask[4];
    GLint depth_func, blend_src_rgb, blend_dst_rgb, blend_src_alpha, blend_dst_alpha, blend_eq_rgb, blend_eq_alpha, front_face;
    GLint program, vao, active_texture;
    GLfloat clear_color[4];
    GLdouble clear_depth;
};

static SavedState save_state() {
    SavedState s = {};
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &s.fbo);
    glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &s.read_fbo);
    glGetIntegerv(GL_VIEWPORT, s.viewport);
    for (int i = 0; i < 8; ++i) {
        glGetIntegerv(GL_DRAW_BUFFER0 + i, &s.draw_buffer[i]);
        if (s.draw_buffer[i] != GL_NONE) s.draw_buffer_count = i + 1;
    }
    s.depth_test   = glIsEnabled(GL_DEPTH_TEST);
    s.blend        = glIsEnabled(GL_BLEND);
    s.cull_face    = glIsEnabled(GL_CULL_FACE);
    s.scissor_test = glIsEnabled(GL_SCISSOR_TEST);
    glGetBooleanv(GL_DEPTH_WRITEMASK, &s.depth_mask);
    glGetBooleanv(GL_COLOR_WRITEMASK, s.color_mask);
    glGetIntegerv(GL_DEPTH_FUNC, &s.depth_func);
    glGetIntegerv(GL_BLEND_SRC_RGB, &s.blend_src_rgb);
    glGetIntegerv(GL_BLEND_DST_RGB, &s.blend_dst_rgb);
    glGetIntegerv(GL_BLEND_SRC_ALPHA, &s.blend_src_alpha);
    glGetIntegerv(GL_BLEND_DST_ALPHA, &s.blend_dst_alpha);
    glGetIntegerv(GL_BLEND_EQUATION_RGB, &s.blend_eq_rgb);
    glGetIntegerv(GL_BLEND_EQUATION_ALPHA, &s.blend_eq_alpha);
    glGetIntegerv(GL_FRONT_FACE, &s.front_face);
    glGetIntegerv(GL_CURRENT_PROGRAM, &s.program);
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &s.vao);
    glGetIntegerv(GL_ACTIVE_TEXTURE, &s.active_texture);
    glGetFloatv(GL_COLOR_CLEAR_VALUE, s.clear_color);
    glGetDoublev(GL_DEPTH_CLEAR_VALUE, &s.clear_depth);
    return s;
}

static void restore_state(const SavedState& s) {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, s.fbo);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, s.read_fbo);
    glViewport(s.viewport[0], s.viewport[1], s.viewport[2], s.viewport[3]);
    // glDrawBuffers rejects the FRONT/BACK tokens of the default framebuffer (see volumerender_utils.cpp)
    if (s.fbo == 0) {
        glDrawBuffer(s.draw_buffer_count > 0 ? (GLenum)s.draw_buffer[0] : GL_NONE);
    } else if (s.draw_buffer_count > 0) {
        glDrawBuffers(s.draw_buffer_count, (const GLenum*)s.draw_buffer);
    } else {
        glDrawBuffer(GL_NONE);
    }
    if (s.depth_test)   glEnable(GL_DEPTH_TEST);   else glDisable(GL_DEPTH_TEST);
    if (s.blend)        glEnable(GL_BLEND);        else glDisable(GL_BLEND);
    if (s.cull_face)    glEnable(GL_CULL_FACE);    else glDisable(GL_CULL_FACE);
    if (s.scissor_test) glEnable(GL_SCISSOR_TEST); else glDisable(GL_SCISSOR_TEST);
    glDepthMask(s.depth_mask);
    glColorMask(s.color_mask[0], s.color_mask[1], s.color_mask[2], s.color_mask[3]);
    glDepthFunc(s.depth_func);
    glBlendFuncSeparate(s.blend_src_rgb, s.blend_dst_rgb, s.blend_src_alpha, s.blend_dst_alpha);
    glBlendEquationSeparate(s.blend_eq_rgb, s.blend_eq_alpha);
    glFrontFace(s.front_face);
    glUseProgram(s.program);
    glBindVertexArray(s.vao);
    glActiveTexture(s.active_texture);
    glClearColor(s.clear_color[0], s.clear_color[1], s.clear_color[2], s.clear_color[3]);
    glClearDepth(s.clear_depth);
}

static void bind_texture(int unit, GLenum target, GLuint tex) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(target, tex);
}

static void attach(GLuint fbo, GLuint tex0, GLuint tex1) {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex0, 0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, tex1, 0);
    const GLenum bufs[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
    glDrawBuffers(2, bufs);
}

static uint64_t desc_hash(const RenderDesc& desc) {
    uint64_t h = md_hash64(&desc.render_target.width, sizeof(desc.render_target.width), desc.scene_hash);
    h = md_hash64(&g_experiment, sizeof(g_experiment), h);
    h = md_hash64(&desc.render_target.height, sizeof(desc.render_target.height), h);
    h = md_hash64(&desc.render_target.depth, sizeof(desc.render_target.depth), h);
    h = md_hash64(&desc.matrix.view, sizeof(desc.matrix.view), h);
    h = md_hash64(&desc.matrix.proj_no_jitter, sizeof(desc.matrix.proj_no_jitter), h);
    h = md_hash64(&desc.shading, sizeof(desc.shading), h);
    for (size_t i = 0; i < desc.num_surfaces; ++i) {
        const Surface& s = desc.surfaces[i];
        const uint64_t mesh_id[3] = {(uint64_t)(uintptr_t)s.mesh, s.mesh ? s.mesh->version : 0, s.mesh ? s.mesh->vao : 0};
        h = md_hash64(mesh_id, sizeof(mesh_id), h);
        h = md_hash64(&s.index_offset, sizeof(s.index_offset), h);
        h = md_hash64(&s.index_count, sizeof(s.index_count), h);
        h = md_hash64(&s.model, sizeof(s.model), h);
        h = md_hash64(&s.color, sizeof(s.color), h);
        h = md_hash64(&s.optical_density, sizeof(s.optical_density), h);
        h = md_hash64(&s.coloring, sizeof(s.coloring), h);
        if (s.coloring != Coloring_Uniform) {
            h = md_hash64(&s.color_volume, sizeof(s.color_volume), h);
            h = md_hash64(&s.field_volume, sizeof(s.field_volume), h);
            h = md_hash64(&s.field_colormap, sizeof(s.field_colormap), h);
            h = md_hash64(&s.world_to_volume, sizeof(s.world_to_volume), h);
            h = md_hash64(&s.field_range_beg, sizeof(s.field_range_beg), h);
            h = md_hash64(&s.field_range_end, sizeof(s.field_range_end), h);
            // The lines change the image like the colours do, the transient one too: without it in the hash the
            // mean would keep the image from before it, and stop drawing once converged
            if (isolines_active(s.isolines)) {
                h = md_hash64(&s.isolines.count, sizeof(s.isolines.count), h);
                h = md_hash64(s.isolines.values, sizeof(float) * MIN(s.isolines.count, (uint32_t)ISOLINES_MAX), h);
                h = md_hash64(&s.isolines.hover, sizeof(s.isolines.hover), h);
                if (s.isolines.hover) h = md_hash64(&s.isolines.hover_value, sizeof(s.isolines.hover_value), h);
                h = md_hash64(&s.isolines.width, sizeof(s.isolines.width), h);
            }
        }
    }
    return h ? h : 1;
}

static bool drawable(const Surface& s) {
    return s.mesh && s.mesh->vao && s.index_count > 0 && s.index_offset + s.index_count <= s.mesh->index_count;
}

static void draw_surface(const Surface& s) {
    glBindVertexArray(s.mesh->vao);
    glDrawElements(GL_TRIANGLES, (GLsizei)s.index_count, GL_UNSIGNED_INT, (const void*)(sizeof(uint32_t) * s.index_offset));
}

bool render(Target* t, const RenderDesc& desc, bool* converged) {
    ASSERT(t);
    if (converged) *converged = true;
    const uint32_t w = desc.render_target.width;
    const uint32_t h = desc.render_target.height;
    if (!gl.ok || !desc.render_target.color || w == 0 || h == 0) return false;

    const SavedState saved = save_state();
    if (!target_init(t, w, h)) {
        restore_state(saved);
        return false;
    }

    const uint64_t hash = desc_hash(desc);
    if (hash != t->hash) {
        t->hash  = hash;
        t->frame = 0;
    }

    bool any = false;
    uint64_t tris = 0;
    vec3_t lo = {FLT_MAX, FLT_MAX, FLT_MAX}, hi = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (size_t i = 0; i < desc.num_surfaces; ++i) {
        const Surface& s = desc.surfaces[i];
        if (!drawable(s)) continue;
        any = true;
        tris += s.index_count / 3;
        // world space bounds, for the reference distance of the optical depth
        for (int c = 0; c < 8; ++c) {
            const vec3_t p = {(c & 1) ? s.mesh->bounds_max.x : s.mesh->bounds_min.x,
                              (c & 2) ? s.mesh->bounds_max.y : s.mesh->bounds_min.y,
                              (c & 4) ? s.mesh->bounds_max.z : s.mesh->bounds_min.z};
            const vec3_t q = vec3_from_vec4(s.model * vec4_from_vec3(p, 1.0f));
            lo = vec3_min(lo, q);
            hi = vec3_max(hi, q);
        }
    }

    t->triangles = tris;

    PUSH_GPU_SECTION("ISOSURFACE MESHES")
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_CULL_FACE);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glFrontFace(GL_CCW);
    glBlendEquation(GL_FUNC_ADD);

    if (any && (t->frame < MAX_FRAMES || desc.keep_drawing)) {
        const mat4_t& V = desc.matrix.view;
        const mat4_t& P = desc.matrix.proj;
        const bool use_depth = desc.render_target.depth != 0;
        glViewport(0, 0, (GLsizei)w, (GLsizei)h);
        glEnable(GL_MULTISAMPLE);

        // Distances along the rays are taken from the nearest point of the bounds (or the eye, inside them)
        const bool  ortho = P.elem[3][3] == 1.0f;
        const vec3_t center = vec3_lerp(lo, hi, 0.5f);
        const float radius = vec3_length(hi - lo) * 0.5f;
        const vec3_t center_view = vec3_from_vec4(V * vec4_from_vec3(center, 1.0f));
        const float od_ref = MAX(0.0f, (ortho ? -center_view.z : vec3_length(center_view)) - radius);

        const postprocess_pipeline::Lighting light = postprocess_pipeline::lighting(desc.shading.mode);
        const float env = MAX(desc.shading.environment, 0.0f);
        const vec3_t sky    = light.sky_radiance * env;
        const vec3_t ground = light.ground_radiance * env;
        const float n2 = desc.shading.ior;
        const float F0 = powf((1.0f - n2) / (1.0f + n2), 2.0f);
        const mat4_t inv_proj = mat4_inverse(P);
        const float inv_res[2] = {1.0f / (float)w, 1.0f / (float)h};

        bind_texture(0, GL_TEXTURE_2D, desc.render_target.depth);
        bind_texture(5, GL_TEXTURE_2D, gl.tex_blue_noise);

        // What the draws of the frame share, in both programs
        auto set_frame_uniforms = [&](const MeshProgram& m) {
            glUseProgram(m.program);
            glUniform1i(m.tex_depth, 0);
            glUniform1i(m.tex_tod, 1);
            glUniform1i(m.tex_color_volume, 2);
            glUniform1i(m.tex_field, 3);
            glUniform1i(m.tex_field_colormap, 4);
            glUniform1i(m.tex_blue_noise, 5);
            glUniform1i(m.tex_peel_depth0, 6);
            glUniform1i(m.tex_peel_depth1, 7);
            glUniform1i(m.tex_peel_opacity, 8);
            glUniform1i(m.use_depth, use_depth ? 1 : 0);
            glUniform1ui(m.frame, t->frame);
            glUniform1i(m.samples, (GLint)t->samples);
            glUniformMatrix4fv(m.view, 1, GL_FALSE, V.elem[0]);
            glUniformMatrix4fv(m.proj, 1, GL_FALSE, P.elem[0]);
            glUniformMatrix4fv(m.inv_proj, 1, GL_FALSE, inv_proj.elem[0]);
            glUniform2fv(m.inv_res, 1, inv_res);
            glUniform1i(m.ortho, ortho ? 1 : 0);
            glUniform1f(m.od_ref, od_ref);
            glUniform3fv(m.light_dir, 1, light.light_dir.elem);
            glUniform3fv(m.dir_radiance, 1, light.light_irradiance.elem);
            glUniform3fv(m.sky_radiance, 1, sky.elem);
            glUniform3fv(m.ground_radiance, 1, ground.elem);
            glUniform1f(m.roughness, desc.shading.roughness);
            glUniform1f(m.F0, F0);
            glUniform1f(m.basic, desc.shading.mode == postprocess_pipeline::ShadingMode_Basic ? 1.0f : 0.0f);
            glUniform1f(m.basic_ambient, CLAMP(desc.shading.ambient, 0.0f, 1.0f));
            glUniform1f(m.basic_specular, MAX(desc.shading.specular, 0.0f));
            glUniform1f(m.linear_colors, desc.shading.colors == postprocess_pipeline::ColorSpace_Linear ? 1.0f : 0.0f);
            glUniform1f(m.light_strength, MAX(desc.shading.light, 0.0f));
        };
        // The program of the passes per sample: the experiment's, if it has one of its own
        const bool no_fragdepth = g_experiment == Experiment_NoFragDepth || g_experiment == Experiment_OneSampleNoFragDepth;
        const MeshProgram& mesh = no_fragdepth ? gl.mesh_nfd : gl.mesh;
        set_frame_uniforms(gl.mesh_tod);
        set_frame_uniforms(mesh);

        auto draw_surfaces = [&](const MeshProgram& m, Pass pass) {
            glUseProgram(m.program);
            glUniform1i(m.pass, (GLint)pass);
            for (size_t i = 0; i < desc.num_surfaces; ++i) {
                const Surface& s = desc.surfaces[i];
                if (!drawable(s)) continue;
                int coloring = 0;
                if (s.coloring == Coloring_ColorVolume && s.color_volume) {
                    coloring = 1;
                    bind_texture(2, GL_TEXTURE_3D, s.color_volume);
                } else if (s.coloring == Coloring_VertexField && s.mesh->vbo_value && s.field_colormap) {
                    coloring = 3;
                    bind_texture(4, GL_TEXTURE_2D, s.field_colormap);
                } else if (s.coloring == Coloring_Field && s.field_volume && s.field_colormap) {
                    coloring = 2;
                    bind_texture(3, GL_TEXTURE_3D, s.field_volume);
                    bind_texture(4, GL_TEXTURE_2D, s.field_colormap);
                }
                const float ext = s.field_range_end - s.field_range_beg;
                glUniformMatrix4fv(m.model, 1, GL_FALSE, s.model.elem[0]);
                glUniform4fv(m.color, 1, s.color.elem);
                glUniform1i(m.coloring, coloring);
                glUniformMatrix4fv(m.world_to_volume, 1, GL_FALSE, s.world_to_volume.elem[0]);
                glUniform1f(m.field_beg, s.field_range_beg);
                glUniform1f(m.field_inv_ext, ext != 0.0f ? 1.0f / ext : 0.0f);
                glUniform1f(m.tau, MAX(s.optical_density, 0.0f) * OPTICAL_SCALE);
                glUniform1ui(m.surface_id, (GLuint)i);
                // Lines only on the values of the vertices: they carry the gradient the width is taken from
                const bool lines = coloring == 3 && isolines_active(s.isolines);
                const GLint line_count = lines ? (GLint)MIN(s.isolines.count, (uint32_t)ISOLINES_MAX) : 0;
                glUniform1i(m.isoline_count, line_count);
                if (line_count > 0) glUniform1fv(m.isoline_values, line_count, s.isolines.values);
                glUniform1i(m.isoline_hover, lines && s.isolines.hover ? 1 : 0);
                glUniform1f(m.isoline_hover_value, s.isolines.hover_value);
                glUniform1f(m.isoline_width, MAX(s.isolines.width, 0.0f));
                if (g_experiment == Experiment_GeometryOnly) glEnable(GL_RASTERIZER_DISCARD);
                draw_surface(s);
                if (g_experiment == Experiment_GeometryOnly) glDisable(GL_RASTERIZER_DISCARD);
            }
        };

        // The coverage of all the surfaces together and the optical depth of what they enclose, at pixel rate (its
        // own program): it compares with nothing per sample, so the coverage of the samples is all it needs of them
        PUSH_GPU_SECTION("Transmittance")
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, t->fbo_tod);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        {
            const GLfloat zero[4] = {0, 0, 0, 0};
            glClearBufferfv(GL_COLOR, 0, zero);
        }
        glDisable(GL_DEPTH_TEST);
        glDepthMask(GL_FALSE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);
        draw_surfaces(gl.mesh_tod, Pass_Transmittance);
        POP_GPU_SECTION()

        // The two nearest layers, exactly: their depths and opacities in each sample
        PUSH_GPU_SECTION("Peel")
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glDepthMask(GL_TRUE);
        glClearDepth(1.0);
        for (int k = 0; k < 2; ++k) {
            PUSH_GPU_SECTION(k == 0 ? "Peel 0" : "Peel 1")
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, t->fbo_peel[k]);
            glDrawBuffer(GL_COLOR_ATTACHMENT0);
            if (k == 0) {
                const GLfloat zero[4] = {0, 0, 0, 0};
                glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
                glClearBufferfv(GL_COLOR, 0, zero);     // opacity 0: no layer
            }
            glClear(GL_DEPTH_BUFFER_BIT);
            glColorMask(k == 0, k == 1, GL_FALSE, GL_FALSE);
            if (k == 1) bind_texture(6, GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_peel_depth[0]);
            draw_surfaces(mesh, k == 0 ? Pass_Peel0 : Pass_Peel1);
            POP_GPU_SECTION()
        }
        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
        POP_GPU_SECTION()

        // The stochastic depth of the layers behind them
        PUSH_GPU_SECTION("Stochastic depth")
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, t->fbo_color);
        glDrawBuffer(GL_NONE);
        glClear(GL_DEPTH_BUFFER_BIT);
        bind_texture(7, GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_peel_depth[1]);
        draw_surfaces(mesh, Pass_Depth);
        POP_GPU_SECTION()

        // Every surface times its visibility
        PUSH_GPU_SECTION("Accumulate")
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        {
            const GLfloat zero[4] = {0, 0, 0, 0};
            glClearBufferfv(GL_COLOR, 0, zero);
        }
        glDepthMask(GL_FALSE);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);
        bind_texture(1, GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_tod);
        bind_texture(8, GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_peel_opacity);
        draw_surfaces(mesh, Pass_Accumulate);
        const int ms_units[4] = {1, 6, 7, 8};
        for (int unit : ms_units) bind_texture(unit, GL_TEXTURE_2D_MULTISAMPLE, 0);
        glDisable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glDisable(GL_BLEND);
        POP_GPU_SECTION()

        // Resolved into the running mean
        PUSH_GPU_SECTION("Mean")
        const uint32_t next = 1 - t->cur;
        attach(t->fbo, t->tex_history_color[next], t->tex_history_t[next]);
        auto& a = gl.accumulate;
        glUseProgram(a.program);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_color);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, t->tex_ms_tod);
        bind_texture(2, GL_TEXTURE_2D, t->tex_history_color[t->cur]);
        bind_texture(3, GL_TEXTURE_2D, t->tex_history_t[t->cur]);
        glUniform1i(a.tex_color, 0);
        glUniform1i(a.tex_tod, 1);
        glUniform1i(a.tex_history_color, 2);
        glUniform1i(a.tex_history_t, 3);
        glUniform1i(a.samples, (GLint)t->samples);
        glUniform1f(a.weight, 1.0f / (float)(MIN(t->frame, MAX_FRAMES - 1) + 1));   // past MAX_FRAMES: keep_drawing
        glBindVertexArray(gl.vao);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);
        POP_GPU_SECTION()

        t->cur = next;
        t->frame += 1;
    }
    if (!any) {
        t->frame = 0;
    }

    // Over the caller's target
    bool done = false;
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, desc.render_target.color, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    const GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
    if (status == GL_FRAMEBUFFER_COMPLETE) {
        glViewport(0, 0, (GLsizei)w, (GLsizei)h);
        if (desc.render_target.clear_color) {
            glClearColor(0, 0, 0, 0);
            glClear(GL_COLOR_BUFFER_BIT);
        }
        if (any && t->frame > 0) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
            glUseProgram(gl.composite.program);
            bind_texture(0, GL_TEXTURE_2D, t->tex_history_color[t->cur]);
            bind_texture(1, GL_TEXTURE_2D, t->tex_history_t[t->cur]);
            glUniform1i(gl.composite.tex_color, 0);
            glUniform1i(gl.composite.tex_t, 1);
            glBindVertexArray(gl.vao);
            glDrawArrays(GL_TRIANGLES, 0, 3);
        }
        done = true;
    } else {
        MD_LOG_ERROR("Isosurface meshes: the colour target framebuffer is incomplete (0x%04X)", (unsigned int)status);
    }
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);

    for (int unit = 0; unit < 6; ++unit) {
        bind_texture(unit, GL_TEXTURE_2D, 0);
        bind_texture(unit, GL_TEXTURE_3D, 0);
    }
    POP_GPU_SECTION()
    restore_state(saved);

    if (converged) *converged = !any || (t->frame >= MAX_FRAMES && !desc.keep_drawing);
    return done;
}

}  // namespace isomesh
