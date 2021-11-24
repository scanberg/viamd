#include "ramachandran.h"

#include "ramachandran/density_gen.inl"
#include "ramachandran/density_gly.inl"
#include "ramachandran/density_pro.inl"
#include "ramachandran/density_pre.inl"

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_array.inl>
#include <md_molecule.h>

#include "gfx/gl.h"
#include "gfx/gl_utils.h"
#include "image.h"
#include "task_system.h"

#include <string.h>

static const uint32_t tex_dim = 512;

typedef double density_map_t[180][180];

namespace ramachandran {

static GLuint fbo = 0;
static GLuint vao = 0;

namespace map {
    static GLuint program = 0;
    static GLint uniform_loc_tex_den = -1;
    static GLint uniform_loc_viewport = -1;
    static GLint uniform_loc_inv_res = -1;
    static GLint uniform_loc_map_colors = -1;
    static GLint uniform_loc_map_offset = -1;
    static GLint uniform_loc_map_length = -1;
    static GLint uniform_loc_map_range = -1;
};

namespace iso {
    static GLuint program = 0;
    static GLint uniform_loc_tex_den = -1;
    static GLint uniform_loc_viewport = -1;
    static GLint uniform_loc_inv_res = -1;
    static GLint uniform_loc_iso_values = -1;
    static GLint uniform_loc_iso_colors = -1;
    static GLint uniform_loc_iso_offset = -1;
    static GLint uniform_loc_iso_length = -1;
}

namespace blur {
    static GLuint program = 0;
    static GLuint tex = 0;
    static GLint uniform_loc_tex = -1;
    static GLint uniform_loc_weights = -1;
    static GLint uniform_loc_num_weights = -1;
    static GLint uniform_loc_inv_res = -1;
    static GLint uniform_loc_step = -1;
}

constexpr const char* v_fs_quad_src = R"(
#version 150 core

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
}
)";

constexpr const char* f_shader_gaussian_blur_src = R"(
#version 330 core

layout(location = 0) out vec4 out_frag;

uniform sampler2D u_tex;

uniform float u_weights[16];
uniform int   u_num_weights;
uniform vec2  u_inv_res;
uniform vec2  u_step;

void main() {
    vec2 coord = vec2(gl_FragCoord.xy) * u_inv_res;

    out_frag = vec4(0,0,0,0);
    for (int i = -u_num_weights + 1; i < u_num_weights; ++i) {
        out_frag += texture(u_tex, coord + u_step * i) * u_weights[abs(i)];
    }
}
)";

constexpr const char* f_shader_map_src = R"(
#version 330 core

layout(location = 0) out vec4 out_frag0;
layout(location = 1) out vec4 out_frag1;
layout(location = 2) out vec4 out_frag2;
layout(location = 3) out vec4 out_frag3;

uniform sampler2D u_tex_den;

uniform vec4 u_viewport;
uniform vec2 u_inv_res;

uniform vec4 u_map_colors[64];
uniform uint u_map_offset[4];
uniform uint u_map_length[4];
uniform vec2 u_map_range[4];

vec4 map_density(float val, uint map_idx) {
    uint offset = u_map_offset[map_idx];
    uint length = u_map_length[map_idx];
    vec2 range  = u_map_range[map_idx];
    
    val = clamp((val - range.x) / (range.y - range.x), 0, 1);
    float s = val * float(length);
    float t = fract(s);
    uint i0 = clamp(uint(s) + 0U, 0U, length - 1U);
    uint i1 = clamp(uint(s) + 1U, 0U, length - 1U);
    return mix(u_map_colors[offset + i0], u_map_colors[offset + i1], t);
}

void main() {
    vec2 coords = vec2(gl_FragCoord.xy) * u_inv_res;
    vec2 uv = u_viewport.xy + coords * u_viewport.zw;
    vec4 d = texture(u_tex_den, uv);
    out_frag0 = map_density(d.x, 0U);
    out_frag1 = map_density(d.y, 1U);
    out_frag2 = map_density(d.z, 2U);
    out_frag3 = map_density(d.w, 3U);
}
)";

constexpr const char* f_shader_iso_src = R"(
#version 330 core

layout(location = 0) out vec4 out_frag0;
layout(location = 1) out vec4 out_frag1;
layout(location = 2) out vec4 out_frag2;
layout(location = 3) out vec4 out_frag3;

uniform sampler2D u_tex_den;
uniform vec4 u_viewport;
uniform vec2 u_inv_res;

uniform float u_iso_values[64];
uniform vec4  u_iso_colors[64];
uniform uint  u_iso_offset[4];
uniform uint  u_iso_length[4];

vec4 map_density(float val, uint idx) {
    vec4 color = vec4(0,0,0,0);
    uint offset = u_iso_offset[idx];
    uint length = u_iso_length[idx];

    for (uint i = 0U; i < length; ++i) {
        if (val >= u_iso_values[offset + i]) {
            color = u_iso_colors[offset + i];
        }
    }

    return color;
}

void main() {
    vec2 coords = vec2(gl_FragCoord.xy) * u_inv_res;
    vec2 uv = u_viewport.xy + coords * u_viewport.zw;
    vec4 d = texture(u_tex_den, uv);
    
    out_frag0 = map_density(d.x, 0U);
    out_frag1 = map_density(d.y, 1U);
    out_frag2 = map_density(d.z, 2U);
    out_frag3 = map_density(d.w, 3U);
}
)";

double b3(double t) {
    double at = fabs(t);
    double t1 = -at + 1.0;
    double t2 = -at + 2.0;
    t1 = t1 * t1 * t1 * (2.0 / 3.0);
    t2 = t2 * t2 * t2 / 6.0;
    return lerp(t2 - t1, t2, step(1.0, at));
}

double bspline(double p0, double p1, double p2, double p3, double s) {
    return p0 * b3(s + 1.0) + p1 * b3(s) + p2 * b3(s - 1.0) + p3 * b3(s - 2.0);
}

static inline double lerp_sample(const double map[180][180], double x, double y) {
    int i_x[2] = { (int)x, (int)x + 1};
    int i_y[2] = { (int)y, (int)y + 1};

    if (i_x[0] < 0)   i_x[0] += 180;
    if (i_x[0] > 179) i_x[0] -= 180;

    if (i_x[1] < 0)   i_x[1] += 180;
    if (i_x[1] > 179) i_x[1] -= 180;

    if (i_y[0] < 0)   i_y[0] += 180;
    if (i_y[0] > 179) i_y[0] -= 180;

    if (i_y[1] < 0)   i_y[1] += 180;
    if (i_y[1] > 179) i_y[1] -= 180;

    double t_x = x - (int)x;
    double t_y = y - (int)y;

    double dy[2] = {
        lerp(map[i_x[0]][i_y[0]], map[i_x[1]][i_y[0]], t_x),
        lerp(map[i_x[0]][i_y[1]], map[i_x[1]][i_y[1]], t_x)
    };

    return lerp (dy[0], dy[1], t_y);
}

static inline double gauss_sample(const double map[180][180], int x, int y) {
    double d = 0.0;
    const double w[] = {70 / 256.0, 56 / 256.0, 28 / 256.0, 8 / 256.0, 1 / 256.0};
    const int k_size = 9;
    for (int xx = 0; xx < k_size; ++xx) {
        int ix = x - k_size/2 + xx;
        if (ix < 0) ix += 180;
        else if (ix >= 180) ix -= 180;
        double k_x = w[abs(xx - k_size / 2)];
        for (int yy = 0; yy < k_size; ++yy) {
            double k_y = w[abs(yy - k_size / 2)];
            int iy = y - k_size/2 + yy;
            if (iy < 0) iy += 180;
            else if (iy >= 180) iy -= 180;
            d += map[ix][iy] * k_x * k_y;
        }
    }
    return d;
}

static inline double gauss_sample_map(const double map[180][180], double x, double y) {
    double cx = x * 180;
    double cy = y * 180;
    int    ix = (int)cx;
    int    iy = (int)cy;
    double tx = cx - ix;
    double ty = cy - iy;

    double y0 = lerp(gauss_sample(map, ix + 0, iy + 0), gauss_sample(map, ix + 1, iy + 0), tx);
    double y1 = lerp(gauss_sample(map, ix + 0, iy + 1), gauss_sample(map, ix + 1, iy + 1), tx);
    return lerp(y0, y1, ty);
}

static inline double linear_sample_map(const double map[180][180], double x, double y) {
    double cx = x * 180;
    double cy = y * 180;
    return lerp_sample(map, cx, cy);
}

void initialize() {
    if (!map::program) {
        GLuint v_shader = gl::compile_shader_from_source(v_fs_quad_src, GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_source(f_shader_map_src, GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(f_shader);
        };

        map::program = glCreateProgram();
        const GLuint shaders[] = {v_shader, f_shader};
        gl::attach_link_detach(map::program, shaders, ARRAY_SIZE(shaders));

        map::uniform_loc_tex_den    = glGetUniformLocation(map::program, "u_tex_den");
        map::uniform_loc_viewport   = glGetUniformLocation(map::program, "u_viewport");
        map::uniform_loc_inv_res    = glGetUniformLocation(map::program, "u_inv_res");
        map::uniform_loc_map_colors = glGetUniformLocation(map::program, "u_map_colors");
        map::uniform_loc_map_offset = glGetUniformLocation(map::program, "u_map_offset");
        map::uniform_loc_map_length = glGetUniformLocation(map::program, "u_map_length");
        map::uniform_loc_map_range  = glGetUniformLocation(map::program, "u_map_range");
    }

    if (!iso::program) {
        GLuint v_shader = gl::compile_shader_from_source(v_fs_quad_src, GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_source(f_shader_iso_src, GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
        glDeleteShader(f_shader);
        };

        iso::program = glCreateProgram();
        const GLuint shaders[] = {v_shader, f_shader};
        gl::attach_link_detach(iso::program, shaders, ARRAY_SIZE(shaders));

        iso::uniform_loc_tex_den    = glGetUniformLocation(iso::program, "u_tex_den");
        iso::uniform_loc_viewport   = glGetUniformLocation(iso::program, "u_viewport");
        iso::uniform_loc_inv_res    = glGetUniformLocation(iso::program, "u_inv_res");
        iso::uniform_loc_iso_values = glGetUniformLocation(iso::program, "u_iso_values");
        iso::uniform_loc_iso_colors = glGetUniformLocation(iso::program, "u_iso_colors");
        iso::uniform_loc_iso_offset = glGetUniformLocation(iso::program, "u_iso_offset");
        iso::uniform_loc_iso_length = glGetUniformLocation(iso::program, "u_iso_length");
    }

    if (!blur::program) {
        GLuint v_shader = gl::compile_shader_from_source(v_fs_quad_src, GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_source(f_shader_gaussian_blur_src, GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(f_shader);
        };

        blur::program = glCreateProgram();
        const GLuint shaders[] = {v_shader, f_shader};
        gl::attach_link_detach(blur::program, shaders, ARRAY_SIZE(shaders));

        blur::uniform_loc_tex         = glGetUniformLocation(blur::program, "u_tex");
        blur::uniform_loc_weights     = glGetUniformLocation(blur::program, "u_weights");
        blur::uniform_loc_num_weights = glGetUniformLocation(blur::program, "u_num_weights");
        blur::uniform_loc_inv_res     = glGetUniformLocation(blur::program, "u_inv_res");
        blur::uniform_loc_step        = glGetUniformLocation(blur::program, "u_step");
    }

    if (!blur::tex) {
        glGenTextures(1, &blur::tex);
        glBindTexture(GL_TEXTURE_2D, blur::tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, tex_dim, tex_dim);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!fbo) {
        glGenFramebuffers(1, &fbo);
    }

    if (!vao) {
        glGenVertexArrays(1, &vao);
    }
}

void shutdown() {
    if (fbo) glDeleteFramebuffers(1, &fbo);
    if (vao) glDeleteVertexArrays(1, &vao);
}

}  // namespace ramachandran

static void blur_density(uint32_t density_tex, uint32_t num_passes) {
    // Backup GL state
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    GLboolean last_enable_blend = glIsEnabled(GL_BLEND);
    GLboolean last_enable_cull_face = glIsEnabled(GL_CULL_FACE);
    GLboolean last_enable_depth_test = glIsEnabled(GL_DEPTH_TEST);
    GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);

    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    const float weights[] = { 252 / 1024.0f, 210 / 1024.0f, 120 / 1024.0f, 45 / 1024.0f, 10 / 1024.0f, 1 / 1024.0f };

    glViewport(0, 0, tex_dim, tex_dim);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ramachandran::fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glBindVertexArray(ramachandran::vao);

    glUseProgram(ramachandran::blur::program);

    const float one_over_dim = 1.0f / tex_dim;

    glUniform1i(ramachandran::blur::uniform_loc_tex, 0);
    glUniform1fv(ramachandran::blur::uniform_loc_weights, ARRAY_SIZE(weights), weights);
    glUniform1i(ramachandran::blur::uniform_loc_num_weights, ARRAY_SIZE(weights));
    glUniform2f(ramachandran::blur::uniform_loc_inv_res, one_over_dim, one_over_dim);

    glActiveTexture(GL_TEXTURE0);
    for (uint32_t i = 0; i < num_passes; ++i) {
        glUniform2f(ramachandran::blur::uniform_loc_step, one_over_dim, 0);
        glBindTexture(GL_TEXTURE_2D, density_tex);
        glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, ramachandran::blur::tex, 0);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glUniform2f(ramachandran::blur::uniform_loc_step, 0, one_over_dim);
        glBindTexture(GL_TEXTURE_2D, ramachandran::blur::tex);
        glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, density_tex, 0);
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }

    glBindVertexArray(0);
    glUseProgram(0);

    // Restore modified GL state
    if (last_enable_blend)
        glEnable(GL_BLEND);
    else
        glDisable(GL_BLEND);
    if (last_enable_cull_face)
        glEnable(GL_CULL_FACE);
    else
        glDisable(GL_CULL_FACE);
    if (last_enable_depth_test)
        glEnable(GL_DEPTH_TEST);
    else
        glDisable(GL_DEPTH_TEST);
    if (last_enable_scissor_test)
        glEnable(GL_SCISSOR_TEST);
    else
        glDisable(GL_SCISSOR_TEST);
    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

static void init_rama_rep(rama_rep_t* rep) {
    glGenTextures(1, &rep->den_tex);
    glGenTextures(4, rep->map_tex);
    glGenTextures(4, rep->iso_tex);

    glBindTexture(GL_TEXTURE_2D, rep->den_tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, tex_dim, tex_dim);

    for (int i = 0; i < 4; ++i) {
        glBindTexture(GL_TEXTURE_2D, rep->map_tex[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, tex_dim, tex_dim);
    }

    for (int i = 0; i < 4; ++i) {
        glBindTexture(GL_TEXTURE_2D, rep->iso_tex[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, tex_dim, tex_dim);
    }
}

static void free_rama_rep(rama_rep_t* rep) {
    glDeleteTextures(1, &rep->den_tex);
    glDeleteTextures(4, rep->map_tex);
    glDeleteTextures(4, rep->iso_tex);
}

bool rama_free(rama_data_t* data) {
    if (!data) {
        return false;
    }
    free_rama_rep(&data->ref);
    free_rama_rep(&data->full);
    free_rama_rep(&data->filt);

    return true;
}

bool rama_init(rama_data_t* data) {
    ASSERT(data);

    rama_free(data);
    memset(data, 0, sizeof(rama_data_t));
    
    const density_map_t* densities[4] = {
        &density_gen,
        &density_gly,
        &density_pro,
        &density_pre
    };

    init_rama_rep(&data->ref);
    init_rama_rep(&data->full);
    init_rama_rep(&data->filt);

    const int64_t mem_size = sizeof(float) * tex_dim * tex_dim * 4;
    float* density_map = (float*)md_alloc(default_allocator, mem_size);
    defer { md_free(default_allocator, density_map, mem_size); };

    // Create reference densities since these never change
    // Resample reference textures into a nicer power of two texture format using cubic bspline upsampling

    for (int y = 0; y < tex_dim; ++y) {
        double v = 1.0 - (y / (double)(tex_dim - 1));
        for (int x = 0; x < tex_dim; ++x) {
            double u = (x / (double)(tex_dim - 1));
            int idx = 4 * (y * tex_dim + x);
            density_map[idx + 0] = (float)ramachandran::gauss_sample_map(*densities[0], u, v);
            density_map[idx + 1] = (float)ramachandran::gauss_sample_map(*densities[1], u, v);
            density_map[idx + 2] = (float)ramachandran::gauss_sample_map(*densities[2], u, v);
            density_map[idx + 3] = (float)ramachandran::gauss_sample_map(*densities[3], u, v);
        }
    }

    glBindTexture(GL_TEXTURE_2D, data->ref.den_tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, tex_dim, tex_dim, GL_RGBA, GL_FLOAT, density_map);
    glBindTexture(GL_TEXTURE_2D, 0);

    return true;
}

struct UserData {
    float* density_tex;
    rama_rep_t* rep;
    const md_backbone_angles_t* angles;
    const uint32_t* type_indices[4];
    uint32_t frame_beg;
    uint32_t frame_end;
    uint32_t frame_stride;
};

task_system::ID rama_rep_compute_density(rama_rep_t* rep, const md_backbone_angles_t* angles, const uint32_t* rama_type_indices[4], uint32_t frame_beg, uint32_t frame_end, uint32_t frame_stride) {
    const int64_t tex_size = sizeof(float) * tex_dim * tex_dim * 4;
    float* density_tex = (float*)md_alloc(default_allocator, tex_size);
    memset(density_tex, 0, tex_size);

    UserData* user_data = (UserData*)md_alloc(default_allocator, sizeof(UserData));
    user_data->density_tex = density_tex;
    user_data->rep = rep;
    user_data->angles = angles;
    memcpy(user_data->type_indices, rama_type_indices, 4 * sizeof(uint32_t*));
    user_data->frame_beg = frame_beg;
    user_data->frame_end = frame_end;
    user_data->frame_stride = frame_stride;

    task_system::ID id = task_system::pool_enqueue("Compute rama density", user_data, [](void* user_data) {
        UserData* data = (UserData*)user_data;

        const float scl = 1.0 / (2.0 * 3.14159265357989);

        const uint32_t frame_beg = data->frame_beg;
        const uint32_t frame_end = data->frame_end;
        const uint32_t frame_stride = data->frame_stride;
        const md_backbone_angles_t* angles = data->angles;

        for (uint32_t f = frame_beg; f < frame_end; ++f) {
            for (uint32_t c = 0; c < 4; ++c) {
                const uint32_t* indices = data->type_indices[c];
                const uint32_t num_indices = (uint32_t)md_array_size(data->type_indices[c]);
                if (num_indices) {
                    for (uint32_t i = 0; i < num_indices; ++i) {
                        uint32_t idx = f * frame_stride + indices[i];
                        if ((angles[idx].phi == 0 && angles[idx].psi == 0)) continue;
                        float u = angles[idx].phi * scl + 0.5f;
                        float v = 1.0f - (angles[idx].psi * scl + 0.5f);
                        uint32_t x = (uint32_t)(u * tex_dim) & (tex_dim - 1);
                        uint32_t y = (uint32_t)(v * tex_dim) & (tex_dim - 1);
                        data->density_tex[4 * (y * tex_dim + x) + c] += 1.0f;
                    }
                }
            }
        }
    },
    [](void* user_data) {
        task_system::main_enqueue("Update rama texture", user_data, [](void* user_data) {
            UserData* data = (UserData*)user_data;
            glBindTexture(GL_TEXTURE_2D, data->rep->den_tex);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, tex_dim, tex_dim, GL_RGBA, GL_FLOAT, data->density_tex);
            glBindTexture(GL_TEXTURE_2D, 0);

            blur_density(data->rep->den_tex, 16);

            md_free(default_allocator, data->density_tex, sizeof(float) * tex_dim * tex_dim * 4);
            md_free(default_allocator, data, sizeof(UserData));
        });
    });

    return id;
}

void rama_rep_render_map(rama_rep_t* rep, const float viewport[4], const rama_colormap_t colormap[4]) {
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glViewport(0, 0, tex_dim, tex_dim);

    const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ramachandran::fbo);
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);

    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rep->map_tex[0], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, rep->map_tex[1], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, rep->map_tex[2], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, rep->map_tex[3], 0);

    vec4_t vp = {viewport[0], viewport[1], viewport[2] - viewport[0], viewport[3] - viewport[1]};
    vec4_t   colors[64] = {0};
    uint32_t offset[4] = {0};
    uint32_t length[4] = {0};
    vec2_t   range[4] = {0};

    uint32_t idx = 0;
    for (uint32_t i = 0; i < 4; ++i) {
        offset[i] = idx;
        length[i] = colormap[i].count;
        range[i] = { colormap[i].min_value, colormap[i].max_value };
        for (uint32_t j = 0; j < colormap[i].count; ++j) {
            colors[idx++] = vec4_from_u32(colormap[i].colors[j]);
        }
    }

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, rep->den_tex);

    glUseProgram(ramachandran::map::program);
    glUniform1i(ramachandran::map::uniform_loc_tex_den, 0);
    glUniform2f(ramachandran::map::uniform_loc_inv_res, 1.0f / tex_dim, 1.0f / tex_dim);
    glUniform4fv(ramachandran::map::uniform_loc_viewport, 1, vp.elem);
    glUniform4fv(ramachandran::map::uniform_loc_map_colors, ARRAY_SIZE(colors), colors[0].elem);
    glUniform1uiv(ramachandran::map::uniform_loc_map_offset, ARRAY_SIZE(offset), offset);
    glUniform1uiv(ramachandran::map::uniform_loc_map_length, ARRAY_SIZE(length), length);
    glUniform2fv(ramachandran::map::uniform_loc_map_range, ARRAY_SIZE(range), range[0].elem);

    glBindVertexArray(ramachandran::vao);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    glBindVertexArray(0);
    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}

void rama_rep_render_iso(rama_rep_t* rep, const float viewport[4], const rama_isomap_t isomap[4]) {
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glViewport(0, 0, tex_dim, tex_dim);

    const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ramachandran::fbo);
    glDrawBuffers(ARRAY_SIZE(draw_buffers), draw_buffers);

    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rep->iso_tex[0], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, rep->iso_tex[1], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, rep->iso_tex[2], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, rep->iso_tex[3], 0);

    vec4_t vp = {viewport[0], viewport[1], viewport[2] - viewport[0], viewport[3] - viewport[1]};

    uint32_t offset[4]  = {0};
    uint32_t length[4]  = {0};
    vec4_t   colors[64] = {0};
    float    values[64] = {0};

    uint32_t idx = 0;
    for (uint32_t i = 0; i < 4; ++i) {
        offset[i] = idx;
        length[i] = isomap[i].count;
        for (uint32_t j = 0; j < isomap[i].count; ++j) {
            colors[idx] = vec4_from_u32(isomap[i].colors[j]);
            values[idx] = isomap[i].values[j];
            idx += 1;
        }
    }

    glUseProgram(ramachandran::iso::program);

    glUniform1i (ramachandran::iso::uniform_loc_tex_den, 0);
    glUniform4fv(ramachandran::iso::uniform_loc_viewport, 1, viewport);
    glUniform2f (ramachandran::iso::uniform_loc_inv_res, 1.0f / tex_dim, 1.0f / tex_dim);
    glUniform4fv(ramachandran::iso::uniform_loc_viewport, 1, vp.elem);
    glUniform1fv(ramachandran::iso::uniform_loc_iso_values, ARRAY_SIZE(values), values);
    glUniform4fv(ramachandran::iso::uniform_loc_iso_colors, ARRAY_SIZE(colors), colors[0].elem);
    glUniform1uiv(ramachandran::iso::uniform_loc_iso_offset, ARRAY_SIZE(offset), offset);
    glUniform1uiv(ramachandran::iso::uniform_loc_iso_length, ARRAY_SIZE(length), length);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, rep->den_tex);

    glBindVertexArray(ramachandran::vao);
    
    // DRAW!
    glDrawArrays(GL_TRIANGLES, 0, 3);

    glBindVertexArray(0);
    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}