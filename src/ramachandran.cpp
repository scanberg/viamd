#include "ramachandran.h"

#include "ramachandran/density_gen.inl"
#include "ramachandran/density_gly.inl"
#include "ramachandran/density_pro.inl"
#include "ramachandran/density_pre.inl"

#include "ramachandran/angles_gen.inl"
#include "ramachandran/angles_gly.inl"
#include "ramachandran/angles_pro.inl"
#include "ramachandran/angles_pre.inl"

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
#include <stdlib.h>

static const uint32_t density_tex_dim = 512;
static const uint32_t tex_dim = 1024;

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
    static GLint uniform_loc_iso_level_colors = -1;
    static GLint uniform_loc_iso_contour_colors = -1;
    static GLint uniform_loc_iso_offset = -1;
    static GLint uniform_loc_iso_length = -1;
    static GLint uniform_loc_iso_contour_line_scale = -1;
}

namespace blur {
    static GLuint program = 0;
    static GLuint tex = 0;
    static GLint uniform_loc_tex = -1;
    static GLint uniform_loc_inv_res = -1;
    static GLint uniform_loc_step = -1;
}

constexpr str_t v_fs_quad_src = MAKE_STR(R"(
#version 150 core

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
}
)");

constexpr str_t f_shader_box_blur_src = MAKE_STR(R"(
#version 330 core

layout(location = 0) out vec4 out_frag;

uniform sampler2D u_tex;

uniform vec2  u_inv_res;
uniform vec2  u_step;

#define KERNEL_RAD 8

void main() {
    vec2 coord = vec2(gl_FragCoord.xy) * u_inv_res;

    out_frag = vec4(0,0,0,0);
    for (int i = -KERNEL_RAD; i <= KERNEL_RAD; ++i) {
        out_frag += texture(u_tex, coord + u_step * i);
    }

    out_frag /= float(KERNEL_RAD * 2);
}
)");

constexpr str_t f_shader_map_src = MAKE_STR(R"(
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
)");

constexpr str_t f_shader_iso_src = MAKE_STR(R"(
#version 330 core

layout(location = 0) out vec4 out_frag0;
layout(location = 1) out vec4 out_frag1;
layout(location = 2) out vec4 out_frag2;
layout(location = 3) out vec4 out_frag3;

uniform sampler2D u_tex_den;
uniform vec4 u_viewport;
uniform vec2 u_inv_res;

uniform float u_iso_values[32];
uniform vec4  u_iso_level_colors[32];
uniform vec4  u_iso_contour_colors[32];
uniform uint  u_iso_offset[4];
uniform uint  u_iso_length[4];
uniform float u_iso_contour_line_scale;

vec4 map_density(float val, uint idx) {
    vec4 base = vec4(0,0,0,0);
    vec4 contour = vec4(0,0,0,0);
    uint offset = u_iso_offset[idx];
    uint length = u_iso_length[idx];

    uint i = offset;
    for (; i < offset + length - 1U; ++i) {
        if (u_iso_values[i] <= val && val < u_iso_values[i + 1U]) break;
    }

    uint i0 = i;
    uint i1 = min(i + 1U, offset + length - 1U);

    // We interpolate the colors between index i and i + 1
    float v0 = u_iso_values[i0];
    float v1 = u_iso_values[i1];

    vec4 b0 = u_iso_level_colors[i0];
    vec4 b1 = u_iso_level_colors[i1];

    float dv = fwidth(val);
    float band = dv * 2.0 * u_iso_contour_line_scale;
    base = mix(b0, b1, smoothstep(v1 - band, v1 + band, val));

    for (uint i = offset; i < offset + length; ++i) {
        float  v = u_iso_values[i]; 
        contour += u_iso_contour_colors[i] * smoothstep(v - band, v, val) * (1.0 - smoothstep(v, v + band, val));
    }

    return contour + base * (1.0 - contour.a);
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
)");

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

static inline double sample_map(const double map[180][180], int x, int y) {
    return map[x][y];
}

static inline double lerp_sample(const double map[180][180], double x, double y) {
    const int i_x[2] = { ((int)x % 180), ((int)x + 1) % 180};
    const int i_y[2] = { (int)y % 180, ((int)y + 1) % 180};

    const double t_x = x - (int)x;
    const double t_y = y - (int)y;

    return lerp(
        lerp(sample_map(map, i_x[0], i_y[0]), sample_map(map, i_x[1], i_y[0]), t_x),
        lerp(sample_map(map, i_x[0], i_y[1]), sample_map(map, i_x[1], i_y[1]), t_x),
        t_y);
}

static inline double linear_sample_map(const double map[180][180], double x, double y) {
    double cx = x * 179;
    double cy = y * 179;
    return lerp_sample(map, cx, cy);
}

/*
static inline double gauss_sample(const double map[180][180], int x, int y) {
    double d = 0.0;
    const double w[] = {70 / 256.0, 56 / 256.0, 28 / 256.0, 8 / 256.0, 1 / 256.0};
    const int k_size = 9;
    for (int xx = 0; xx < k_size; ++xx) {
        int ix = x - k_size/2 + xx;
        double k_x = w[abs(xx - k_size / 2)];
        for (int yy = 0; yy < k_size; ++yy) {
            double k_y = w[abs(yy - k_size / 2)];
            int iy = y - k_size/2 + yy;
            d += sample_map(map, x, y) * k_x * k_y;
        }
    }
    return d;
}

static inline double box_blur(double map[180][180]) {
    const int KERNEL_RADIUS = 3;
    double tmp[180] = {0};
    for (int y = 0; y < 180; ++y) {
        for (int x = 0; x < 180; ++x) {
            tmp[x] = sample_map(map, x - 2, y) + sample_map(map, x - 1, y) + sample_map(map, x, y) + sample_map(map, x + 1, y) + sample_map(map, x + 2, y);
        }
        for (int x = 0; x < 180; ++x) {
            map[x][y] = tmp[x] * 0.2;
        }
    }

    for (int x = 0; x < 180; ++x) {
        for (int y = 0; y < 180; ++y) {
            tmp[x] = sample_map(map, x, y - 2) + sample_map(map, x, y - 1) + sample_map(map, x, y) + sample_map(map, x, y + 1) + sample_map(map, x, y + 2);
        }
        for (int y = 0; y < 180; ++y) {
            map[x][y] = tmp[y] * 0.2;
        }
    }
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
*/


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
        gl::attach_link_detach(map::program, shaders, (int)ARRAY_SIZE(shaders));

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
        gl::attach_link_detach(iso::program, shaders, (int)ARRAY_SIZE(shaders));

        iso::uniform_loc_tex_den    = glGetUniformLocation(iso::program, "u_tex_den");
        iso::uniform_loc_viewport   = glGetUniformLocation(iso::program, "u_viewport");
        iso::uniform_loc_inv_res    = glGetUniformLocation(iso::program, "u_inv_res");
        iso::uniform_loc_iso_values = glGetUniformLocation(iso::program, "u_iso_values");
        iso::uniform_loc_iso_level_colors = glGetUniformLocation(iso::program, "u_iso_level_colors");
        iso::uniform_loc_iso_contour_colors = glGetUniformLocation(iso::program, "u_iso_contour_colors");
        iso::uniform_loc_iso_offset = glGetUniformLocation(iso::program, "u_iso_offset");
        iso::uniform_loc_iso_length = glGetUniformLocation(iso::program, "u_iso_length");
        iso::uniform_loc_iso_contour_line_scale = glGetUniformLocation(iso::program, "u_iso_contour_line_scale");
    }

    if (!blur::program) {
        GLuint v_shader = gl::compile_shader_from_source(v_fs_quad_src, GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_source(f_shader_box_blur_src, GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(f_shader);
        };

        blur::program = glCreateProgram();
        const GLuint shaders[] = {v_shader, f_shader};
        gl::attach_link_detach(blur::program, shaders, (int)ARRAY_SIZE(shaders));

        blur::uniform_loc_tex         = glGetUniformLocation(blur::program, "u_tex");
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
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, density_tex_dim, density_tex_dim);
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

static inline void blur_rows_acc(vec4_t* out, const vec4_t* in, int dim, int kernel_width) {
    const int mod = dim - 1;
    const float scl = 1.0f / (2 * kernel_width + 1);

    for (int row = 0; row < dim; ++row) {
        const vec4_t* src_row = in + dim * row;
        vec4_t* dst_row = out + dim * row;

        vec4_t acc = vec4_zero();
        for (int x = -(kernel_width+1); x < kernel_width; ++x) {
            acc = acc + src_row[x & mod];
        }

        int x = 0;
        for (; x < kernel_width + 1; ++x) {
            acc = vec4_max(vec4_zero(), acc - src_row[(x -(kernel_width+1)) & mod] + src_row[x + kernel_width]);
            dst_row[x] = acc * scl;
        }

        for (; x < dim - kernel_width; ++x) {
            acc = vec4_max(vec4_zero(), acc - src_row[x -(kernel_width+1)] + src_row[x + kernel_width]);
            dst_row[x] = acc * scl;
        }

        for (; x < dim; ++x) {
            acc = vec4_max(vec4_zero(), acc - src_row[x -(kernel_width+1)] + src_row[(x + kernel_width) & mod]);
            dst_row[x] = acc * scl;
        }
    }
}

static inline void transpose(vec4_t* dst, const vec4_t* src, int n) {
#if 0
    // Na�ve version
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < n; ++y) {
            dst[y * n + x] = src[x * n + y];
        }
    }
#endif
    // Block transpose
    const int block = 8;
    ASSERT(n % block == 0);

    for (int i = 0; i < n; i += block) {
        for (int j = 0; j < n; j += block) {
            for (int k = i; k < i + block; ++k) {
                for (int l = j; l < j + block; ++l) {
                    dst[k + l*n] = src[l + k*n];
                }
            }
        }
    }
}

static void boxes_for_gauss(int* box_w, int n, float sigma) {  // Number of boxes, standard deviation
    ASSERT(box_w);
    float wIdeal = sqrtf((12 * sigma * sigma / n) + 1);  // Ideal averaging filter width
    int wl = (int)wIdeal;
    if (wl % 2 == 0) wl--;
    int wu = wl + 2;

    float mIdeal = (12 * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4);
    int m = (int)(mIdeal + 0.5f);

    for (int i = 0; i < n; i++) box_w[i] = (i < m ? wl : wu);
}

static void blur_density_box(vec4_t* data, int dim, int num_passes) {
    ASSERT(dim > 0 && (dim & (dim - 1)) == 0); // Ensure dimension is power of two

    const md_allocator_i* alloc = default_allocator;    // Thread safe allocator (And also, temp allocator may not accomodate such large allocation!)
    vec4_t* tmp_data = (vec4_t*)md_alloc(alloc, dim * dim * sizeof(vec4_t));
    defer { md_free(alloc, tmp_data, dim * dim * sizeof(vec4_t)); };

    const int kernel_width = 4;

    for (int rep = 0; rep < num_passes; ++rep) {
        blur_rows_acc(tmp_data, data, dim, kernel_width);
        blur_rows_acc(data, tmp_data, dim, kernel_width);
    }
    transpose(tmp_data, data, dim);

    for (int rep = 0; rep < num_passes; ++rep) {
        blur_rows_acc(data, tmp_data, dim, kernel_width);
        blur_rows_acc(tmp_data, data, dim, kernel_width);
    }
    transpose(data, tmp_data, dim);
}

static void blur_density_gaussian(vec4_t* data, int dim, float sigma) {
    ASSERT(dim > 0 && (dim & (dim - 1)) == 0); // Ensure dimension is power of two

    const md_allocator_i* alloc = default_allocator;    // Thread safe allocator (And also, temp allocator may not accomodate such large allocation!)
    vec4_t* tmp_data = (vec4_t*)md_alloc(alloc, dim * dim * sizeof(vec4_t));
    defer { md_free(alloc, tmp_data, dim * dim * sizeof(vec4_t)); };

    int box_w[3];
    boxes_for_gauss(box_w, 3, sigma);

    blur_rows_acc(tmp_data, data, dim, box_w[0]);
    blur_rows_acc(data, tmp_data, dim, box_w[1]);
    blur_rows_acc(tmp_data, data, dim, box_w[2]);
    transpose(data, tmp_data, dim);

    blur_rows_acc(tmp_data, data, dim, box_w[0]);
    blur_rows_acc(data, tmp_data, dim, box_w[1]);
    blur_rows_acc(tmp_data, data, dim, box_w[2]);
    transpose(data, tmp_data, dim);
}

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

    glViewport(0, 0, density_tex_dim, density_tex_dim);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ramachandran::fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glBindVertexArray(ramachandran::vao);

    glUseProgram(ramachandran::blur::program);

    const float one_over_dim = 1.0f / density_tex_dim;

    glUniform1i(ramachandran::blur::uniform_loc_tex, 0);
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
    ASSERT(rep);

    glGenTextures(1, &rep->den_tex);
    glGenTextures(4, rep->map_tex);
    glGenTextures(4, rep->iso_tex);

    glBindTexture(GL_TEXTURE_2D, rep->den_tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA32F, density_tex_dim, density_tex_dim);


    for (int i = 0; i < 4; ++i) {
        glBindTexture(GL_TEXTURE_2D, rep->map_tex[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, tex_dim, tex_dim);
    }

    for (int i = 0; i < 4; ++i) {
        glBindTexture(GL_TEXTURE_2D, rep->iso_tex[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, tex_dim, tex_dim);
    }

    glBindTexture(GL_TEXTURE_2D, 0);
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
    
    const density_map_t* densities[4] = {
        &density_gen,
        &density_gly,
        &density_pro,
        &density_pre
    };

    init_rama_rep(&data->ref);
    init_rama_rep(&data->full);
    init_rama_rep(&data->filt);

    const int64_t mem_size = sizeof(float) * density_tex_dim * density_tex_dim * 4;
    float* density_map = (float*)md_alloc(default_allocator, mem_size);
    memset(density_map, 0, mem_size);
    defer { md_free(default_allocator, density_map, mem_size); };

    // Create reference densities since these never change
    // Resample reference textures into a nicer power of two texture format using some upsampling scheme

    double ref_sum[4] = {0};

    const uint32_t dim = 180;
    for (uint32_t y = 0; y < dim; ++y) {
        for (uint32_t x = 0; x < dim; ++x) {
            ref_sum[0] += (*densities[0])[y][x];
            ref_sum[1] += (*densities[1])[y][x];
            ref_sum[2] += (*densities[2])[y][x];
            ref_sum[3] += (*densities[3])[y][x];
        }
    }

    double density_sum[4] = {0};

    /*
    const float angle_to_coord_scale = 1.0 / 360.0;
    const float angle_to_coord_offset = 0.5f;

    for (uint32_t i = 0; i < ARRAY_SIZE(ref_rama_gen); ++i) {
        float u = ref_rama_gen[i][1] * angle_to_coord_scale + angle_to_coord_offset;
        float v = ref_rama_gen[i][0] * angle_to_coord_scale + angle_to_coord_offset;
        uint32_t x = (uint32_t)(u * density_tex_dim) & (density_tex_dim - 1);
        uint32_t y = (uint32_t)(v * density_tex_dim) & (density_tex_dim - 1);
        density_map[4 * (y * density_tex_dim + x) + 0] += 1.0f;
    }

    for (uint32_t i = 0; i < ARRAY_SIZE(ref_rama_gly); ++i) {
        float u = ref_rama_gly[i][1] * angle_to_coord_scale + angle_to_coord_offset;
        float v = ref_rama_gly[i][0] * angle_to_coord_scale + angle_to_coord_offset;
        uint32_t x = (uint32_t)(u * density_tex_dim) & (density_tex_dim - 1);
        uint32_t y = (uint32_t)(v * density_tex_dim) & (density_tex_dim - 1);
        density_map[4 * (y * density_tex_dim + x) + 1] += 1.0f;
    }

    for (uint32_t i = 0; i < ARRAY_SIZE(ref_rama_pro); ++i) {
        float u = ref_rama_pro[i][1] * angle_to_coord_scale + angle_to_coord_offset;
        float v = ref_rama_pro[i][0] * angle_to_coord_scale + angle_to_coord_offset;
        uint32_t x = (uint32_t)(u * density_tex_dim) & (density_tex_dim - 1);
        uint32_t y = (uint32_t)(v * density_tex_dim) & (density_tex_dim - 1);
        density_map[4 * (y * density_tex_dim + x) + 2] += 1.0f;
    }

    for (uint32_t i = 0; i < ARRAY_SIZE(ref_rama_pre); ++i) {
        float u = ref_rama_pre[i][1] * angle_to_coord_scale + angle_to_coord_offset;
        float v = ref_rama_pre[i][0] * angle_to_coord_scale + angle_to_coord_offset;
        uint32_t x = (uint32_t)(u * density_tex_dim) & (density_tex_dim - 1);
        uint32_t y = (uint32_t)(v * density_tex_dim) & (density_tex_dim - 1);
        density_map[4 * (y * density_tex_dim + x) + 3] += 1.0f;
    }
    */

    density_sum[0] = (double)ARRAY_SIZE(ref_rama_gen);
    density_sum[1] = (double)ARRAY_SIZE(ref_rama_gly);
    density_sum[2] = (double)ARRAY_SIZE(ref_rama_pro);
    density_sum[3] = (double)ARRAY_SIZE(ref_rama_pre);

    for (uint32_t y = 0; y < density_tex_dim; ++y) {
        double v = (y / (double)(density_tex_dim - 1));
        for (uint32_t x = 0; x < density_tex_dim; ++x) {
            double u = (x / (double)(density_tex_dim - 1));
            uint32_t idx = 4 * (y * density_tex_dim + x);
            density_map[idx + 0] = (float)ramachandran::linear_sample_map(*densities[0], u, v);
            density_map[idx + 1] = (float)ramachandran::linear_sample_map(*densities[1], u, v);
            density_map[idx + 2] = (float)ramachandran::linear_sample_map(*densities[2], u, v);
            density_map[idx + 3] = (float)ramachandran::linear_sample_map(*densities[3], u, v);

            density_sum[0] += density_map[idx + 0];
            density_sum[1] += density_map[idx + 1];
            density_sum[2] += density_map[idx + 2];
            density_sum[3] += density_map[idx + 3];
        }
    }

    data->ref.den_sum[0] = (float)density_sum[0];
    data->ref.den_sum[1] = (float)density_sum[1];
    data->ref.den_sum[2] = (float)density_sum[2];
    data->ref.den_sum[3] = (float)density_sum[3];

    blur_density_box((vec4_t*)density_map, density_tex_dim, 1);

    glBindTexture(GL_TEXTURE_2D, data->ref.den_tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, density_tex_dim, density_tex_dim, GL_RGBA, GL_FLOAT, density_map);
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
    float sigma;
};

task_system::ID rama_rep_compute_density(rama_rep_t* rep, const md_backbone_angles_t* angles, const uint32_t* rama_type_indices[4], uint32_t frame_beg, uint32_t frame_end, uint32_t frame_stride, float sigma) {

    const int64_t tex_size = sizeof(float) * density_tex_dim * density_tex_dim * 4;
    float* density_tex = (float*)md_alloc(default_allocator, tex_size);
    memset(density_tex, 0, tex_size);

    UserData user_data = {
        .density_tex = density_tex,
        .rep = rep,
        .angles = angles,
        .frame_beg = frame_beg,
        .frame_end = frame_end,
        .frame_stride = frame_stride,
        .sigma = sigma,
    };

    memcpy(user_data.type_indices, rama_type_indices, 4 * sizeof(uint32_t*));

    task_system::ID id = task_system::pool_enqueue("Compute rama density", [data = user_data]() {
        const float angle_to_coord_scale = 1.0f / (2.0f * PI);
        const float angle_to_coord_offset = 0.5f;

        const uint32_t frame_beg = data.frame_beg;
        const uint32_t frame_end = data.frame_end;
        const uint32_t frame_stride = data.frame_stride;
        const md_backbone_angles_t* angles = data.angles;

        double sum[4] = {0,0,0,0};

        for (uint32_t f = frame_beg; f < frame_end; ++f) {
            for (uint32_t c = 0; c < 4; ++c) {
                const uint32_t* indices = data.type_indices[c];
                const uint32_t num_indices = (uint32_t)md_array_size(data.type_indices[c]);
                if (num_indices) {
                    for (uint32_t i = 0; i < num_indices; ++i) {
                        uint32_t idx = f * frame_stride + indices[i];
                        if ((angles[idx].phi == 0 && angles[idx].psi == 0)) continue;
                        float u = angles[idx].phi * angle_to_coord_scale + angle_to_coord_offset;
                        float v = angles[idx].psi * angle_to_coord_scale + angle_to_coord_offset;
                        uint32_t x = (uint32_t)(u * density_tex_dim) & (density_tex_dim - 1);
                        uint32_t y = (uint32_t)(v * density_tex_dim) & (density_tex_dim - 1);
                        data.density_tex[4 * (y * density_tex_dim + x) + c] += 1.0f;
                        sum[c] += 1.0;
                    }
                }
            }
        }

        //blur_density_cpu((vec4_t*)data.density_tex, density_tex_dim, 8);
        blur_density_gaussian((vec4_t*)data.density_tex, density_tex_dim, data.sigma);

        data.rep->den_sum[0] = (float)sum[0];
        data.rep->den_sum[1] = (float)sum[1];
        data.rep->den_sum[2] = (float)sum[2];
        data.rep->den_sum[3] = (float)sum[3];
    });

    task_system::main_enqueue("Update rama texture", [data = user_data]() {
       
        glBindTexture(GL_TEXTURE_2D, data.rep->den_tex);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, density_tex_dim, density_tex_dim, GL_RGBA, GL_FLOAT, data.density_tex);
        glBindTexture(GL_TEXTURE_2D, 0);
        
        md_free(default_allocator, data.density_tex, sizeof(float) * density_tex_dim * density_tex_dim * 4);

    }, id);

    return id;
}

void rama_rep_render_map(rama_rep_t* rep, const float viewport[4], const rama_colormap_t colormap[4], uint32_t display_res) {
    (void)display_res;

    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glViewport(0, 0, tex_dim, tex_dim);

    const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ramachandran::fbo);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);

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
    glUniform4fv(ramachandran::map::uniform_loc_map_colors,  (int)ARRAY_SIZE(colors), colors[0].elem);
    glUniform1uiv(ramachandran::map::uniform_loc_map_offset, (int)ARRAY_SIZE(offset), offset);
    glUniform1uiv(ramachandran::map::uniform_loc_map_length, (int)ARRAY_SIZE(length), length);
    glUniform2fv(ramachandran::map::uniform_loc_map_range,   (int)ARRAY_SIZE(range), range[0].elem);

    glBindVertexArray(ramachandran::vao);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    glBindVertexArray(0);
    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}

void rama_rep_render_iso(rama_rep_t* rep, const float viewport[4], const rama_isomap_t isomap[4], uint32_t display_res) {
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glViewport(0, 0, tex_dim, tex_dim);

    const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ramachandran::fbo);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);

    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rep->iso_tex[0], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, rep->iso_tex[1], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, rep->iso_tex[2], 0);
    glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, rep->iso_tex[3], 0);

    vec4_t vp = {viewport[0], viewport[1], viewport[2] - viewport[0], viewport[3] - viewport[1]};

    const uint32_t cap = 32;
    uint32_t offset[4]  = {0};
    uint32_t length[4]  = {0};
    vec4_t   level_colors[cap] = {0};
    vec4_t   contour_colors[cap] = {0};
    float    values[cap] = {0};

    uint32_t idx = 0;
    for (uint32_t i = 0; i < 4; ++i) {
        offset[i] = idx;
        length[i] = isomap[i].count;
        for (uint32_t j = 0; j < MIN(isomap[i].count, cap); ++j) {
            if (isomap[i].level_colors) {
                level_colors[idx] = vec4_from_u32(isomap[i].level_colors[j]);
            }
            if (isomap[i].contour_colors) {
                contour_colors[idx] = vec4_from_u32(isomap[i].contour_colors[j]);
            }
            values[idx] = isomap[i].values[j];
            idx += 1;
        }
    }

    if (display_res == 0) {
        display_res = tex_dim;
    }

    if (display_res > tex_dim) {
        display_res = tex_dim;
    }

    float contour_line_scale = (float)tex_dim / (float)display_res;

    glUseProgram(ramachandran::iso::program);

    glUniform1i (ramachandran::iso::uniform_loc_tex_den, 0);
    glUniform4fv(ramachandran::iso::uniform_loc_viewport, 1, viewport);
    glUniform2f (ramachandran::iso::uniform_loc_inv_res, 1.0f / tex_dim, 1.0f / tex_dim);
    glUniform4fv(ramachandran::iso::uniform_loc_viewport, 1, vp.elem);
    glUniform1fv(ramachandran::iso::uniform_loc_iso_values, (int)ARRAY_SIZE(values), values);
    glUniform4fv(ramachandran::iso::uniform_loc_iso_level_colors, (int)ARRAY_SIZE(level_colors), level_colors[0].elem);
    glUniform4fv(ramachandran::iso::uniform_loc_iso_contour_colors, (int)ARRAY_SIZE(contour_colors), contour_colors[0].elem);
    glUniform1uiv(ramachandran::iso::uniform_loc_iso_offset, (int)ARRAY_SIZE(offset), offset);
    glUniform1uiv(ramachandran::iso::uniform_loc_iso_length, (int)ARRAY_SIZE(length), length);
    glUniform1f(ramachandran::iso::uniform_loc_iso_contour_line_scale, contour_line_scale);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, rep->den_tex);

    glBindVertexArray(ramachandran::vao);
    
    // DRAW!
    glDrawArrays(GL_TRIANGLES, 0, 3);

    glBindVertexArray(0);
    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
}