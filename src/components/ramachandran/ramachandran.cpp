#include <event.h>

#include "density_gen.inl"
#include "density_gly.inl"
#include "density_pro.inl"
#include "density_pre.inl"

#if 0
#include "angles_gen.inl"
#include "angles_gly.inl"
#include "angles_pro.inl"
#include "angles_pre.inl"
#endif

#include <viamd.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <md_molecule.h>

#include "gfx/gl.h"
#include "gfx/gl_utils.h"
#include "image.h"
#include "task_system.h"

#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <imgui_internal.h>
#include <implot_internal.h>

#include <string.h>
#include <stdlib.h>

// @TODO:
// - Move Ramachandran indices FROM molecule.backbone into this, it is only relevant here
// - Merge the Compute tasks into a single task which iterates over the set of frames once


static const uint32_t density_tex_dim = 512;
static const uint32_t tex_dim = 1024;

typedef double density_map_t[180][180];

namespace ramachandran {

enum RamachandranDisplayMode {
    IsoLevels,
    IsoLines,
    Colormap,
};

// Representation (density, colormap, iso)
struct rama_rep_t {
    uint32_t map_tex[4] = {};
    uint32_t iso_tex[4] = {};
    float    den_sum[4] = {};	// Density sum for each ramachandran type (General, Glycine, Proline and PreProline)
    uint32_t den_tex = 0;		// This uses 4 channels, one for each ramachandran type (General, Glycine, Proline and PreProline)
};

struct rama_colormap_t {
    const uint32_t* colors;
    uint32_t count;
    float min_value;
    float max_value;
};

struct rama_isomap_t {
    const float* values;
    const uint32_t* level_colors;	// Optional
    const uint32_t* contour_colors; // Optional
    uint32_t count;
};

constexpr str_t v_fs_quad_src = STR_LIT(R"(
#version 150 core

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
}
)");

constexpr str_t f_shader_map_src = STR_LIT(R"(
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

constexpr str_t f_shader_iso_src = STR_LIT(R"(
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

static inline double b3(double t) {
    double at = fabs(t);
    double t1 = -at + 1.0;
    double t2 = -at + 2.0;
    t1 = t1 * t1 * t1 * (2.0 / 3.0);
    t2 = t2 * t2 * t2 / 6.0;
    return lerp(t2 - t1, t2, step(1.0, at));
}

static inline double bspline(double p0, double p1, double p2, double p3, double s) {
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

static inline void transpose(vec4_t* dst, const vec4_t* src, int dim) {
    // Block transpose
    const int block = 8;
    const int n = dim;
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

    const md_allocator_i* alloc = md_get_temp_allocator();    // Thread safe allocator!
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

    const md_allocator_i* alloc = md_get_temp_allocator();    // Thread safe allocator!
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

static void rama_rep_init(rama_rep_t* rep) {
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

static void rama_rep_free(rama_rep_t* rep) {
    glDeleteTextures(1, &rep->den_tex);
    glDeleteTextures(4, rep->map_tex);
    glDeleteTextures(4, rep->iso_tex);
}

struct Ramachandran : viamd::EventHandler {
    struct {
        GLuint program = 0;
        GLint uniform_loc_tex_den = -1;
        GLint uniform_loc_viewport = -1;
        GLint uniform_loc_inv_res = -1;
        GLint uniform_loc_map_colors = -1;
        GLint uniform_loc_map_offset = -1;
        GLint uniform_loc_map_length = -1;
        GLint uniform_loc_map_range = -1;
    } map;

    struct {
        GLuint program = 0;
        GLint uniform_loc_tex_den = -1;
        GLint uniform_loc_viewport = -1;
        GLint uniform_loc_inv_res = -1;
        GLint uniform_loc_iso_values = -1;
        GLint uniform_loc_iso_level_colors = -1;
        GLint uniform_loc_iso_contour_colors = -1;
        GLint uniform_loc_iso_offset = -1;
        GLint uniform_loc_iso_length = -1;
        GLint uniform_loc_iso_contour_line_scale = -1;
    } iso;

    GLuint fbo = 0;
    GLuint vao = 0;

    char input[256] = "all";
    char error[256] = "";

    bool evaluate = true;
    bool input_valid = false;
    bool show_window = false;

    uint64_t backbone_fingerprint = 0;
    uint64_t full_fingerprint = 0;
    uint64_t filt_fingerprint = 0;

    struct {
        rama_rep_t ref  = {};
        rama_rep_t full = {};
        rama_rep_t filt = {};
    } rama_data;

    md_array(uint32_t) rama_type_indices[4] = {};

    struct {
        ImVec4 border_color     = {0.0f, 0.0f, 0.0f, 1.0f};
        ImVec4 base_color       = {0.9f, 0.9f, 0.9f, 0.9f};
        ImVec4 selection_color  = {0.5f, 0.5f, 0.8f, 0.9f};
        ImVec4 highlight_color  = {0.8f, 0.8f, 0.5f, 0.9f};
        float base_radius       = 4.5f;
        float selected_radius   = 6.5f;
    } style;

    float blur_sigma = 5.0f;

    // DRAW STATE
    ImPlotRect selection_rect;
    SelectionOperator op = SelectionOperator::Or;
    bool is_selecting[4] = {false, false, false, false};

    float ref_alpha  = 0.85f;
    float full_alpha = 0.85f;
    float filt_alpha = 0.85f;

    RamachandranDisplayMode display_mode[3] = { IsoLevels, IsoLines, IsoLines };
    ImPlotColormap colormap[3] = { ImPlotColormap_Hot, ImPlotColormap_Plasma, ImPlotColormap_Viridis };
    ImVec4 isoline_colors[3]   = { {1,1,1,1}, {1,1,1,1}, {1,1,1,1} };

    int layout_mode = 0;
    bool show_layer[4] = {true, true, true, true};

    ImPlotRect viewrect = ImPlotRect(-180, 180, -180, 180);

    task_system::ID compute_density_full = 0;
    task_system::ID compute_density_filt = 0;

    md_allocator_i* arena = 0;

    Ramachandran() { viamd::event_system_register_handler(*this); }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                initialize(state);
                break;
            }
            case viamd::EventType_ViamdShutdown:
                shutdown();
                break;
            case viamd::EventType_ViamdFrameTick: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                update(state);
                draw(state);
                break;
            }
            case viamd::EventType_ViamdDrawMenu:
                ImGui::Checkbox("Ramachandran", &show_window);
                break;
            default:
                break;
            }
        }
    }

    void initialize(ApplicationState& state) {
        arena = md_arena_allocator_create(state.allocator.persistent, MEGABYTES(1));

        if (!map.program) {
            GLuint v_shader = gl::compile_shader_from_source(v_fs_quad_src, GL_VERTEX_SHADER);
            GLuint f_shader = gl::compile_shader_from_source(f_shader_map_src, GL_FRAGMENT_SHADER);
            defer {
                glDeleteShader(v_shader);
                glDeleteShader(f_shader);
            };

            map.program = glCreateProgram();
            const GLuint shaders[] = {v_shader, f_shader};
            gl::attach_link_detach(map.program, shaders, (int)ARRAY_SIZE(shaders));

            map.uniform_loc_tex_den    = glGetUniformLocation(map.program, "u_tex_den");
            map.uniform_loc_viewport   = glGetUniformLocation(map.program, "u_viewport");
            map.uniform_loc_inv_res    = glGetUniformLocation(map.program, "u_inv_res");
            map.uniform_loc_map_colors = glGetUniformLocation(map.program, "u_map_colors");
            map.uniform_loc_map_offset = glGetUniformLocation(map.program, "u_map_offset");
            map.uniform_loc_map_length = glGetUniformLocation(map.program, "u_map_length");
            map.uniform_loc_map_range  = glGetUniformLocation(map.program, "u_map_range");
        }

        if (!iso.program) {
            GLuint v_shader = gl::compile_shader_from_source(v_fs_quad_src, GL_VERTEX_SHADER);
            GLuint f_shader = gl::compile_shader_from_source(f_shader_iso_src, GL_FRAGMENT_SHADER);
            defer {
                glDeleteShader(v_shader);
            glDeleteShader(f_shader);
            };

            iso.program = glCreateProgram();
            const GLuint shaders[] = {v_shader, f_shader};
            gl::attach_link_detach(iso.program, shaders, (int)ARRAY_SIZE(shaders));

            iso.uniform_loc_tex_den    = glGetUniformLocation(iso.program, "u_tex_den");
            iso.uniform_loc_viewport   = glGetUniformLocation(iso.program, "u_viewport");
            iso.uniform_loc_inv_res    = glGetUniformLocation(iso.program, "u_inv_res");
            iso.uniform_loc_iso_values = glGetUniformLocation(iso.program, "u_iso_values");
            iso.uniform_loc_iso_level_colors = glGetUniformLocation(iso.program, "u_iso_level_colors");
            iso.uniform_loc_iso_contour_colors = glGetUniformLocation(iso.program, "u_iso_contour_colors");
            iso.uniform_loc_iso_offset = glGetUniformLocation(iso.program, "u_iso_offset");
            iso.uniform_loc_iso_length = glGetUniformLocation(iso.program, "u_iso_length");
            iso.uniform_loc_iso_contour_line_scale = glGetUniformLocation(iso.program, "u_iso_contour_line_scale");
        }

        if (!fbo) {
            glGenFramebuffers(1, &fbo);
        }

        if (!vao) {
            glGenVertexArrays(1, &vao);
        }

        rama_rep_init(&rama_data.ref);
        rama_rep_init(&rama_data.full);
        rama_rep_init(&rama_data.filt);

        rama_init_ref(&rama_data.ref);
    }

    void shutdown() {
        rama_rep_free(&rama_data.ref);
        rama_rep_free(&rama_data.full);
        rama_rep_free(&rama_data.filt);

        if (fbo) glDeleteFramebuffers(1, &fbo);
        if (vao) glDeleteVertexArrays(1, &vao);
    }

    void update(ApplicationState& state) {
        const size_t num_frames = md_trajectory_num_frames(state.mold.traj);
        if (show_window && state.mold.mol.protein_backbone.count > 0 && num_frames > 0) {
            if (backbone_fingerprint != state.trajectory_data.backbone_angles.fingerprint) {

                if (task_system::task_is_running(compute_density_full)) {
                    task_system::task_interrupt(compute_density_full);
                }
                if (task_system::task_is_running(compute_density_filt)) {
                    task_system::task_interrupt(compute_density_filt);
                }

                if (!task_system::task_is_running(compute_density_full) &&
                    !task_system::task_is_running(compute_density_filt))
                {
                    backbone_fingerprint = state.trajectory_data.backbone_angles.fingerprint;

                    md_array_shrink(rama_type_indices[0], 0);
                    md_array_shrink(rama_type_indices[1], 0);
                    md_array_shrink(rama_type_indices[2], 0);
                    md_array_shrink(rama_type_indices[3], 0);

                    for (uint32_t i = 0; i < (uint32_t)md_array_size(state.mold.mol.protein_backbone.ramachandran_type); ++i) {
                        switch (state.mold.mol.protein_backbone.ramachandran_type[i]) {
                        case MD_RAMACHANDRAN_TYPE_GENERAL: md_array_push(rama_type_indices[0], i, arena); break;
                        case MD_RAMACHANDRAN_TYPE_GLYCINE: md_array_push(rama_type_indices[1], i, arena); break;
                        case MD_RAMACHANDRAN_TYPE_PROLINE: md_array_push(rama_type_indices[2], i, arena); break;
                        case MD_RAMACHANDRAN_TYPE_PREPROL: md_array_push(rama_type_indices[3], i, arena); break;
                        default: break;
                        }
                    }

                    full_fingerprint = 0;
                    filt_fingerprint = 0;
                }
            }

            if (full_fingerprint != state.trajectory_data.backbone_angles.fingerprint) {
                if (!task_system::task_is_running(compute_density_full)) {
                    full_fingerprint = state.trajectory_data.backbone_angles.fingerprint;

                    const uint32_t frame_beg = 0;
                    const uint32_t frame_end = (uint32_t)num_frames;
                    const uint32_t frame_stride = (uint32_t)state.trajectory_data.backbone_angles.stride;

                    compute_density_full = rama_rep_compute_density(&rama_data.full, state.trajectory_data.backbone_angles.data, frame_beg, frame_end, frame_stride);
                } else {
                    task_system::task_interrupt(compute_density_full);
                }
            }

            if (filt_fingerprint != state.timeline.filter.fingerprint) {
                if (!task_system::task_is_running(compute_density_filt)) {
                    filt_fingerprint = state.timeline.filter.fingerprint;

                    const uint32_t frame_beg = MIN((uint32_t)state.timeline.filter.beg_frame, (uint32_t)num_frames - 1);
                    const uint32_t frame_end = MIN((uint32_t)state.timeline.filter.end_frame + 1, (uint32_t)num_frames);
                    const uint32_t frame_stride = (uint32_t)state.trajectory_data.backbone_angles.stride;

                    compute_density_filt = rama_rep_compute_density(&rama_data.filt, state.trajectory_data.backbone_angles.data, frame_beg, frame_end, frame_stride);
                }
                else {
                    task_system::task_interrupt(compute_density_filt);
                }
            }
        }
    }

    void draw(ApplicationState& state) {
        if (!show_window) return;

        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(2, 2));
        defer { ImGui::PopStyleVar(1); };

        ImGui::SetNextWindowSize({300,350}, ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Ramachandran", &show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_MenuBar)) {
            constexpr const char* plot_labels[4] = {"##General", "##Glycine", "##Proline", "##Pre-Proline"};
            constexpr const char* layer_labels[3] = { "Reference", "Full Trajectory", "Filtered Trajectory" };
            constexpr const char* option_labels[3] = { "IsoLevels", "IsoLines", "Colormap" };

            constexpr const float min_ext = -180.0f;
            constexpr const float max_ext = 180.0f;
            constexpr const float reset_coords[2] = { min_ext, max_ext };
            constexpr const char* x_lbl = "\xc2\xb0\xcf\x86";   // utf8 Degree Phi
            constexpr const char* y_lbl = "\xc2\xb0\xcf\x88";   // utf8 Degree Psi

            constexpr const ImPlotFlags flags = ImPlotFlags_Equal | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect; // | ImPlotFlags_AntiAliased;
            constexpr const ImPlotFlags subplotflags = ImPlotSubplotFlags_NoResize | ImPlotSubplotFlags_NoMenus;
            constexpr const ImPlotAxisFlags axis_flags = ImPlotAxisFlags_Foreground | ImPlotAxisFlags_NoLabel | ImPlotAxisFlags_NoTickLabels;

            if (viewrect.Size().x < 1) {
                viewrect.X.Max = viewrect.X.Min + 1;
            }
            if (viewrect.Size().y < 1) {
                viewrect.Y.Max = viewrect.Y.Min + 1;
            }

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Layers")) {

                    ImGui::Separator();
                    ImGui::Text("Layers");
                    for (int i = 0; i < 3; ++i) {
                        ImGui::Checkbox(layer_labels[i], &show_layer[i]);
                        if (show_layer[i]) {
                            ImGui::PushID(i);
                            if (ImGui::BeginCombo(layer_labels[i], option_labels[display_mode[i]])) {
                                if (ImGui::Selectable(option_labels[IsoLevels], display_mode[i] == IsoLevels)) display_mode[i] = IsoLevels;
                                if (ImGui::Selectable(option_labels[IsoLines],  display_mode[i] == IsoLines))  display_mode[i] = IsoLines;
                                if (ImGui::Selectable(option_labels[Colormap],  display_mode[i] == Colormap))  display_mode[i] = Colormap;
                                ImGui::EndCombo();
                            }
                            if (display_mode[i] == Colormap) {
                                ImPlot::ColormapSelection("Color Map", &colormap[i]);
                            } else if (display_mode[i] == IsoLines) {
                                ImGui::ColorEdit4Minimal("Line Color", &isoline_colors[i].x);
                            }
                            ImGui::PopID();
                        }
                        ImGui::Separator();
                    }
                    ImGui::Checkbox("Current", &show_layer[3]);
                    if (show_layer[3]) {
                        ImGui::SliderFloat("Point Size", &style.base_radius, 1.0f, 10.0f);
                        ImGui::ColorEdit4Minimal("Point Color", &style.base_color.x);
                    }
                    ImGui::Separator();
                    if (ImGui::SliderFloat("Density Blur Sigma", &blur_sigma, 0.1f, 10.0f)) {
                        full_fingerprint = 0;
                        filt_fingerprint = 0;
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("View")) {
                    if (ImGui::Selectable("Side-By-Side", layout_mode == 0)) layout_mode = 0;
                    if (ImGui::Selectable("General",      layout_mode == 1)) layout_mode = 1;
                    if (ImGui::Selectable("Glycine",      layout_mode == 2)) layout_mode = 2;
                    if (ImGui::Selectable("Proline",      layout_mode == 3)) layout_mode = 3;
                    if (ImGui::Selectable("Preproline",   layout_mode == 4)) layout_mode = 4;
                    ImGui::EndMenu();
                }
                ImGui::EndMenuBar();
            }

            const auto& mol = state.mold.mol;
            md_bitfield_t* selection_mask = &state.selection.selection_mask;
            md_bitfield_t* highlight_mask = &state.selection.highlight_mask;

            const int plot_offset = MAX(0, layout_mode - 1);
            const int plot_cols = (layout_mode == 0) ? 2 : 1;
            const int plot_rows = (layout_mode == 0) ? 2 : 1;
            const int num_plots = plot_cols * plot_rows;

            ImPlotInputMap& input_map = ImPlot::GetInputMap();
            input_map.OverrideMod = ImGuiMod_Shift;

            auto formatter = [](double value, char* buff, int size, void* user_data) -> int {
                const char* suffix = (const char*)user_data;
                value = deperiodize(value, 0, 360.0);
                return snprintf(buff, size, "%g%s", value, suffix);
            };

            ImVec2 plot_size = {0,0};

            ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(2,2));
            defer { ImPlot::PopStyleVar(); };

            if (ImPlot::BeginSubplots("##Ramaplots", plot_rows, plot_cols, ImVec2(-1, -1), subplotflags | ImPlotSubplotFlags_ShareItems)) {
                for (int plot_idx = plot_offset; plot_idx < plot_offset + num_plots; ++plot_idx) {
                    if (ImPlot::BeginPlot(plot_labels[plot_idx], ImVec2(), flags)) {
                        ImPlotPoint view_mid = { (viewrect.X.Min + viewrect.X.Max) * 0.5, (viewrect.Y.Min + viewrect.Y.Max) * 0.5 };

                        ImPlot::SetupAxesLimits(min_ext, max_ext, min_ext, max_ext, ImPlotCond_Once);
                        ImPlot::SetupAxisLinks(ImAxis_X1, &viewrect.X.Min, &viewrect.X.Max);
                        ImPlot::SetupAxisLinks(ImAxis_Y1, &viewrect.Y.Min, &viewrect.Y.Max);
                        // @NOTE(Robin): This wont work out of the box due to the periodic domain.
                        //ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, -720, +720);
                        //ImPlot::SetupAxisLimitsConstraints(ImAxis_Y1, -720, +720);
                        ImPlot::SetupAxes(x_lbl, y_lbl, axis_flags, axis_flags);

                        ImPlot::SetupAxisFormat(ImAxis_X1, formatter, (void*)x_lbl);
                        ImPlot::SetupAxisFormat(ImAxis_Y1, formatter, (void*)y_lbl);

                        ImPlot::SetupFinish();

                        viewrect = ImPlot::GetPlotLimits();

                        ImPlot::PushStyleVar(ImPlotStyleVar_Marker, ImPlotMarker_Square);
                        ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0,0,0,0));
                        ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0,0,0,0));
                        ImPlot::PlotScatter("##Hidden reset helper", reset_coords, reset_coords, 2);
                        ImPlot::PopStyleColor(2);
                        ImPlot::PopStyleVar();

                        //ImPlot::PlotDummy("Reference");
                        //ImPlot::PlotDummy("Full");
                        //ImPlot::PlotDummy("Filtered");
                        //ImPlot::PlotDummy("Current");

                        //ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(3)->Color = 0xFFFFFFFF;

                        //bool show_ref  = ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(0)->Show;
                        //bool show_full = ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(1)->Show;
                        //bool show_filt = ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(2)->Show;
                        //bool show_curr = ImPlot::GetCurrentContext()->CurrentSubplot->Items.GetLegendItem(3)->Show;
                        const bool show_ref  = show_layer[0];
                        const bool show_full = show_layer[1];
                        const bool show_filt = show_layer[2];
                        const bool show_curr = show_layer[3];

                        ImVec2 plot_min = ImPlot::PlotToPixels(viewrect.Min());
                        ImVec2 plot_max = ImPlot::PlotToPixels(viewrect.Max());
                        plot_size = { fabsf(plot_max.x - plot_min.x), fabsf(plot_max.y - plot_min.y) };

                        if (show_ref || show_full || show_filt) {
                            ImPlot::PushPlotClipRect();
                            ImDrawList* dl = ImPlot::GetPlotDrawList();

                            const uint32_t ref_tex  = display_mode[0] == Colormap ? rama_data.ref.map_tex[plot_idx]  : rama_data.ref.iso_tex[plot_idx];
                            const uint32_t full_tex = display_mode[1] == Colormap ? rama_data.full.map_tex[plot_idx] : rama_data.full.iso_tex[plot_idx];
                            const uint32_t filt_tex = display_mode[2] == Colormap ? rama_data.filt.map_tex[plot_idx] : rama_data.filt.iso_tex[plot_idx];

                            if (show_ref)  dl->AddImage((ImTextureID)(intptr_t)ref_tex,  plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, ref_alpha));
                            if (show_full) dl->AddImage((ImTextureID)(intptr_t)full_tex, plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, full_alpha));
                            if (show_filt) dl->AddImage((ImTextureID)(intptr_t)filt_tex, plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, filt_alpha));

                            ImPlot::PopPlotClipRect();
                        }

                        const bool hovered = ImPlot::IsPlotHovered();
                        const bool shift   = ImGui::GetIO().KeyShift;
                        ImPlotPoint mouse_coord = ImPlot::GetPlotMousePos();

                        if (is_selecting[plot_idx]) {
                            selection_rect.X.Max = mouse_coord.x;
                            selection_rect.Y.Max = mouse_coord.y;
                            if (selection_rect.Size().x != 0 && selection_rect.Size().y != 0) {
                                md_bitfield_copy(highlight_mask, selection_mask);
                                ImPlot::DragRect(0, &selection_rect.X.Min, &selection_rect.Y.Min, &selection_rect.X.Max, &selection_rect.Y.Max, ImVec4(1,1,1,0.5f), ImPlotDragToolFlags_NoInputs);
                            }

                            if (!shift) {
                                is_selecting[plot_idx] = false;
                            }
                        } else if (hovered && shift && (ImGui::IsMouseClicked(0) || ImGui::IsMouseClicked(1))) {
                            is_selecting[plot_idx] = true;
                            selection_rect = ImPlotRect(mouse_coord.x, mouse_coord.x, mouse_coord.y, mouse_coord.y);
                            op = ImGui::IsMouseClicked(0) ? SelectionOperator::Or : SelectionOperator::AndNot;
                        }

                        const double cut_d2 = fabs(ImPlot::PixelsToPlot(9,0).x - ImPlot::PixelsToPlot(0,0).x);
                        double min_d2 = DBL_MAX;
                        int64_t mouse_hover_idx = -1;

                        if (hovered) {
                            md_bitfield_clear(highlight_mask);
                        }

                        if (show_curr && mol.protein_backbone.angle) {
                            const uint32_t* indices = rama_type_indices[plot_idx];

                            double min_x = MIN(selection_rect.X.Min, selection_rect.X.Max);
                            double max_x = MAX(selection_rect.X.Min, selection_rect.X.Max);
                            double min_y = MIN(selection_rect.Y.Min, selection_rect.Y.Max);
                            double max_y = MAX(selection_rect.Y.Min, selection_rect.Y.Max);

                            double ref_x = (min_x + max_x) * 0.5;
                            double ref_y = (min_y + max_y) * 0.5;

                            mouse_coord.x = deperiodize(mouse_coord.x, ref_x, 360.0);
                            mouse_coord.y = deperiodize(mouse_coord.y, ref_y, 360.0);

                            for (size_t i = 0; i < md_array_size(indices); ++i) {
                                uint32_t idx = indices[i];

                                if (mol.protein_backbone.angle[idx].phi == 0 && mol.protein_backbone.angle[idx].psi == 0) continue;

                                ImPlotPoint coord = ImPlotPoint(RAD_TO_DEG(mol.protein_backbone.angle[idx].phi), RAD_TO_DEG(mol.protein_backbone.angle[idx].psi));
                                coord.x = deperiodize(coord.x, ref_x, 360.0);
                                coord.y = deperiodize(coord.y, ref_y, 360.0);

                                if (is_selecting[plot_idx]) {
                                    if (min_x <= coord.x && coord.x <= max_x && min_y <= coord.y && coord.y <= max_y) {
                                        md_residue_idx_t res_idx = mol.protein_backbone.residue_idx[idx];
                                        if (res_idx < (int)mol.residue.count) {
                                            md_range_t range = md_residue_atom_range(mol.residue, res_idx);
                                            modify_field(highlight_mask, range, op);
                                        }
                                    }
                                }
                                if (hovered) {
                                    double d2 = (coord.x - mouse_coord.x) * (coord.x - mouse_coord.x) + (coord.y - mouse_coord.y) * (coord.y - mouse_coord.y);
                                    if (d2 < cut_d2 && d2 < min_d2) {
                                        min_d2 = d2;
                                        mouse_hover_idx = idx;
                                    }
                                }
                            }

                            if (is_selecting[plot_idx]) {
                                grow_mask_by_selection_granularity(highlight_mask, state.selection.granularity, mol);
                            }

                            if (mouse_hover_idx != -1) {
                                if (mouse_hover_idx < (int64_t)mol.protein_backbone.count) {
                                    md_residue_idx_t res_idx = mol.protein_backbone.residue_idx[mouse_hover_idx];
                                    if (res_idx < (int)mol.residue.count) {
                                        md_range_t range = md_residue_atom_range(mol.residue, res_idx);
                                        modify_field(highlight_mask, range, SelectionOperator::Or);
                                        grow_mask_by_selection_granularity(highlight_mask, state.selection.granularity, mol);
                                    }
                                    str_t lbl = LBL_TO_STR(mol.residue.name[res_idx]);
                                    ImGui::SetTooltip("res[%d]: %.*s %d", res_idx + 1, (int)lbl.len, lbl.ptr, mol.residue.id[res_idx]);
                                }
                            }

                            uint32_t* selection_indices = 0;
                            uint32_t* highlight_indices = 0;

                            for (uint32_t i = 0; i < (uint32_t)md_array_size(indices); ++i) {
                                uint32_t idx = indices[i];
                                if (mol.protein_backbone.angle[idx].phi == 0 && mol.protein_backbone.angle[idx].psi == 0) continue;

                                int64_t atom_idx = mol.protein_backbone.atoms[idx].ca;
                                if (md_bitfield_test_bit(highlight_mask, atom_idx)) {
                                    md_array_push(highlight_indices, idx, state.allocator.frame);
                                }

                                if (md_bitfield_test_bit(selection_mask, atom_idx)) {
                                    md_array_push(selection_indices, idx, state.allocator.frame);
                                }
                            }

                            struct UserData {
                                const vec2_t*   coords;
                                const uint32_t* indices;
                                ImPlotPoint     view_center; // We use this to deperiodize the coordinates
                            };

                            auto index_getter = [](int i, void* user_data) -> ImPlotPoint {
                                const UserData* data = (UserData*)user_data;
                                uint32_t idx = data->indices[i];

                                if (data->coords[idx].x == 0 && data->coords[idx].y == 0) {
                                    return { INFINITY, INFINITY }; // Hide by INF!
                                }

                                double x = deperiodize(RAD_TO_DEG(data->coords[idx].x), data->view_center.x, 360.0);
                                double y = deperiodize(RAD_TO_DEG(data->coords[idx].y), data->view_center.y, 360.0);
                                return { x, y };
                            };

                            ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, style.base_radius, style.base_color, 1.2f, style.border_color);
                            ImPlot::SetNextLineStyle(ImVec4(1, 1, 1, 1));
                            if (md_array_size(indices) > 0) {
                                UserData user_data = { (const vec2_t*)(mol.protein_backbone.angle), indices, view_mid };
                                ImPlot::PlotScatterG("##Current", index_getter, &user_data, (int)md_array_size(indices));
                            }

                            if (md_array_size(selection_indices) > 0) {
                                UserData user_data = { (const vec2_t*)(mol.protein_backbone.angle), selection_indices, view_mid };
                                ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, style.base_radius + 1, style.selection_color, 2.0f, style.border_color);
                                ImPlot::PlotScatterG("##Selection", index_getter, &user_data, (int)md_array_size(selection_indices));
                            }

                            if (md_array_size(highlight_indices) > 0) {
                                UserData user_data = { (const vec2_t*)(mol.protein_backbone.angle), highlight_indices, view_mid };
                                ImPlot::SetNextMarkerStyle(ImPlotMarker_Square, style.base_radius + 1, style.highlight_color, 2.0f, style.border_color);
                                ImPlot::PlotScatterG("##Highlight", index_getter, &user_data, (int)md_array_size(highlight_indices));
                            }
                        }

                        if (ImGui::IsMouseReleased(0) || ImGui::IsMouseReleased(1)) {
                            if (is_selecting[plot_idx]) {
                                is_selecting[plot_idx] = false;
                                auto size = selection_rect.Size();
                                if (size.x == 0 && size.y == 0) {
                                    if (mouse_hover_idx != -1) {
                                        modify_field(selection_mask, highlight_mask, op);
                                    } else {
                                        if (ImGui::IsMouseReleased(1)) {
                                            md_bitfield_clear(selection_mask);
                                        }
                                    }
                                }
                                else {
                                    // Commit whatever is in the highlight mask
                                    modify_field(selection_mask, highlight_mask, SelectionOperator::Set);
                                }
                            } 
                        }

                        ImPlot::EndPlot();
                    }
                }

                if (viewrect.X.Max < viewrect.X.Min) {
                    double tmp = viewrect.X.Min;
                    viewrect.X.Min = viewrect.X.Max;
                    viewrect.X.Max = tmp;
                }
                if (viewrect.Y.Max < viewrect.Y.Min) {
                    double tmp = viewrect.Y.Min;
                    viewrect.Y.Min = viewrect.Y.Max;
                    viewrect.Y.Max = tmp;
                }

                const double ratio = viewrect.Size().x / viewrect.Size().y;
                const double mid_x = (viewrect.X.Min + viewrect.X.Max) * 0.5;
                const double mid_y = (viewrect.Y.Min + viewrect.Y.Max) * 0.5;

                if (viewrect.Size().x > 360) {
                    viewrect.X.Min = mid_x - 180;
                    viewrect.X.Max = mid_x + 180;
                    viewrect.Y.Min = mid_y - 0.5 * viewrect.Size().x / ratio;
                    viewrect.Y.Max = mid_y + 0.5 * viewrect.Size().x / ratio;
                } else if (viewrect.Size().x < 10) {
                    viewrect.X.Min = mid_x - 5;
                    viewrect.X.Max = mid_x + 5;
                    viewrect.Y.Min = mid_y - 0.5 * viewrect.Size().x / ratio;
                    viewrect.Y.Max = mid_y + 0.5 * viewrect.Size().x / ratio;
                }

                if (viewrect.Size().y > 360) {
                    viewrect.Y.Min = mid_y - 180;
                    viewrect.Y.Max = mid_y + 180;
                    viewrect.X.Min = mid_x - 0.5 * viewrect.Size().y * ratio;
                    viewrect.X.Max = mid_x + 0.5 * viewrect.Size().y * ratio;
                } else if (viewrect.Size().y < 10) {
                    viewrect.Y.Min = mid_y - 5;
                    viewrect.Y.Max = mid_y + 5;
                    viewrect.X.Min = mid_x - 0.5 * viewrect.Size().y * ratio;
                    viewrect.X.Max = mid_x + 0.5 * viewrect.Size().y * ratio;
                }

                ImPlot::EndSubplots();
            }

            vec4_t viewport = { (float)viewrect.X.Min / 180.0f, (float)viewrect.Y.Min / 180.0f, (float)viewrect.X.Max / 180.0f, (float)viewrect.Y.Max / 180.0f };
            viewport = viewport * 0.5f + 0.5f;

            const float* ref_sum  = rama_data.ref.den_sum;
            const float* full_sum = rama_data.full.den_sum;
            const float* filt_sum = rama_data.filt.den_sum;

            const float ref_iso_values[4][3] = {
                {0, 0.0005f, 0.02f},  // 99.95%, 98% for General
                {0, 0.0020f, 0.02f},  // 99.80%, 98% for Others
                {0, 0.0020f, 0.02f},
                {0, 0.0020f, 0.02f},
            };

            const uint32_t ref_iso_level_colors[4][3] = {
                {0x00000000, 0xFFFFE8B3, 0xFFFFD97F},
                {0x00000000, 0xFFC5E8FF, 0xFF7FCCFF},
                {0x00000000, 0xFFC5FFD0, 0xFF8CFF7F},
                {0x00000000, 0xFFFFE8B3, 0xFFFFD97F}
            };

            PUSH_GPU_SECTION("RENDER RAMA");

            uint32_t max_dim = (uint32_t)MAX(plot_size.x, plot_size.y);

            rama_rep_t* reps[3] = {
                &rama_data.ref,
                &rama_data.full,
                &rama_data.filt
            };

            const float scl[3][4] = {
                { 1.0f, 1.0f, 1.0f, 1.0f },
                { full_sum[0] / ref_sum[0], full_sum[1] / ref_sum[1], full_sum[2] / ref_sum[2], full_sum[3] / ref_sum[3] },
                { filt_sum[0] / ref_sum[0], filt_sum[1] / ref_sum[1], filt_sum[2] / ref_sum[2], filt_sum[3] / ref_sum[3] },
            };

            for (uint32_t i = 0; i < 3; ++i) {
                if (display_mode[i] == IsoLevels || display_mode[i] == IsoLines) {
                    const uint32_t count = 3;
                    uint32_t colors[4][count] = {0};
                    uint32_t  lines[4][count] = {0};
                    float    values[4][count] = {0};

                    MEMCPY(values, ref_iso_values, sizeof(ref_iso_values));
                    for (uint32_t j = 0; j < count; ++j) {
                        values[0][j] *= MAX(scl[i][0], FLT_EPSILON);
                        values[1][j] *= MAX(scl[i][1], FLT_EPSILON);
                        values[2][j] *= MAX(scl[i][2], FLT_EPSILON);
                        values[3][j] *= MAX(scl[i][3], FLT_EPSILON);
                    }

                    if (display_mode[i] == IsoLevels) {
                        MEMCPY(colors, ref_iso_level_colors, sizeof(colors));
                        MEMCPY(lines,  ref_iso_level_colors, sizeof(lines));
                    }
                    else {
                        //
                        ImVec4 rgba_base  = isoline_colors[i];
                        ImVec4 rgba_light = isoline_colors[i];
                        float h, s, v;
                        ImGui::ColorConvertRGBtoHSV(rgba_base.x, rgba_base.y, rgba_base.z, h, s, v);
                        ImGui::ColorConvertHSVtoRGB(h, s * 0.5f, v, rgba_light.x, rgba_light.y, rgba_light.z);

                        uint32_t line_colors[3] = {
                            0,
                            ImGui::ColorConvertFloat4ToU32(rgba_base),
                            ImGui::ColorConvertFloat4ToU32(rgba_light),
                        };
                        for (int j = 0; j < 4; ++j) {
                            lines[j][0] = line_colors[0];
                            lines[j][1] = line_colors[1];
                            lines[j][2] = line_colors[2];
                        }
                    }

                    rama_isomap_t maps[4] = {
                        {
                            .values = values[0],
                            .level_colors = colors[0],
                            .contour_colors = lines[0],
                            .count = count,
                        },
                        {
                            .values = values[1],
                            .level_colors = colors[1],
                            .contour_colors = lines[1],
                            .count = count,
                        },
                        {
                            .values = values[2],
                            .level_colors = colors[2],
                            .contour_colors = lines[2],
                            .count = count,
                        },
                        {
                            .values = values[3],
                            .level_colors = colors[3],
                            .contour_colors = lines[3],
                            .count = count,
                        },
                    };

                    rama_rep_render_iso(reps[i], viewport.elem, maps, max_dim);
                } else if (display_mode[i] == Colormap){
                    // Fill in colors based on active color_map
                    uint32_t colors[32] = {0};
                    uint32_t num_colors = ImPlot::GetColormapSize(colormap[i]);
                    for (uint32_t j = 0; j < num_colors; ++j) {
                        colors[j] = ImPlot::GetColormapColorU32(j, colormap[i]);
                    }
                    // Mask off the first color to be transparent
                    colors[0] = colors[0] & 0x00FFFFFFU;

                    rama_colormap_t maps[4] = {
                        {
                            .colors = colors,
                            .count = num_colors,
                            .min_value = 0,
                            .max_value = 0.5f * scl[i][0],
                        },
                    {
                        .colors = colors,
                        .count = num_colors,
                        .min_value = 0,
                        .max_value = 0.5f * scl[i][1],
                    },
                    {
                        .colors = colors,
                        .count = num_colors,
                        .min_value = 0,
                        .max_value = 0.5f * scl[i][2],
                    },
                    {
                        .colors = colors,
                        .count = num_colors,
                        .min_value = 0,
                        .max_value = 0.5f * scl[i][3],
                    },
                    };
                    rama_rep_render_map(reps[i], viewport.elem, maps, max_dim);
                }
            }

            POP_GPU_SECTION();

            ImPlot::MapInputDefault();
        }
        ImGui::End();
    }

    bool rama_init_ref(rama_rep_t* rep) {
        const density_map_t* densities[4] = {
            &density_gen,
            &density_gly,
            &density_pro,
            &density_pre
        };

        const size_t mem_size = sizeof(float) * density_tex_dim * density_tex_dim * 4;
        float* density_map = (float*)md_alloc(md_get_heap_allocator(), mem_size);
        defer { md_free(md_get_heap_allocator(), density_map, mem_size); };
        MEMSET(density_map, 0, mem_size);

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

        for (uint32_t y = 0; y < density_tex_dim; ++y) {
            double v = (y / (double)(density_tex_dim - 1));
            for (uint32_t x = 0; x < density_tex_dim; ++x) {
                double u = (x / (double)(density_tex_dim - 1));
                uint32_t idx = 4 * (y * density_tex_dim + x);
                density_map[idx + 0] = (float)linear_sample_map(*densities[0], u, v);
                density_map[idx + 1] = (float)linear_sample_map(*densities[1], u, v);
                density_map[idx + 2] = (float)linear_sample_map(*densities[2], u, v);
                density_map[idx + 3] = (float)linear_sample_map(*densities[3], u, v);

                density_sum[0] += density_map[idx + 0];
                density_sum[1] += density_map[idx + 1];
                density_sum[2] += density_map[idx + 2];
                density_sum[3] += density_map[idx + 3];
            }
        }

        rep->den_sum[0] = (float)density_sum[0];
        rep->den_sum[1] = (float)density_sum[1];
        rep->den_sum[2] = (float)density_sum[2];
        rep->den_sum[3] = (float)density_sum[3];

        blur_density_box((vec4_t*)density_map, density_tex_dim, 1);

        glBindTexture(GL_TEXTURE_2D, rep->den_tex);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, density_tex_dim, density_tex_dim, GL_RGBA, GL_FLOAT, density_map);
        glBindTexture(GL_TEXTURE_2D, 0);

        return true;
    }

    task_system::ID rama_rep_compute_density(rama_rep_t* rep, const md_backbone_angles_t* angles, uint32_t frame_beg, uint32_t frame_end, uint32_t frame_stride) {
        struct UserData {
            uint64_t alloc_size;
            vec4_t* density_tex;
            rama_rep_t* rep;
            const md_backbone_angles_t* angles;
            const uint32_t* type_indices[4];
            uint32_t frame_beg;
            uint32_t frame_end;
            uint32_t frame_stride;
            float sigma;
            md_allocator_i* alloc;
        };

        md_allocator_i* alloc = md_get_heap_allocator();

        uint64_t tex_size = sizeof(vec4_t) * density_tex_dim * density_tex_dim;
        uint64_t alloc_size = sizeof(UserData) + tex_size + alignof(vec4_t);
        UserData* user_data = (UserData*)md_alloc(alloc, alloc_size);
        vec4_t* density_tex = (vec4_t*)NEXT_ALIGNED_ADDRESS(user_data + 1, alignof(vec4_t));
        memset(density_tex, 0, tex_size);

        user_data->alloc_size = alloc_size;
        user_data->density_tex = density_tex;
        user_data->rep = rep;
        user_data->angles = angles;
        user_data->type_indices[0] = rama_type_indices[0];
        user_data->type_indices[1] = rama_type_indices[1];
        user_data->type_indices[2] = rama_type_indices[2];
        user_data->type_indices[3] = rama_type_indices[3];
        user_data->frame_beg = frame_beg;
        user_data->frame_end = frame_end;
        user_data->frame_stride = frame_stride;
        user_data->sigma = blur_sigma;
        user_data->alloc = alloc;

        task_system::ID async_task = task_system::create_pool_task(STR_LIT("Rama density"), [data = user_data]() {
            const float angle_to_coord_scale = 1.0f / (2.0f * PI);
            const float angle_to_coord_offset = 0.5f;

            const uint32_t frame_beg = data->frame_beg;
            const uint32_t frame_end = data->frame_end;
            const uint32_t frame_stride = data->frame_stride;
            const md_backbone_angles_t* angles = data->angles;

            double sum[4] = {0,0,0,0};

            for (uint32_t f = frame_beg; f < frame_end; ++f) {
                for (uint32_t c = 0; c < 4; ++c) {
                    const uint32_t* indices = data->type_indices[c];
                    const uint32_t num_indices = (uint32_t)md_array_size(data->type_indices[c]);
                    vec4_t val = {0, 0, 0, 0};
                    val.elem[c] = 1.0f;
                    if (num_indices) {
                        for (uint32_t i = 0; i < num_indices; ++i) {
                            uint32_t idx = f * frame_stride + indices[i];
                            if ((angles[idx].phi == 0 && angles[idx].psi == 0)) continue;
                            float u = angles[idx].phi * angle_to_coord_scale + angle_to_coord_offset;
                            float v = angles[idx].psi * angle_to_coord_scale + angle_to_coord_offset;
                            uint32_t x = (uint32_t)(u * density_tex_dim) & (density_tex_dim - 1);
                            uint32_t y = (uint32_t)(v * density_tex_dim) & (density_tex_dim - 1);
                            ASSERT(x < density_tex_dim);
                            ASSERT(y < density_tex_dim);
                            data->density_tex[y * density_tex_dim + x] += val;
                            sum[c] += 1.0;
                        }
                    }
                }
            }

            blur_density_gaussian(data->density_tex, density_tex_dim, data->sigma);

            data->rep->den_sum[0] = (float)sum[0];
            data->rep->den_sum[1] = (float)sum[1];
            data->rep->den_sum[2] = (float)sum[2];
            data->rep->den_sum[3] = (float)sum[3];
        });

        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update rama texture"), [data = user_data]() {
            gl::set_texture_2D_data(data->rep->den_tex, data->density_tex, GL_RGBA32F);
            md_free(data->alloc, data, data->alloc_size);
        });

        task_system::set_task_dependency(main_task, async_task);
        task_system::enqueue_task(async_task);

        return async_task;
    }

    void rama_rep_render_map(rama_rep_t* rep, const float viewport[4], const rama_colormap_t colmap[4], uint32_t display_res) {
        (void)display_res;

        glDisable(GL_BLEND);
        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);

        glViewport(0, 0, tex_dim, tex_dim);

        const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
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
            length[i] = colmap[i].count;
            range[i]  = { colmap[i].min_value, colmap[i].max_value };
            for (uint32_t j = 0; j < colmap[i].count; ++j) {
                colors[idx++] = vec4_from_u32(colmap[i].colors[j]);
            }
        }

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, rep->den_tex);

        glUseProgram(map.program);
        glUniform1i(map.uniform_loc_tex_den, 0);
        glUniform2f(map.uniform_loc_inv_res, 1.0f / tex_dim, 1.0f / tex_dim);
        glUniform4fv(map.uniform_loc_viewport, 1, vp.elem);
        glUniform4fv(map.uniform_loc_map_colors,  (int)ARRAY_SIZE(colors), colors[0].elem);
        glUniform1uiv(map.uniform_loc_map_offset, (int)ARRAY_SIZE(offset), offset);
        glUniform1uiv(map.uniform_loc_map_length, (int)ARRAY_SIZE(length), length);
        glUniform2fv(map.uniform_loc_map_range,   (int)ARRAY_SIZE(range), range[0].elem);

        glBindVertexArray(vao);

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

        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
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

        glUseProgram(iso.program);

        glUniform1i (iso.uniform_loc_tex_den, 0);
        glUniform4fv(iso.uniform_loc_viewport, 1, viewport);
        glUniform2f (iso.uniform_loc_inv_res, 1.0f / tex_dim, 1.0f / tex_dim);
        glUniform4fv(iso.uniform_loc_viewport, 1, vp.elem);
        glUniform1fv(iso.uniform_loc_iso_values, (int)ARRAY_SIZE(values), values);
        glUniform4fv(iso.uniform_loc_iso_level_colors, (int)ARRAY_SIZE(level_colors), level_colors[0].elem);
        glUniform4fv(iso.uniform_loc_iso_contour_colors, (int)ARRAY_SIZE(contour_colors), contour_colors[0].elem);
        glUniform1uiv(iso.uniform_loc_iso_offset, (int)ARRAY_SIZE(offset), offset);
        glUniform1uiv(iso.uniform_loc_iso_length, (int)ARRAY_SIZE(length), length);
        glUniform1f(iso.uniform_loc_iso_contour_line_scale, contour_line_scale);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, rep->den_tex);

        glBindVertexArray(vao);

        // DRAW!
        glDrawArrays(GL_TRIANGLES, 0, 3);

        glBindVertexArray(0);
        glUseProgram(0);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

};

static Ramachandran instance = {};

}  // namespace ramachandran
