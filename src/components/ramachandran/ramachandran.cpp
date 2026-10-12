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

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#include <core/md_array.h>
#include <core/md_hash.h>
#include <core/md_bitfield.h>
#include <md_system.h>
#include <md_util.h>
#include <md_filter.h>

#include <viamd_event.h>
#include <viamd.h>
#include <serialization_utils.h>
#include <event.h>
#include <components/analysis_density_shaders.h>
#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <image.h>
#include <task_system.h>

#include <imgui_widgets.h>
#include <implot_widgets.h>
#include <imgui_internal.h>
#include <implot_internal.h>

#include <string.h>
#include <stdlib.h>
#include <memory>
#include <vector>
#include <cmath>

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
    float    den_max[4] = {};	// The peak density of each type
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

constexpr str_t v_fs_quad_src = analysis_density::vertex_shader;
constexpr str_t f_shader_map_src = analysis_density::map_shader;
constexpr str_t f_shader_iso_src = analysis_density::iso_shader;

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

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    vec4_t* tmp_data = (vec4_t*)md_temp_alloc(temp, dim * dim * sizeof(vec4_t));

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

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    vec4_t* tmp_data = (vec4_t*)md_temp_alloc(temp, dim * dim * sizeof(vec4_t));

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

// What is kept for one series of the plot: the rows of its property as residues of the backbone, and
// the densities of its trajectory. Kept for as long as the application runs, as a task may still be on its
// way to deliver into it when the series leaves the plot.
struct RamaSeries {
    uint32_t uid = 0;
    char path[SERIES_PATH_CAP] = "";
    rama_rep_t rep = {};
    task_system::ID compute = 0;
    uint64_t generation = 1;
    uint64_t rep_generation = 0;
    uint64_t fingerprint = 0;
    uint64_t requested = 0;
    uint64_t pending_hash = 0;          // What the density should be made of, and since when
    double pending_time = 0.0;
    uint64_t map_hash = 0;
    bool mapped = false;                // The rows are known to be residues of the backbone
    bool used = false;
    // The rows of each ramachandran type (general, glycine, proline, pre-proline)
    md_array(uint32_t) type_indices[4] = {};
    // The backbone segment of each row, when mapped
    md_array(uint32_t) segments = 0;
};

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

    // The shortcut: a residue selection which becomes a ramachandran() property of the script
    char input[256] = "all";
    char error[256] = "";
    char property_name[64] = "rama";
    char property_error[256] = "";
    md_bitfield_t residue_mask = {};
    bool input_valid = false;
    uint64_t input_hash = 0;

    // A property just created, added to the plot once the script has evaluated it
    char pending_name[64] = "";
    double pending_time = 0.0;

    bool show_window = false;

    // What the plot shows. Only the first subplot is used, the layout is the one of the timeline's
    PlotSubplot plot[PLOT_MAX_SUBPLOTS];
    std::vector<std::unique_ptr<RamaSeries>> series_states;

    struct {
        rama_rep_t ref  = {};
    } rama_data;

    struct {
        ImVec4 base_outline         = {1.0f, 1.0f, 1.0f, 1.0f};
        ImVec4 base_fill            = {1.0f, 1.0f, 1.0f, 1.0f};
        ImVec4 selection_outline    = {0.3f, 0.3f, 0.8f, 0.9f};
        ImVec4 selection_fill       = {0.3f, 0.3f, 0.8f, 0.9f};
        ImVec4 highlight_outline    = {0.8f, 0.8f, 0.3f, 0.9f};
        ImVec4 highlight_fill       = {0.8f, 0.8f, 0.3f, 0.9f};
        float interaction_size      = 5.0f; // Size of points that are subjected to selection and highlight
        float interaction_weight    = 2.0f;
    } style;


    // DRAW STATE
    ImPlotRect selection_rect;
    SelectionOperator op = SelectionOperator::Or;
    bool is_selecting[4] = {false, false, false, false};

    // The reference layer, the one which is not a series
    float ref_alpha  = 0.85f;

    RamachandranDisplayMode ref_mode = IsoLevels;
    ImPlotColormap ref_colormap = ImPlotColormap_Hot;
    ImVec4 ref_line_color = {1,1,1,1};

    int layout_mode = 0;
    bool show_reference = true;

    ImPlotRect viewrect = ImPlotRect(-180, 180, -180, 180);

    ApplicationState* app_state = 0;

    Ramachandran() { viamd::event_system_register_handler(*this); }

    // ## Workspace: how the layers are drawn

    void reset_workspace_settings() {
        str_copy_to_char_buf(input, sizeof(input), STR_LIT("all"));
        str_copy_to_char_buf(property_name, sizeof(property_name), STR_LIT("rama"));
        property_error[0] = '\0';
        pending_name[0] = '\0';
        input_hash = 0;
        input_valid = false;
        plot_clear(plot, PLOT_MAX_SUBPLOTS);
        ref_alpha = 0.85f;
        ref_mode = IsoLevels;
        ref_colormap = ImPlotColormap_Hot;
        ref_line_color = {1,1,1,1};
        show_reference = true;
        layout_mode = 0;
        invalidate_densities();
    }

    void serialize(viamd::serialization_state_t& state) {
        viamd::write_section_header(state, STR_LIT("Ramachandran"));
        viamd::write_str(state, STR_LIT("Filter"), str_from_cstr(input));
        viamd::write_str(state, STR_LIT("PropertyName"), str_from_cstr(property_name));
        viamd::write_flt(state, STR_LIT("ReferenceAlpha"), ref_alpha);
        viamd::write_int(state, STR_LIT("ReferenceMode"), (int)ref_mode);
        viamd::write_int(state, STR_LIT("ReferenceColormap"), (int)ref_colormap);
        viamd::write_vec4(state, STR_LIT("ReferenceLineColor"), vec4_set(ref_line_color.x, ref_line_color.y, ref_line_color.z, ref_line_color.w));
        viamd::write_bool(state, STR_LIT("ShowReference"), show_reference);
        viamd::write_int(state, STR_LIT("LayoutMode"), layout_mode);
        plot_layout_serialize(state, STR_LIT("RamachandranLayout"), STR_LIT("RamachandranSeries"), plot, 1);
    }

    void deserialize(viamd::deserialization_state_t& state) {
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "Filter")) {
                str_t filter;
                if (viamd::extract_str(filter, arg)) str_copy_to_char_buf(input, sizeof(input), filter);
            } else if (str_eq_cstr(ident, "PropertyName")) {
                viamd::extract_to_char_buf(property_name, sizeof(property_name), arg);
            } else if (str_eq_cstr(ident, "ReferenceAlpha")) {
                viamd::extract_flt(ref_alpha, arg);
            } else if (str_eq_cstr(ident, "ReferenceMode")) {
                int m;
                if (viamd::extract_int(m, arg)) ref_mode = (RamachandranDisplayMode)CLAMP(m, (int)IsoLevels, (int)Colormap);
            } else if (str_eq_cstr(ident, "ReferenceColormap")) {
                int c;
                if (viamd::extract_int(c, arg)) ref_colormap = (ImPlotColormap)c;
            } else if (str_eq_cstr(ident, "ReferenceLineColor")) {
                vec4_t c;
                if (viamd::extract_vec4(c, arg)) ref_line_color = ImVec4(c.x, c.y, c.z, c.w);
            } else if (str_eq_cstr(ident, "ShowReference")) {
                viamd::extract_bool(show_reference, arg);
            } else if (str_eq_cstr(ident, "LayoutMode")) {
                viamd::extract_int(layout_mode, arg);
            }
        }
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                initialize(state);
                workspace_register_window("Ramachandran", &show_window);
                break;
            }
            case viamd::EventType_ViamdSerialize:
                serialize(*(viamd::serialization_state_t*)e.payload);
                break;
            case viamd::EventType_ViamdDeserializeBegin:
                reset_workspace_settings();
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                if (str_eq(viamd::section_header(state), STR_LIT("Ramachandran"))) {
                    deserialize(state);
                } else {
                    int num_subplots = 1;
                    plot_layout_deserialize(state, STR_LIT("RamachandranLayout"), STR_LIT("RamachandranSeries"), app_state, plot, &num_subplots);
                }
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
            case viamd::EventType_ViamdSystemInit: {
                ApplicationState& state = *(ApplicationState*)e.payload;
                on_topology_init(state);
                break;
            }
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Ramachandran", &show_window);
                break;
            default:
                break;
            }
        }
    }

    void initialize(ApplicationState& state) {
        app_state = &state;
        md_bitfield_init(&residue_mask, state.allocator.persistent);

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
        rama_init_ref(&rama_data.ref);
    }

    void shutdown() {
        for (auto& st : series_states) {
            task_system::task_interrupt_and_wait_for(st->compute);
            ++st->generation;
            rama_rep_free(&st->rep);
            for (int i = 0; i < 4; ++i) md_array_free(st->type_indices[i], md_get_heap_allocator());
            md_array_free(st->segments, md_get_heap_allocator());
        }
        series_states.clear();
        rama_rep_free(&rama_data.ref);

        if (fbo) glDeleteFramebuffers(1, &fbo);
        if (vao) glDeleteVertexArrays(1, &vao);
    }

    void on_topology_init(ApplicationState& state) {
        (void)state;
        for (auto& st : series_states) {
            task_system::task_interrupt_and_wait_for(st->compute);
            st->map_hash = 0;
        }
        invalidate_densities();
        input_hash = 0;
        input_valid = false;
    }

    // The densities are rebuilt: a setting which enters them changed
    void invalidate_densities() {
        for (auto& st : series_states) {
            ++st->generation;
            st->fingerprint = 0;
        }
    }

    RamaSeries* find_state(uint32_t uid) {
        for (auto& st : series_states) {
            if (st->uid == uid) return st.get();
        }
        return nullptr;
    }

    // The state of a series. One which no series uses is taken over before a new one is made: its
    // textures are large, and a task still on its way to it delivers nothing, as its generation moves on.
    RamaSeries& acquire_state(uint32_t uid, const char* path) {
        RamaSeries* st = find_state(uid);
        if (!st) {
            for (auto& candidate : series_states) {
                if (!candidate->used) { st = candidate.get(); break; }
            }
        }
        if (!st) {
            series_states.push_back(std::make_unique<RamaSeries>());
            st = series_states.back().get();
            rama_rep_init(&st->rep);
            st->uid = 0;
        }
        if (st->uid != uid || strcmp(st->path, path) != 0) {
            st->uid = uid;
            str_copy_to_char_buf(st->path, sizeof(st->path), str_from_cstr(path));
            st->map_hash = 0;
            st->mapped = false;
            ++st->generation;
            st->rep_generation = 0;
            st->fingerprint = 0;
        }
        return *st;
    }

    // The residues of the backbone which the rows of a property are, in the order of the system: what the
    // script's ramachandran() lists for the atoms its visualization marks. Without them the rows are all drawn
    // as general residues, and cannot be selected.
    void map_series(ApplicationState& state, RamaSeries& st, const SeriesVectorView& v) {
        const md_system_t& sys = state.mold.sys;
        struct { uint64_t ir; uint64_t topology; uint32_t count; uint32_t pad; } key = { state.script.ir_fingerprint, sys.topology_version, v.count, 0 };
        const uint64_t hash = md_hash64(&key, sizeof(key), 0) | 1;
        if (st.map_hash == hash) return;
        st.map_hash = hash;
        st.mapped = false;
        ++st.generation;
        st.rep_generation = 0;
        st.fingerprint = 0;

        md_allocator_i* alloc = md_get_heap_allocator();
        for (int i = 0; i < 4; ++i) md_array_shrink(st.type_indices[i], 0);
        md_array_shrink(st.segments, 0);

        const md_protein_backbone_data_t& bb = sys.protein_backbone;
        if (v.vis_payload && bb.segment.count > 0) {
            md_bitfield_t mask = {};
            md_bitfield_init(&mask, alloc);
            md_script_vis_t vis = {};
            md_script_vis_init(&vis, alloc);
            const md_script_vis_ctx_t ctx = { state.script.eval_ir, &sys, &state.mold.state };
            if (md_script_vis_eval_payload(&vis, v.vis_payload, -1, &ctx, MD_SCRIPT_VISUALIZE_ATOMS)) {
                md_bitfield_copy(&mask, &vis.atom_mask);
                for (size_t r = 0; r < bb.range.count; ++r) {
                    const uint32_t beg = bb.range.offset[r];
                    const uint32_t end = bb.range.offset[r + 1];
                    for (uint32_t i = beg + 1; i + 1 < end; ++i) {
                        const md_urange_t atoms = md_component_atom_range(&sys.component, bb.segment.comp_idx[i]);
                        if (md_bitfield_scan(&mask, atoms.beg, atoms.end)) md_array_push(st.segments, i, alloc);
                    }
                }
            }
            md_script_vis_free(&vis);
            md_bitfield_free(&mask);
        }

        if (md_array_size(st.segments) == v.count) {
            st.mapped = true;
            md_temp_scope_t temp = md_temp_begin();
            md_ramachandran_type_t* types = md_temp_alloc_array(temp, md_ramachandran_type_t, bb.segment.count);
            md_util_backbone_ramachandran_classify(types, bb.segment.count, &sys);
            for (uint32_t row = 0; row < v.count; ++row) {
                switch (types[st.segments[row]]) {
                case MD_RAMACHANDRAN_TYPE_GLYCINE: md_array_push(st.type_indices[1], row, alloc); break;
                case MD_RAMACHANDRAN_TYPE_PROLINE: md_array_push(st.type_indices[2], row, alloc); break;
                case MD_RAMACHANDRAN_TYPE_PREPROL: md_array_push(st.type_indices[3], row, alloc); break;
                default:                           md_array_push(st.type_indices[0], row, alloc); break;
                }
            }
            md_temp_end(temp);
        } else {
            md_array_shrink(st.segments, 0);
            for (uint32_t row = 0; row < v.count; ++row) md_array_push(st.type_indices[0], row, alloc);
        }
    }

    // The residue selection of the shortcut, evaluated when it changes
    void update_input(const ApplicationState& state) {
        const uint64_t hash = md_hash64_str(str_from_cstr(input), state.mold.sys.topology_version ^ state.script.ir_fingerprint) | 1;
        if (hash == input_hash) return;
        input_hash = hash;
        MEMSET(error, 0, sizeof(error));
        bool dynamic = false;
        input_valid = md_filter(&residue_mask, str_from_cstr(input), &state.mold.sys, &state.mold.state, state.script.ir, &dynamic, error, sizeof(error));
        if (input_valid && dynamic) {
            input_valid = false;
            str_copy_to_char_buf(error, sizeof(error), STR_LIT("Use a frame-independent residue selection."));
        }
    }

    void create_property(ApplicationState& state) {
        property_error[0] = '\0';
        if (!md_script_identifier_name_valid(str_from_cstr(property_name))) {
            snprintf(property_error, sizeof(property_error), "Enter a valid script identifier as the property name.");
            return;
        }
        const bool whole_system = str_eq(str_trim(str_from_cstr(input)), STR_LIT("all"));
        const str_t statement = str_printf(state.allocator.frame, "%s = ramachandran(%s);", property_name, whole_system ? "" : input);
        if (!script_append_statement(&state, statement)) {
            snprintf(property_error, sizeof(property_error), "Could not create the property: check the name and selection for script errors or a duplicate identifier.");
            return;
        }
        snprintf(pending_name, sizeof(pending_name), "%s", property_name);
        pending_time = ImGui::GetTime();
    }

    // Adds a property created from the shortcut once the script has evaluated it
    void update_pending(ApplicationState& state) {
        if (!pending_name[0]) return;
        char path[SERIES_PATH_CAP];
        snprintf(path, sizeof(path), "script/%s", pending_name);
        const SeriesKey key = series_key(SeriesSource_Script, str_from_cstr(path));
        SeriesVectorView v;
        if (series_resolve_vector(&v, &state, key, 2)) {
            plot_add_vector_series(&state, plot[0], key, false);
            pending_name[0] = '\0';
        } else if (ImGui::GetTime() - pending_time > 20.0) {
            pending_name[0] = '\0';
        }
    }

    static bool shows_density(const PlotSeries& s) {
        return s.iso_lines || s.iso_levels || s.density_map;
    }

    // The pixels of blur of a density, from the smoothing of the series
    static float blur_of(const PlotSeries& s) { return s.smoothing * 2.5f; }

    void update(ApplicationState& state) {
        app_state = &state;
        update_input(state);
        update_pending(state);

        const PlotSubplot& sp = plot[0];
        SeriesVectorView views[PLOT_MAX_SERIES_PER_SUBPLOT];
        bool resolved[PLOT_MAX_SERIES_PER_SUBPLOT] = {};

        // The states which are in use stay with their series, before any is taken over
        for (auto& st : series_states) st->used = false;
        for (int j = 0; j < sp.count; ++j) {
            resolved[j] = series_resolve_vector(&views[j], &state, sp.series[j].key, 2);
            if (!resolved[j]) continue;
            if (RamaSeries* st = find_state(sp.series[j].uid)) st->used = true;
        }

        for (int j = 0; j < sp.count; ++j) {
            if (!resolved[j]) continue;
            const PlotSeries& s = sp.series[j];
            const SeriesVectorView& v = views[j];
            RamaSeries& st = acquire_state(s.uid, s.key.path);
            st.used = true;
            map_series(state, st, v);

            const size_t num_frames = v.num_frames;
            if (!show_window || num_frames == 0 || !shows_density(s)) continue;

            // The density is derived from the frames of the property, so it is rebuilt exactly when
            // the property's version, or the blur, moves. A filtered property holds zeros outside
            // its frames, which the density leaves out.
            const md_backbone_angles_t* angles = (const md_backbone_angles_t*)v.data;
            const uint32_t stride = v.count;
            const float sigma = blur_of(s);

            struct { uint64_t version; float sigma; uint32_t pad; } key = { v.version, sigma, 0 };
            const uint64_t hash = md_hash64(&key, sizeof(key), st.generation);
            st.requested = hash;

            // A property changes with every frame range the script evaluates: the settled one is turned into a
            // density, and the one before it stays on screen until then
            const double now = ImGui::GetTime();
            if (hash != st.pending_hash) {
                st.pending_hash = hash;
                st.pending_time = now;
            }
            const bool settled = st.rep_generation == 0 || now - st.pending_time > 0.4;
            if (st.fingerprint != hash && settled) {
                if (!task_system::task_is_running(st.compute)) {
                    st.fingerprint = hash;
                    st.compute = rama_rep_compute_density(&st, angles, (uint32_t)num_frames, stride, sigma);
                } else {
                    task_system::task_interrupt(st.compute);
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
            constexpr const char* option_labels[3] = { "IsoLevels", "IsoLines", "Colormap" };

            constexpr const float min_ext = -180.0f;
            constexpr const float max_ext = 180.0f;
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

            PlotSubplot& sp = plot[0];

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Properties")) {
                    series_vector_properties_menu(&state, sp, 2);
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Layers")) {

                    ImGui::Checkbox("Reference", &show_reference);
                    if (show_reference) {
                        if (ImGui::BeginCombo("Display", option_labels[ref_mode])) {
                            if (ImGui::Selectable(option_labels[IsoLevels], ref_mode == IsoLevels)) ref_mode = IsoLevels;
                            if (ImGui::Selectable(option_labels[IsoLines],  ref_mode == IsoLines))  ref_mode = IsoLines;
                            if (ImGui::Selectable(option_labels[Colormap],  ref_mode == Colormap))  ref_mode = Colormap;
                            ImGui::EndCombo();
                        }
                        if (ref_mode == Colormap) {
                            ImPlot::ColormapSelection("Color Map", &ref_colormap);
                        } else if (ref_mode == IsoLines) {
                            ImGui::ColorEdit4Minimal("Line Color", &ref_line_color.x);
                        }
                        ImGui::SliderFloat("Opacity", &ref_alpha, 0.0f, 1.0f);
                    }
                    ImGui::SeparatorText("Series");
                    series_vector_layers_menu(&state, sp, 2);
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
                if (ImGui::MenuItem("Export", nullptr, false, sp.count > 0)) {
                    property_export_open(&state, PropertyExportKind_Ramachandran, plot_first_series_key(&sp, 1));
                }
                ImGui::EndMenuBar();
            }

            if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
                state.animation.mode = state.animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
            }

            // The shortcut
            ImGui::SetNextItemWidth(-1);
            ImGui::InputQuery("##residue-filter", input, sizeof(input), input_valid);
            if (ImGui::IsItemHovered()) {
                if (input_valid) md_bitfield_copy(&state.selection.highlight_mask, &residue_mask);
                else if (error[0]) ImGui::SetTooltip("%s", error);
            }
            ImGui::SetNextItemWidth(MAX(60.0f, ImGui::GetContentRegionAvail().x - 145.0f));
            ImGui::InputText("##Property Name", property_name, sizeof(property_name));
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Name of the property.");
            ImGui::SameLine();
            ImGui::BeginDisabled(!input_valid);
            if (ImGui::Button("Create Property")) create_property(state);
            ImGui::EndDisabled();
            if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                ImGui::SetTooltip("Add '%s = ramachandran(...)' for this selection to the script and show it here.", property_name);
            }
            if (property_error[0]) ImGui::TextWrapped("%s", property_error);
            if (!input_valid && error[0]) ImGui::TextWrapped("%s", error);
            if (pending_name[0]) ImGui::TextDisabled("Waiting for '%s' to be evaluated...", pending_name);

            // The series, which this may add to or remove from
            series_vector_bar(&state, sp, 2, false);

            SeriesVectorView views[PLOT_MAX_SERIES_PER_SUBPLOT];
            bool resolved[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
            RamaSeries* states[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
            char plot_ids[PLOT_MAX_SERIES_PER_SUBPLOT][SERIES_PATH_CAP + 72];
            for (int j = 0; j < sp.count; ++j) {
                resolved[j] = series_resolve_vector(&views[j], &state, sp.series[j].key, 2);
                series_label(views[j].label, sizeof(views[j].label), &state, sp.series[j].key);
                snprintf(plot_ids[j], sizeof(plot_ids[j]), "%s##%u", views[j].label, sp.series[j].uid);
                if (resolved[j]) states[j] = find_state(sp.series[j].uid);
            }

            const auto& sys = state.mold.sys;
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
                value = deperiodize_ortho(value, 0, 360.0);
                return snprintf(buff, size, "%g%s", value, suffix);
            };

            ImVec2 plot_size = {0,0};
            SeriesAction action = SeriesAction_None;
            int action_idx = -1;

            ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(2,2));
            defer { ImPlot::PopStyleVar(); };

            if (ImPlot::BeginSubplots("##Ramaplots", plot_rows, plot_cols, ImVec2(-1, -1), subplotflags | ImPlotSubplotFlags_ShareItems)) {
                for (int plot_idx = plot_offset; plot_idx < plot_offset + num_plots; ++plot_idx) {
                    if (ImPlot::BeginPlot(plot_labels[plot_idx], ImVec2(), flags)) {
                        ImPlotPoint view_mid = { (viewrect.X.Min + viewrect.X.Max) * 0.5, (viewrect.Y.Min + viewrect.Y.Max) * 0.5 };

                        ImPlot::SetupAxesLimits(min_ext, max_ext, min_ext, max_ext, ImPlotCond_Once);
                        ImPlot::SetupAxisLinks(ImAxis_X1, &viewrect.X.Min, &viewrect.X.Max);
                        ImPlot::SetupAxisLinks(ImAxis_Y1, &viewrect.Y.Min, &viewrect.Y.Max);
                        ImPlot::SetupAxes(x_lbl, y_lbl, axis_flags, axis_flags);

                        ImPlot::SetupAxisFormat(ImAxis_X1, formatter, (void*)x_lbl);
                        ImPlot::SetupAxisFormat(ImAxis_Y1, formatter, (void*)y_lbl);

                        ImPlot::SetupFinish();

                        viewrect = ImPlot::GetPlotLimits();

                        ImVec2 plot_min = ImPlot::PlotToPixels(viewrect.Min());
                        ImVec2 plot_max = ImPlot::PlotToPixels(viewrect.Max());
                        plot_size = { fabsf(plot_max.x - plot_min.x), fabsf(plot_max.y - plot_min.y) };

                        // The layers, below the points: the reference, then the densities of each series
                        {
                            ImPlot::PushPlotClipRect();
                            ImDrawList* dl = ImPlot::GetPlotDrawList();
                            if (show_reference) {
                                const uint32_t ref_tex = ref_mode == Colormap ? rama_data.ref.map_tex[plot_idx] : rama_data.ref.iso_tex[plot_idx];
                                dl->AddImage((ImTextureID)(intptr_t)ref_tex, plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, ref_alpha));
                            }
                            for (int j = 0; j < sp.count; ++j) {
                                if (!states[j]) continue;
                                const RamaSeries& st = *states[j];
                                const PlotSeries& s = sp.series[j];
                                if (!shows_density(s)) continue;
                                // The legend shows and hides a series as a whole
                                const ImPlotItem* item = ImPlot::GetItem(plot_ids[j]);
                                if (item && !item->Show) continue;
                                const rama_rep_t& rep = st.rep;
                                if (st.rep_generation != st.generation || rep.den_sum[plot_idx] <= 0) continue;
                                if (s.density_map) {
                                    dl->AddImage((ImTextureID)(intptr_t)rep.map_tex[plot_idx], plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, s.density_alpha));
                                }
                                if (s.iso_levels || s.iso_lines) {
                                    dl->AddImage((ImTextureID)(intptr_t)rep.iso_tex[plot_idx], plot_min, plot_max, {0,0}, {1,1}, ImColor(1.0f, 1.0f, 1.0f, s.density_alpha));
                                }
                            }
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
                        }

                        const double cut_d2 = fabs(ImPlot::PixelsToPlot(9,0).x - ImPlot::PixelsToPlot(0,0).x);
                        double min_d2 = DBL_MAX;
                        int hover_series = -1;
                        int64_t hover_row = -1;

                        if (hovered) {
                            md_bitfield_clear(highlight_mask);
                        }

                        double min_x = MIN(selection_rect.X.Min, selection_rect.X.Max);
                        double max_x = MAX(selection_rect.X.Min, selection_rect.X.Max);
                        double min_y = MIN(selection_rect.Y.Min, selection_rect.Y.Max);
                        double max_y = MAX(selection_rect.Y.Min, selection_rect.Y.Max);

                        double ref_x = (min_x + max_x) * 0.5;
                        double ref_y = (min_y + max_y) * 0.5;

                        mouse_coord.x = deperiodize_ortho(mouse_coord.x, ref_x, 360.0);
                        mouse_coord.y = deperiodize_ortho(mouse_coord.y, ref_y, 360.0);

                        // The points of the current frame, a series each
                        const md_backbone_angles_t* curr_angles[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                        for (int j = 0; j < sp.count; ++j) {
                            if (!states[j] || views[j].num_frames == 0 || !sp.series[j].show_points) continue;
                            const uint32_t frame = MIN((uint32_t)MAX(0.0, state.animation.frame + 0.5), views[j].num_frames - 1);
                            curr_angles[j] = (const md_backbone_angles_t*)views[j].data + (size_t)frame * views[j].count;
                            const RamaSeries& st = *states[j];
                            const uint32_t* indices = st.type_indices[plot_idx];

                            for (size_t i = 0; i < md_array_size(indices); ++i) {
                                const uint32_t row = indices[i];
                                const md_backbone_angles_t bb_angles = curr_angles[j][row];
                                if (bb_angles.phi == 0 && bb_angles.psi == 0) continue;

                                ImPlotPoint coord = ImPlotPoint(RAD_TO_DEG(bb_angles.phi), RAD_TO_DEG(bb_angles.psi));
                                coord.x = deperiodize_ortho(coord.x, ref_x, 360.0);
                                coord.y = deperiodize_ortho(coord.y, ref_y, 360.0);

                                if (is_selecting[plot_idx] && st.mapped) {
                                    if (min_x <= coord.x && coord.x <= max_x && min_y <= coord.y && coord.y <= max_y) {
                                        md_component_idx_t comp_idx = sys.protein_backbone.segment.comp_idx[st.segments[row]];
                                        md_urange_t range = md_component_atom_range(&sys.component, comp_idx);
                                        if (range.beg != range.end) {
                                            modify_field(highlight_mask, range, SelectionOperator::Or);
                                        }
                                    }
                                }
                                if (hovered) {
                                    double d2 = (coord.x - mouse_coord.x) * (coord.x - mouse_coord.x) + (coord.y - mouse_coord.y) * (coord.y - mouse_coord.y);
                                    if (d2 < cut_d2 && d2 < min_d2) {
                                        min_d2 = d2;
                                        hover_series = j;
                                        hover_row = row;
                                    }
                                }
                            }
                        }

                        if (is_selecting[plot_idx]) {
                            grow_mask_by_selection_granularity(highlight_mask, state.selection.granularity, sys);
                        }

                        if (hover_series != -1) {
                            const RamaSeries& st = *states[hover_series];
                            if (st.mapped) {
                                md_component_idx_t comp_idx = sys.protein_backbone.segment.comp_idx[st.segments[hover_row]];
                                md_urange_t range = md_component_atom_range(&sys.component, comp_idx);
                                if (range.beg != range.end) {
                                    modify_field(highlight_mask, range, SelectionOperator::Or);
                                    grow_mask_by_selection_granularity(highlight_mask, state.selection.granularity, sys);
                                }
                                str_t lbl = md_component_name(&sys.component, comp_idx);
                                md_sequence_id_t seq_id = md_component_seq_id(&sys.component, comp_idx);
                                ImGui::SetTooltip("res[%d]: %.*s %d", comp_idx + 1, (int)lbl.len, lbl.ptr, seq_id);
                            } else {
                                ImGui::SetTooltip("%s [%d]", views[hover_series].label, (int)hover_row + 1);
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

                            double x = deperiodize_ortho(RAD_TO_DEG(data->coords[idx].x), data->view_center.x, 360.0);
                            double y = deperiodize_ortho(RAD_TO_DEG(data->coords[idx].y), data->view_center.y, 360.0);
                            return { x, y };
                        };

                        for (int j = 0; j < sp.count; ++j) {
                            series_sync_visibility(sp.series[j], plot_ids[j]);
                            if (!curr_angles[j]) {
                                // Still an entry of the legend
                                if (states[j]) ImPlot::PlotDummy(plot_ids[j]);
                                continue;
                            }
                            const RamaSeries& st = *states[j];
                            const PlotSeries& s = sp.series[j];
                            const uint32_t* indices = st.type_indices[plot_idx];

                            uint32_t* selection_indices = 0;
                            uint32_t* highlight_indices = 0;
                            if (st.mapped) {
                                for (uint32_t i = 0; i < (uint32_t)md_array_size(indices); ++i) {
                                    const uint32_t row = indices[i];
                                    const md_backbone_angles_t bb_angles = curr_angles[j][row];
                                    if (bb_angles.phi == 0 && bb_angles.psi == 0) continue;

                                    const int64_t atom_idx = sys.protein_backbone.segment.atoms[st.segments[row]].ca;
                                    if (md_bitfield_test_bit(highlight_mask, atom_idx)) {
                                        md_array_push(highlight_indices, row, state.allocator.frame);
                                    }
                                    if (md_bitfield_test_bit(selection_mask, atom_idx)) {
                                        md_array_push(selection_indices, row, state.allocator.frame);
                                    }
                                }
                            }

                            const char* id = plot_ids[j];

                            if (md_array_size(indices) > 0) {
                                UserData user_data = { (const vec2_t*)(curr_angles[j]), indices, view_mid };
                                ImPlot::SetNextMarkerStyle(s.marker, s.marker_size, s.color, 1.0f, s.color);
                                ImPlot::PlotScatterG(id, index_getter, &user_data, (int)md_array_size(indices));
                            }

                            if (md_array_size(selection_indices) > 0) {
                                UserData user_data = { (const vec2_t*)(curr_angles[j]), selection_indices, view_mid };
                                ImPlot::SetNextMarkerStyle(s.marker, style.interaction_size, style.selection_fill, style.interaction_weight, style.selection_outline);
                                ImPlot::PlotScatterG("##Selection", index_getter, &user_data, (int)md_array_size(selection_indices));
                            }

                            if (md_array_size(highlight_indices) > 0) {
                                UserData user_data = { (const vec2_t*)(curr_angles[j]), highlight_indices, view_mid };
                                ImPlot::SetNextMarkerStyle(s.marker, style.interaction_size, style.highlight_fill, style.interaction_weight, style.highlight_outline);
                                ImPlot::PlotScatterG("##Highlight", index_getter, &user_data, (int)md_array_size(highlight_indices));
                            }
                        }

                        // The legend entries edit their series like the chips do
                        for (int j = 0; j < sp.count; ++j) {
                            if (ImPlot::BeginLegendPopup(plot_ids[j])) {
                                const SeriesAction a = series_vector_editor(&state, sp.series[j], 2, false);
                                if (a != SeriesAction_None) { action = a; action_idx = j; }
                                ImPlot::EndLegendPopup();
                            }
                        }

                        if (ImGui::IsMouseReleased(0) || ImGui::IsMouseReleased(1)) {
                            if (is_selecting[plot_idx]) {
                                op = ImGui::IsMouseReleased(0) ? SelectionOperator::Or : SelectionOperator::AndNot;
                                is_selecting[plot_idx] = false;
                                auto size = selection_rect.Size();
                                if (size.x == 0 && size.y == 0) {
                                    if (hover_series != -1) {
                                        modify_field(selection_mask, highlight_mask, op);
                                    } else {
                                        if (ImGui::IsMouseReleased(1)) {
                                            md_bitfield_clear(selection_mask);
                                        }
                                    }
                                }
                                else {
                                    // Commit whatever is in the highlight mask
                                    modify_field(selection_mask, highlight_mask, op);
                                }
                            } 
                        }

                        if (ImPlot::FitThisFrame()) {
                            ImPlotPlot* plot = ImPlot::GetCurrentPlot();
                            plot->Axes[ImAxis_X1].FitExtents = ImPlotRange(-180, 180);
                            plot->Axes[ImAxis_Y1].FitExtents = ImPlotRange(-180, 180);
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

            // Renders the maps of one density into its representation: contours and the bands between them, and a colormap.
            // The sums relative to the reference scale the iso values.
            auto render_layer = [&](rama_rep_t* rep, bool fill, bool lines, bool map, ImPlotColormap cmap, const float scl[4], ImVec4 color, bool reference_levels, int levels, const float* peak) {
                if (fill || lines) {
                    // The reference has its fixed probabilities. A series has contours at evenly spaced fractions of
                    // the peak of each type, the first value being the empty level.
                    const uint32_t count = reference_levels ? 3 : (uint32_t)CLAMP(levels, 1, 7) + 1;
                    uint32_t fill_cols[4][8] = {0};
                    uint32_t line_cols[4][8] = {0};
                    float    values[4][8] = {0};

                    for (uint32_t j = 0; j < 4; ++j) {
                        for (uint32_t i = 0; i < count; ++i) {
                            if (reference_levels) {
                                values[j][i] = ref_iso_values[j][i] * MAX(scl[j], FLT_EPSILON);
                            } else {
                                // A type without samples has no contours, and still needs values which differ
                                values[j][i] = peak[j] > 0.0f ? peak[j] * (float)i / (float)count : (float)i;
                            }
                        }
                    }

                    // The contours and the bands are tones of one colour, from light to full with the level
                    ImVec4 base  = color;
                    ImVec4 light = color;
                    float h, sat, val;
                    ImGui::ColorConvertRGBtoHSV(base.x, base.y, base.z, h, sat, val);
                    ImGui::ColorConvertHSVtoRGB(h, sat * 0.5f, val, light.x, light.y, light.z);
                    const uint32_t ref_line_tones[3] = {
                        0,
                        ImGui::ColorConvertFloat4ToU32(base),
                        ImGui::ColorConvertFloat4ToU32(light),
                    };

                    for (uint32_t i = 1; i < count; ++i) {
                        const float t = reference_levels ? 1.0f : (float)i / (float)(count - 1);
                        const ImVec4 tone = ImVec4(light.x + (base.x - light.x) * t, light.y + (base.y - light.y) * t,
                                                   light.z + (base.z - light.z) * t, base.w);
                        const uint32_t line_tone = reference_levels ? ref_line_tones[i] : ImGui::ColorConvertFloat4ToU32(tone);
                        const uint32_t fill_tone = ImGui::ColorConvertFloat4ToU32(ImVec4(tone.x, tone.y, tone.z, tone.w * (0.15f + 0.45f * t)));
                        for (int j = 0; j < 4; ++j) {
                            if (!reference_levels && !(peak[j] > 0.0f)) continue;
                            if (fill) fill_cols[j][i] = reference_levels ? ref_iso_level_colors[j][i] : fill_tone;
                            // Bands alone get contours of their own colour, which keeps their edges clean
                            line_cols[j][i] = lines ? line_tone : fill_cols[j][i];
                        }
                    }

                    rama_isomap_t maps[4] = {
                        { .values = values[0], .level_colors = fill_cols[0], .contour_colors = line_cols[0], .count = count },
                        { .values = values[1], .level_colors = fill_cols[1], .contour_colors = line_cols[1], .count = count },
                        { .values = values[2], .level_colors = fill_cols[2], .contour_colors = line_cols[2], .count = count },
                        { .values = values[3], .level_colors = fill_cols[3], .contour_colors = line_cols[3], .count = count },
                    };

                    rama_rep_render_iso(rep, viewport.elem, maps, max_dim);
                }
                if (map) {
                    // Fill in colors based on active color_map
                    uint32_t colors[32] = {0};
                    uint32_t num_colors = ImPlot::GetColormapSize(cmap);
                    for (uint32_t j = 0; j < num_colors; ++j) {
                        colors[j] = ImPlot::GetColormapColorU32(j, cmap);
                    }
                    // Mask off the first color to be transparent
                    colors[0] = colors[0] & 0x00FFFFFFU;

                    rama_colormap_t maps[4] = {
                        { .colors = colors, .count = num_colors, .min_value = 0, .max_value = 0.5f * scl[0] },
                        { .colors = colors, .count = num_colors, .min_value = 0, .max_value = 0.5f * scl[1] },
                        { .colors = colors, .count = num_colors, .min_value = 0, .max_value = 0.5f * scl[2] },
                        { .colors = colors, .count = num_colors, .min_value = 0, .max_value = 0.5f * scl[3] },
                    };
                    rama_rep_render_map(rep, viewport.elem, maps, max_dim);
                }
            };

            const float ref_scl[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
            if (show_reference) {
                render_layer(&rama_data.ref, ref_mode == IsoLevels, ref_mode == IsoLines, ref_mode == Colormap, ref_colormap, ref_scl, ref_line_color, true, 2, nullptr);
            }

            const float* ref_sum = rama_data.ref.den_sum;
            for (int j = 0; j < sp.count; ++j) {
                if (!states[j]) continue;
                RamaSeries& st = *states[j];
                const PlotSeries& s = sp.series[j];
                if (!shows_density(s)) continue;
                rama_rep_t* rep = &st.rep;
                if (st.rep_generation != st.generation) continue;
                const float scl[4] = {
                    rep->den_sum[0] / ref_sum[0], rep->den_sum[1] / ref_sum[1], rep->den_sum[2] / ref_sum[2], rep->den_sum[3] / ref_sum[3],
                };
                render_layer(rep, s.iso_levels, s.iso_lines, s.density_map, s.density_colormap, scl, s.color, false, s.iso_count, rep->den_max);
            }

            POP_GPU_SECTION();

            ImPlot::MapInputDefault();

            series_visibility_done(sp);
            series_apply_action(&state, sp, action_idx, action);
            series_vector_drag_out(sp);
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
        md_temp_scope_t temp_scope = md_temp_begin();
        defer { md_temp_end(temp_scope); };
        float* density_map = (float*)md_temp_alloc(temp_scope, mem_size);
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

    task_system::ID rama_rep_compute_density(RamaSeries* st, const md_backbone_angles_t* angles, uint32_t num_frames, uint32_t frame_stride, float sigma) {
        if (num_frames == 0 || frame_stride == 0) return 0;
        const uint32_t frame_beg = 0;
        const uint32_t frame_end = num_frames;

        struct UserData {
            uint64_t alloc_size;
            vec4_t* density_tex;
            RamaSeries* st;
            const md_backbone_angles_t* angles;     // The frames [frame_beg, frame_end) of the property, copied: it outlives nothing
            const uint32_t* type_indices[4];
            uint32_t frame_beg;
            uint32_t frame_end;
            uint32_t frame_stride;
            float sigma;
            md_allocator_i* alloc;
            uint32_t complete;
            float density_sum[4];
            float density_max[4];
            uint64_t generation;
            uint64_t fingerprint;
        };

        md_allocator_i* alloc = md_get_heap_allocator();

        const uint64_t tex_size = sizeof(vec4_t) * density_tex_dim * density_tex_dim;
        const uint64_t angles_size = sizeof(md_backbone_angles_t) * (uint64_t)(frame_end - frame_beg) * frame_stride;
        const uint64_t alloc_size = sizeof(UserData) + tex_size + angles_size + 2 * alignof(vec4_t);
        UserData* user_data = (UserData*)md_alloc(alloc, alloc_size);
        vec4_t* density_tex = (vec4_t*)NEXT_ALIGNED_ADDRESS(user_data + 1, alignof(vec4_t));
        memset(density_tex, 0, tex_size);
        md_backbone_angles_t* angle_copy = (md_backbone_angles_t*)NEXT_ALIGNED_ADDRESS((char*)density_tex + tex_size, alignof(vec4_t));
        // The property lives in the evaluation of the script, which a recompilation frees
        memcpy(angle_copy, angles + (size_t)frame_beg * frame_stride, angles_size);

        user_data->alloc_size = alloc_size;
        user_data->density_tex = density_tex;
        user_data->st = st;
        user_data->angles = angle_copy;
        // The worker owns its snapshot of the rows of each type
        for (int i = 0; i < 4; ++i) {
            user_data->type_indices[i] = md_array_create(uint32_t, md_array_size(st->type_indices[i]), alloc);
            if (md_array_size(st->type_indices[i])) {
                MEMCPY((void*)user_data->type_indices[i], st->type_indices[i], md_array_bytes(st->type_indices[i]));
            }
        }
        user_data->frame_beg = frame_beg;
        user_data->frame_end = frame_end;
        user_data->frame_stride = frame_stride;
        user_data->sigma = sigma;
        user_data->alloc = alloc;
        user_data->complete = 0;
        user_data->generation = st->generation;
        user_data->fingerprint = st->fingerprint;

        task_system::ID async_task = task_system::create_pool_task(STR_LIT("Rama density"), [data = user_data]() {
            const float angle_to_coord_scale = 1.0f / (2.0f * PI);
            const float angle_to_coord_offset = 0.5f;

            const uint32_t num_frames = data->frame_end - data->frame_beg;
            const uint32_t frame_stride = data->frame_stride;
            const md_backbone_angles_t* angles = data->angles;

            double sum[4] = {0,0,0,0};

            for (uint32_t f = 0; f < num_frames; ++f) {
                for (uint32_t c = 0; c < 4; ++c) {
                    const uint32_t* indices = data->type_indices[c];
                    const uint32_t num_indices = (uint32_t)md_array_size(data->type_indices[c]);
                    vec4_t val = {0, 0, 0, 0};
                    val.elem[c] = 1.0f;
                    if (num_indices) {
                        for (uint32_t i = 0; i < num_indices; ++i) {
                            const size_t idx = (size_t)f * frame_stride + indices[i];
                            if (!std::isfinite(angles[idx].phi) || !std::isfinite(angles[idx].psi)) continue;
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

            for (int i = 0; i < 4; ++i) data->density_sum[i] = (float)sum[i];
            float peak[4] = {0, 0, 0, 0};
            for (size_t i = 0; i < (size_t)density_tex_dim * density_tex_dim; ++i) {
                for (int c = 0; c < 4; ++c) peak[c] = MAX(peak[c], data->density_tex[i].elem[c]);
            }
            for (int i = 0; i < 4; ++i) data->density_max[i] = peak[i];
            data->complete = 1;
        });

        task_system::ID main_task = task_system::create_main_task(STR_LIT("##Update rama texture"), [data = user_data]() {
            RamaSeries* st = data->st;
            const uint64_t fingerprint = st->requested;
            if (data->complete && data->generation == st->generation && data->fingerprint == fingerprint) {
                rama_rep_t* rep = &st->rep;
                for (int i = 0; i < 4; ++i) {
                    rep->den_sum[i] = data->density_sum[i];
                    rep->den_max[i] = data->density_max[i];
                }
                gl::set_texture_2D_data(rep->den_tex, 0, data->density_tex, GL_RGBA32F);
                st->rep_generation = data->generation;
            }
            for (int i = 0; i < 4; ++i) md_array_free(data->type_indices[i], data->alloc);
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
