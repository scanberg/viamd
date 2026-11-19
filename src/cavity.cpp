#include <cavity.h>
#include <gostshyp.h>

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/gbuffer.h>
#include <gfx/view_param.h>

#include <md_system.h>
#include <md_system_state.h>
#include <core/md_log.h>
#include <core/md_hash.h>
#include <core/md_vec_math.h>

#include <imgui.h>
#include <implot.h>

#include <shaders.inl>

#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>

static const float GLYPH_MIN_PX = 1.5f;     // pixels: the radius a glyph is never thinner than, to be seen and hovered
static const double REF_PERCENTILE = 0.99;  // the value a quantity that is not a length maps to its scale at

struct CavityData {
    uint64_t source_hash = 0;
    GostshypTesserae t;

    // Where the cavity is: the rigid motion from the file's frame to the world's (gostshyp_place), and a hash of the
    // positions it was made for
    GostshypPlacement placement;
    uint64_t placement_hash = 0;

    // The sizes and the colours, and what they were made of
    uint64_t glyph_hash = 0;
    ColorScaleSpan span = {};

    // GPU
    GLuint buf[3] = {};             // tess_a (point), tess_b (normal, colour value), glyph (radius, height)
    GLuint tex[3] = {};
    GLuint colormap_tex = 0;
    int    colormap = -1;
    uint32_t version = 0;
};

enum { BUF_TESS_A, BUF_TESS_B, BUF_GLYPH, BUF_COUNT };

bool cavity_available(const md_system_t& sys) {
    return gostshyp_available(sys);
}

// ---------------------------------------------------------------------------
// The quantities
// ---------------------------------------------------------------------------

static double attribute_value(const GostshypTesserae& t, CavityAttribute a, uint32_t j) {
    switch (a) {
    case CavityAttribute::Area:              return t.area[j];
    case CavityAttribute::DecayLength:       return t.lambda[j];
    case CavityAttribute::WallEnergyPerArea: return gostshyp_energy_per_area(t, j);
    case CavityAttribute::WallEnergy:        return gostshyp_energy(t, j);
    case CavityAttribute::WallDensity:       return gostshyp_wall_density(t, j);
    case CavityAttribute::Amplitude:         return t.amp[j];
    case CavityAttribute::AmplitudePerArea:  return t.area[j] > 0 ? t.amp[j] / t.area[j] : 0.0;
    case CavityAttribute::Switching:         return t.switching[j];
    case CavityAttribute::LebedevWeight:     return t.lebedev[j];
    case CavityAttribute::Constant:
    default:                                 return 1.0;
    }
}

double cavity_attribute_value(const CavityRepresentation& cav, CavityAttribute a, uint32_t tessera) {
    const CavityData* d = cav.data;
    if (!d || tessera >= d->t.num) return 0.0;
    return attribute_value(d->t, a, tessera);
}

// The q-quantile of the quantity's magnitudes over the tesserae
static double quantile_magnitude(const GostshypTesserae& t, CavityAttribute a, double q) {
    std::vector<double> v(t.num);
    for (uint32_t j = 0; j < t.num; ++j) v[j] = fabs(attribute_value(t, a, j));
    if (v.empty()) return 0.0;
    const size_t k = MD_MIN(v.size() - 1, (size_t)(q * (double)v.size()));
    std::nth_element(v.begin(), v.begin() + k, v.end());
    return v[k];
}

// A size (Angstrom) of each tessera from a mapping: a length as it is, times the scale; another quantity over its 99th
// percentile, at the scale in Angstrom - by its square root for a radius (the cross section goes with the value),
// linearly for a height; a constant, the scale
static void sizes(std::vector<float>* out, const GostshypTesserae& t, const CavityMapping& m, bool radius) {
    out->assign(t.num, 0.0f);
    if (m.attribute == CavityAttribute::Constant) {
        for (uint32_t j = 0; j < t.num; ++j) (*out)[j] = m.scale;
        return;
    }
    if (cavity_attribute_is_length(m.attribute)) {
        for (uint32_t j = 0; j < t.num; ++j) {
            const double v = fmax(attribute_value(t, m.attribute, j), 0.0);
            (*out)[j] = (float)(m.scale * (m.attribute == CavityAttribute::Area ? sqrt(v / M_PI) : v));
        }
        return;
    }
    const double ref = quantile_magnitude(t, m.attribute, REF_PERCENTILE);
    for (uint32_t j = 0; j < t.num; ++j) {
        const double x = ref > 0 ? fabs(attribute_value(t, m.attribute, j)) / ref : 0.0;
        (*out)[j] = (float)(m.scale * (radius ? sqrt(x) : x));
    }
}

// ---------------------------------------------------------------------------
// GPU buffers
// ---------------------------------------------------------------------------

static void upload_buffer(CavityData* d, int i, GLenum format, const void* data, size_t bytes) {
    if (!d->buf[i]) glGenBuffers(1, &d->buf[i]);
    if (!d->tex[i]) glGenTextures(1, &d->tex[i]);
    glBindBuffer(GL_TEXTURE_BUFFER, d->buf[i]);
    // A buffer texture of nothing is not allowed to be bound: one element at least
    static const uint32_t zero[4] = {};
    glBufferData(GL_TEXTURE_BUFFER, bytes ? bytes : sizeof(zero), bytes ? data : zero, GL_STATIC_DRAW);
    glBindTexture(GL_TEXTURE_BUFFER, d->tex[i]);
    glTexBuffer(GL_TEXTURE_BUFFER, format, d->buf[i]);
    glBindTexture(GL_TEXTURE_BUFFER, 0);
    glBindBuffer(GL_TEXTURE_BUFFER, 0);
}

static void upload_glyphs(CavityData* d, const CavityRepresentation& cav) {
    const GostshypTesserae& t = d->t;
    const uint32_t n = t.num;
    std::vector<float> r, h;
    sizes(&r, t, cav.radius, true);
    sizes(&h, t, cav.height, false);
    std::vector<float> a(4 * n), b(4 * n), g(4 * n, 0.0f);
    for (uint32_t j = 0; j < n; ++j) {
        for (int k = 0; k < 3; ++k) {
            a[4 * j + k] = (float)t.pos[3 * j + k];
            b[4 * j + k] = (float)t.nrm[3 * j + k];
        }
        a[4 * j + 3] = 1.0f;
        b[4 * j + 3] = (float)attribute_value(t, cav.color, j);
        g[4 * j + 0] = r[j];
        g[4 * j + 1] = h[j];
    }
    upload_buffer(d, BUF_TESS_A, GL_RGBA32F, a.data(), a.size() * sizeof(float));
    upload_buffer(d, BUF_TESS_B, GL_RGBA32F, b.data(), b.size() * sizeof(float));
    upload_buffer(d, BUF_GLYPH,  GL_RGBA32F, g.data(), g.size() * sizeof(float));

    // The spread of the colour's quantity, for its scale
    d->span = {};
    if (cav.color != CavityAttribute::Constant && n > 0) {
        std::vector<double> v(n);
        for (uint32_t j = 0; j < n; ++j) v[j] = attribute_value(t, cav.color, j);
        std::sort(v.begin(), v.end());
        d->span.valid = true;
        d->span.min = (float)v.front();
        d->span.max = (float)v.back();
        d->span.lo  = (float)v[(size_t)(0.01 * (n - 1))];
        d->span.hi  = (float)v[(size_t)(0.99 * (n - 1))];
    }
}

static void free_gpu(CavityData* d) {
    for (int i = 0; i < BUF_COUNT; ++i) {
        if (d->tex[i]) glDeleteTextures(1, &d->tex[i]);
        if (d->buf[i]) glDeleteBuffers(1, &d->buf[i]);
        d->tex[i] = d->buf[i] = 0;
    }
    if (d->colormap_tex) glDeleteTextures(1, &d->colormap_tex);
    d->colormap_tex = 0;
    d->colormap = -1;
}

static GLuint colormap_texture(CavityData* d, int colormap) {
    if (d->colormap_tex && d->colormap == colormap) return d->colormap_tex;
    uint32_t pixels[256];
    for (int i = 0; i < 256; ++i) {
        ImVec4 col = ImPlot::SampleColormap((float)i / 255.0f, colormap);
        col.w = 1.0f;
        pixels[i] = ImGui::ColorConvertFloat4ToU32(col);
    }
    if (!d->colormap_tex) {
        gl::init_texture_2D(&d->colormap_tex, 256, 1, GL_RGBA8);
        glBindTexture(GL_TEXTURE_2D, d->colormap_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    gl::set_texture_2D_data(d->colormap_tex, 0, pixels, GL_RGBA8);
    d->colormap = colormap;
    return d->colormap_tex;
}

// The frame of the file to the world (Angstrom): the placement
static md_f32x4x4_t file_to_world(const CavityData& d) {
    md_f32x4x4_t M = {};
    for (int c = 0; c < 3; ++c)
        for (int r = 0; r < 3; ++r) M.elem[c][r] = (float)d.placement.R[3 * c + r];
    for (int r = 0; r < 3; ++r) M.elem[3][r] = (float)d.placement.t[r];
    M.elem[3][3] = 1.0f;
    return M;
}

// ---------------------------------------------------------------------------
// Update
// ---------------------------------------------------------------------------

static uint64_t glyph_hash(const CavityRepresentation& cav) {
    uint64_t h = md_hash64(&cav.radius, sizeof(cav.radius), 0);
    h = md_hash64(&cav.height, sizeof(cav.height), h);
    h = md_hash64(&cav.color, sizeof(cav.color), h);
    return h;
}

bool cavity_update(CavityRepresentation* cav, const md_system_t& sys, const md_system_state_t& state) {
    if (!cavity_available(sys)) return false;
    if (!cav->data) cav->data = new CavityData();
    CavityData* d = cav->data;

    const uint64_t h = gostshyp_source_hash(sys);
    if (h != d->source_hash) {
        free_gpu(d);
        *d = CavityData{};
        d->source_hash = h;
        if (!gostshyp_read(&d->t, sys)) return false;
    }
    if (d->t.num == 0) return false;

    // Placed on the QM atoms as they are
    if (state.xyz && state.num_atoms) {
        uint64_t ph = 0;
        for (uint32_t i = 0; i < d->t.num_atoms; ++i) {
            if (d->t.sys_idx[i] < state.num_atoms) ph = md_hash64(&state.xyz[d->t.sys_idx[i]], sizeof(md_f32x3_t), ph);
        }
        if (ph != d->placement_hash) {
            d->placement_hash = ph;
            gostshyp_place(&d->placement, d->t, state.xyz, state.num_atoms, d->t.sys_idx.data(), 1.0);
            d->version++;
        }
    }

    const uint64_t gh = glyph_hash(*cav);
    if (gh != d->glyph_hash || !d->tex[BUF_GLYPH]) {
        d->glyph_hash = gh;
        upload_glyphs(d, *cav);
        d->version++;
    }
    color_scale_update_range(&cav->scale, d->span);
    return true;
}

void cavity_invalidate(CavityRepresentation* cav) {
    if (cav->data) cav->data->source_hash = 0;
}

void cavity_free(CavityRepresentation* cav) {
    if (!cav->data) return;
    free_gpu(cav->data);
    delete cav->data;
    cav->data = nullptr;
}

ColorScaleSpan cavity_span(const CavityRepresentation& cav) {
    return cav.data ? cav.data->span : ColorScaleSpan{};
}

uint64_t cavity_hash(const CavityRepresentation& cav) {
    uint64_t h = glyph_hash(cav);
    h = md_hash64(&cav.shape, sizeof(cav.shape), h);
    h = md_hash64(&cav.scale, sizeof(cav.scale), h);
    h = md_hash64(&cav.unlit, sizeof(cav.unlit), h);
    if (cav.data) h = md_hash64(&cav.data->version, sizeof(cav.data->version), h);
    return h;
}

uint32_t cavity_tessera_count(const CavityRepresentation& cav) { return cav.data ? cav.data->t.num : 0; }

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------

static struct {
    GLuint glyph = 0;
    GLuint vao = 0;
    bool   failed = false;
} gl_state;

static bool ensure_program() {
    if (gl_state.failed) return false;
    if (!gl_state.glyph) {
        GLuint vs = gl::compile_shader_from_source({(const char*)cavity_glyph_vert, cavity_glyph_vert_size}, GL_VERTEX_SHADER);
        GLuint fs = gl::compile_shader_from_source({(const char*)cavity_glyph_frag, cavity_glyph_frag_size}, GL_FRAGMENT_SHADER);
        if (vs && fs) {
            gl_state.glyph = glCreateProgram();
            const GLuint shaders[] = {vs, fs};
            if (!gl::attach_link_detach(gl_state.glyph, shaders, 2)) {
                MD_LOG_ERROR("Cavity: the glyph program did not link");
                glDeleteProgram(gl_state.glyph);
                gl_state.glyph = 0;
            }
        }
        if (vs) glDeleteShader(vs);
        if (fs) glDeleteShader(fs);
        if (!gl_state.glyph) {
            gl_state.failed = true;
            return false;
        }
    }
    if (!gl_state.vao) glGenVertexArrays(1, &gl_state.vao);
    return true;
}

static void jitter_uv(float out[2], const md_f32x4x4_t& P) {
    // As md_gl reads it off the projection
    if (P.elem[2][3] == 0.0f) {
        out[0] = -P.elem[3][0] * 0.5f;
        out[1] = -P.elem[3][1] * 0.5f;
    } else {
        out[0] = P.elem[2][0] * 0.5f;
        out[1] = P.elem[2][1] * 0.5f;
    }
}

// The program of the glyphs with its buffers and uniforms, for the G-buffer (proj, jittered) or the marks (unjittered)
static void glyph_program(const CavityRepresentation& cav, const CavityData& d, const ViewParam& view, const md_f32x4x4_t& proj, bool marks) {
    const GLuint prog = gl_state.glyph;
    glUseProgram(prog);
    const char* names[BUF_COUNT] = {"u_tess_a", "u_tess_b", "u_glyph"};
    for (int i = 0; i < BUF_COUNT; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_BUFFER, d.tex[i]);
        glUniform1i(glGetUniformLocation(prog, names[i]), i);
    }
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, colormap_texture((CavityData*)&d, cav.scale.colormap));
    glUniform1i(glGetUniformLocation(prog, "u_colormap"), 3);
    // The depth of the marks on a unit of its own even when unused: two types of sampler on one unit fail the draw
    glUniform1i(glGetUniformLocation(prog, "u_tex_depth"), 8);

    const md_f32x4x4_t file_to_view = view.matrix.curr.view * file_to_world(d);
    const md_f32x4x4_t curr_view_to_prev_clip = view.matrix.prev.proj * view.matrix.prev.view * view.matrix.inv.view;
    float jit[4];
    jitter_uv(jit + 0, view.matrix.curr.proj);
    jitter_uv(jit + 2, view.matrix.prev.proj);
    const bool ortho = proj.elem[2][3] == 0.0f;
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    const float p11 = proj.elem[1][1];

    glUniformMatrix4fv(glGetUniformLocation(prog, "u_cavity_to_view"), 1, GL_FALSE, &file_to_view.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_view_to_clip"), 1, GL_FALSE, &proj.elem[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_curr_view_to_prev_clip"), 1, GL_FALSE, &curr_view_to_prev_clip.elem[0][0]);
    glUniform4fv(glGetUniformLocation(prog, "u_jitter_uv"), 1, jit);
    glUniform1i(glGetUniformLocation(prog, "u_ortho"), ortho ? 1 : 0);
    glUniform1i(glGetUniformLocation(prog, "u_only"), -1);
    glUniform1i(glGetUniformLocation(prog, "u_sphere"), cav.shape == CavityGlyphShape::Sphere ? 1 : 0);
    glUniform2f(glGetUniformLocation(prog, "u_range"), cav.scale.range_beg, cav.scale.range_end);
    glUniform1f(glGetUniformLocation(prog, "u_min_px"), GLYPH_MIN_PX);
    glUniform1f(glGetUniformLocation(prog, "u_pixel"), 2.0f / (p11 * (float)MD_MAX(vp[3], 1)));
    glUniform1i(glGetUniformLocation(prog, "u_color_mode"), cav.color == CavityAttribute::Constant ? 1 : 0);
    glUniform1i(glGetUniformLocation(prog, "u_unlit"), cav.unlit ? 1 : 0);
    glUniform1i(glGetUniformLocation(prog, "u_marks"), marks ? 1 : 0);
    glUniform1f(glGetUniformLocation(prog, "u_proj_11"), p11);
}

static void glyph_unbind() {
    for (int i = 0; i < BUF_COUNT; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_BUFFER, 0);
    }
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0);
    glUseProgram(0);
}

void cavity_draw_glyphs(const CavityRepresentation& cav, const ViewParam& view, uint32_t tessera_base) {
    const CavityData* d = cav.data;
    if (!d || d->t.num == 0 || !d->tex[BUF_GLYPH] || !ensure_program()) return;
    glyph_program(cav, *d, view, view.matrix.curr.proj, false);
    glUniform1ui(glGetUniformLocation(gl_state.glyph, "u_tessera_base"), tessera_base);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDisable(GL_CULL_FACE);
    glBindVertexArray(gl_state.vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, (GLsizei)d->t.num);
    glBindVertexArray(0);
    glEnable(GL_CULL_FACE);
    glyph_unbind();
}

void cavity_draw_marks(const CavityRepresentation& cav, uint32_t tessera, const GBuffer& gbuf, const ViewParam& view) {
    const CavityData* d = cav.data;
    if (!d || tessera >= d->t.num || !d->tex[BUF_GLYPH] || !ensure_program()) return;
    // The marks are of the image shown: unjittered, at its size (gbuffer_marks_begin bound the target and blending)
    glViewport(0, 0, gbuf.out_width, gbuf.out_height);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glyph_program(cav, *d, view, view.matrix.curr.proj_no_jitter, true);
    const GLuint prog = gl_state.glyph;
    glUniform1i(glGetUniformLocation(prog, "u_only"), (int)tessera);
    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, gbuf.tex.depth);
    const md_f32x4x4_t clip_to_view = md_f32x4x4_inverse(view.matrix.curr.proj);
    glUniformMatrix4fv(glGetUniformLocation(prog, "u_clip_to_view"), 1, GL_FALSE, &clip_to_view.elem[0][0]);
    glUniform2f(glGetUniformLocation(prog, "u_inv_size"), 1.0f / gbuf.out_width, 1.0f / gbuf.out_height);
    glBindVertexArray(gl_state.vao);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
    glBindVertexArray(0);
    glActiveTexture(GL_TEXTURE8);
    glBindTexture(GL_TEXTURE_2D, 0);
    glyph_unbind();
    glEnable(GL_CULL_FACE);
}

// ---------------------------------------------------------------------------
// What it is
// ---------------------------------------------------------------------------

bool cavity_tessera_info(CavityTesseraInfo* info, const CavityRepresentation& cav, uint32_t tessera) {
    const CavityData* d = cav.data;
    if (!d || tessera >= d->t.num) return false;
    const GostshypTesserae& t = d->t;
    *info = {};
    info->atom = t.atom[tessera];
    info->sys_atom = t.sys_idx[info->atom];
    info->area = (float)t.area[tessera];
    info->decay_length = (float)t.lambda[tessera];
    info->wall_energy = (float)gostshyp_energy_per_area(t, tessera);
    info->wall_density = (float)gostshyp_wall_density(t, tessera);
    info->amplitude = (float)t.amp[tessera];
    info->switching = (float)t.switching[tessera];
    return true;
}

bool cavity_summary(CavitySummary* out, const CavityRepresentation& cav) {
    const CavityData* d = cav.data;
    if (!d || d->t.num == 0) return false;
    const GostshypTesserae& t = d->t;
    *out = {};
    out->pressure_gpa = (float)t.pressure_gpa;
    out->num_tesserae = t.num;
    for (uint32_t i = 0; i < t.num_atoms; ++i) out->num_spheres += t.radius[i] > 0;
    double A = 0, E = 0, V = 0, buried = 0;
    for (uint32_t j = 0; j < t.num; ++j) {
        A += t.area[j];
        E += gostshyp_energy(t, j);
        V += t.area[j] * t.lambda[j];
        for (uint32_t k = 0; k < t.num_atoms; ++k) {
            const double Rk = t.radius[k];
            if (k == t.atom[j] || !(Rk > 0)) continue;
            const double dx = t.pos[3 * j] - t.qm_xyz[3 * k], dy = t.pos[3 * j + 1] - t.qm_xyz[3 * k + 1], dz = t.pos[3 * j + 2] - t.qm_xyz[3 * k + 2];
            if (dx * dx + dy * dy + dz * dz < Rk * Rk * 0.9999) {
                out->num_buried += 1;
                buried += t.area[j];
                break;
            }
        }
    }
    out->tessera_area = (float)A;
    out->wall_energy = (float)E;
    out->effective_volume = (float)V;
    out->buried_area = A > 0 ? (float)(buried / A) : 0.0f;
    return true;
}
