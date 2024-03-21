#include "volumerender_utils.h"

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/postprocessing_utils.h>
#include <color_utils.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <shaders.inl>

#define PUSH_GPU_SECTION(lbl)                                                                       \
    {                                                                                               \
        if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); \
    }
#define POP_GPU_SECTION()                       \
    {                                           \
        if (glPopDebugGroup) glPopDebugGroup(); \
    }

static constexpr str_t v_shader_src_fs_quad = STR_LIT(
    R"(
#version 150 core

out vec2 tc;

uniform vec2 u_tc_scl = vec2(1,1);

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
	tc = (gl_Position.xy * 0.5 + 0.5) * u_tc_scl;
}
)");

static constexpr str_t f_shader_src_median = STR_LIT(
R"(
#version 150 core

uniform sampler2D T;

out vec4 out_frag;

// Change these 2 defines to change precision
#define vec vec3

#define s2(a, b)				temp = a; a = min(a, b); b = max(temp, b);
#define mn3(a, b, c)			s2(a, b); s2(a, c);
#define mx3(a, b, c)			s2(b, c); s2(a, c);

#define mnmx3(a, b, c)			mx3(a, b, c); s2(a, b);                                   // 3 exchanges
#define mnmx4(a, b, c, d)		s2(a, b); s2(c, d); s2(a, c); s2(b, d);                   // 4 exchanges
#define mnmx5(a, b, c, d, e)	s2(a, b); s2(c, d); mn3(a, c, e); mx3(b, d, e);           // 6 exchanges
#define mnmx6(a, b, c, d, e, f) s2(a, d); s2(b, e); s2(c, f); mn3(a, b, c); mx3(d, e, f); // 7 exchanges

void main() {

  vec4 col = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2( 0,  0), 0);

  vec v[6];

  v[0] = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2(-1, -1), 0).rgb;
  v[1] = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2( 0, -1), 0).rgb;
  v[2] = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2( 1, -1), 0).rgb;
  v[3] = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2(-1,  0), 0).rgb;
  v[4] = col.rgb;
  v[5] = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2( 1,  0), 0).rgb;

  // Starting with a subset of size 6, remove the min and max each time
  vec temp;
  mnmx6(v[0], v[1], v[2], v[3], v[4], v[5]);

  v[5] = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2(-1,  1), 0).rgb;

  mnmx5(v[1], v[2], v[3], v[4], v[5]);

  v[5] = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2( 0,  1), 0).rgb;

  mnmx4(v[2], v[3], v[4], v[5]);

  v[5] = texelFetch(T, ivec2(gl_FragCoord.xy) + ivec2( 1,  1), 0).rgb;

  mnmx3(v[3], v[4], v[5]);
  out_frag = vec4(v[4], col.a);
}
)");

namespace volume {

static struct {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ubo = 0;
    GLuint fbo = 0;

    GLuint tex_entry  = 0;
    GLuint tex_exit   = 0;
    GLuint tex_result = 0;

    uint32_t width  = 0;
    uint32_t height = 0;

    struct {
        GLuint entry_exit = 0;
        GLuint entry_exit_depth = 0;
        GLuint dvr_only = 0;
        GLuint iso_only = 0;
        GLuint dvr_and_iso = 0;
        GLuint median = 0;
    } program;
} gl;

struct UniformData {
    mat4_t view_to_model_mat;
    mat4_t model_to_view_mat;
    mat4_t inv_proj_mat;
    mat4_t model_view_proj_mat;

    vec2_t inv_res;
    float density_scale;
    float _pad0;

    vec3_t clip_volume_min;
    float _pad1;
    vec3_t clip_volume_max;
    float time;

    vec3_t gradient_spacing_world_space;
    float _pad2;
    mat4_t gradient_spacing_tex_space;
};

void initialize() {
    GLuint v_shader_entry_exit          = gl::compile_shader_from_source({(const char*)entryexit_vert, entryexit_vert_size}, GL_VERTEX_SHADER);
    GLuint f_shader_entry_exit          = gl::compile_shader_from_source({(const char*)entryexit_frag, entryexit_frag_size}, GL_FRAGMENT_SHADER);
    GLuint f_shader_entry_exit_depth    = gl::compile_shader_from_source({(const char*)entryexit_frag, entryexit_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define SAMPLE_DEPTH"));

    GLuint v_shader_vol                 = gl::compile_shader_from_source(v_shader_src_fs_quad, GL_VERTEX_SHADER);
    GLuint f_shader_dvr_only            = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_DVR"));
    GLuint f_shader_iso_only            = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_ISO"));
    GLuint f_shader_dvr_and_iso         = gl::compile_shader_from_source({(const char*)raycaster_frag, raycaster_frag_size}, GL_FRAGMENT_SHADER, STR_LIT("#define INCLUDE_DVR\n#define INCLUDE_ISO"));

    defer {
        glDeleteShader(v_shader_vol);
        glDeleteShader(f_shader_entry_exit);
        glDeleteShader(f_shader_dvr_only);
        glDeleteShader(f_shader_iso_only);
        glDeleteShader(f_shader_dvr_and_iso);
    };

    if (v_shader_entry_exit == 0 || v_shader_vol == 0 || f_shader_entry_exit == 0|| f_shader_dvr_only == 0 || f_shader_iso_only == 0 || f_shader_dvr_and_iso == 0) {
        MD_LOG_ERROR("shader compilation failed, shader program for raycasting will not be updated");
        return;
    }

    if (!gl.program.entry_exit) gl.program.entry_exit = glCreateProgram();
    if (!gl.program.entry_exit_depth) gl.program.entry_exit_depth = glCreateProgram();
    if (!gl.program.dvr_only) gl.program.dvr_only = glCreateProgram();
    if (!gl.program.iso_only) gl.program.iso_only = glCreateProgram();
    if (!gl.program.dvr_and_iso) gl.program.dvr_and_iso = glCreateProgram();

    {
        const GLuint shaders[] = {v_shader_entry_exit, f_shader_entry_exit};
        gl::attach_link_detach(gl.program.entry_exit, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_entry_exit, f_shader_entry_exit_depth};
        gl::attach_link_detach(gl.program.entry_exit_depth, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_dvr_only};
        gl::attach_link_detach(gl.program.dvr_only, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_iso_only};
        gl::attach_link_detach(gl.program.iso_only, shaders, (int)ARRAY_SIZE(shaders));
    }
    {
        const GLuint shaders[] = {v_shader_vol, f_shader_dvr_and_iso};
        gl::attach_link_detach(gl.program.dvr_and_iso, shaders, (int)ARRAY_SIZE(shaders));
    }

    if (!gl.vbo) {
        // https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip
        constexpr uint8_t cube_strip[42] = {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1,
                                            0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
        glGenBuffers(1, &gl.vbo);
        glBindBuffer(GL_ARRAY_BUFFER, gl.vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cube_strip), cube_strip, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (!gl.vao) {
        glGenVertexArrays(1, &gl.vao);
        glBindVertexArray(gl.vao);
        glBindBuffer(GL_ARRAY_BUFFER, gl.vbo);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_UNSIGNED_BYTE, GL_FALSE, 0, (const GLvoid*)0);
        glBindVertexArray(0);
    }

    if (!gl.ubo) {
        glGenBuffers(1, &gl.ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(UniformData), 0, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
    }

    if (!gl.tex_entry) {
        glGenTextures(1, &gl.tex_entry);
    }

    if (!gl.tex_exit) {
        glGenTextures(1, &gl.tex_exit);
    }

    if (!gl.tex_result) {
        glGenTextures(1, &gl.tex_result);
    }

    if (!gl.fbo) {
        glGenFramebuffers(1, &gl.fbo);
    }
}

void shutdown() {}

mat4_t compute_model_to_world_matrix(vec3_t min_world_aabb, vec3_t max_world_aabb) {
    vec3_t ext = max_world_aabb - min_world_aabb;
    vec3_t off = min_world_aabb;
    return {
            ext.x, 0, 0, 0,
            0, ext.y, 0, 0,
            0, 0, ext.z, 0,
            off.x, off.y, off.z, 1,
        };
}

mat4_t compute_world_to_model_matrix(vec3_t min_world_aabb, vec3_t max_world_aabb) {
    vec3_t ext = max_world_aabb - min_world_aabb;
    vec3_t off = min_world_aabb;
    return {
            1.0f / ext.x, 0, 0, 0,
            0, 1.0f / ext.y, 0, 0,
            0, 0, 1.0f / ext.z, 0,
            -off.x, -off.y, -off.z, 1,
        };
}

mat4_t compute_model_to_texture_matrix(int dim_x, int dim_y, int dim_z) {
    (void)dim_x;
    (void)dim_y;
    (void)dim_z;
    return mat4_ident();
}

mat4_t compute_texture_to_model_matrix(int dim_x, int dim_y, int dim_z) {
    (void)dim_x;
    (void)dim_y;
    (void)dim_z;
    return mat4_ident();
}

void render_volume(const RenderDesc& desc) {
    if (!desc.direct_volume_rendering_enabled && !desc.isosurface_enabled) return;

    int    iso_count = CLAMP((int)desc.iso_surface.count, 0, 8);
    float  iso_values[8];
    vec4_t iso_colors[8];

    MEMCPY(iso_values, desc.iso_surface.values, iso_count * sizeof(float));
    MEMCPY(iso_colors, desc.iso_surface.colors, iso_count * sizeof(vec4_t));

    // Sort on iso value
    for (int i = 0; i < iso_count - 1; ++i) {
        for (int j = i + 1; j < iso_count; ++j) {
            if (iso_values[j] < iso_values[i]) {
                float  val_tmp = iso_values[i];
                vec4_t col_tmp = iso_colors[i];
                iso_values[i] = iso_values[j];
                iso_colors[i] = iso_colors[j];
                iso_values[j] = val_tmp;
                iso_colors[j] = col_tmp;
            }
        }
    }

    GLint bound_fbo;
    GLint bound_viewport[4];
    GLint bound_draw_buffer[8] = {0};
    GLint bound_draw_buffer_count = 0;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &bound_fbo);
    glGetIntegerv(GL_VIEWPORT, bound_viewport);
    for (int i = 0; i < 8; ++i) {
        glGetIntegerv(GL_DRAW_BUFFER0 + i, &bound_draw_buffer[i]);
        // @NOTE: Assume that its tightly packed and if we stumple upon a zero draw buffer index, we enterpret that as the 'end'
        if (bound_draw_buffer[i] == GL_NONE) {
            bound_draw_buffer_count = i;
            break;
        }
    }

    if (gl.width < desc.render_target.width ||
        gl.height < desc.render_target.height)
    {
        gl.width = desc.render_target.width;
        gl.height = desc.render_target.height;
        gl::init_texture_2D(&gl.tex_entry,  gl.width, gl.height, GL_RGB16);
        gl::init_texture_2D(&gl.tex_exit,   gl.width, gl.height, GL_RGB16);
        gl::init_texture_2D(&gl.tex_result, gl.width, gl.height, GL_RGBA8);
    }

    const mat4_t model_to_view_matrix = mat4_mul(desc.matrix.view, desc.matrix.model);

    static float time = 0.0f;
    time += 1.0f / 100.0f;
    if (time > 100.0) time -= 100.0f;
    if (!desc.temporal.enabled) {
        time = 0.0f;
    }

    UniformData data;
    data.view_to_model_mat = mat4_inverse(model_to_view_matrix);
    data.model_to_view_mat = model_to_view_matrix;
    data.inv_proj_mat      = desc.matrix.inv_proj;
    data.model_view_proj_mat = desc.matrix.proj * model_to_view_matrix;
    data.inv_res = {1.f / (float)(desc.render_target.width), 1.f / (float)(desc.render_target.height)};
    data.density_scale = desc.global_scaling.density;
    data.clip_volume_min = desc.clip_volume.min;
    data.clip_volume_max = desc.clip_volume.max;
    data.time = time;
    data.gradient_spacing_world_space = desc.voxel_spacing;
    data.gradient_spacing_tex_space = data.view_to_model_mat * mat4_scale(desc.voxel_spacing.x, desc.voxel_spacing.y, desc.voxel_spacing.z);

    glBindBuffer(GL_UNIFORM_BUFFER, gl.ubo);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UniformData), &data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    bool use_depth = desc.render_target.depth;

    if (use_depth) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, desc.render_target.depth);
    }

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, gl.ubo);
    glBindVertexArray(gl.vao);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tex_entry, 0);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gl.tex_exit,  0);
    const GLuint draw_bufs[] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
    glDrawBuffers(2, draw_bufs);

    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);

    glViewport(0, 0, desc.render_target.width, desc.render_target.height);

    glEnable(GL_CULL_FACE);
    {
        PUSH_GPU_SECTION("VOLUME ENTRY / EXIT");
        
        const GLuint prog = use_depth ? gl.program.entry_exit_depth : gl.program.entry_exit;
        const GLint uniform_block_index = glGetUniformBlockIndex(prog, "UniformData");
        const GLint uniform_loc_tex_depth = glGetUniformLocation(prog, "u_tex_depth");

        glUseProgram(prog);
        glUniform1i(uniform_loc_tex_depth, 0);
        glUniformBlockBinding(prog, uniform_block_index, 0);

        glCullFace(GL_FRONT);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 42);
        POP_GPU_SECTION()
    }

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl.tex_entry);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, gl.tex_exit);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, desc.texture.volume);

    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, desc.texture.transfer_function);

    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.tex_result, 0);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    if (desc.render_target.color) {
        ASSERT(glIsTexture(desc.render_target.color));
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, desc.render_target.color, 0);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
    } else {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
        glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
        glDrawBuffers(bound_draw_buffer_count, (GLenum*)bound_draw_buffer);
    }

    PUSH_GPU_SECTION("VOLUME RAYCASTING")
    {
        const GLuint vol_prog = desc.direct_volume_rendering_enabled ? (desc.isosurface_enabled ? gl.program.dvr_and_iso : gl.program.dvr_only) : gl.program.iso_only;

        const GLint uniform_block_index     = glGetUniformBlockIndex(vol_prog, "UniformData");
        const GLint uniform_loc_tex_entry   = glGetUniformLocation(vol_prog, "u_tex_entry");
        const GLint uniform_loc_tex_exit    = glGetUniformLocation(vol_prog, "u_tex_exit");
        const GLint uniform_loc_tex_tf      = glGetUniformLocation(vol_prog, "u_tex_tf");
        const GLint uniform_loc_tex_volume  = glGetUniformLocation(vol_prog, "u_tex_volume");
        const GLint uniform_loc_iso_values  = glGetUniformLocation(vol_prog, "u_iso.values");
        const GLint uniform_loc_iso_colors  = glGetUniformLocation(vol_prog, "u_iso.colors");
        const GLint uniform_loc_iso_count   = glGetUniformLocation(vol_prog, "u_iso.count");

        glUseProgram(vol_prog);

        glUniform1i(uniform_loc_tex_entry, 0);
        glUniform1i(uniform_loc_tex_exit,  1);
        glUniform1i(uniform_loc_tex_volume, 2);
        glUniform1i(uniform_loc_tex_tf, 3);
        glUniform1fv(uniform_loc_iso_values, (GLsizei)iso_count, (const float*)iso_values);
        glUniform4fv(uniform_loc_iso_colors, (GLsizei)iso_count, (const float*)iso_colors);
        glUniform1i(uniform_loc_iso_count, (int)iso_count);
        glUniformBlockBinding(vol_prog, uniform_block_index, 0);

        glDrawArrays(GL_TRIANGLES, 0, 3);

        glBindVertexArray(0);
        glUseProgram(0);
    }
    POP_GPU_SECTION()

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    if (!desc.render_target.color) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, bound_fbo);
        glViewport(bound_viewport[0], bound_viewport[1], bound_viewport[2], bound_viewport[3]);
        glDrawBuffers(bound_draw_buffer_count, (GLenum*)bound_draw_buffer);
    }

}

}  // namespace volume
