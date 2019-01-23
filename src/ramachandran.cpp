#include "ramachandran.h"
#include <core/gl.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <gfx/gl_utils.h>
#include "image.h"

namespace ramachandran {

// Accumulation texture data
constexpr int acc_width = 1024;
constexpr int acc_height = 1024;

static GLuint gui_tex = 0;
static GLuint seg_tex = 0;
static GLuint acc_tex = 0;
static GLuint col_tex = 0;

static GLuint coord_tex = 0;
static GLuint coord_buf = 0;
static GLuint fbo = 0;
static GLuint program = 0;

static GLuint vao = 0;
static GLuint ibo = 0;
static GLuint vbo = 0;

static GLint uniform_loc_coord_tex = -1;
static GLint uniform_loc_instance_offset = -1;
static GLint uniform_loc_radius = -1;
static GLint uniform_loc_color = -1;
static GLint uniform_loc_outline = -1;

static Image src_img = {};  // This is the unmodified source image (Ramachandran plot)
static Image gui_img = {};  // Visual representation which is used in gui-elements
static Image seg_img = {};  // Segmentation version of image (Blurred, to hide low res artifacts and for smoother transitions)
static Image col_img = {};  // Color version of image (Even more blurred, for smoother transitions of secondary structure colors)

GLuint get_accumulation_texture() { return acc_tex; }
GLuint get_gui_texture() { return gui_tex; }
GLuint get_segmentation_texture() { return seg_tex; }
GLuint get_color_texture() { return col_tex; }

const Image& get_gui_image() { return gui_img; }
const Image& get_segmentation_image() { return seg_img; }
const Image& get_color_image() { return col_img; }

// @NOTE: This should generate a quad with a certain size in texture coordinates
constexpr const char* v_shader_src = R"(
#version 150 core

uniform int u_instance_offset = 0;
uniform samplerBuffer u_tex_coord;
uniform float u_radius;
out vec2 uv;

void main() {
	int VID = gl_VertexID;
	int IID = gl_InstanceID + u_instance_offset;

	vec2 coord = texelFetch(u_tex_coord, IID).xy;
	uv = vec2(VID / 2, VID % 2) * 2.0 - 1.0; 

	gl_Position = vec4(coord * 2.0 - 1.0 + uv * u_radius, 0, 1);
}
)";

// @NOTE: Do some radial falloff based on uv coordinate
constexpr const char* f_shader_src = R"(
#version 150 core

uniform vec4 u_color;
uniform float u_outline;
in vec2 uv;
out vec4 out_frag;

float step_dist(float edge, float dist) {
	float factor = 1.0; // <-- value can be played around with a bit
	float mask = step(edge, dist);
	float step_w = factor * length(vec2(dFdx(mask), dFdy(mask)));
	return smoothstep(-step_w/2.0, step_w/2.0, mask);
}

void main() {
	float dist = sqrt(dot(uv, uv));
	if (u_outline > 0) {
		vec4 inner_color = u_color.rgba;
		vec4 rim_color = vec4(0,0,0,1);
		vec4 outer_color = vec4(0,0,0,0);
		out_frag = mix(inner_color, mix(rim_color, outer_color, step_dist(1.0, dist)), step_dist(1.0 - u_outline, dist));
	} else {
		float falloff = max(0, 1.0 - dist);
		out_frag = vec4(u_color.rgb, u_color.a * falloff);	
	}
}
)";

void init_map(Image* img, GLuint tex, const ColorMap& color_map, int blur_level) {
    ASSERT(src_img);
    ASSERT(img);
    ASSERT(tex);

    init_image(img, src_img);

    constexpr uint32 IN_BACKGROUND = 0xFFFFFFFF;
    constexpr uint32 IN_ALPHA_HIGH = 0xFF0000FF;
    constexpr uint32 IN_ALPHA_MID = 0xFF7F7FFF;
    constexpr uint32 IN_BETA_HIGH = 0xFFFF0000;
    constexpr uint32 IN_BETA_MID = 0xFFFF7F7F;
    constexpr uint32 IN_LEFT_ALPHA_HIGH = 0xFF00FF00;
    constexpr uint32 IN_LEFT_ALPHA_MID = 0xFF7FFF7F;
    constexpr uint32 IN_P_MID = 0xFF7FFFFF;

    for (int i = 0; i < src_img.width * src_img.height; i++) {
        const uint32 in_color = src_img.data[i];
        uint32& out_color = img->data[i];
        switch (in_color) {
            case IN_BACKGROUND:
                out_color = math::convert_color(color_map.region_color[Region_None]);
                break;
            case IN_ALPHA_HIGH:
                out_color = math::convert_color(color_map.region_color[Region_AlphaHigh]);
                break;
            case IN_ALPHA_MID:
                out_color = math::convert_color(color_map.region_color[Region_AlphaMid]);
                break;
            case IN_BETA_HIGH:
                out_color = math::convert_color(color_map.region_color[Region_BetaHigh]);
                break;
            case IN_BETA_MID:
                out_color = math::convert_color(color_map.region_color[Region_BetaMid]);
                break;
            case IN_LEFT_ALPHA_HIGH:
                out_color = math::convert_color(color_map.region_color[Region_LeftAlphaHigh]);
                break;
            case IN_LEFT_ALPHA_MID:
                out_color = math::convert_color(color_map.region_color[Region_LeftAlphaMid]);
                break;
            case IN_P_MID:
                out_color = math::convert_color(color_map.region_color[Region_PMid]);
                break;
            default:
                break;
        }
    }

    gaussian_blur(img, blur_level);

    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, img->width, img->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img->data);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void init_gui_map(const ColorMap& color_map, int blur_level) { init_map(&gui_img, gui_tex, color_map, blur_level); }
void init_segmentation_map(const ColorMap& color_map, int blur_level) { init_map(&seg_img, seg_tex, color_map, blur_level); }
void init_color_map(const ColorMap& color_map, int blur_level) { init_map(&col_img, col_tex, color_map, blur_level); }

void initialize() {
    if (!read_image(&src_img, VIAMD_IMAGE_DIR "/ramachandran.bmp")) {
        LOG_ERROR("Could not read ramachandran map!");
        return;
    }

    if (!gui_tex) {
        glGenTextures(1, &gui_tex);
        glBindTexture(GL_TEXTURE_2D, gui_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!seg_tex) {
        glGenTextures(1, &seg_tex);
        glBindTexture(GL_TEXTURE_2D, seg_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!col_tex) {
        glGenTextures(1, &col_tex);
        glBindTexture(GL_TEXTURE_2D, col_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    ColorMap seg_color_map{};
    seg_color_map.region_color[Region_None] = vec4(0, 0, 0, 0);
    seg_color_map.region_color[Region_AlphaHigh] = vec4(1, 0, 0, 1);
    seg_color_map.region_color[Region_AlphaMid] = vec4(0, 0, 0, 0);
    seg_color_map.region_color[Region_BetaHigh] = vec4(0, 0, 1, 1);
    seg_color_map.region_color[Region_BetaMid] = vec4(0, 0, 0, 0);
    seg_color_map.region_color[Region_LeftAlphaHigh] = vec4(1, 0, 0, 1);
    seg_color_map.region_color[Region_LeftAlphaMid] = vec4(0, 0, 0, 0);
    seg_color_map.region_color[Region_PMid] = vec4(0, 0, 0, 0);

    init_segmentation_map(seg_color_map, 2);
    init_color_map(ColorMap(), 4);
    init_gui_map(ColorMap(), 2);

    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    GLuint v_shader = glCreateShader(GL_VERTEX_SHADER);
    GLuint f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    defer {
        glDeleteShader(v_shader);
        glDeleteShader(f_shader);
    };
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        LOG_ERROR("Compiling ramachandran vertex shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Compiling ramachandran fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        LOG_ERROR("Linking ramachandran program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    uniform_loc_coord_tex = glGetUniformLocation(program, "u_coord_tex");
    uniform_loc_instance_offset = glGetUniformLocation(program, "u_instance_offset");
    uniform_loc_radius = glGetUniformLocation(program, "u_radius");
    uniform_loc_color = glGetUniformLocation(program, "u_color");
    uniform_loc_outline = glGetUniformLocation(program, "u_outline");

    if (!acc_tex) {
        glGenTextures(1, &acc_tex);
        glBindTexture(GL_TEXTURE_2D, acc_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, acc_width, acc_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!coord_buf) {
        glGenBuffers(1, &coord_buf);
    }

    if (!coord_tex) {
        glGenTextures(1, &coord_tex);
    }

    if (!fbo) {
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, acc_tex, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    if (!vao) {
        glGenVertexArrays(1, &vao);
    }

    if (!vbo) {
        glGenBuffers(1, &vbo);
    }

    if (!ibo) {
        constexpr unsigned char data[4] = {0, 1, 2, 3};
        glGenBuffers(1, &ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 4, data, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
}

void shutdown() {
    if (seg_tex) glDeleteTextures(1, &seg_tex);
    if (acc_tex) glDeleteTextures(1, &acc_tex);
    if (col_tex) glDeleteTextures(1, &col_tex);
    if (coord_buf) glDeleteBuffers(1, &coord_buf);
    if (coord_tex) glDeleteTextures(1, &coord_tex);
    if (fbo) glDeleteFramebuffers(1, &fbo);
    if (vbo) glDeleteBuffers(1, &vbo);
    if (vbo) glDeleteVertexArrays(1, &vao);
}

void clear_accumulation_texture() {
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);

    glViewport(0, 0, acc_width, acc_height);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

void compute_accumulation_texture(Array<const vec2> angles, vec4 color, float radius, float outline) {
    constexpr float ONE_OVER_TWO_PI = 1.f / (2.f * math::PI);

    struct Coord {
        unsigned short x, y;
    };

    // Use fast scratch memory here
    Coord* coords = (Coord*)TMP_MALLOC((angles.count) * sizeof(Coord));
    defer { TMP_FREE(coords); };

    int32 count = 0;
    for (const auto& angle : angles) {
        if (angle.x == 0 || angle.y == 0) continue;
        vec2 coord = vec2(angle.x, angle.y) * ONE_OVER_TWO_PI + 0.5f;  // [-PI, PI] -> [0, 1]
        coord.y = 1.f - coord.y;
        coords[count].x = (unsigned short)(coord.x * 0xffff);
        coords[count].y = (unsigned short)(coord.y * 0xffff);
        count++;
    }

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, count * sizeof(Coord), coords, GL_STREAM_DRAW);

    // Backup GL state
    GLint last_polygon_mode[2];
    glGetIntegerv(GL_POLYGON_MODE, last_polygon_mode);
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    GLenum last_blend_src_rgb;
    glGetIntegerv(GL_BLEND_SRC_RGB, (GLint*)&last_blend_src_rgb);
    GLenum last_blend_dst_rgb;
    glGetIntegerv(GL_BLEND_DST_RGB, (GLint*)&last_blend_dst_rgb);
    GLenum last_blend_src_alpha;
    glGetIntegerv(GL_BLEND_SRC_ALPHA, (GLint*)&last_blend_src_alpha);
    GLenum last_blend_dst_alpha;
    glGetIntegerv(GL_BLEND_DST_ALPHA, (GLint*)&last_blend_dst_alpha);
    GLenum last_blend_equation_rgb;
    glGetIntegerv(GL_BLEND_EQUATION_RGB, (GLint*)&last_blend_equation_rgb);
    GLenum last_blend_equation_alpha;
    glGetIntegerv(GL_BLEND_EQUATION_ALPHA, (GLint*)&last_blend_equation_alpha);
    GLboolean last_enable_blend = glIsEnabled(GL_BLEND);
    GLboolean last_enable_cull_face = glIsEnabled(GL_CULL_FACE);
    GLboolean last_enable_depth_test = glIsEnabled(GL_DEPTH_TEST);
    GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);

    // RENDER TO ACCUMULATION TEXTURE

    glViewport(0, 0, acc_width, acc_height);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Texture 0
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_BUFFER, coord_tex);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RG16, vbo);

    glUseProgram(program);
    glUniform1i(uniform_loc_coord_tex, 0);

    // Draw
    glUniform1f(uniform_loc_radius, radius * 0.01f);
    glUniform1i(uniform_loc_instance_offset, 0);
    glUniform4fv(uniform_loc_color, 1, &color[0]);
    glUniform1f(uniform_loc_outline, outline);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glDrawElementsInstanced(GL_TRIANGLE_STRIP, 4, GL_UNSIGNED_BYTE, 0, count);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    // Restore modified GL state
    glBlendEquationSeparate(last_blend_equation_rgb, last_blend_equation_alpha);
    glBlendFuncSeparate(last_blend_src_rgb, last_blend_dst_rgb, last_blend_src_alpha, last_blend_dst_alpha);
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
    glPolygonMode(GL_FRONT_AND_BACK, last_polygon_mode[0]);
    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

}  // namespace ramachandran
