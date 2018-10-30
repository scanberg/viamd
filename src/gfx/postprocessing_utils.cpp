/*-----------------------------------------------------------------------
  Copyright (c) 2014, NVIDIA. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.
   * Neither the name of its contributors may be used to endorse
     or promote products derived from this software without specific
     prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
-----------------------------------------------------------------------*/

// Shaders for HBAO are taken from nVidias examples and are copyright protected as stated above

#include "postprocessing_utils.h"
#include <core/types.h>
#include <core/common.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <gfx/gl_utils.h>
#include <stdio.h>

namespace postprocessing {

// @TODO: Use half-res render targets for SSAO
// @TODO: Use shared textures for all postprocessing operations
// @TODO: Use some kind of unified pipeline for all post processing operations

static struct {
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint v_shader_fs_quad = 0;
    GLuint tex_width = 0;
    GLuint tex_height = 0;

    struct {
        GLuint fbo = 0;
        GLuint program = 0;
        struct {
            GLuint color_coc = 0;
        } tex;
        struct {
            GLint tex_depth = -1;
            GLint tex_color = -1;
            GLint focus_point = -1;
            GLint focus_scale = -1;
        } uniform_loc;
    } half_res;

    struct {
        GLuint fbo = 0;
        GLuint texture = 0;
        GLuint program_persp = 0;
        GLuint program_ortho = 0;
        struct {
            GLint clip_info = -1;
            GLint tex_depth = -1;
        } uniform_loc;
    } linearize_depth;

    struct {
        GLuint tex_random = 0;
        GLuint ubo_hbao_data = 0;

        struct {
            GLuint fbo = 0;
            GLuint texture = 0;
            GLuint program = 0;

            struct {
                GLint control_buffer = -1;
                GLint tex_linear_depth = -1;
                GLint tex_normal = -1;
                GLint tex_random = -1;
            } uniform_loc;
        } hbao;

        struct {
            GLuint fbo = 0;
            GLuint texture = 0;
            GLuint program_first = 0;
            GLuint program_second = 0;
            struct {
                GLint sharpness = -1;
                GLint inv_res_dir = -1;
                GLint texture = -1;
            } uniform_loc;
        } blur;
    } ssao;

    struct {
        GLuint program = 0;
        struct {
            GLint tex_half_res = -1;
            GLint tex_color = -1;
            GLint tex_depth = -1;
            GLint pixel_size = -1;
            GLint focus_point = -1;
            GLint focus_scale = -1;
        } uniform_loc;
    } bokeh_dof;

    struct {
        GLuint program = 0;
    } bloom;

    struct {
        GLuint program = 0;
        struct {
            GLint mode = -1;
            GLint tex_color = -1;
        } uniform_loc;
    } tonemapping;

} gl;

static const char* v_shader_src_fs_quad = R"(
#version 150 core

out vec2 tc;

void main() {
	uint idx = uint(gl_VertexID) % 3U;
	gl_Position = vec4(
		(float( idx     &1U)) * 4.0 - 1.0,
		(float((idx>>1U)&1U)) * 4.0 - 1.0,
		0, 1.0);
	tc = gl_Position.xy * 0.5 + 0.5;
}
)";

static const char* f_shader_src_linearize_depth = R"(
#ifndef PERSPECTIVE
#define PERSPECTIVE 1
#endif

// z_n * z_f,  z_n - z_f,  z_f, *not used*
uniform vec4 u_clip_info;
uniform sampler2D u_tex_depth;

float ReconstructCSZ(float d, vec4 clip_info) {
#ifdef PERSPECTIVE
    return (clip_info[0] / (clip_info[1] * d + clip_info[2]));
#else
    return (clipInfo[1]+clipInfo[2] - d * clipInfo[1]);
#endif
}

out vec4 out_frag;

void main() {
  float d = texelFetch(u_tex_depth, ivec2(gl_FragCoord.xy), 0).x;
  out_frag = vec4(ReconstructCSZ(d, u_clip_info));
}
)";

static bool setup_program(GLuint* program, const char* name, const char* f_shader_src, const char* defines = nullptr) {
    ASSERT(program);
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    auto f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    if (defines) {
        const char* sources[2] = {defines, f_shader_src};
        glShaderSource(f_shader, 2, sources, 0);
    } else {
        glShaderSource(f_shader, 1, &f_shader_src, 0);
    }

    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        LOG_ERROR("Error while compiling %s shader:\n%s", name, buffer);
        return false;
    }

    if (!*program) {
        *program = glCreateProgram();
    } else {
        // TODO: DETATCH ANY SHADERS?
    }

    glAttachShader(*program, gl.v_shader_fs_quad);
    glAttachShader(*program, f_shader);
    glLinkProgram(*program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, *program)) {
        LOG_ERROR("Error while linking %s program:\n%s", name, buffer);
        return false;
    }

    glDetachShader(*program, gl.v_shader_fs_quad);
    glDetachShader(*program, f_shader);
    glDeleteShader(f_shader);

    return true;
}

static bool is_orthographic_proj_matrix(const mat4& proj_mat) { return math::length2(vec2(proj_mat[3])) > 0.f; }

namespace ssao {
#ifndef AO_RANDOM_TEX_SIZE
#define AO_RANDOM_TEX_SIZE 4
#endif

#ifndef AO_MAX_SAMPLES
#define AO_MAX_SAMPLES 1
#endif

#ifndef AO_DIRS
#define AO_DIRS 8
#endif

#ifndef AO_SAMPLES
#define AO_SAMPLES 4
#endif

#ifndef AO_BLUR
#define AO_BLUR 1
#endif

static GLuint fbo_hbao = 0;
static GLuint fbo_blur = 0;

static GLuint tex_random = 0;
static GLuint tex_hbao = 0;
static GLuint tex_blur = 0;

static GLuint prog_hbao = 0;
static GLuint prog_blur_first = 0;
static GLuint prog_blur_second = 0;

static GLint uniform_block_index_hbao_control_buffer = -1;
static GLint uniform_loc_hbao_tex_linear_depth = -1;
static GLint uniform_loc_hbao_tex_normal = -1;
static GLint uniform_loc_hbao_tex_random = -1;

static GLint uniform_loc_blur_sharpness = -1;
static GLint uniform_loc_blur_inv_res_dir = -1;
static GLint uniform_loc_blur_texture = -1;

static GLuint ubo_hbao_data = 0;

static GLuint tex_width;
static GLuint tex_height;

struct HBAOData {
    float radius_to_screen;
    float r2;
    float neg_inv_r2;
    float n_dot_v_bias;

    vec2 inv_full_res;
    vec2 inv_quarter_res;

    float ao_multiplier;
    float pow_exponent;
    vec2 _pad0;

    vec4 proj_info;

    vec2 proj_scale;
    int proj_ortho;
    int _pad1;
};

static const char* f_shader_src_hbao = R"(
#pragma optionNV(unroll all)

#ifndef AO_BLUR
#define AO_BLUR 0
#endif

#define M_PI 3.14159265

#ifndef AO_STEPS
#define AO_STEPS 4
#endif

#ifndef AO_DIRS
#define AO_DIRS 8
#endif

#ifndef AO_USE_NORMAL
#define AO_USE_NORMAL 1
#endif

#ifndef AO_RANDOM_TEX_SIZE
#define AO_RANDOM_TEX_SIZE 4
#endif

struct HBAOData {
  float   radius_to_screen;
  float   r2;
  float   neg_inv_r2;
  float   n_dot_v_bias;
 
  vec2    inv_full_res;
  vec2    inv_quarter_res;
  
  float   ao_multiplier;
  float   pow_exponent;
  vec2    _pad0;
  
  vec4    proj_info;

  vec2    proj_scale;
  int     proj_ortho;
  int     _pad1;
  
  //vec4    offsets[AO_RANDOM_TEX_SIZE*AO_RANDOM_TEX_SIZE];
  //vec4    jitters[AO_RANDOM_TEX_SIZE*AO_RANDOM_TEX_SIZE];
};

// tweakables
const float NUM_STEPS = AO_STEPS;
const float NUM_DIRECTIONS = AO_DIRS; // tex_random/jitter initialization depends on this

layout(std140) uniform u_control_buffer {
  HBAOData control;
};

uniform sampler2D u_tex_linear_depth;
uniform sampler2D u_tex_normal;
uniform sampler2D u_tex_random;

in vec2 tc;
out vec4 out_frag;

void OutputColor(vec4 color) {
  out_frag = color;
}

//----------------------------------------------------------------------------------

vec3 UVToView(vec2 uv, float eye_z) {
  return vec3((uv * control.proj_info.xy + control.proj_info.zw) * (control.proj_ortho != 0 ? 1. : eye_z), eye_z);
}

vec3 FetchViewPos(vec2 uv) {
  float ViewDepth = textureLod(u_tex_linear_depth, uv, 0).x;
  return UVToView(uv, ViewDepth);
}

vec3 MinDiff(vec3 P, vec3 Pr, vec3 Pl) {
  vec3 V1 = Pr - P;
  vec3 V2 = P - Pl;
  return (dot(V1,V1) < dot(V2,V2)) ? V1 : V2;
}

vec3 ReconstructNormal(vec2 UV, vec3 P) {
  vec3 Pr = FetchViewPos(UV + vec2(control.inv_full_res.x, 0));
  vec3 Pl = FetchViewPos(UV + vec2(-control.inv_full_res.x, 0));
  vec3 Pt = FetchViewPos(UV + vec2(0, control.inv_full_res.y));
  vec3 Pb = FetchViewPos(UV + vec2(0, -control.inv_full_res.y));
  return normalize(cross(MinDiff(P, Pr, Pl), MinDiff(P, Pt, Pb)));
}

vec3 DecodeNormal(vec2 enc) {
    vec2 fenc = enc*4-2;
    float f = dot(fenc,fenc);
    float g = sqrt(1-f/4.0);
    vec3 n;
    n.xy = fenc*g;
    n.z = 1-f/2.0;
    return n;
}

vec3 FetchViewNormal(vec2 uv) {
	vec2 enc = textureLod(u_tex_normal, uv, 0).xy;
	vec3 n = DecodeNormal(enc);
	return n * vec3(1,1,-1);
}

//----------------------------------------------------------------------------------
float Falloff(float DistanceSquare) {
  // 1 scalar mad instruction
  return DistanceSquare * control.neg_inv_r2 + 1.0;
}

//----------------------------------------------------------------------------------
// P = view-space position at the kernel center
// N = view-space normal at the kernel center
// S = view-space position of the current sample
//----------------------------------------------------------------------------------
float ComputeAO(vec3 P, vec3 N, vec3 S) {
  vec3 V = S - P;
  float VdotV = dot(V, V);
  float NdotV = dot(N, V) * 1.0/sqrt(VdotV);

  // Use saturate(x) instead of max(x,0.f) because that is faster on Kepler
  return clamp(NdotV - control.n_dot_v_bias,0,1) * clamp(Falloff(VdotV),0,1);
}

//----------------------------------------------------------------------------------
vec2 RotateDirection(vec2 Dir, vec2 CosSin) {
  return vec2(Dir.x*CosSin.x - Dir.y*CosSin.y,
              Dir.x*CosSin.y + Dir.y*CosSin.x);
}

//----------------------------------------------------------------------------------
vec4 GetJitter() {
  // (cos(Alpha),sin(Alpha),rand1,rand2)
  return textureLod(u_tex_random, (gl_FragCoord.xy / AO_RANDOM_TEX_SIZE), 0);
}

//----------------------------------------------------------------------------------
float ComputeCoarseAO(vec2 FullResUV, float RadiusPixels, vec4 Rand, vec3 ViewPosition, vec3 ViewNormal) {
  // Divide by NUM_STEPS+1 so that the farthest samples are not fully attenuated
  float StepSizePixels = RadiusPixels / (NUM_STEPS + 1);

  const float Alpha = 2.0 * M_PI / NUM_DIRECTIONS;
  float AO = 0;

  for (float DirectionIndex = 0; DirectionIndex < NUM_DIRECTIONS; ++DirectionIndex)
  {
    float Angle = Alpha * DirectionIndex;

    // Compute normalized 2D direction
    vec2 Direction = RotateDirection(vec2(cos(Angle), sin(Angle)), Rand.xy);

    // Jitter starting sample within the first step
    float RayPixels = (Rand.z * StepSizePixels + 1.0);

    for (float StepIndex = 0; StepIndex < NUM_STEPS; ++StepIndex)
    {
      vec2 SnappedUV = round(RayPixels * Direction) * control.inv_full_res + FullResUV;
      vec3 S = FetchViewPos(SnappedUV);

      RayPixels += StepSizePixels;

      AO += ComputeAO(ViewPosition, ViewNormal, S);
    }
  }

  AO *= control.ao_multiplier / (NUM_DIRECTIONS * NUM_STEPS);

  return clamp(1.0 - AO, 0, 1);
}

//----------------------------------------------------------------------------------
void main() {
  vec2 uv = tc;
  vec3 ViewPosition = FetchViewPos(uv);

  // Reconstruct view-space normal from nearest neighbors
#if AO_USE_NORMAL
  vec3 ViewNormal = FetchViewNormal(uv);
#else
  vec3 ViewNormal = vec3(0,0,-1);
#endif

  // Compute projection of disk of radius control.R into screen space
  float RadiusPixels = control.radius_to_screen / (control.proj_ortho != 0 ? 1.0 : ViewPosition.z);

  // Get jitter vector for the current full-res pixel
  vec4 Rand = GetJitter();

  float AO = ComputeCoarseAO(uv, RadiusPixels, Rand, ViewPosition, ViewNormal);

#if AO_BLUR
  OutputColor(vec4(pow(AO, control.pow_exponent), ViewPosition.z, 0, 1));
#else
  OutputColor(vec4(vec3(pow(AO, control.pow_exponent)), 1));
#endif 
}
)";

static const char* f_shader_src_hbao_blur = R"(
#pragma optionNV(unroll all)

const float KERNEL_RADIUS = 3;
  
uniform float u_sharpness;
uniform vec2  u_inv_res_dir; // either set x to 1/width or y to 1/height
uniform sampler2D u_texture;

in vec2 tc;
out vec4 out_frag;

#ifndef AO_BLUR_PRESENT
#define AO_BLUR_PRESENT 1
#endif

//-------------------------------------------------------------------------

float BlurFunction(vec2 uv, float r, float center_c, float center_d, inout float w_total)
{
  vec2  aoz = texture( u_texture, uv ).xy;
  float c = aoz.x;
  float d = aoz.y;
  
  const float BlurSigma = float(KERNEL_RADIUS) * 0.5;
  const float BlurFalloff = 1.0 / (2.0*BlurSigma*BlurSigma);
  
  float ddiff = (d - center_d) * u_sharpness;
  float w = exp2(-r*r*BlurFalloff - ddiff*ddiff);
  w_total += w;

  return c*w;
}

void main()
{
  vec2  aoz = texture( u_texture, tc ).xy;
  float center_c = aoz.x;
  float center_d = aoz.y;
  
  float c_total = center_c;
  float w_total = 1.0;
  
  for (float r = 1; r <= KERNEL_RADIUS; ++r)
  {
    vec2 uv = tc + u_inv_res_dir * r;
    c_total += BlurFunction(uv, r, center_c, center_d, w_total);  
  }
  
  for (float r = 1; r <= KERNEL_RADIUS; ++r)
  {
    vec2 uv = tc - u_inv_res_dir * r;
    c_total += BlurFunction(uv, r, center_c, center_d, w_total);  
  }
  
#if AO_BLUR_PRESENT
  out_frag = vec4(vec3(c_total/w_total), 1);
#else
  out_frag = vec4(c_total/w_total, center_d, 0, 1);
#endif
}
)";

void setup_ubo_hbao_data(GLuint ubo, int width, int height, const mat4& proj_mat, float intensity, float radius, float bias) {
    ASSERT(ubo);
    constexpr float METERS_TO_VIEWSPACE = 1.f;
    const float* proj_data = &proj_mat[0][0];

    bool is_ortho = is_orthographic_proj_matrix(proj_mat);

    vec4 proj_info;
    float proj_scl;
    if (!is_ortho) {
        proj_info = vec4(2.0f / (proj_data[4 * 0 + 0]),                          // (x) * (R - L)/N
                         2.0f / (proj_data[4 * 1 + 1]),                          // (y) * (T - B)/N
                         -(1.0f - proj_data[4 * 2 + 0]) / proj_data[4 * 0 + 0],  // L/N
                         -(1.0f + proj_data[4 * 2 + 1]) / proj_data[4 * 1 + 1]   // B/N
        );

        // proj_scl = float(height) / (math::tan(fovy * 0.5f) * 2.0f);
        proj_scl = float(height) * proj_data[4 * 1 + 1] * 0.5f;
        proj_scl = float(height) / (math::tan(math::PI / 8.f) * 2.0f);
    } else {
        proj_info = vec4(2.0f / (proj_data[4 * 0 + 0]),                          // ((x) * R - L)
                         2.0f / (proj_data[4 * 1 + 1]),                          // ((y) * T - B)
                         -(1.0f + proj_data[4 * 3 + 0]) / proj_data[4 * 0 + 0],  // L
                         -(1.0f - proj_data[4 * 3 + 1]) / proj_data[4 * 1 + 1]   // B
        );
        proj_scl = float(height) / proj_info[1];
    }

    float r = radius * METERS_TO_VIEWSPACE;

    HBAOData data;
    data.radius_to_screen = r * 0.5f * proj_scl;
    data.r2 = r * r;
    data.neg_inv_r2 = -1.f / (r * r);
    data.n_dot_v_bias = math::clamp(bias, 0.f, 1.f);

    data.inv_full_res = vec2(1.f / float(width), 1.f / float(height));
    data.inv_quarter_res = vec2(1.f / float((width + 3) / 4), 1.f / float((height + 3) / 4));

    data.ao_multiplier = 1.f / (1.f - data.n_dot_v_bias);
    data.pow_exponent = math::max(intensity, 0.f);

    data.proj_info = proj_info;
    // @NOTE: Is one needed?
    data.proj_scale = vec2(proj_scl);
    data.proj_ortho = is_ortho ? 1 : 0;

    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(HBAOData), &data, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}

void initialize_rnd_tex(GLuint rnd_tex, int num_direction) {
    ASSERT(AO_MAX_SAMPLES == 1);
    constexpr int BUFFER_SIZE = AO_RANDOM_TEX_SIZE * AO_RANDOM_TEX_SIZE * AO_MAX_SAMPLES;
    signed short buffer[BUFFER_SIZE * 4];

    for (int i = 0; i < BUFFER_SIZE; i++) {
#define SCALE ((1 << 15))
        float rand1 = math::rnd();
        float rand2 = math::rnd();
        float angle = 2.f * math::PI * rand1 / (float)num_direction;

        buffer[i * 4 + 0] = (signed short)(SCALE * math::cos(angle));
        buffer[i * 4 + 1] = (signed short)(SCALE * math::sin(angle));
        buffer[i * 4 + 2] = (signed short)(SCALE * rand2);
        buffer[i * 4 + 3] = (signed short)(SCALE * 0);
#undef SCALE
    }

    // @TODO: If MSAA and AO_MAX_SAMPLES > 1, then this probably has to go into a texture array
    glBindTexture(GL_TEXTURE_2D, rnd_tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16_SNORM, AO_RANDOM_TEX_SIZE, AO_RANDOM_TEX_SIZE, 0, GL_RGBA, GL_SHORT, buffer);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

float compute_sharpness(float radius) { return 30.f / radius; }

void initialize(int width, int height) {
    // @TODO: dynamically generate this
    const char* defines = R"(
        #version 150 core
        #define AO_RANDOM_TEX_SIZE 4
        #define PERSPECTIVE 1
        #define AO_BLUR 0
        #define AO_STEPS 4
        #define AO_DIRS 8
        #define AO_USE_NORMAL 1
    )";

    const char* define_blur_first = R"(
        #version 150 core
        #define AO_BLUR_PRESENT 0
    )";

    const char* define_blur_second = R"(
        #version 150 core
        #define AO_BLUR_PRESENT 1
    )";

    setup_program(&prog_hbao, "hbao", f_shader_src_hbao, defines);
    setup_program(&prog_blur_first, "hbao first blur", f_shader_src_hbao_blur, define_blur_first);
    setup_program(&prog_blur_second, "hbao second blur", f_shader_src_hbao_blur, define_blur_second);

    uniform_block_index_hbao_control_buffer = glGetUniformBlockIndex(prog_hbao, "u_control_buffer");
    uniform_loc_hbao_tex_linear_depth = glGetUniformLocation(prog_hbao, "u_tex_linear_depth");
    uniform_loc_hbao_tex_normal = glGetUniformLocation(prog_hbao, "u_tex_normal");
    uniform_loc_hbao_tex_random = glGetUniformLocation(prog_hbao, "u_tex_random");

    ASSERT(uniform_block_index_hbao_control_buffer != -1);
    ASSERT(uniform_loc_hbao_tex_linear_depth != -1);
    ASSERT(uniform_loc_hbao_tex_normal != -1);
    ASSERT(uniform_loc_hbao_tex_random != -1);

    uniform_loc_blur_sharpness = glGetUniformLocation(prog_blur_second, "u_sharpness");
    uniform_loc_blur_inv_res_dir = glGetUniformLocation(prog_blur_second, "u_inv_res_dir");
    uniform_loc_blur_texture = glGetUniformLocation(prog_blur_second, "u_texture");

    ASSERT(uniform_loc_blur_sharpness != -1);
    ASSERT(uniform_loc_blur_inv_res_dir != -1);
    ASSERT(uniform_loc_blur_texture != -1);

    if (!fbo_hbao) glGenFramebuffers(1, &fbo_hbao);
    if (!fbo_blur) glGenFramebuffers(1, &fbo_blur);

    if (!tex_random) glGenTextures(1, &tex_random);
    if (!tex_hbao) glGenTextures(1, &tex_hbao);
    if (!tex_blur) glGenTextures(1, &tex_blur);

    initialize_rnd_tex(tex_random, AO_DIRS);

    glBindTexture(GL_TEXTURE_2D, tex_hbao);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, width, height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, tex_blur);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG16F, width, height, 0, GL_RG, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glBindTexture(GL_TEXTURE_2D, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, fbo_hbao);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_hbao, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, fbo_blur);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_blur, 0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    if (!ubo_hbao_data) glGenBuffers(1, &ubo_hbao_data);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo_hbao_data);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(HBAOData), nullptr, GL_DYNAMIC_DRAW);

    tex_width = width;
    tex_height = height;
}

void shutdown() {
    if (fbo_hbao) glDeleteFramebuffers(1, &fbo_hbao);
    if (fbo_blur) glDeleteFramebuffers(1, &fbo_blur);

    if (tex_random) glDeleteTextures(1, &tex_random);
    if (tex_blur) glDeleteTextures(1, &tex_blur);

    if (ubo_hbao_data) glDeleteBuffers(1, &ubo_hbao_data);

    if (prog_hbao) glDeleteProgram(prog_hbao);
    if (prog_blur_second) glDeleteProgram(prog_blur_second);
    if (prog_blur_first) glDeleteProgram(prog_blur_first);
}

}  // namespace ssao

namespace deferred {

static GLuint prog_deferred = 0;
static GLint uniform_loc_texture_depth = -1;
static GLint uniform_loc_texture_color = -1;
static GLint uniform_loc_texture_normal = -1;
static GLint uniform_loc_inv_proj_mat = -1;
static const char* f_shader_src_deferred = R"(
#version 150 core

uniform sampler2D u_texture_depth;
uniform sampler2D u_texture_color;
uniform sampler2D u_texture_normal;

uniform mat4 u_inv_proj_mat;

in vec2 tc;
out vec4 out_frag;

vec4 depth_to_view_coord(vec2 tex_coord, float depth) {
    vec4 clip_coord = vec4(vec3(tex_coord, depth) * 2.0 - 1.0, 1.0);
    vec4 view_coord = u_inv_proj_mat * clip_coord;
    return view_coord / view_coord.w;
}

float fresnel(float H_dot_V) {
    const float n1 = 1.0;
    const float n2 = 1.5;
    const float R0 = pow((n1-n2)/(n1+n2), 2);

    return R0 + (1.0 - R0)*pow(1.0 - H_dot_V, 5);
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

vec3 shade(vec3 color, vec3 V, vec3 N) {
    const vec3 env_radiance = vec3(0.5);
    const vec3 dir_radiance = vec3(0.5);
    const vec3 L = normalize(vec3(1));
    const float spec_exp = 10.0;

    vec3 H = normalize(L + V);
    float H_dot_V = max(0.0, dot(H, V));
    float N_dot_H = max(0.0, dot(N, H));
    float N_dot_L = max(0.0, dot(N, L));
    float fr = fresnel(H_dot_V);

    vec3 diffuse = color.rgb * (env_radiance + N_dot_L * dir_radiance);
    vec3 specular = dir_radiance * pow(N_dot_H, spec_exp);

    return mix(diffuse, specular, fr);
}

void main() {
	float depth = texelFetch(u_texture_depth, ivec2(gl_FragCoord.xy), 0).x;
	if (depth == 1.0) discard;
	vec4 color = texelFetch(u_texture_color, ivec2(gl_FragCoord.xy), 0);
	vec3 normal = decode_normal(texelFetch(u_texture_normal, ivec2(gl_FragCoord.xy), 0).xy);
	vec4 view_coord = depth_to_view_coord(tc, depth);

	vec3 N = normal;
	vec3 V = -normalize(view_coord.xyz);
	vec3 result = shade(color.rgb, V, N);

	out_frag = vec4(result, color.a);
}
)";

void initialize() {
    if (!prog_deferred) setup_program(&prog_deferred, "deferred", f_shader_src_deferred);
    uniform_loc_texture_depth = glGetUniformLocation(prog_deferred, "u_texture_depth");
    uniform_loc_texture_color = glGetUniformLocation(prog_deferred, "u_texture_color");
    uniform_loc_texture_normal = glGetUniformLocation(prog_deferred, "u_texture_normal");
    uniform_loc_inv_proj_mat = glGetUniformLocation(prog_deferred, "u_inv_proj_mat");
}

void shutdown() {
    if (prog_deferred) glDeleteProgram(prog_deferred);
}
}  // namespace deferred

namespace tonemapping {
static GLuint prog_tonemap = 0;
static GLint uniform_loc_texture = -1;
static const char* f_shader_src_tonemap = R"(
#version 150 core

uniform sampler2D u_texture;
out vec4 out_frag;

vec3 passthrough(vec3 c) {
    return c;
}

vec3 reinhard(vec3 c) {
    c *= 1;  // Hardcoded Exposure Adjustment
    c = c / (c + vec3(1.0));
    return pow(c, vec3(1.0 / 2.2));
}

vec3 uncharted2_tonemap(vec3 x) {
    const float A = 0.15;
    const float B = 0.50;
    const float C = 0.10;
    const float D = 0.20;
    const float E = 0.02;
    const float F = 0.30;
    return ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}

vec3 uncharted2(vec3 c) {
    const float W = 11.2;
    //c *= 2;  // Hardcoded Exposure Adjustment

    float exp_bias = 0.5;
    vec3 curr = uncharted2_tonemap(exp_bias * c);

    vec3 white_scale = vec3(1.0) / uncharted2_tonemap(vec3(W));
    vec3 color = curr * white_scale;
      
    return pow(color, vec3(1.0/2.2));
}

vec3 hejl_dawsson(vec3 c) {
   c *= 1;  // Hardcoded Exposure Adjustment
   vec3 x = max(vec3(0), c-vec3(0.004));
   return (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
}

void main() {
	vec4 color = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
	out_frag = vec4(color.rgb, color.a);
}
)";

void initialize() {
    if (!prog_tonemap) setup_program(&prog_tonemap, "tonemap", f_shader_src_tonemap);
    uniform_loc_texture = glGetUniformLocation(prog_tonemap, "u_texture");
}

void shutdown() {
    if (prog_tonemap) glDeleteProgram(prog_tonemap);
}

}  // namespace tonemapping

namespace bokeh_dof {

const char* f_shader_half_res_src = R"(
#version 150 core

uniform sampler2D u_tex_depth; // Linear depth
uniform sampler2D u_tex_color;

uniform float u_focus_point;
uniform float u_focus_scale;

const float MAX_BLUR_SIZE = 20.0; 

float getBlurSize(float d, float fp, float fs) {
	float coc = clamp((1.0 / fp - 1.0 / d) * fs, -1.0, 1.0);
	return abs(coc);
}

in vec2 tc;
out vec4 out_frag;

void main() {
	float depth = texture(u_tex_depth, tc).r;
	vec3  color = texture(u_tex_color, tc).rgb;
	float coc   = getBlurSize(depth, u_focus_point, u_focus_scale);

	out_frag = vec4(color, coc);
}
)";

// From http://blog.tuxedolabs.com/2018/05/04/bokeh-depth-of-field-in-single-pass.html
const char* f_shader_src = R"(
#version 150 core

uniform sampler2D uHalfRes; // Half res color (rgb) and coc (a)
uniform sampler2D uColor;   // Image to be processed 
uniform sampler2D uDepth;   // Linear depth, where 1.0 == far plane 

uniform vec2 uPixelSize;    // The size of a pixel: vec2(1.0/width, 1.0/height) 
uniform float uFocusPoint; 
uniform float uFocusScale;

const float GOLDEN_ANGLE = 2.39996323; 
const float MAX_BLUR_SIZE = 20.0; 
const float RAD_SCALE = 2.0; // Smaller = nicer blur, larger = faster

float getBlurSize(float depth, float focusPoint, float focusScale)
{
	float coc = clamp((1.0 / focusPoint - 1.0 / depth)*focusScale, -1.0, 1.0);
	return abs(coc) * MAX_BLUR_SIZE;
}

vec3 depthOfField(vec2 tex_coord, float focus_point, float focus_scale)
{
	float center_depth  = texture(uDepth, tex_coord).r;
	vec3  center_color  = texture(uColor, tex_coord).rgb;
	float center_coc    = getBlurSize(center_depth, focus_point, focus_scale);
	vec4  color_coc_sum = vec4(center_color, center_coc);

	float contrib_sum   = 1.0;
	float radius        = RAD_SCALE;
	int	  offset_idx    = 0;
	float ang           = 0.0;

	for (; radius < MAX_BLUR_SIZE; ang += GOLDEN_ANGLE)
	{
		vec2 tc = tex_coord + vec2(cos(ang), sin(ang)) * uPixelSize * radius;
		float sample_depth = texture(uDepth, tc).r;
		vec3  sample_color = texture(uColor, tc).rgb;
		float sample_coc   = getBlurSize(sample_depth, focus_point, focus_scale);
		vec4 sample_color_coc = vec4(sample_color, sample_coc);

/*
		if (sample_depth > center_depth)
			sample_coc = clamp(sample_coc, 0.0, center_coc*2.0);
*/

		color_coc_sum     += mix(color_coc_sum / contrib_sum, sample_color_coc, smoothstep(radius-0.5, radius+0.5, sample_color_coc.a));
		contrib_sum       += 1.0;
		radius            += RAD_SCALE/radius;

		if (color_coc_sum.a / contrib_sum > 10.0) break;
	}

	for (; radius < MAX_BLUR_SIZE; ang += GOLDEN_ANGLE) {
		vec2 tc = tex_coord + vec2(cos(ang), sin(ang)) * uPixelSize * radius;
		vec4  sample_color_coc = texture(uHalfRes, tc) * vec4(1, 1, 1, MAX_BLUR_SIZE);
		
/*
		float sample_depth = texture(uDepth, tc).r;
		if (sample_depth > center_depth) sample_color_coc.a = clamp(sample_color_coc.a, 0.0, center_coc*2.0);
*/

		color_coc_sum     += mix(color_coc_sum / contrib_sum, sample_color_coc, smoothstep(radius-0.5, radius+0.5, sample_color_coc.a));
		contrib_sum       += 1.0;
		radius            += RAD_SCALE/radius;
	}

	return color_coc_sum.rgb /= contrib_sum;
}

in vec2 tc;
out vec4 out_frag;

void main() {
	const float focus_point = 0.5;
	const float focus_scale = 0.5;

	vec3 dof = depthOfField(tc, uFocusPoint, uFocusScale);

	out_frag = vec4(dof, 1);
}

)";

}  // namespace bokeh_dof

namespace fxaa {}

void initialize(int width, int height) {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    if (!gl.vao) glGenVertexArrays(1, &gl.vao);
    if (!gl.vbo) glGenBuffers(1, &gl.vbo);

    glBindVertexArray(gl.vao);
    glBindBuffer(GL_ARRAY_BUFFER, gl.vbo);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (const GLvoid*)0);
    glBindVertexArray(0);

    if (!gl.v_shader_fs_quad) {
        gl.v_shader_fs_quad = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(gl.v_shader_fs_quad, 1, &v_shader_src_fs_quad, 0);
        glCompileShader(gl.v_shader_fs_quad);
        if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, gl.v_shader_fs_quad)) {
            LOG_ERROR("Error while compiling postprocessing fs-quad vertex shader:\n%s", buffer);
        }
    }

    // LINEARIZE DEPTH
    if (!gl.linearize_depth.program_persp)
        setup_program(&gl.linearize_depth.program_persp, "linearize depth persp", f_shader_src_linearize_depth,
                      "#version 150 core\n#define PERSPECTIVE 1");
    if (!gl.linearize_depth.program_ortho)
        setup_program(&gl.linearize_depth.program_ortho, "linearize depth ortho", f_shader_src_linearize_depth,
                      "#version 150 core\n#define PERSPECTIVE 0");
    if (!gl.linearize_depth.fbo) glGenFramebuffers(1, &gl.linearize_depth.fbo);
    if (!gl.linearize_depth.texture) glGenTextures(1, &gl.linearize_depth.texture);

    glBindTexture(GL_TEXTURE_2D, gl.linearize_depth.texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linearize_depth.fbo);
    glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.linearize_depth.texture, 0);
    GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE) {
        LOG_ERROR("Something went wrong");
    }
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    gl.linearize_depth.uniform_loc.clip_info = glGetUniformLocation(gl.linearize_depth.program_persp, "u_clip_info");
    gl.linearize_depth.uniform_loc.tex_depth = glGetUniformLocation(gl.linearize_depth.program_persp, "u_tex_depth");

    // HALF RES
    if (!gl.half_res.program) {
        setup_program(&gl.half_res.program, "half-res", bokeh_dof::f_shader_half_res_src);
        if (gl.half_res.program) {
            gl.half_res.uniform_loc.tex_depth = glGetUniformLocation(gl.half_res.program, "u_tex_depth");
            gl.half_res.uniform_loc.tex_color = glGetUniformLocation(gl.half_res.program, "u_tex_color");
            gl.half_res.uniform_loc.focus_point = glGetUniformLocation(gl.half_res.program, "u_focus_point");
            gl.half_res.uniform_loc.focus_scale = glGetUniformLocation(gl.half_res.program, "u_focus_scale");
        }
    }

    if (!gl.half_res.tex.color_coc) {
        glGenTextures(1, &gl.half_res.tex.color_coc);
    }
    glBindTexture(GL_TEXTURE_2D, gl.half_res.tex.color_coc);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width / 2, height / 2, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);

    if (!gl.half_res.fbo) {
        glGenFramebuffers(1, &gl.half_res.fbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.half_res.fbo);
        glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gl.half_res.tex.color_coc, 0);
        GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            LOG_ERROR("Something went wrong");
        }
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    }

    // DOF
    if (!gl.bokeh_dof.program) {
        if (setup_program(&gl.bokeh_dof.program, "bokeh dof", bokeh_dof::f_shader_src)) {
            gl.bokeh_dof.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.program, "uHalfRes");
            gl.bokeh_dof.uniform_loc.tex_color = glGetUniformLocation(gl.bokeh_dof.program, "uColor");
            gl.bokeh_dof.uniform_loc.tex_depth = glGetUniformLocation(gl.bokeh_dof.program, "uDepth");
            gl.bokeh_dof.uniform_loc.pixel_size = glGetUniformLocation(gl.bokeh_dof.program, "uPixelSize");
            gl.bokeh_dof.uniform_loc.focus_point = glGetUniformLocation(gl.bokeh_dof.program, "uFocusPoint");
            gl.bokeh_dof.uniform_loc.focus_scale = glGetUniformLocation(gl.bokeh_dof.program, "uFocusScale");
        }
    }

    gl.tex_width = width;
    gl.tex_height = height;

    ssao::initialize(width, height);
    deferred::initialize();
    tonemapping::initialize();
}

void shutdown() {
    ssao::shutdown();
    deferred::shutdown();
    tonemapping::shutdown();

    if (gl.vao) glDeleteVertexArrays(1, &gl.vao);
    if (gl.vbo) glDeleteBuffers(1, &gl.vbo);
    if (gl.v_shader_fs_quad) glDeleteShader(gl.v_shader_fs_quad);
}

void linearize_depth(GLuint depth_tex, float near_plane, float far_plane, bool orthographic = false) {
    const vec4 clip_info(near_plane * far_plane, near_plane - far_plane, far_plane, 0);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.linearize_depth.fbo);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);

    if (orthographic)
        glUseProgram(gl.linearize_depth.program_ortho);
    else
        glUseProgram(gl.linearize_depth.program_persp);
    glUniform1i(gl.linearize_depth.uniform_loc.tex_depth, 0);
    glUniform4fv(gl.linearize_depth.uniform_loc.clip_info, 1, &clip_info[0]);

    // ASSUME THAT THE APPROPRIATE FS_QUAD VAO IS BOUND
    glDrawArrays(GL_TRIANGLES, 0, 3);
}

void apply_ssao(GLuint depth_tex, GLuint normal_tex, const mat4& proj_matrix, float intensity, float radius, float bias) {
    ASSERT(glIsTexture(depth_tex));
    ASSERT(glIsTexture(normal_tex));

    ssao::setup_ubo_hbao_data(ssao::ubo_hbao_data, ssao::tex_width, ssao::tex_height, proj_matrix, intensity, radius, bias);

    const float n = proj_matrix[3][2] / (proj_matrix[2][2] - 1.f);
    const float f = (proj_matrix[2][2] - 1.f) * n / (proj_matrix[2][2] + 1.f);
    const vec2 inv_res = vec2(1.f) / vec2(ssao::tex_width, ssao::tex_height);
    const float sharpness = ssao::compute_sharpness(radius);

    glBindVertexArray(gl.vao);

    GLint last_fbo;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    // GLint draw_buffers[8];
    // for (int i = 0; i < 8; i++) glGetIntegerv(GL_DRAW_BUFFER0 + i, draw_buffers + i);

    glViewport(0, 0, ssao::tex_width, ssao::tex_height);

    // LINEARIZE_DEPTH
    linearize_depth(depth_tex, n, f, is_orthographic_proj_matrix(proj_matrix));

    // RENDER HBAO
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ssao::fbo_hbao);
    glUseProgram(ssao::prog_hbao);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl.linearize_depth.texture);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, normal_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, ssao::tex_random);

    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ssao::ubo_hbao_data);
    glUniformBlockBinding(ssao::prog_hbao, ssao::uniform_block_index_hbao_control_buffer, 0);
    glUniform1i(ssao::uniform_loc_hbao_tex_linear_depth, 0);
    glUniform1i(ssao::uniform_loc_hbao_tex_normal, 1);
    glUniform1i(ssao::uniform_loc_hbao_tex_random, 2);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    // BLUR FIRST
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, ssao::fbo_blur);
    glUseProgram(ssao::prog_blur_first);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ssao::tex_hbao);
    glUniform1i(ssao::uniform_loc_blur_texture, 0);
    glUniform1f(ssao::uniform_loc_blur_sharpness, sharpness);
    glUniform2f(ssao::uniform_loc_blur_inv_res_dir, inv_res.x, 0);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    // BLEND RESULTS WITH PREVIOUSLY BOUND FBO AND ITS COLOR TEXTURE
    glEnable(GL_BLEND);
    glBlendFunc(GL_ZERO, GL_SRC_COLOR);
    glColorMask(1, 1, 1, 0);

    // BLUR SECOND
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);
    // glDrawBuffers(8, (GLenum*)draw_buffers);
    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
    glUseProgram(ssao::prog_blur_second);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, ssao::tex_blur);
    glUniform1i(ssao::uniform_loc_blur_texture, 0);
    glUniform1f(ssao::uniform_loc_blur_sharpness, sharpness);
    glUniform2f(ssao::uniform_loc_blur_inv_res_dir, 0, inv_res.y);

    glDrawArrays(GL_TRIANGLES, 0, 3);

    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glColorMask(1, 1, 1, 1);

    glBindVertexArray(0);
}

void render_deferred(GLuint depth_tex, GLuint color_tex, GLuint normal_tex, const mat4& inv_proj_matrix) {
    ASSERT(glIsTexture(depth_tex));
    ASSERT(glIsTexture(color_tex));
    ASSERT(glIsTexture(normal_tex));

    glUseProgram(deferred::prog_deferred);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, normal_tex);
    glUniform1i(deferred::uniform_loc_texture_depth, 0);
    glUniform1i(deferred::uniform_loc_texture_color, 1);
    glUniform1i(deferred::uniform_loc_texture_normal, 2);
    glUniformMatrix4fv(deferred::uniform_loc_inv_proj_mat, 1, GL_FALSE, &inv_proj_matrix[0][0]);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

void half_res_color_coc(GLuint linear_depth_tex, GLuint color_tex, float focus_point, float focus_scale) {
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    glViewport(0, 0, gl.tex_width / 2, gl.tex_height / 2);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gl.half_res.fbo);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, linear_depth_tex);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glUseProgram(gl.half_res.program);

    glUniform1i(gl.half_res.uniform_loc.tex_depth, 0);
    glUniform1i(gl.half_res.uniform_loc.tex_color, 1);
    glUniform1f(gl.half_res.uniform_loc.focus_point, focus_point);
    glUniform1f(gl.half_res.uniform_loc.focus_scale, focus_scale);

    // ASSUME THAT THE APPROPRIATE FS_QUAD VAO IS BOUND
    glDrawArrays(GL_TRIANGLES, 0, 3);

    glViewport(last_viewport[0], last_viewport[1], last_viewport[2], last_viewport[3]);
}

void apply_dof(GLuint depth_tex, GLuint color_tex, const mat4& proj_matrix, float focus_point, float focus_scale) {
    ASSERT(glIsTexture(depth_tex));
    ASSERT(glIsTexture(color_tex));

    const float n = proj_matrix[3][2] / (proj_matrix[2][2] - 1.f);
    const float f = (proj_matrix[2][2] - 1.f) * n / (proj_matrix[2][2] + 1.f);
    const bool ortho = is_orthographic_proj_matrix(proj_matrix);
    const vec2 pixel_size = vec2(1.f / gl.tex_width, 1.f / gl.tex_height);

    GLint last_fbo;
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &last_fbo);

    glBindVertexArray(gl.vao);

    linearize_depth(depth_tex, n, f, ortho);

    half_res_color_coc(gl.linearize_depth.texture, color_tex, focus_point, focus_scale);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gl.half_res.tex.color_coc);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, gl.linearize_depth.texture);

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, color_tex);

    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, last_fbo);

    glUseProgram(gl.bokeh_dof.program);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_half_res, 0);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_depth, 1);
    glUniform1i(gl.bokeh_dof.uniform_loc.tex_color, 2);
    glUniform2f(gl.bokeh_dof.uniform_loc.pixel_size, pixel_size.x, pixel_size.y);
    glUniform1f(gl.bokeh_dof.uniform_loc.focus_point, focus_point);
    glUniform1f(gl.bokeh_dof.uniform_loc.focus_scale, focus_scale);

    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
    glActiveTexture(GL_TEXTURE0);
}

void apply_tonemapping(GLuint color_tex) {
    ASSERT(glIsTexture(color_tex));
    glUseProgram(tonemapping::prog_tonemap);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, color_tex);
    glUniform1i(tonemapping::uniform_loc_texture, 0);
    glBindVertexArray(gl.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glUseProgram(0);
}

}  // namespace postprocessing
