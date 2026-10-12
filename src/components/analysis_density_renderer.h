#pragma once
#include <components/analysis_density.h>
#include <components/analysis_density_shaders.h>
#include <imgui.h>
#include <implot.h>
#include <gfx/gl.h>
#include <gfx/gl_utils.h>

namespace analysis_density {

// Scalar adapter for the same map/contour shaders used by the four-channel Rama
// renderer. Density preparation is domain-specific; color and isoline rendering
// are deliberately shared.
struct Renderer {
    GLuint program[2] = {};
    GLuint density_texture = 0;
    GLuint image_texture = 0;
    GLuint framebuffer = 0;
    GLuint vao = 0;
    static constexpr int image_size = 512;

    void destroy() {
        for (GLuint& p : program) { if (p) glDeleteProgram(p); p = 0; }
        glDeleteTextures(1, &density_texture);
        glDeleteTextures(1, &image_texture);
        glDeleteFramebuffers(1, &framebuffer);
        glDeleteVertexArrays(1, &vao);
        density_texture = image_texture = framebuffer = vao = 0;
    }

    bool initialize() {
        if (program[0] && program[1]) return true;
        const GLuint vertex = gl::compile_shader_from_source(vertex_shader, GL_VERTEX_SHADER);
        const str_t source[] = {map_shader, iso_shader};
        for (int i = 0; i < 2; ++i) {
            const GLuint fragment = gl::compile_shader_from_source(source[i], GL_FRAGMENT_SHADER);
            program[i] = glCreateProgram();
            const GLuint shaders[] = {vertex, fragment};
            gl::attach_link_detach(program[i], shaders, 2);
            glDeleteShader(fragment);
            GLint linked = 0;
            glGetProgramiv(program[i], GL_LINK_STATUS, &linked);
            if (!linked) {
                glDeleteShader(vertex);
                destroy();
                return false;
            }
        }
        glDeleteShader(vertex);
        glGenTextures(1, &density_texture);
        glGenTextures(1, &image_texture);
        glGenFramebuffers(1, &framebuffer);
        glGenVertexArrays(1, &vao);
        glBindTexture(GL_TEXTURE_2D, image_texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, image_size, image_size, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        for (GLuint tex : {density_texture, image_texture}) {
            glBindTexture(GL_TEXTURE_2D, tex);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        }
        return true;
    }

    // line_color is the colour of the isolines and of the filled levels, the colormap is for Colormap
    bool render(const Triangle& density, Display display, ImPlotColormap colormap, ImVec4 line_color = ImVec4(1,1,1,1)) {
        GLint old_program, old_vao, old_fbo, old_texture, old_active, viewport[4];
        glGetIntegerv(GL_CURRENT_PROGRAM, &old_program);
        glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &old_vao);
        glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &old_fbo);
        glGetIntegerv(GL_ACTIVE_TEXTURE, &old_active);
        glGetIntegerv(GL_VIEWPORT, viewport);
        glActiveTexture(GL_TEXTURE0);
        glGetIntegerv(GL_TEXTURE_BINDING_2D, &old_texture);
        const GLboolean blend = glIsEnabled(GL_BLEND), cull = glIsEnabled(GL_CULL_FACE);
        const GLboolean depth = glIsEnabled(GL_DEPTH_TEST), scissor = glIsEnabled(GL_SCISSOR_TEST);
        auto restore = [&]() {
            glBindTexture(GL_TEXTURE_2D, old_texture);
            glActiveTexture(old_active);
            glUseProgram(old_program);
            glBindVertexArray(old_vao);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, old_fbo);
            glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
            if (blend) glEnable(GL_BLEND);
            if (cull) glEnable(GL_CULL_FACE);
            if (depth) glEnable(GL_DEPTH_TEST);
            if (scissor) glEnable(GL_SCISSOR_TEST);
        };
        if (!initialize()) { restore(); return false; }
        const auto data = density.texture_data();
        const int dim = Triangle::resolution + 1;
        glBindTexture(GL_TEXTURE_2D, density_texture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, dim, dim, 0, GL_RED, GL_FLOAT, data.data());
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, framebuffer);
        glFramebufferTexture(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, image_texture, 0);
        const GLenum attachment = GL_COLOR_ATTACHMENT0;
        glDrawBuffers(1, &attachment);
        if (glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
            restore();
            return false;
        }
        glDisable(GL_BLEND);
        glDisable(GL_CULL_FACE);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_SCISSOR_TEST);
        glViewport(0, 0, image_size, image_size);
        const GLuint p = program[display == Colormap ? 0 : 1];
        glUseProgram(p);
        auto loc = [p](const char* name) { return glGetUniformLocation(p, name); };
        glUniform1i(loc("u_tex_den"), 0);
        glUniform1i(loc("u_triangle_domain"), 1);
        glUniform2f(loc("u_inv_res"), 1.0f / image_size, 1.0f / image_size);
        glUniform4f(loc("u_viewport"), 0, 0, 1, 1);
        const GLuint offsets[4] = {};
        if (display == Colormap) {
            ImVec4 colors[16];
            for (int i = 0; i < 16; ++i) {
                colors[i] = ImPlot::SampleColormap(i / 15.0f, colormap);
                if (i == 0) colors[i].w = 0;
            }
            const GLuint lengths[4] = {16,16,16,16};
            const float ranges[8] = {0,1,0,1,0,1,0,1};
            glUniform4fv(loc("u_map_colors"), 16, &colors[0].x);
            glUniform1uiv(loc("u_map_offset"), 4, offsets);
            glUniform1uiv(loc("u_map_length"), 4, lengths);
            glUniform2fv(loc("u_map_range"), 4, ranges);
        } else {
            ImVec4 levels[8], lines[8];
            float values[8];
            for (int i = 0; i < 8; ++i) {
                values[i] = i / 8.0f;
                // Contours alone, or the bands between them: the other is a layer of its own
                lines[i] = display == Isolines ? line_color : ImVec4(0,0,0,0);
                levels[i] = display == FilledLevels ? ImVec4(line_color.x, line_color.y, line_color.z, line_color.w * (0.1f + 0.5f * i / 7.0f)) : ImVec4(0,0,0,0);
                if (i == 0) levels[i].w = lines[i].w = 0;
            }
            const GLuint lengths[4] = {8,8,8,8};
            glUniform1fv(loc("u_iso_values"), 8, values);
            glUniform4fv(loc("u_iso_level_colors"), 8, &levels[0].x);
            glUniform4fv(loc("u_iso_contour_colors"), 8, &lines[0].x);
            glUniform1uiv(loc("u_iso_offset"), 4, offsets);
            glUniform1uiv(loc("u_iso_length"), 4, lengths);
            glUniform1f(loc("u_iso_contour_line_scale"), 1.0f);
        }
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        restore();
        return true;
    }

    void draw(float opacity) const {
        if (!image_texture) return;
        ImPlot::PlotImage("##Density", (ImTextureID)(intptr_t)image_texture,
            ImPlotPoint(0,0), ImPlotPoint(1,Triangle::height),
            ImVec2(0,1), ImVec2(1,0), ImVec4(1,1,1,opacity));
    }
};
}
