// Renderer - Abstracts all rendering operations for VIAMD
// This module encapsulates OpenGL rendering calls and provides a clean interface
// for rendering molecular representations, immediate graphics, and post-processing

#pragma once

#include <core/md_vec_math.h>
#include <md_gl.h>

#include "gl.h"
#include "view_param.h"
#include "postprocessing_utils.h"

// Forward declarations
struct ApplicationState;
struct Representation;

namespace renderer {

// Initialization and cleanup
void initialize();
void shutdown();

// Main rendering entry points
void begin_frame(ApplicationState* state);
void end_frame(ApplicationState* state);

// Rendering passes
void render_scene(ApplicationState* state);
void render_molecular_representations(ApplicationState* state);
void render_volume_representations(ApplicationState* state);
void render_immediate_graphics(ApplicationState* state);
void render_ui_overlay(ApplicationState* state);

// Viewport and framebuffer management
void set_viewport(int x, int y, int width, int height);
void clear_framebuffer(const vec4_t& color, float depth = 1.0f, int stencil = 0);
void bind_default_framebuffer();

// GBuffer operations
void clear_gbuffer(GBuffer& gbuffer);
void fill_gbuffer(ApplicationState* state);

// Representation rendering
void draw_representations_opaque(ApplicationState* state);
void draw_representations_transparent(ApplicationState* state);
void draw_representations_opaque_lean_and_mean(ApplicationState* state, uint32_t mask);

// Post-processing
void apply_postprocessing(ApplicationState* state);

// Immediate mode rendering
namespace immediate {
    void begin(const mat4_t& view_matrix, const mat4_t& proj_matrix);
    void end();
    void flush();
    
    void draw_simulation_box(ApplicationState* state);
    void draw_selection_overlay(ApplicationState* state);
    void draw_highlight_overlay(ApplicationState* state);
}

} // namespace renderer
