// Renderer implementation
// This file contains all OpenGL rendering code for VIAMD molecular visualization

#include "renderer.h"
#include "immediate_draw_utils.h"
#include "postprocessing_utils.h"
#include "volumerender_utils.h"
#include "gl_utils.h"

#include <viamd.h>
#include <color_utils.h>
#include <event.h>

#include <imgui.h>

#include <core/md_log.h>
#include <core/md_array.h>

// External allocators
extern md_allocator_i* frame_alloc;
extern md_allocator_i* persistent_alloc;

namespace renderer {

// Helper macros
#define PUSH_GPU_SECTION(lbl) { if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); }
#define POP_GPU_SECTION()     { if (glPopDebugGroup) glPopDebugGroup(); }

#define HIGHLIGHT_PULSE_TIME_SCALE  5.0
#define HIGHLIGHT_PULSE_ALPHA_SCALE 0.1

void initialize() {
    ::immediate::initialize();
}

void shutdown() {
    ::immediate::shutdown();
}

void set_viewport(int x, int y, int width, int height) {
    glViewport(x, y, width, height);
    glScissor(x, y, width, height);
}

void clear_framebuffer(const vec4_t& color, float depth, int stencil) {
    glClearColor(color.x, color.y, color.z, color.w);
    glClearDepth(depth);
    glClearStencil(stencil);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
}

void bind_default_framebuffer() {
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glDrawBuffer(GL_BACK);
}

void clear_gbuffer(GBuffer& gbuffer) {
    ::clear_gbuffer(&gbuffer);
}

void begin_frame(ApplicationState* state) {
    (void)state;
    // Frame initialization if needed
}

void end_frame(ApplicationState* state) {
    (void)state;
    // Frame cleanup if needed
}

void draw_representations_opaque(ApplicationState* data) {
    ASSERT(data);
    if (data->mold.sys.atom.count == 0) return;
    
    const size_t num_representations = md_array_size(data->representation.reps);
    if (num_representations == 0) return;

    md_array(md_gl_draw_op_t) draw_ops = 0;
    for (size_t i = 0; i < num_representations; ++i) {
        const Representation& rep = data->representation.reps[i];
        
        if (rep.type > RepresentationType::Cartoon) continue;
        
        if (rep.enabled && rep.type_is_valid) {
            md_gl_draw_op_t op = {
                .type = (md_gl_rep_type_t)rep.type,
                .args = {},
                .rep = data->representation.reps[i].md_rep,
                .model_matrix = NULL,
            };
            MEMCPY(&op.args, &rep.scale, sizeof(op.args));
            md_array_push(draw_ops, op, frame_alloc);
        }
    }
    
    md_gl_draw_args_t args = {
        .shaders = data->mold.gl_shaders,
        .draw_operations = {
            .count = (uint32_t)md_array_size(draw_ops),
            .ops = draw_ops,
        },
        .view_transform = {
            .view_matrix = &data->view.param.matrix.curr.view.elem[0][0],
            .proj_matrix = &data->view.param.matrix.curr.proj.elem[0][0],
            .prev_view_matrix = &data->view.param.matrix.prev.view.elem[0][0],
            .prev_proj_matrix = &data->view.param.matrix.prev.proj.elem[0][0],
        },
    };
    
    md_gl_draw(&args);
}

void draw_representations_transparent(ApplicationState* state) {
    ASSERT(state);
    if (state->mold.sys.atom.count == 0) return;
    
    const size_t num_representations = md_array_size(state->representation.reps);
    if (num_representations == 0) return;
    
    for (size_t i = 0; i < num_representations; ++i) {
        const Representation& rep = state->representation.reps[i];
        if (!rep.enabled) continue;
        if (rep.type != RepresentationType::ElectronicStructure) continue;
        
        const IsoDesc* iso = nullptr;
        
        switch (rep.electronic_structure.type) {
        case ElectronicStructureType::MolecularOrbital:
        case ElectronicStructureType::NaturalTransitionOrbitalParticle:
        case ElectronicStructureType::NaturalTransitionOrbitalHole:
            iso = &rep.electronic_structure.iso_psi;
            break;
        case ElectronicStructureType::MolecularOrbitalDensity:
        case ElectronicStructureType::NaturalTransitionOrbitalDensityParticle:
        case ElectronicStructureType::NaturalTransitionOrbitalDensityHole:
        case ElectronicStructureType::AttachmentDensity:
        case ElectronicStructureType::DetachmentDensity:
        case ElectronicStructureType::ElectronDensity:
            iso = &rep.electronic_structure.iso_den;
            break;
        default:
            ASSERT(false);
        }
        
        volume::RenderDesc desc = {
            .render_target = {
                .depth = state->gbuffer.tex.depth,
                .color  = state->gbuffer.tex.transparency,
                .width  = state->gbuffer.width,
                .height = state->gbuffer.height,
            },
            .texture = {
                .volume = rep.electronic_structure.vol.tex_id,
                .transfer_function = rep.electronic_structure.dvr.tf_tex,
            },
            .matrix = {
                .model = rep.electronic_structure.vol.texture_to_world,
                .view  = state->view.param.matrix.curr.view,
                .proj  = state->view.param.matrix.curr.proj,
                .inv_proj = state->view.param.matrix.inv.proj,
            },
            .clip_volume = {
                .min = {0,0,0},
                .max = {1,1,1},
            },
            .temporal = {
                .enabled = state->visuals.temporal_aa.enabled,
            },
            .iso = {
                .enabled = iso->enabled,
                .count   = iso->count,
                .values  = iso->values,
                .colors  = iso->colors,
            },
            .dvr = {
                .enabled = rep.electronic_structure.dvr.enabled,
                .min_tf_value = -1.0f,
                .max_tf_value =  1.0f,
            },
            .shading = {
                .env_radiance = state->visuals.background.color * state->visuals.background.intensity * 0.25f,
                .roughness = 0.3f,
                .dir_radiance = {10,10,10},
                .ior = 1.5f,
            },
            .voxel_spacing = rep.electronic_structure.vol.voxel_size,
        };
        
        volume::render_volume(desc);
    }
}

void draw_representations_opaque_lean_and_mean(ApplicationState* data, uint32_t mask) {
    md_gl_draw_op_t* draw_ops = 0;
    for (size_t i = 0; i < md_array_size(data->representation.reps); ++i) {
        const Representation& rep = data->representation.reps[i];
        
        if (rep.type > RepresentationType::Cartoon) continue;
        
        if (rep.enabled && rep.type_is_valid) {
            md_gl_draw_op_t op = {
                .type = (md_gl_rep_type_t)rep.type,
                .args = {},
                .rep = data->representation.reps[i].md_rep,
                .model_matrix = NULL,
            };
            MEMCPY(&op.args, &rep.scale, sizeof(op.args));
            md_array_push(draw_ops, op, frame_alloc);
        }
    }
    
    md_gl_draw_args_t args = {
        .shaders = data->mold.gl_shaders_lean_and_mean,
        .draw_operations = {
            .count = md_array_size(draw_ops),
            .ops = draw_ops,
        },
        .view_transform = {
            .view_matrix = &data->view.param.matrix.curr.view.elem[0][0],
            .proj_matrix = &data->view.param.matrix.curr.proj.elem[0][0],
        },
        .atom_mask = mask,
    };
    
    md_gl_draw(&args);
}

void fill_gbuffer(ApplicationState* data) {
    const GLenum draw_buffers[] = {
        GL_COLOR_ATTACHMENT_COLOR,
        GL_COLOR_ATTACHMENT_NORMAL,
        GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_PICKING,
        GL_COLOR_ATTACHMENT_TRANSPARENCY
    };
    
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    
    // Enable all draw buffers
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, data->gbuffer.fbo);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
    
    PUSH_GPU_SECTION("G-Buffer fill")
    
    // Draw simulation box if enabled
    if (data->simulation_box.enabled && data->mold.sys.unitcell.flags != 0) {
        PUSH_GPU_SECTION("Draw Simulation Box")
        const mat4_t basis_model_mat = md_unitcell_basis_mat4(&data->mold.sys.unitcell);
        ::immediate::set_model_view_matrix(data->view.param.matrix.curr.view);
        ::immediate::set_proj_matrix(data->view.param.matrix.curr.proj);
        ::immediate::draw_box_wireframe({0,0,0}, {1,1,1}, basis_model_mat, convert_color(data->simulation_box.color));
        ::immediate::render();
        POP_GPU_SECTION()
    }
    
    // Draw velocity of static objects
    PUSH_GPU_SECTION("Blit Static Velocity")
    glDrawBuffer(GL_COLOR_ATTACHMENT_VELOCITY);
    glDepthMask(0);
    postprocessing::blit_static_velocity(data->gbuffer.tex.depth, data->view.param);
    glDepthMask(1);
    POP_GPU_SECTION()
    
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);
    
    // Draw representations
    PUSH_GPU_SECTION("Draw Opaque")
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
    draw_representations_opaque(data);
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderOpaque, viamd::EventPayloadType_ApplicationState, data);
    POP_GPU_SECTION()
    
    glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
    
    // Selection and highlight overlays
    PUSH_GPU_SECTION("Selection")
    const bool atom_selection_empty = md_bitfield_popcount(&data->selection.selection_mask) == 0;
    const bool atom_highlight_empty = md_bitfield_popcount(&data->selection.highlight_mask) == 0;
    
    glDepthMask(0);
    
    if (!atom_selection_empty) {
        glColorMask(0, 0, 0, 0);
        
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_EQUAL);
        
        glEnable(GL_STENCIL_TEST);
        glStencilMask(0xFF);
        
        glClearStencil(1);
        glClear(GL_STENCIL_BUFFER_BIT);
        
        glStencilFunc(GL_GREATER, 0x02, 0xFF);
        glStencilOp(GL_KEEP, GL_ZERO, GL_REPLACE);
        draw_representations_opaque_lean_and_mean(data, AtomBit_Selected | AtomBit_Visible);
        
        glDisable(GL_DEPTH_TEST);
        
        glStencilMask(0x0);
        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
        glColorMask(1, 1, 1, 1);
        
        glStencilFunc(GL_EQUAL, 2, 0xFF);
        postprocessing::blit_color(data->selection.color.selection.visible);
        
        glStencilFunc(GL_EQUAL, 0, 0xFF);
        postprocessing::blit_color(data->selection.color.selection.hidden);
    }
    
    if (!atom_highlight_empty) {
        glColorMask(0, 0, 0, 0);
        
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_EQUAL);
        
        glEnable(GL_STENCIL_TEST);
        glStencilMask(0xFF);
        
        glClearStencil(1);
        glClear(GL_STENCIL_BUFFER_BIT);
        
        glStencilFunc(GL_GREATER, 0x02, 0xFF);
        glStencilOp(GL_KEEP, GL_ZERO, GL_REPLACE);
        draw_representations_opaque_lean_and_mean(data, AtomBit_Highlighted | AtomBit_Visible);
        
        glDisable(GL_DEPTH_TEST);
        
        glStencilMask(0x0);
        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
        glColorMask(1, 1, 1, 1);
        
        glStencilFunc(GL_EQUAL, 2, 0xFF);
        vec4_t col_vis = data->selection.color.highlight.visible;
        col_vis.w += sin(ImGui::GetTime() * HIGHLIGHT_PULSE_TIME_SCALE) * HIGHLIGHT_PULSE_ALPHA_SCALE;
        postprocessing::blit_color(col_vis);
        
        glStencilFunc(GL_EQUAL, 0, 0xFF);
        postprocessing::blit_color(data->selection.color.highlight.hidden);
    }
    
    glDisable(GL_STENCIL_TEST);
    
    if (!atom_selection_empty) {
        PUSH_GPU_SECTION("Desaturate")
        const float saturation = data->selection.color.saturation;
        glDrawBuffer(GL_COLOR_ATTACHMENT_COLOR);
        postprocessing::scale_hsv(data->gbuffer.tex.color, vec3_t{1, saturation, 1});
        POP_GPU_SECTION()
    }
    
    glDepthFunc(GL_LESS);
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);
    POP_GPU_SECTION()
    
    // Draw transparent representations
    PUSH_GPU_SECTION("Draw Transparent")
    draw_representations_transparent(data);
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderTransparent, viamd::EventPayloadType_ApplicationState, data);
    POP_GPU_SECTION()
    
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glDisable(GL_DEPTH_TEST);
    POP_GPU_SECTION()
}

void apply_postprocessing(ApplicationState* state) {
    PUSH_GPU_SECTION("Postprocessing")
    postprocessing::Descriptor desc;

    desc.background.color = state->visuals.background.color * state->visuals.background.intensity;

    desc.ambient_occlusion.enabled = state->visuals.ssao.enabled;
    desc.ambient_occlusion.intensity = state->visuals.ssao.intensity;
    desc.ambient_occlusion.radius = state->visuals.ssao.radius;
    desc.ambient_occlusion.bias = state->visuals.ssao.bias;

    desc.tonemapping.enabled = state->visuals.tonemapping.enabled;
    desc.tonemapping.mode = state->visuals.tonemapping.tonemapper;
    desc.tonemapping.exposure = state->visuals.tonemapping.exposure;
    desc.tonemapping.gamma = state->visuals.tonemapping.gamma;

    desc.depth_of_field.enabled = state->visuals.dof.enabled;
    desc.depth_of_field.focus_depth = state->visuals.dof.focus_depth;
    desc.depth_of_field.focus_scale = state->visuals.dof.focus_scale;

    desc.fxaa.enabled = state->visuals.fxaa.enabled;

    constexpr float MOTION_BLUR_REFERENCE_DT = 1.0f / 60.0f;
    const float dt_compensation = MOTION_BLUR_REFERENCE_DT / (float)state->app.timing.delta_s;
    const float motion_scale = state->visuals.temporal_aa.motion_blur.motion_scale * dt_compensation;
    desc.temporal_aa.enabled = state->visuals.temporal_aa.enabled;
    desc.temporal_aa.feedback_min = state->visuals.temporal_aa.feedback_min;
    desc.temporal_aa.feedback_max = state->visuals.temporal_aa.feedback_max;
    desc.temporal_aa.motion_blur.enabled = state->visuals.temporal_aa.motion_blur.enabled;
    desc.temporal_aa.motion_blur.motion_scale = motion_scale;

    desc.sharpen.enabled = state->visuals.temporal_aa.enabled && state->visuals.sharpen.enabled;
    desc.sharpen.weight  = state->visuals.sharpen.weight;

    desc.input_textures.depth = state->gbuffer.tex.depth;
    desc.input_textures.color = state->gbuffer.tex.color;
    desc.input_textures.normal = state->gbuffer.tex.normal;
    desc.input_textures.velocity = state->gbuffer.tex.velocity;
    desc.input_textures.transparency = state->gbuffer.tex.transparency;

    postprocessing::shade_and_postprocess(desc, state->view.param);
    POP_GPU_SECTION()
}

void render_scene(ApplicationState* state) {
    clear_gbuffer(state->gbuffer);
    fill_gbuffer(state);
}

namespace immediate {
    void begin(const mat4_t& view_matrix, const mat4_t& proj_matrix) {
        ::immediate::set_model_view_matrix(view_matrix);
        ::immediate::set_proj_matrix(proj_matrix);
    }
    
    void end() {
        ::immediate::render();
    }
    
    void flush() {
        ::immediate::render();
    }
}

} // namespace renderer
