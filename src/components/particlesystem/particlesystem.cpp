#include <event.h>
#include <viamd.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>

#include <gfx/gl.h>
#include <gfx/gl_utils.h>

#include <imgui_widgets.h>
#include <shaders.inl>

namespace particlesystem {

struct ParticleParams {
    mat4_t   texture_to_world;
    mat4_t   world_to_texture;
    float    dt;
    float    scalar_min;
    float    grad_magn_min;
	float    min_lifetime_in_frames;
	float	 max_lifetime_in_frames;
    uint32_t num_particles;
    uint32_t seed;
    uint32_t max_attempts;
};

struct RenderParams {
    mat4_t model_view_proj;
    vec4_t particle_color;
    float particle_size;
    float max_lifetime_in_frames;
    uint32_t _pad[2];
};

struct ParticleSystem : viamd::EventHandler {
    bool show_window = false;
    bool enabled = false;
    bool initialized = false;
    
    // Particle parameters
    uint32_t num_particles = 10000;
    float min_lifetime_in_frames = 100.0f;
    float max_lifetime_in_frames = 1000.0f;
    float scalar_min = 0.02f;
    float grad_magn_min = 0.01f;
    float timestep = 10.0f;
    float particle_size = 2.0f;
    vec4_t particle_color = {1.0f, 1.0f, 1.0f, 0.8f};
    uint32_t max_attempts = 64;
    
    // OpenGL resources
    GLuint particle_buffer = 0;
    GLuint params_ubo = 0;
    GLuint render_ubo = 0;
    GLuint render_vao = 0;
    
    GLuint advect_program = 0;
    GLuint render_program = 0;
    
    // Volume texture reference
    GLuint volume_texture = 0;
    int volume_dim[3] = {0, 0, 0};
    mat4_t volume_texture_to_world = mat4_ident();
    mat4_t volume_world_to_texture = mat4_ident();
    
    // Frame counter for seed diversity
    uint32_t frame_counter = 0;
    
    ParticleSystem() {
        viamd::event_system_register_handler(*this);
    }
    
    virtual ~ParticleSystem() {
        cleanup_gl_resources();
    }
    
    void initialize_gl_resources() {
        // Compile compute shader
        GLuint advect_shader = gl::compile_shader_from_source(
            {(const char*)advect_particles_comp, advect_particles_comp_size},
            GL_COMPUTE_SHADER
        );
        
        // Create compute program
        advect_program = glCreateProgram();
        if (!gl::attach_link_detach(advect_program, &advect_shader, 1)) {
            glDeleteProgram(advect_program);
            advect_program = 0;
        }
        glDeleteShader(advect_shader);
        
        // Compile render shaders
        GLuint vert_shader = gl::compile_shader_from_source( {(const char*)render_particles_vert, render_particles_vert_size}, GL_VERTEX_SHADER );
        GLuint frag_shader = gl::compile_shader_from_source( {(const char*)render_particles_frag, render_particles_frag_size}, GL_FRAGMENT_SHADER );
        
        // Create render program
        render_program = glCreateProgram();
        {
            const GLuint shaders[] = { vert_shader, frag_shader };
            if (!gl::attach_link_detach(render_program, shaders, 2)) {
                glDeleteProgram(render_program);
                render_program = 0;
            }
        }
        glDeleteShader(vert_shader);
        glDeleteShader(frag_shader);
        
        // Create uniform buffer (shared by seed and advect)
        glGenBuffers(1, &params_ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, params_ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(ParticleParams), nullptr, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        glGenBuffers(1, &render_ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, render_ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(RenderParams), nullptr, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        // Create empty VAO (required by OpenGL, even though we use SSBO)
        glGenVertexArrays(1, &render_vao);
    }
    
    void cleanup_gl_resources() {
        if (particle_buffer) glDeleteBuffers(1, &particle_buffer);
        if (params_ubo) glDeleteBuffers(1, &params_ubo);
        if (render_ubo) glDeleteBuffers(1, &render_ubo);
        if (render_vao) glDeleteVertexArrays(1, &render_vao);
        if (advect_program) glDeleteProgram(advect_program);
        if (render_program) glDeleteProgram(render_program);
    }
    
    void initialize_particles() {
        if (initialized && particle_buffer) {
            glDeleteBuffers(1, &particle_buffer);
        }
        
        // Create particle buffer with zero-initialized data
        // All particles start with life = 0, so they'll be reseeded on first update
        std::vector<vec4_t> initial_particles(num_particles, vec4_t{0, 0, 0, 0});
        
        glGenBuffers(1, &particle_buffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, particle_buffer);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(vec4_t) * num_particles, initial_particles.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        
        initialized = true;
    }
    
    void update_particles(float dt) {
        if (!initialized || !volume_texture || !advect_program) return;

        PUSH_GPU_SECTION("Update Particle System")
        
        // Increment frame counter for seed diversity
        frame_counter++;
        
        ParticleParams params = {};
        params.texture_to_world = volume_texture_to_world;
        params.world_to_texture = volume_world_to_texture;
        params.dt = dt;
        params.scalar_min = scalar_min;
        params.grad_magn_min = grad_magn_min;
        params.min_lifetime_in_frames = min_lifetime_in_frames;
        params.max_lifetime_in_frames = max_lifetime_in_frames;
        params.num_particles = num_particles;
        params.seed = (frame_counter * 1664525u + 1013904223u);
        params.max_attempts = max_attempts;
        
        glBindBuffer(GL_UNIFORM_BUFFER, params_ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(ParticleParams), &params);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        glUseProgram(advect_program);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particle_buffer);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, params_ubo);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_3D, volume_texture);
        
        uint32_t group_count = (num_particles + 255) / 256;
        glDispatchCompute(group_count, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, 0);
        glUseProgram(0);

        POP_GPU_SECTION()
    }
    
    void render_particles(const mat4_t& view_proj_matrix) {
        if (!initialized || !render_program) return;
        
        // Debug: Check if we have particles
        if (num_particles == 0) {
            MD_LOG_DEBUG("Particle system: num_particles is 0");
            return;
        }
        
        PUSH_GPU_SECTION("Draw Particle System")

        RenderParams params = {};
        params.model_view_proj = view_proj_matrix;
        params.particle_color = particle_color;
        params.particle_size = particle_size;
        params.max_lifetime_in_frames = max_lifetime_in_frames;
        
        glBindBuffer(GL_UNIFORM_BUFFER, render_ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(RenderParams), &params);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        // Save and set OpenGL state
        GLboolean depth_test_enabled = glIsEnabled(GL_DEPTH_TEST);
        GLboolean depth_mask;
        glGetBooleanv(GL_DEPTH_WRITEMASK, &depth_mask);
        
        glDisable(GL_DEPTH_TEST);  // Don't depth test particles
        glDepthMask(GL_FALSE);      // Don't write to depth buffer
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_PROGRAM_POINT_SIZE);
        
        glUseProgram(render_program);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particle_buffer);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, render_ubo);
        
        glBindVertexArray(render_vao);
        glDrawArrays(GL_POINTS, 0, num_particles);
        glBindVertexArray(0);
        
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, 0);
        glUseProgram(0);
        
        // Restore OpenGL state
        glDisable(GL_PROGRAM_POINT_SIZE);
        glDisable(GL_BLEND);
        if (depth_test_enabled) glEnable(GL_DEPTH_TEST);
        glDepthMask(depth_mask);

        POP_GPU_SECTION()
    }
    
    void draw_ui(const ApplicationState& app_state) {
        if (!show_window) return;
        
        if (ImGui::Begin("Particle System", &show_window)) {
            if (ImGui::Checkbox("Enabled", &enabled)) {
                if (enabled && !initialized) {
                    initialize_particles();
                }
            }
            
            ImGui::Separator();
            ImGui::Text("Particle Parameters");

            // @TODO: Enlist all representations which have a volume texture and allow user to select one
            size_t num_representations = md_array_size(app_state.representation.reps);
            
            const int cap = 16;
            const char* volume_names[cap] = { "None" };
			int volume_rep_indices[cap] = { -1 };
            int len = 1;

            for (size_t i = 0; i < num_representations && i < cap - 1; ++i) {
				const Representation& rep = app_state.representation.reps[i];
                if (rep.electronic_structure.vol.format == VolumeFormat::R16G16B16A16_FLOAT || 
                    rep.electronic_structure.vol.format == VolumeFormat::R32G32B32A32_FLOAT) {
                    volume_names[len] = app_state.representation.reps[i].name;
					volume_rep_indices[len] = (int)i;
					len++;
                }
            }

			static int vol_index = 0;
            if (ImGui::Combo("Volume Texture", &vol_index, volume_names, len)) {
                if (vol_index != 0) {
					int rep_idx = volume_rep_indices[vol_index];
                    const Representation& rep = app_state.representation.reps[rep_idx];
                    volume_texture = rep.electronic_structure.vol.tex_id;
                    volume_dim[0]  = rep.electronic_structure.vol.dim[0];
                    volume_dim[1]  = rep.electronic_structure.vol.dim[1];
                    volume_dim[2]  = rep.electronic_structure.vol.dim[2];
					volume_texture_to_world = rep.electronic_structure.vol.texture_to_world;
                    volume_world_to_texture = mat4_inverse(volume_texture_to_world);
                } else {
                    volume_texture = 0;
                    volume_dim[0] = volume_dim[1] = volume_dim[2] = 0;
                    volume_texture_to_world = mat4_ident();
                    volume_world_to_texture = mat4_ident();
                }
                // Initialize particles when volume changes
                initialize_particles();
            }

            bool params_changed = false;
            params_changed |= ImGui::SliderInt("Num Particles", (int*)&num_particles, 1000, 1000000);
            ImGui::RangeSliderFloat("Lifetime Range", &min_lifetime_in_frames, &max_lifetime_in_frames, 100.0f, 10000.0f);
            ImGui::SliderFloat("Scalar Min", &scalar_min, 0.0f, 0.2f, "%.7f");
            ImGui::SliderFloat("Gradient Magnitude Min", &grad_magn_min, 0.0f, 0.2f, "%.7f");
            ImGui::SliderInt("Max Attempts", (int*)&max_attempts, 1, 256);

            if (params_changed) {
                initialize_particles();
            }
            
            ImGui::Separator();
            ImGui::Text("Rendering");
            ImGui::SliderFloat("Particle Size", &particle_size, 1.0f, 20.0f);
            ImGui::ColorEdit4("Particle Color", &particle_color.x);
            
            ImGui::Separator();
            ImGui::Text("Simulation");
            ImGui::SliderFloat("Timestep", &timestep, -1.0f, 1.0f);
            if (ImGui::Button("Reinitialize")) {
                initialize_particles();
            }
        }
        ImGui::End();
    }
    
    void process_events(const viamd::Event* events, size_t num_events) override {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event& event = events[i];
            
            switch (event.type) {
                case viamd::EventType_ViamdInitialize: {
                    // Component initialization
                    initialize_gl_resources();
                    break;
                }
                
                case viamd::EventType_ViamdShutdown: {
                    cleanup_gl_resources();
                    break;
                }
                
                case viamd::EventType_ViamdFrameTick: {
                    ASSERT(event.payload_type == viamd::EventPayloadType_ApplicationState);
                    ApplicationState* app_state = (ApplicationState*)event.payload;
                    if (enabled && initialized) {
                        update_particles(timestep);
                    }
                    if (show_window) {
                        draw_ui(*app_state);
                    }
                    break;
                }
                
                case viamd::EventType_ViamdRenderTransparent: {
                    if (enabled && initialized && event.payload_type == viamd::EventPayloadType_ApplicationState) {
                        auto* app_state = (ApplicationState*)event.payload;
                        if (app_state) {
                            // Compute view and projection matrices from camera
                            mat4_t view_proj = app_state->view.param.matrix.curr.proj * app_state->view.param.matrix.curr.view;
                            render_particles(view_proj);
                        }
                    }
                    break;
                }
                
                case viamd::EventType_ViamdWindowDrawMenu: {
                    if (ImGui::BeginMenu("Components")) {
                        ImGui::MenuItem("Particle System", nullptr, &show_window);
                        ImGui::EndMenu();
                    }
                    break;
                }
                
                default:
                    break;
            }
        }
    }
};

static ParticleSystem instance = {};

}  // namespace particlesystem
