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

struct Particle {
    vec4_t position;       // xyz: position, w: age
    vec4_t velocity;       // xyz: velocity, w: lifetime
    vec4_t trail_pos[8];   // Trail positions
};

struct SeedParams {
    uint32_t volume_dim[3];
    uint32_t num_particles;
    vec3_t volume_min;
    float scalar_min;
    vec3_t volume_max;
    float scalar_max;
    float min_lifetime;
    float max_lifetime;
    uint32_t random_seed;
    float _pad0;
};

struct AdvectParams {
    uint32_t volume_dim[3];
    uint32_t num_particles;
    vec3_t volume_min;
    float dt;
    vec3_t volume_max;
    float scalar_min;
    float scalar_max;
    float min_lifetime;
    float max_lifetime;
    uint32_t random_seed;
};

struct RenderParams {
    mat4_t model_view_proj;
    vec3_t volume_min;
    float trail_width;
    vec3_t volume_max;
    float _pad0;
    vec4_t particle_color;
    uint32_t num_particles;
    uint32_t num_trail_segments;
    uint32_t _pad1;
    uint32_t _pad2;
};

struct ParticleSystem : viamd::EventHandler {
    bool show_window = false;
    bool enabled = false;
    bool initialized = false;
    
    // Particle parameters
    uint32_t num_particles = 1000;
    float min_lifetime = 1.0f;
    float max_lifetime = 5.0f;
    float scalar_min = 0.2f;
    float scalar_max = 0.8f;
    float timestep = 0.016f;
    float trail_width = 5.0f;
    vec4_t particle_color = {1.0f, 1.0f, 1.0f, 0.8f};
    uint32_t num_trail_segments = 8;
    
    // OpenGL resources
    GLuint particle_buffer = 0;
    GLuint seed_ubo = 0;
    GLuint advect_ubo = 0;
    GLuint render_ubo = 0;
    GLuint render_vao = 0;
    GLuint render_vbo = 0;
    
    GLuint seed_program = 0;
    GLuint advect_program = 0;
    GLuint render_program = 0;
    
    // Volume texture reference
    GLuint volume_texture = 0;
    int volume_dim[3] = {0, 0, 0};
    vec3_t volume_min = {0, 0, 0};
    vec3_t volume_max = {1, 1, 1};
    
    // Frame counter for seed diversity
    uint32_t frame_counter = 0;
    
    ParticleSystem() {
        viamd::event_system_register_handler(*this);
        initialize_gl_resources();
    }
    
    virtual ~ParticleSystem() {
        cleanup_gl_resources();
    }
    
    void initialize_gl_resources() {
        // Compile compute shaders
        GLuint seed_shader = gl::compile_shader_from_source(
            {(const char*)seed_particles_comp, seed_particles_comp_size},
            GL_COMPUTE_SHADER
        );
        GLuint advect_shader = gl::compile_shader_from_source(
            {(const char*)advect_particles_comp, advect_particles_comp_size},
            GL_COMPUTE_SHADER
        );
        
        // Create compute programs
        seed_program = glCreateProgram();
        glAttachShader(seed_program, seed_shader);
        glLinkProgram(seed_program);
        glDetachShader(seed_program, seed_shader);
        glDeleteShader(seed_shader);
        
        advect_program = glCreateProgram();
        glAttachShader(advect_program, advect_shader);
        glLinkProgram(advect_program);
        glDetachShader(advect_program, advect_shader);
        glDeleteShader(advect_shader);
        
        // Compile render shaders
        GLuint vert_shader = gl::compile_shader_from_source(
            {(const char*)render_particles_vert, render_particles_vert_size},
            GL_VERTEX_SHADER
        );
        GLuint frag_shader = gl::compile_shader_from_source(
            {(const char*)render_particles_frag, render_particles_frag_size},
            GL_FRAGMENT_SHADER
        );
        
        // Create render program
        render_program = glCreateProgram();
        glAttachShader(render_program, vert_shader);
        glAttachShader(render_program, frag_shader);
        glLinkProgram(render_program);
        glDetachShader(render_program, vert_shader);
        glDetachShader(render_program, frag_shader);
        glDeleteShader(vert_shader);
        glDeleteShader(frag_shader);
        
        // Create uniform buffers
        glGenBuffers(1, &seed_ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, seed_ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(SeedParams), nullptr, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        glGenBuffers(1, &advect_ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, advect_ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(AdvectParams), nullptr, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        glGenBuffers(1, &render_ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, render_ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(RenderParams), nullptr, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        // Create VAO for rendering
        glGenVertexArrays(1, &render_vao);
        glBindVertexArray(render_vao);
        
        glGenBuffers(1, &render_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, render_vbo);
        
        glEnableVertexAttribArray(0);
        glVertexAttribIPointer(0, 1, GL_UNSIGNED_INT, sizeof(uint32_t), 0);
        
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void cleanup_gl_resources() {
        if (particle_buffer) glDeleteBuffers(1, &particle_buffer);
        if (seed_ubo) glDeleteBuffers(1, &seed_ubo);
        if (advect_ubo) glDeleteBuffers(1, &advect_ubo);
        if (render_ubo) glDeleteBuffers(1, &render_ubo);
        if (render_vbo) glDeleteBuffers(1, &render_vbo);
        if (render_vao) glDeleteVertexArrays(1, &render_vao);
        if (seed_program) glDeleteProgram(seed_program);
        if (advect_program) glDeleteProgram(advect_program);
        if (render_program) glDeleteProgram(render_program);
    }
    
    void initialize_particles() {
        if (initialized && particle_buffer) {
            glDeleteBuffers(1, &particle_buffer);
        }
        
        // Create particle buffer
        glGenBuffers(1, &particle_buffer);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, particle_buffer);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(Particle) * num_particles, nullptr, GL_DYNAMIC_DRAW);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        
        // Seed particles
        if (volume_texture && seed_program) {
            SeedParams params;
            params.volume_dim[0] = volume_dim[0];
            params.volume_dim[1] = volume_dim[1];
            params.volume_dim[2] = volume_dim[2];
            params.num_particles = num_particles;
            params.volume_min = volume_min;
            params.volume_max = volume_max;
            params.scalar_min = scalar_min;
            params.scalar_max = scalar_max;
            params.min_lifetime = min_lifetime;
            params.max_lifetime = max_lifetime;
            // Use frame counter for better seed diversity across rapid calls
            params.random_seed = (frame_counter * 1664525u + 1013904223u);
            
            glBindBuffer(GL_UNIFORM_BUFFER, seed_ubo);
            glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(SeedParams), &params);
            glBindBuffer(GL_UNIFORM_BUFFER, 0);
            
            glUseProgram(seed_program);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particle_buffer);
            glBindBufferBase(GL_UNIFORM_BUFFER, 0, seed_ubo);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_3D, volume_texture);
            
            uint32_t group_count = (num_particles + 255) / 256;
            glDispatchCompute(group_count, 1, 1);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
            glBindBufferBase(GL_UNIFORM_BUFFER, 0, 0);
            glUseProgram(0);
        }
        
        // Prepare render VBO with particle indices (sequential 0, 1, 2, ...)
        // Each particle has (num_trail_segments + 1) vertices
        uint32_t vertex_count = num_particles * (num_trail_segments + 1);
        uint32_t* indices = new uint32_t[vertex_count];
        for (uint32_t i = 0; i < vertex_count; ++i) {
            indices[i] = i;
        }
        
        glBindBuffer(GL_ARRAY_BUFFER, render_vbo);
        glBufferData(GL_ARRAY_BUFFER, vertex_count * sizeof(uint32_t), indices, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        
        delete[] indices;
        initialized = true;
    }
    
    void update_particles(float dt) {
        if (!initialized || !volume_texture || !advect_program) return;
        
        // Increment frame counter for seed diversity
        frame_counter++;
        
        AdvectParams params;
        params.volume_dim[0] = volume_dim[0];
        params.volume_dim[1] = volume_dim[1];
        params.volume_dim[2] = volume_dim[2];
        params.num_particles = num_particles;
        params.volume_min = volume_min;
        params.volume_max = volume_max;
        params.dt = dt;
        params.scalar_min = scalar_min;
        params.scalar_max = scalar_max;
        params.min_lifetime = min_lifetime;
        params.max_lifetime = max_lifetime;
        // Use frame counter for better seed diversity
        params.random_seed = (frame_counter * 1664525u + 1013904223u);
        
        glBindBuffer(GL_UNIFORM_BUFFER, advect_ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(AdvectParams), &params);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        glUseProgram(advect_program);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particle_buffer);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, advect_ubo);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_3D, volume_texture);
        
        uint32_t group_count = (num_particles + 255) / 256;
        glDispatchCompute(group_count, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, 0);
        glUseProgram(0);
    }
    
    void render_particles(const mat4_t& view_proj_matrix) {
        if (!initialized || !render_program) return;
        
        RenderParams params;
        params.model_view_proj = view_proj_matrix;
        params.volume_min = volume_min;
        params.volume_max = volume_max;
        params.trail_width = trail_width;
        params.particle_color = particle_color;
        params.num_particles = num_particles;
        params.num_trail_segments = num_trail_segments;
        
        glBindBuffer(GL_UNIFORM_BUFFER, render_ubo);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(RenderParams), &params);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_PROGRAM_POINT_SIZE);
        
        glUseProgram(render_program);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particle_buffer);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, render_ubo);
        
        glBindVertexArray(render_vao);
        uint32_t vertex_count = num_particles * (num_trail_segments + 1);
        glDrawArrays(GL_POINTS, 0, vertex_count);
        glBindVertexArray(0);
        
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, 0);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, 0);
        glUseProgram(0);
        
        glDisable(GL_PROGRAM_POINT_SIZE);
        glDisable(GL_BLEND);
    }
    
    void draw_ui() {
        if (!show_window) return;
        
        if (ImGui::Begin("Particle System", &show_window)) {
            ImGui::Checkbox("Enabled", &enabled);
            
            ImGui::Separator();
            ImGui::Text("Particle Parameters");
            
            bool params_changed = false;
            params_changed |= ImGui::SliderInt("Num Particles", (int*)&num_particles, 100, 10000);
            params_changed |= ImGui::SliderFloat("Min Lifetime", &min_lifetime, 0.1f, 10.0f);
            params_changed |= ImGui::SliderFloat("Max Lifetime", &max_lifetime, 0.1f, 10.0f);
            params_changed |= ImGui::SliderFloat("Scalar Min", &scalar_min, 0.0f, 1.0f);
            params_changed |= ImGui::SliderFloat("Scalar Max", &scalar_max, 0.0f, 1.0f);
            
            if (params_changed && initialized) {
                initialize_particles();
            }
            
            ImGui::Separator();
            ImGui::Text("Rendering");
            ImGui::SliderFloat("Trail Width", &trail_width, 1.0f, 20.0f);
            ImGui::ColorEdit4("Particle Color", &particle_color.x);
            
            ImGui::Separator();
            ImGui::Text("Simulation");
            ImGui::SliderFloat("Timestep", &timestep, 0.001f, 0.1f);
            
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
                    break;
                }
                
                case viamd::EventType_ViamdShutdown: {
                    cleanup_gl_resources();
                    break;
                }
                
                case viamd::EventType_ViamdFrameTick: {
                    if (enabled && initialized) {
                        update_particles(timestep);
                    }
                    break;
                }
                
                case viamd::EventType_ViamdRenderTransparent: {
                    if (enabled && initialized && event.payload_type == viamd::EventPayloadType_ApplicationState) {
                        auto* app_state = (ApplicationState*)event.payload;
                        if (app_state) {
                            // Compute view and projection matrices from camera
                            float aspect = (float)app_state->app.framebuffer.width / (float)app_state->app.framebuffer.height;
                            mat4_t view = camera_world_to_view_matrix(app_state->view.camera);
                            mat4_t proj = camera_perspective_projection_matrix(app_state->view.camera, aspect);
                            mat4_t view_proj = mat4_mul(proj, view);
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

static ParticleSystem* g_particle_system = nullptr;

void initialize() {
    if (!g_particle_system) {
        g_particle_system = new ParticleSystem();
    }
}

void shutdown() {
    if (g_particle_system) {
        delete g_particle_system;
        g_particle_system = nullptr;
    }
}

void set_volume_texture(GLuint texture, int dim_x, int dim_y, int dim_z, vec3_t min, vec3_t max) {
    if (g_particle_system) {
        g_particle_system->volume_texture = texture;
        g_particle_system->volume_dim[0] = dim_x;
        g_particle_system->volume_dim[1] = dim_y;
        g_particle_system->volume_dim[2] = dim_z;
        g_particle_system->volume_min = min;
        g_particle_system->volume_max = max;
        
        if (texture && !g_particle_system->initialized) {
            g_particle_system->initialize_particles();
        }
    }
}

void draw_ui() {
    if (g_particle_system) {
        g_particle_system->draw_ui();
    }
}

}  // namespace particlesystem
