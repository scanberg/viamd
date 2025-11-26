#version 430 core

layout (location = 0) in uint particle_id;

struct Particle {
    vec4 position;      // xyz: position, w: age
    vec4 velocity;      // xyz: velocity, w: lifetime
    vec4 trail_pos[8];  // Trail positions for rendering
};

layout (std430, binding = 0) readonly buffer particle_buffer {
    Particle particles[];
};

layout (std140, binding = 0) uniform render_params {
    mat4 model_view_proj;
    vec3 volume_min;
    float trail_width;
    vec3 volume_max;
    float _pad0;
    vec4 particle_color;
    uint num_particles;
    uint num_trail_segments;
    uint _pad1;
    uint _pad2;
};

out VS_OUT {
    vec4 color;
    float age_norm;
} vs_out;

void main() {
    uint particle_idx = particle_id / (num_trail_segments + 1);
    uint segment_idx = particle_id % (num_trail_segments + 1);
    
    if (particle_idx >= num_particles) {
        gl_Position = vec4(0.0);
        return;
    }
    
    vec3 pos;
    float alpha;
    
    if (segment_idx == 0) {
        // Current particle position
        pos = particles[particle_idx].position.xyz;
        alpha = 1.0;
    } else {
        // Trail position
        uint trail_idx = segment_idx - 1;
        if (trail_idx >= 8) {
            gl_Position = vec4(0.0);
            return;
        }
        pos = particles[particle_idx].trail_pos[trail_idx].xyz;
        // Fade out along the trail (trail_idx goes from 0 to min(7, num_trail_segments-1))
        alpha = 1.0 - float(trail_idx) / 8.0;
    }
    
    // Compute normalized age for color variation
    float age = particles[particle_idx].position.w;
    float lifetime = particles[particle_idx].velocity.w;
    float age_norm = clamp(age / lifetime, 0.0, 1.0);
    
    gl_Position = model_view_proj * vec4(pos, 1.0);
    gl_PointSize = trail_width;
    
    vs_out.color = vec4(particle_color.rgb, particle_color.a * alpha);
    vs_out.age_norm = age_norm;
}
