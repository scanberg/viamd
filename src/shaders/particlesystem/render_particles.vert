#version 430 core

layout (std430, binding = 0) readonly buffer particle_buffer {
    vec4 particles[];
};

layout (std140, binding = 0) uniform render_params {
    mat4 texture_to_world;
    mat4 model_view_proj;
    vec4 particle_color;
    float particle_size;
    float max_lifetime;
    uint _pad[2];
};

out VS_OUT {
    vec4 color;
} vs_out;

void main() {
    vec4 particle = particles[gl_VertexID];
    
    // Normalize lifetime to [0, 1] based on max_lifetime
    // Particle.w stores remaining lifetime, so higher value = younger
    float alpha = clamp(particle.w / max_lifetime, 0.0, 1.0);
    
    gl_Position = model_view_proj * vec4(particle.xyz, 1.0);
    gl_PointSize = particle_size;
    vs_out.color = vec4(particle_color.rgb, particle_color.a * alpha);
}
