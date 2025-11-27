#version 430 core

in VS_OUT {
    vec4 color;
} fs_in;

layout (location = 0) out vec4 out_color;

void main() {
    // Create a smooth circular point sprite
    //vec2 coord = gl_PointCoord - vec2(0.5);
    //float dist = length(coord);
    
    //if (dist > 0.5) {
    //    discard;
    //}
    
    // Smooth falloff at edges
    //float alpha_falloff = 1.0 - smoothstep(0.3, 0.5, dist);
    float alpha = fs_in.color.a; // * alpha_falloff;
    
    out_color = vec4(fs_in.color.rgb, alpha);
}
