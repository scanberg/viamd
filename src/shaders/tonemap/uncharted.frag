#version 150 core

uniform sampler2D u_texture;
uniform float u_exposure = 1.0;
uniform float u_gamma = 2.2;

out vec4 out_frag;

// Source from here
// http://filmicworlds.com/blog/filmic-tonemapping-operators/

vec3 uncharted_tonemap(vec3 x) {
    const float A = 0.15;
    const float B = 0.50;
    const float C = 0.10;
    const float D = 0.20;
    const float E = 0.02;
    const float F = 0.30;
    return ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}

vec3 uncharted(vec3 c) {
    const float W = 11.2;
    c *= u_exposure;

    const float exposure_bias = 0.5;
    vec3 curr = uncharted_tonemap(exposure_bias * c);

    vec3 white_scale = vec3(1.0) / uncharted_tonemap(vec3(W));
    vec3 color = curr * white_scale;

    return pow(color, 1.0 / vec3(u_gamma));
}

void main() {
    vec4 color = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
    color.rgb = uncharted(color.rgb);
    out_frag = color;
}