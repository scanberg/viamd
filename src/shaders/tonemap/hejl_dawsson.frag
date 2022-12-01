#version 150 core

uniform sampler2D u_texture;
uniform float u_exposure = 1.0;
uniform float u_gamma = 2.2;

out vec4 out_frag;

// Source from
// https://de.slideshare.net/hpduiker/filmic-tonemapping-for-realtime-rendering-siggraph-2010-color-course

vec3 hejl_dawsson(vec3 c) {
   c *= u_exposure;
   vec3 x = max(vec3(0), c-vec3(0.004));
   return (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
}

void main() {
    vec4 color = texelFetch(u_texture, ivec2(gl_FragCoord.xy), 0);
    color.rgb = hejl_dawsson(color.rgb);
    out_frag = color;
}