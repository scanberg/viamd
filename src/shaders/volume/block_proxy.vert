#version 410 core

// Proxy geometry for the isosurface rays: one box per block of the volume, instanced, collapsed when
// the block cannot hold a surface for the current isovalues. Rasterized into a depth buffer twice, once
// keeping the nearest and once the farthest depth, it gives every pixel the span of its ray that
// reaches a block that matters. Testing the isovalues here means a new isovalue needs no rebuild.
//
// A block matters when an isovalue lies within its [min, max] (a surface may cross it) or when it lies
// entirely inside a surface that absorbs (tau > 0), which matters where the clip box cuts a lobe open.

#ifndef MAX_ISO
#define MAX_ISO 8
#endif

uniform sampler3D u_tex_minmax;
uniform mat4  u_model_to_clip;
uniform ivec3 u_block_dim;
uniform vec3  u_block_ext;      // extent of a block in model space
uniform vec3  u_clip_min;
uniform vec3  u_clip_max;

uniform float u_iso_values[MAX_ISO];
uniform float u_iso_tau[MAX_ISO];
uniform int   u_iso_count;

void main() {
    int id = gl_InstanceID;
    ivec3 b = ivec3(id % u_block_dim.x, (id / u_block_dim.x) % u_block_dim.y, id / (u_block_dim.x * u_block_dim.y));
    vec2 mm = texelFetch(u_tex_minmax, b, 0).xy;

    bool relevant = false;
    for (int i = 0; i < u_iso_count; ++i) {
        float v = u_iso_values[i];
        bool crossing = (mm.x <= v && v <= mm.y);
        bool interior = (u_iso_tau[i] > 0.0) && ((v >= 0.0) ? (mm.x >= v) : (mm.y <= v));
        relevant = relevant || crossing || interior;
    }

    vec3 lo = max(vec3(b) * u_block_ext, u_clip_min);
    vec3 hi = min(vec3(b + 1) * u_block_ext, u_clip_max);

    if (!relevant || any(greaterThanEqual(lo, hi))) {
        // Every vertex of the instance at the same point outside the clip volume: nothing is drawn
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
        return;
    }

    // A unit cube as a 14 vertex triangle strip, from the vertex index
    int  i = gl_VertexID;
    vec3 c = vec3((0x287A >> i) & 1, (0x02AF >> i) & 1, (0x31E3 >> i) & 1);
    gl_Position = u_model_to_clip * vec4(mix(lo, hi, c), 1.0);
}
