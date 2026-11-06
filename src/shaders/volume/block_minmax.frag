#version 410 core

// One texel per block of the volume: the min and max of the voxels that the trilinear interpolation
// anywhere inside the block can reach, which is the block plus a one voxel apron. Any value sampled
// inside the block lies in [min, max], so a block whose range holds no isovalue has no surface in it.
// Rendered one layer (z) of blocks at a time.

uniform sampler3D u_tex_volume;
uniform int u_layer;
uniform int u_block_size;

layout(location = 0) out vec2 out_minmax;

void main() {
    ivec3 dim = textureSize(u_tex_volume, 0);
    ivec3 b   = ivec3(ivec2(gl_FragCoord.xy), u_layer);
    ivec3 lo  = max(b * u_block_size - 1, ivec3(0));
    ivec3 hi  = min((b + 1) * u_block_size, dim - 1);

    float mn =  3.402823e38;
    float mx = -3.402823e38;
    for (int z = lo.z; z <= hi.z; ++z) {
        for (int y = lo.y; y <= hi.y; ++y) {
            for (int x = lo.x; x <= hi.x; ++x) {
                float v = texelFetch(u_tex_volume, ivec3(x, y, z), 0).r;
                mn = min(mn, v);
                mx = max(mx, v);
            }
        }
    }
    out_minmax = vec2(mn, mx);
}
