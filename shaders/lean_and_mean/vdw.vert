#version 150 core
#extension GL_ARB_explicit_attrib_location : enable

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;
uniform vec4 u_jitter_uv;

uniform float u_radius_scale = 1.0;
uniform vec4 u_color;
uniform uint u_mask;

layout (location = 0) in vec3  in_position;
layout (location = 1) in float in_radius;
layout (location = 2) in vec4  in_color;
layout (location = 3) in uint  in_mask;

out VS_GS {
    flat vec4 view_sphere;
    flat vec2 axis_a;
    flat vec2 axis_b;
    flat vec2 center;
    flat float z;
} out_geom;

// From Inigo Quilez!
void proj_sphere(in vec4 sphere, 
                 in float fle,
                 out vec2 axis_a,
                 out vec2 axis_b,
                 out vec2 center) {
    vec3  o = sphere.xyz;
    float r2 = sphere.w*sphere.w;
    float z2 = o.z*o.z; 
    float l2 = dot(o,o);
    
    // axis
    axis_a = fle*sqrt(-r2*(r2-l2)/((l2-z2)*(r2-z2)*(r2-z2)))*vec2( o.x,o.y);
    axis_b = fle*sqrt(-r2*(r2-l2)/((l2-z2)*(r2-z2)*(r2-l2)))*vec2(-o.y,o.x);
    center = -fle*o.z*o.xy/(z2-r2);
}

void main() {
    uint ref_mask = u_mask;
    vec4 color = in_color;
    if ((in_mask & ref_mask) == 0U || color.a == 0.0) {
        out_geom.view_sphere = vec4(0,0,0,0);
    } else {
        vec3 pos = in_position;
        float rad = in_radius * u_radius_scale;
        vec4 view_coord = u_view_mat * vec4(pos, 1.0);
        vec4 view_sphere = vec4(view_coord.xyz, rad);

		float fle = u_proj_mat[1][1]; // Focal length
		float inv_ar = u_proj_mat[0][0] / u_proj_mat[1][1]; // 1.0 / aspect_ratio
		float z = -u_proj_mat[2][2] - u_proj_mat[3][2] / (view_coord.z + rad); // Compute view depth (with bias of radius) 

		vec2 axis_a;
		vec2 axis_b;
		vec2 center;
		proj_sphere(view_sphere, fle, axis_a, axis_b, center);
		vec2 scl = vec2(inv_ar, 1.0);
        axis_a *= scl;
        axis_b *= scl;
        center *= scl;
        center += u_jitter_uv.xy * 0.5;

        // Use center from projection matrix instead to get jitter
        //vec4 cs = u_proj_mat * (view_coord + vec4(0,0,0,0));
        //cs /= cs.w;
        //center = cs.xy;
        //z = cs.z;

        out_geom.view_sphere = view_sphere;
		out_geom.axis_a = axis_a;
		out_geom.axis_b = axis_b;
		out_geom.center = center;
		out_geom.z = z;
	}
}