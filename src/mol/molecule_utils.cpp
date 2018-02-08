#include "molecule_utils.h"
#include <gfx/gl_utils.h>
#include <core/common.h>
#include <imgui.h>

namespace molecule {

void transform_positions(Array<vec3> positions, const mat4& transformation) {
    for (auto& p : positions) {
        p = vec3(transformation * vec4(p, 1));
    }
}

DynamicArray<Bond> compute_bonds(const Array<vec3> atom_pos, const Array<Element> atom_elem, const Array<Residue> residues, Allocator& alloc) {
	ASSERT(atom_pos.count == atom_elem.count);

	DynamicArray<Bond> bonds(alloc);

	auto try_and_create_atom_bond =
		[&pos = atom_pos, &elem = atom_elem, &bonds = bonds](int i, int j) {
		auto d = element::covalent_radius(elem[i]) + element::covalent_radius(elem[j]);
		auto d1 = d + 0.3f;
		auto d2 = d - 0.5f;
		auto v = pos[i] - pos[j];
		auto dist2 = glm::dot(v, v);

		if (dist2 < (d1 * d1) && dist2 >(d2 * d2)) {
			bonds.push_back({ i, j });
		}
	};

	if (residues) {
		// Create connections within residues
		for (int64 i = 0; i < residues.count; i++) {
			const Residue& res = residues[i];
			for (int atom_i = res.beg_atom_idx; atom_i < res.end_atom_idx - 1; atom_i++) {
				for (int atom_j = atom_i + 1; atom_j < res.end_atom_idx; atom_j++) {
					try_and_create_atom_bond(atom_i, atom_j);
				}
			}
		}

		// @TODO:
		// Create connections between consecutive residues (Is it enough and correct to only check
		// consecutive residues?)

		for (int res_i = 0; res_i < residues.count - 1; res_i++) {
			auto res_j = res_i + 1;
			auto atom_i_beg = residues[res_i].beg_atom_idx;
			auto atom_i_end = residues[res_i].end_atom_idx;
			auto atom_j_beg = residues[res_j].beg_atom_idx;
			auto atom_j_end = residues[res_j].end_atom_idx;

			for (int atom_i = atom_i_beg; atom_i < atom_i_end; atom_i++) {
				for (int atom_j = atom_j_beg; atom_j < atom_j_end; atom_j++) {
					try_and_create_atom_bond(atom_i, atom_j);
				}
			}
		}
	}
	else {
		// Brute force N^2 check
		// @TODO: Use spatial hash
		auto atom_count = atom_pos.count;
		for (int atom_i = 0; atom_i < atom_count - 1; ++atom_i) {
			for (int atom_j = 1; atom_j < atom_count; ++atom_j) {
				try_and_create_atom_bond(atom_i, atom_j);
			}
		}
	}

	return bonds;
}

/*
DynamicArray<Backbone> compute_backbones(const Array<Residue> residues, const Array<Bond> bonds, Allocator& alloc) {

}

void compute_backbones(Array<Backbone> backbone_dst, const Array<Residue> residues, const Array<Bond> bonds) {

}
*/

DynamicArray<float>	compute_atom_radii(const Array<Element> elements, Allocator& alloc) {
	DynamicArray<float> radii(elements.count, 0, alloc);
	compute_atom_radii(radii, elements);
	return radii;
}

void compute_atom_radii(Array<float> radii_dst, const Array<Element> elements) {
	ASSERT(radii_dst.count <= elements.count);
	for (int64 i = 0; i < radii_dst.count; i++) {
		radii_dst[i] = element::vdw_radius(elements[i]);
	}
}

DynamicArray<uint32> compute_atom_colors(const MoleculeStructure& mol, ColorMapping mapping, Allocator& alloc) {
	DynamicArray<uint32> colors(mol.atom_elements.count, 0xFFFFFFFF, alloc);
	compute_atom_colors(colors, mol, mapping);
	return colors;
}

void compute_atom_colors(Array<uint32> color_dst, const MoleculeStructure& mol, ColorMapping mapping) {
	// @TODO: Implement different mappings
	(void)mapping;
	
	// CPK
	for (int64 i = 0; i < color_dst.count; i++) {
		color_dst[i] = element::color(mol.atom_elements[i]);
	}
}

namespace draw {

static constexpr int VERTEX_BUFFER_SIZE = MEGABYTES(4);
static GLuint vbo = 0;

namespace vdw {
struct Vertex {
    float32 position[3];
    float32 radius;
    uint32 color;
};

static GLuint vao = 0;

static GLuint v_shader = 0;
static GLuint g_shader = 0;
static GLuint f_shader = 0;
static GLuint program = 0;

static GLint attrib_loc_pos = -1;
static GLint attrib_loc_rad = -1;
static GLint attrib_loc_col = -1;
static GLint uniform_loc_model_mat = -1;
static GLint uniform_loc_view_mat = -1;
static GLint uniform_loc_proj_mat = -1;
static GLint uniform_loc_radius_scl = -1;

static const char* v_shader_src = R"(
#version 150 core

uniform mat4  u_model_mat;
uniform float u_radius_scl;

in vec3	 v_position;
in float v_radius;
in vec4  v_color;

out Vertex {
    vec4 color;
    float radius;
} out_vert;

void main() {
    gl_Position = u_model_mat * vec4(v_position, 1.0);
    out_vert.color = v_color;
    out_vert.radius = v_radius * u_radius_scl;
}
)";

static const char* g_shader_src = R"(
#version 150 core

uniform mat4 u_view_mat;
uniform mat4 u_proj_mat;

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

in Vertex {
    vec4 color;
    float radius;
} in_vert[];

out Fragment {
    flat vec4 color;
    flat vec4 view_sphere;
    smooth vec4 view_pos;
} out_frag;

void main()
{
	vec4 color = in_vert[0].color;
    float radius = in_vert[0].radius;

	if (radius == 0 || color.a == 0) return;

    vec4 pos = gl_in[0].gl_Position;
    vec4 view_pos = u_view_mat * pos;
    float len = length(view_pos.xyz);
    vec3 view_dir = view_pos.xyz / len;

    out_frag.color = color;
    out_frag.view_sphere = vec4(view_pos.xyz, radius);

    view_pos.xyz -= view_dir * radius;
    vec2 uv;

    const float sqrt_two = sqrt(2.0);
    float scl = sqrt_two * radius;

    uv = vec2(-1, -1) * scl;
    out_frag.view_pos = view_pos + vec4(uv, 0, 1);
    gl_Position = u_proj_mat * out_frag.view_pos;
    EmitVertex();

    uv = vec2(-1, 1) * scl;
    out_frag.view_pos = view_pos + vec4(uv, 0, 1);
    gl_Position = u_proj_mat * out_frag.view_pos;
    EmitVertex();

    uv = vec2(1, -1) * scl;
    out_frag.view_pos = view_pos + vec4(uv, 0, 1);
    gl_Position = u_proj_mat * out_frag.view_pos;
    EmitVertex();

    uv = vec2(1, 1) * scl;
    out_frag.view_pos = view_pos + vec4(uv, 0, 1);
    gl_Position = u_proj_mat * out_frag.view_pos;
    EmitVertex();

    EndPrimitive();
}
)";

static const char* f_shader_src = R"(
#version 150 core
#extension GL_ARB_conservative_depth: enable

uniform mat4 u_proj_mat;

in Fragment {
    flat vec4 color;
    flat vec4 view_sphere;
    smooth vec4 view_pos;
} in_frag;

#ifdef GL_EXT_conservative_depth
layout (depth_greater) out float gl_FragDepth;
#endif
out vec4 out_color;

void main() {
    vec3 center = in_frag.view_sphere.xyz;
    float radius = in_frag.view_sphere.w;
    vec3 view_dir = normalize(in_frag.view_pos.xyz);

    vec3 m = -center;
    vec3 d = view_dir;
    float r = radius;
    float b = dot(m, d);
    float c = dot(m, m) - r*r;
    float discr = b*b - c;
    if (discr < 0.0) discard;
    float t = -b -sqrt(discr);

    vec3 view_hit = d * t;
    vec3 view_normal = (view_hit - center) / radius;
    vec4 color = in_frag.color;

    vec4 coord = vec4(0, 0, view_hit.z, 1);
    coord = u_proj_mat * coord;
    coord = coord / coord.w;

    gl_FragDepth = coord.z * 0.5 + 0.5;
    out_color = vec4(color.rgb, color.a);
}
)";

static void initialize() {
    constexpr int BUFFER_SIZE = 1024;
    char buffer[BUFFER_SIZE];

    v_shader = glCreateShader(GL_VERTEX_SHADER);
    g_shader = glCreateShader(GL_GEOMETRY_SHADER);
    f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(g_shader, 1, &g_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);

    glCompileShader(v_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, v_shader)) {
        printf("Error while compiling vdw vertex shader:\n%s\n", buffer);
    }
    glCompileShader(g_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, g_shader)) {
        printf("Error while compiling vdw geometry shader:\n%s\n", buffer);
    }
    glCompileShader(f_shader);
    if (gl::get_shader_compile_error(buffer, BUFFER_SIZE, f_shader)) {
        printf("Error while compiling vdw fragment shader:\n%s\n", buffer);
    }

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, g_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    if (gl::get_program_link_error(buffer, BUFFER_SIZE, program)) {
        printf("Error while linking vdw program:\n%s\n", buffer);
    }

    glDetachShader(program, v_shader);
    glDetachShader(program, g_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(g_shader);
    glDeleteShader(f_shader);

    attrib_loc_pos = glGetAttribLocation(program, "v_position");
    attrib_loc_rad = glGetAttribLocation(program, "v_radius");
    attrib_loc_col = glGetAttribLocation(program, "v_color");
    uniform_loc_model_mat = glGetUniformLocation(program, "u_model_mat");
    uniform_loc_view_mat = glGetUniformLocation(program, "u_view_mat");
    uniform_loc_proj_mat = glGetUniformLocation(program, "u_proj_mat");
    uniform_loc_radius_scl = glGetUniformLocation(program, "u_radius_scl");

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnableVertexAttribArray(attrib_loc_pos);
    glVertexAttribPointer(attrib_loc_pos, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)0);

    glEnableVertexAttribArray(attrib_loc_rad);
    glVertexAttribPointer(attrib_loc_rad, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)12);

    glEnableVertexAttribArray(attrib_loc_col);
    glVertexAttribPointer(attrib_loc_col, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (const GLvoid*)16);

    glBindVertexArray(0);
}

static void shutdown() {
    if (vao) glDeleteVertexArrays(1, &vao);
    if (program) glDeleteProgram(program);
}

}  // namespace vdw

void initialize() {
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, VERTEX_BUFFER_SIZE, nullptr, GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    vdw::initialize();
}

void shutdown() {
    if (vbo) glDeleteBuffers(1, &vbo);

    vdw::shutdown();
}

void draw_vdw(const Array<vec3> atom_positions, const Array<float> atom_radii, const Array<uint32> atom_colors, const mat4& model_mat,
              const mat4& view_mat, const mat4& proj_mat, float radii_scale) {
    int64_t count = atom_positions.count;
    ASSERT(count == atom_radii.count && count == atom_colors.count);
    ASSERT(count * sizeof(vdw::Vertex) < VERTEX_BUFFER_SIZE);

	static bool init = false;
	if (!init) {
		init = true;
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		vdw::Vertex* data = (vdw::Vertex*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
		for (int64_t i = 0; i < count; i++) {
			data[i].position[0] = atom_positions[i][0];
			data[i].position[1] = atom_positions[i][1];
			data[i].position[2] = atom_positions[i][2];
			data[i].radius = atom_radii[i];
			data[i].color = atom_colors[i];
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);
	}

    glBindVertexArray(vdw::vao);
    glUseProgram(vdw::program);
    glUniformMatrix4fv(vdw::uniform_loc_model_mat, 1, GL_FALSE, &model_mat[0][0]);
	glUniformMatrix4fv(vdw::uniform_loc_view_mat, 1, GL_FALSE, &view_mat[0][0]);
	glUniformMatrix4fv(vdw::uniform_loc_proj_mat, 1, GL_FALSE, &proj_mat[0][0]);
    glUniform1f(vdw::uniform_loc_radius_scl, radii_scale);
    glDrawArrays(GL_POINTS, 0, (GLsizei)count);
    glUseProgram(0);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

}  // namespace draw

}  // namespace molecule