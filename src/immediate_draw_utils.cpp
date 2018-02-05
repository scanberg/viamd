#include <immediate_draw_utils.h>
#include <core/gl.h>
#include <core/array.h>

namespace immediate {

struct Vertex {
    float position[3];
    unsigned short tex_coord[2];
    unsigned char color[4];
};

using Index = uint16;

struct DrawCommand {
    Index offset;
    Index count;
    GLenum primitive_type;
};

static DynamicArray<DrawCommand> commands;
static DynamicArray<Vertex> vertices;
static DynamicArray<Index> indices;
// static DynamicArray<Index> point_indices;
// static DynamicArray<Index> line_indices;
// static DynamicArray<Index> triangle_indices;

static GLuint vbo = 0;
static GLuint ibo = 0;
static GLuint vao = 0;

static GLuint v_shader = 0;
static GLuint f_shader = 0;
static GLuint program = 0;

static GLint attrib_loc_pos = -1;
static GLint attrib_loc_tc = -1;
static GLint attrib_loc_col = -1;
static GLint uniform_loc_mvp = -1;

static const char* v_shader_src =
    "#version 150\n"
    "uniform mat4 u_mvp;\n"
    "in vec3 in_pos;\n"
    "in vec2 in_tc;\n"
    "in vec4 in_col;\n"
    "out vec2 tc;\n"
    "out vec4 col;\n"
    "void main(){\n"
    "gl_Position = u_mvp * vec4(in_pos, 1);\n"
    "tc = in_tc;\n"
    "col = in_col;\n"
    "}";

static const char* f_shader_src =
    "#version 150\n"
    "in vec2 tc;\n"
    "in vec4 col;\n"
    "out vec4 out_frag;\n"
    "void main(){\n"
    "out_frag = col;\n"
    "}";

static void log_shader_compile_error(GLuint shader) {
    constexpr int MAX_BUFFER_LENGTH = 1024;
    char buffer[MAX_BUFFER_LENGTH];

    GLint success = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (success == GL_FALSE) {
        int length;
        glGetShaderInfoLog(shader, MAX_BUFFER_LENGTH, &length, buffer);
        printf("%s\n", buffer);
    }
}

static void log_program_link_error(GLuint program) {
    constexpr int MAX_BUFFER_LENGTH = 1024;
    char buffer[MAX_BUFFER_LENGTH];

    GLint success = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (success == GL_FALSE) {
        int length;
        glGetProgramInfoLog(program, MAX_BUFFER_LENGTH, &length, buffer);
        printf("%s\n", buffer);
    }
}

static inline void append_draw_command(GLenum primitive_type) {
    Index count = primitive_type == GL_POINTS ? 1 : (primitive_type == GL_LINES ? 2 : 3);
    if (commands.count > 0) {
        if (commands.back().primitive_type == primitive_type) {
            commands.back().count += count;
        }
    } else {
        DrawCommand cmd{(Index)indices.count, count, primitive_type};
        commands.push_back(cmd);
    }
}

void initialize() {
    v_shader = glCreateShader(GL_VERTEX_SHADER);
    f_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(v_shader, 1, &v_shader_src, 0);
    glShaderSource(f_shader, 1, &f_shader_src, 0);
    glCompileShader(v_shader);
    log_shader_compile_error(v_shader);
    glCompileShader(f_shader);
    log_shader_compile_error(f_shader);

    program = glCreateProgram();
    glAttachShader(program, v_shader);
    glAttachShader(program, f_shader);
    glLinkProgram(program);
    log_program_link_error(program);

    glDetachShader(program, v_shader);
    glDetachShader(program, f_shader);

    glDeleteShader(v_shader);
    glDeleteShader(f_shader);

    attrib_loc_pos = glGetAttribLocation(program, "in_pos");
    attrib_loc_tc = glGetAttribLocation(program, "in_tc");
    attrib_loc_col = glGetAttribLocation(program, "in_col");
    uniform_loc_mvp = glGetUniformLocation(program, "u_mvp");

    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ibo);

    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (const GLvoid*)0);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_UNSIGNED_SHORT, GL_TRUE, sizeof(Vertex), (const GLvoid*)12);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Vertex), (const GLvoid*)16);

    glBindVertexArray(0);

    vertices.reserve(100000);
}

void shutdown() {
    if (vbo) glDeleteBuffers(1, &vbo);
    if (vao) glDeleteVertexArrays(1, &vao);
    if (program) glDeleteProgram(program);
}

void flush(float mvp_matrix[16]) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.count * sizeof(Vertex), vertices.data, GL_STREAM_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.count * sizeof(Index), indices.data, GL_STREAM_DRAW);

    glUseProgram(program);
    glUniformMatrix4fv(uniform_loc_mvp, 1, GL_FALSE, mvp_matrix);
    glBindVertexArray(vao);

    for (const auto& cmd : commands) {
        glDrawElements(cmd.primitive_type, cmd.count, sizeof(Index) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT, reinterpret_cast<const void*>(cmd.offset));
    }

    glBindVertexArray(0);
    glUseProgram(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    vertices.clear();
    indices.clear();
    commands.clear();
}

void draw_point(float pos[3], unsigned char color[4]) {
    Index idx = (Index)vertices.count;
    Vertex vert = {{pos[0], pos[1], pos[2]}, {0, 0}, {color[0], color[1], color[2], color[3]}};
    vertices.push_back(vert);
    indices.push_back(idx);
    append_draw_command(GL_POINTS);
}

void draw_triangle(float v0[3], float v1[3], float v2[3], unsigned char color[4]) {
    Index idx = (Index)vertices.count;
    Vertex vert0 = {{v0[0], v0[1], v0[2]}, {0, 0}, {color[0], color[1], color[2], color[3]}};
    Vertex vert1 = {{v1[0], v1[1], v1[2]}, {0, 0}, {color[0], color[1], color[2], color[3]}};
    Vertex vert2 = {{v2[0], v2[1], v2[2]}, {0, 0}, {color[0], color[1], color[2], color[3]}};

    vertices.push_back(vert0);
    vertices.push_back(vert1);
    vertices.push_back(vert2);

    indices.push_back(idx);
    indices.push_back(idx + 1);
    indices.push_back(idx + 2);

    append_draw_command(GL_TRIANGLES);
}

}
