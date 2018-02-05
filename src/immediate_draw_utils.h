#pragma once

namespace immediate {

void initialize();
void shutdown();
void flush(float mvp_matrix[16]);

void draw_sphere(float pos[3], float radius, unsigned char color[4]);
void draw_point(float pos[3], unsigned char color[4]);
void draw_line(float from[3], float to[3], unsigned char color[4]);
void draw_triangle(float v0[3], float v1[3], float v2[3], unsigned char color[4]);
void draw_screen_space_quad(int x0, int y0, int x1, int y1);

}