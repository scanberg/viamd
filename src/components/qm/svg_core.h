#pragma once

// A scene of rectangles, cubic curves, paths and text, written as an SVG file: what the flow diagram is exported as
// (flow_draw.h).

#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_vec.h>

#include <stdbool.h>
#include <stdint.h>

struct md_allocator_i;
struct md_writer_i;

typedef uint32_t vg_color_t;

#define VG_COLOR_RGBA(r, g, b, a) ((((uint32_t)(r) & 0xFFu) << 24u) | (((uint32_t)(g) & 0xFFu) << 16u) | (((uint32_t)(b) & 0xFFu) << 8u) | ((uint32_t)(a) & 0xFFu))
#define VG_COLOR_RGB(r, g, b)      VG_COLOR_RGBA((r), (g), (b), 255)

typedef enum vg_paint_type_t {
    VG_PAINT_NONE = 0,
    VG_PAINT_COLOR,
    VG_PAINT_LINEAR_GRADIENT,
} vg_paint_type_t;

typedef struct vg_linear_gradient_t {
    md_float2_t p0;
    md_float2_t p1;
    vg_color_t col0;
    vg_color_t col1;
} vg_linear_gradient_t;

typedef struct vg_paint_t {
    vg_paint_type_t type;
    union {
        vg_color_t color;
        vg_linear_gradient_t linear_gradient;
    } data;
} vg_paint_t;

typedef enum vg_line_cap_t {
    VG_LINE_CAP_BUTT = 0,
    VG_LINE_CAP_ROUND,
    VG_LINE_CAP_SQUARE,
} vg_line_cap_t;

typedef enum vg_line_join_t {
    VG_LINE_JOIN_MITER = 0,
    VG_LINE_JOIN_ROUND,
    VG_LINE_JOIN_BEVEL,
} vg_line_join_t;

typedef enum vg_fill_rule_t {
    VG_FILL_RULE_NON_ZERO = 0,
    VG_FILL_RULE_EVEN_ODD,
} vg_fill_rule_t;

// Which point of a line of text its x is: its start, its middle or its end (SVG text-anchor). The
// viewer's own font metrics do the aligning, so text stays aligned in a font the writer cannot measure.
typedef enum vg_text_anchor_t {
    VG_TEXT_ANCHOR_START = 0,
    VG_TEXT_ANCHOR_MIDDLE,
    VG_TEXT_ANCHOR_END,
} vg_text_anchor_t;

typedef struct vg_style_t {
    vg_paint_t fill;
    vg_paint_t stroke;
    float stroke_width;
    float miter_limit;
    vg_line_cap_t line_cap;
    vg_line_join_t line_join;
    vg_fill_rule_t fill_rule;
} vg_style_t;

typedef enum vg_path_verb_t {
    VG_PATH_MOVE_TO = 0,
    VG_PATH_LINE_TO,
    VG_PATH_CUBIC_TO,
    VG_PATH_CLOSE,
} vg_path_verb_t;

typedef enum vg_cmd_type_t {
    VG_CMD_RECT = 0,
    VG_CMD_BEZIER_CUBIC,
    VG_CMD_PATH,
    VG_CMD_TEXT,
} vg_cmd_type_t;

typedef struct vg_cmd_t {
    vg_cmd_type_t type;
    vg_style_t style;
    union {
        struct {
            md_float2_t min;
            md_float2_t max;
            float rounding;
        } rect;
        struct {
            md_float2_t p1;
            md_float2_t p2;
            md_float2_t p3;
            md_float2_t p4;
        } bezier_cubic;
        struct {
            uint32_t verb_offset;
            uint32_t verb_count;
            uint32_t point_offset;
            uint32_t point_count;
        } path;
        struct {
            md_float2_t pos;
            float font_size;
            uint32_t text_offset;
            uint32_t text_len;
            uint32_t anchor;    // vg_text_anchor_t
            float rotation;     // degrees, clockwise, about the middle of the line at pos.x
        } text;
    } data;
} vg_cmd_t;

typedef struct vg_scene_t {
    struct md_allocator_i* alloc;
    md_float2_t size;
    md_array(vg_cmd_t) commands;
    md_array(md_float2_t) points;
    md_array(vg_path_verb_t) verbs;
    md_array(char) strings;
    md_array(md_float2_t) path_points;
    md_array(vg_path_verb_t) path_verbs;
} vg_scene_t;


static inline vg_paint_t vg_paint_none(void) {
    vg_paint_t paint = {};
    paint.type = VG_PAINT_NONE;
    return paint;
}

static inline vg_paint_t vg_paint_color(vg_color_t color) {
    vg_paint_t paint = {};
    paint.type = VG_PAINT_COLOR;
    paint.data.color = color;
    return paint;
}

static inline vg_paint_t vg_paint_linear_gradient(md_float2_t p0, md_float2_t p1, vg_color_t col0, vg_color_t col1) {
    vg_paint_t paint = {};
    paint.type = VG_PAINT_LINEAR_GRADIENT;
    paint.data.linear_gradient.p0 = p0;
    paint.data.linear_gradient.p1 = p1;
    paint.data.linear_gradient.col0 = col0;
    paint.data.linear_gradient.col1 = col1;
    return paint;
}

static inline vg_style_t vg_style_default(void) {
    vg_style_t style = {};
    style.fill.type = VG_PAINT_NONE;
    style.stroke.type = VG_PAINT_NONE;
    style.stroke_width = 1.0f;
    style.miter_limit = 4.0f;
    style.line_cap = VG_LINE_CAP_BUTT;
    style.line_join = VG_LINE_JOIN_MITER;
    style.fill_rule = VG_FILL_RULE_NON_ZERO;
    return style;
}

static inline vg_style_t vg_style_stroke(vg_color_t color, float width) {
    vg_style_t style = vg_style_default();
    style.stroke = vg_paint_color(color);
    style.stroke_width = width;
    return style;
}

vg_color_t vg_color_from_vec4(md_float4_t color);

void vg_scene_init(vg_scene_t* scene, md_float2_t size, struct md_allocator_i* alloc);

bool vg_scene_write_svg_file(const vg_scene_t* scene, md_str_t path);

void vg_add_rect(vg_scene_t* scene, md_float2_t min, md_float2_t max, vg_color_t color, float rounding, float thickness);
void vg_add_rect_filled(vg_scene_t* scene, md_float2_t min, md_float2_t max, vg_color_t color, float rounding);
void vg_add_rect_styled(vg_scene_t* scene, md_float2_t min, md_float2_t max, float rounding, const vg_style_t* style);

void vg_add_bezier_cubic_styled(vg_scene_t* scene, md_float2_t p0, md_float2_t p1, md_float2_t p2, md_float2_t p3, const vg_style_t* style);

void vg_add_text(vg_scene_t* scene, md_float2_t pos, float font_size, vg_color_t color, md_str_t text);
void vg_add_text_styled(vg_scene_t* scene, md_float2_t pos, float font_size, const vg_style_t* style, md_str_t text);
// pos.y is the top of the text, as for vg_add_text; pos.x is the point the anchor names
void vg_add_text_anchored(vg_scene_t* scene, md_float2_t pos, float font_size, vg_color_t color, vg_text_anchor_t anchor, md_str_t text);
// As vg_add_text_anchored, turned by 'rotation' degrees (clockwise, -90 reads upwards) about the middle of the line
// at pos.x: (pos.x, pos.y + font_size / 2)
void vg_add_text_rotated(vg_scene_t* scene, md_float2_t pos, float font_size, vg_color_t color, vg_text_anchor_t anchor, float rotation, md_str_t text);

void vg_path_clear(vg_scene_t* scene);
void vg_path_move_to(vg_scene_t* scene, md_float2_t pos);
void vg_path_line_to(vg_scene_t* scene, md_float2_t pos);
void vg_path_bezier_cubic_curve_to(vg_scene_t* scene, md_float2_t p1, md_float2_t p2, md_float2_t p3);
void vg_path_close(vg_scene_t* scene);
void vg_path_fill_styled(vg_scene_t* scene, const vg_style_t* style);

