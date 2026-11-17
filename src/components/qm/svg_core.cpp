#include "svg_core.h"

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <core/md_writer.h>

#include <math.h>
#include <string.h>

static vg_style_t safe_style(const vg_style_t* style) {
    return style ? *style : vg_style_default();
}

static bool paint_is_none(vg_paint_t paint) {
    return paint.type == VG_PAINT_NONE;
}

static bool paint_is_gradient(vg_paint_t paint) {
    return paint.type == VG_PAINT_LINEAR_GRADIENT;
}

static uint32_t push_points(vg_scene_t* scene, const md_float2_t* points, size_t count) {
    MD_ASSERT(scene);
    MD_ASSERT(scene->alloc);
    const uint32_t offset = (uint32_t)md_array_size(scene->points);
    md_array_push_array(scene->points, points, count, scene->alloc);
    return offset;
}

static uint32_t push_verbs(vg_scene_t* scene, const vg_path_verb_t* verbs, size_t count) {
    MD_ASSERT(scene);
    MD_ASSERT(scene->alloc);
    const uint32_t offset = (uint32_t)md_array_size(scene->verbs);
    md_array_push_array(scene->verbs, verbs, count, scene->alloc);
    return offset;
}

static uint32_t push_string(vg_scene_t* scene, md_str_t str) {
    MD_ASSERT(scene);
    MD_ASSERT(scene->alloc);
    const uint32_t offset = (uint32_t)md_array_size(scene->strings);
    md_array_push_array(scene->strings, str.ptr, str.len, scene->alloc);
    md_array_push(scene->strings, '\0', scene->alloc);
    return offset;
}

static void push_cmd(vg_scene_t* scene, const vg_cmd_t* cmd) {
    MD_ASSERT(scene);
    MD_ASSERT(scene->alloc);
    md_array_push(scene->commands, *cmd, scene->alloc);
}

static void emit_xml_escaped(md_writer_i* w, md_str_t text) {
    for (size_t i = 0; i < text.len; ++i) {
        switch (text.ptr[i]) {
        case '&':
            md_write_cstr(w, "&amp;");
            break;
        case '<':
            md_write_cstr(w, "&lt;");
            break;
        case '>':
            md_write_cstr(w, "&gt;");
            break;
        case '\"':
            md_write_cstr(w, "&quot;");
            break;
        case '\'':
            md_write_cstr(w, "&apos;");
            break;
        default:
            md_write_char(w, text.ptr[i]);
            break;
        }
    }
}

static void color_unpack(vg_color_t color, uint8_t rgba[4]) {
    rgba[0] = (uint8_t)((color >> 24u) & 0xFFu);
    rgba[1] = (uint8_t)((color >> 16u) & 0xFFu);
    rgba[2] = (uint8_t)((color >> 8u) & 0xFFu);
    rgba[3] = (uint8_t)(color & 0xFFu);
}

static void emit_color_attr(md_writer_i* w, const char* attr_name, vg_color_t color) {
    uint8_t rgba[4];
    color_unpack(color, rgba);
    md_writef(w, " %s=\"rgb(%u,%u,%u)\"", attr_name, rgba[0], rgba[1], rgba[2]);
    if (rgba[3] != 255) {
        md_writef(w, " %s-opacity=\"%.9g\"", attr_name, rgba[3] / 255.0);
    }
}

static void emit_paint_attr(md_writer_i* w, const char* attr_name, vg_paint_t paint, size_t cmd_idx) {
    switch (paint.type) {
    default:
    case VG_PAINT_NONE:
        md_writef(w, " %s=\"none\"", attr_name);
        break;
    case VG_PAINT_COLOR:
        emit_color_attr(w, attr_name, paint.data.color);
        break;
    case VG_PAINT_LINEAR_GRADIENT:
        md_writef(w, " %s=\"url(#mdvg-grad-%zu-%s)\"", attr_name, cmd_idx, attr_name);
        break;
    }
}

static const char* line_cap_str(vg_line_cap_t cap) {
    switch (cap) {
    default:
    case VG_LINE_CAP_BUTT:   return "butt";
    case VG_LINE_CAP_ROUND:  return "round";
    case VG_LINE_CAP_SQUARE: return "square";
    }
}

static const char* line_join_str(vg_line_join_t join) {
    switch (join) {
    default:
    case VG_LINE_JOIN_MITER: return "miter";
    case VG_LINE_JOIN_ROUND: return "round";
    case VG_LINE_JOIN_BEVEL: return "bevel";
    }
}

static const char* fill_rule_str(vg_fill_rule_t fill_rule) {
    switch (fill_rule) {
    default:
    case VG_FILL_RULE_NON_ZERO: return "nonzero";
    case VG_FILL_RULE_EVEN_ODD: return "evenodd";
    }
}

static void emit_style(md_writer_i* w, const vg_style_t* style, size_t cmd_idx) {
    emit_paint_attr(w, "fill", style->fill, cmd_idx);
    emit_paint_attr(w, "stroke", style->stroke, cmd_idx);

    if (!paint_is_none(style->stroke)) {
        md_writef(w, " stroke-width=\"%.9g\"", style->stroke_width);
        md_writef(w, " stroke-linecap=\"%s\"", line_cap_str(style->line_cap));
        md_writef(w, " stroke-linejoin=\"%s\"", line_join_str(style->line_join));
        md_writef(w, " stroke-miterlimit=\"%.9g\"", style->miter_limit);
    }
    if (!paint_is_none(style->fill)) {
        md_writef(w, " fill-rule=\"%s\"", fill_rule_str(style->fill_rule));
    }
}

static void emit_gradient_def(md_writer_i* w, size_t cmd_idx, const char* slot, const vg_linear_gradient_t* grad) {
    uint8_t rgba0[4];
    uint8_t rgba1[4];
    color_unpack(grad->col0, rgba0);
    color_unpack(grad->col1, rgba1);

    md_writef(w, "    <linearGradient id=\"mdvg-grad-%zu-%s\" gradientUnits=\"userSpaceOnUse\" x1=\"%.9g\" y1=\"%.9g\" x2=\"%.9g\" y2=\"%.9g\">\n",
        cmd_idx, slot,
        grad->p0.x, grad->p0.y,
        grad->p1.x, grad->p1.y);
    md_writef(w, "      <stop offset=\"0%%\" stop-color=\"rgb(%u,%u,%u)\"", rgba0[0], rgba0[1], rgba0[2]);
    if (rgba0[3] != 255) {
        md_writef(w, " stop-opacity=\"%.9g\"", rgba0[3] / 255.0);
    }
    md_write_cstr(w, "/>\n");
    md_writef(w, "      <stop offset=\"100%%\" stop-color=\"rgb(%u,%u,%u)\"", rgba1[0], rgba1[1], rgba1[2]);
    if (rgba1[3] != 255) {
        md_writef(w, " stop-opacity=\"%.9g\"", rgba1[3] / 255.0);
    }
    md_write_cstr(w, "/>\n");
    md_write_cstr(w, "    </linearGradient>\n");
}

static void emit_defs(md_writer_i* w, const vg_scene_t* scene) {
    bool has_gradients = false;
    for (size_t i = 0; i < md_array_size(scene->commands); ++i) {
        const vg_style_t* style = &scene->commands[i].style;
        has_gradients |= paint_is_gradient(style->fill) || paint_is_gradient(style->stroke);
    }
    if (!has_gradients) return;

    md_write_cstr(w, "  <defs>\n");
    for (size_t i = 0; i < md_array_size(scene->commands); ++i) {
        const vg_style_t* style = &scene->commands[i].style;
        if (paint_is_gradient(style->fill)) {
            emit_gradient_def(w, i, "fill", &style->fill.data.linear_gradient);
        }
        if (paint_is_gradient(style->stroke)) {
            emit_gradient_def(w, i, "stroke", &style->stroke.data.linear_gradient);
        }
    }
    md_write_cstr(w, "  </defs>\n");
}

static void emit_path_data(md_writer_i* w, const vg_scene_t* scene, const vg_cmd_t* cmd) {
    const vg_path_verb_t* verbs = scene->verbs + cmd->data.path.verb_offset;
    const md_float2_t* points = scene->points + cmd->data.path.point_offset;
    size_t point_idx = 0;

    for (size_t i = 0; i < cmd->data.path.verb_count; ++i) {
        switch (verbs[i]) {
        case VG_PATH_MOVE_TO:
            md_writef(w, "M %.9g %.9g ", points[point_idx + 0].x, points[point_idx + 0].y);
            point_idx += 1;
            break;
        case VG_PATH_LINE_TO:
            md_writef(w, "L %.9g %.9g ", points[point_idx + 0].x, points[point_idx + 0].y);
            point_idx += 1;
            break;
        case VG_PATH_CUBIC_TO:
            md_writef(w, "C %.9g %.9g %.9g %.9g %.9g %.9g ",
                points[point_idx + 0].x, points[point_idx + 0].y,
                points[point_idx + 1].x, points[point_idx + 1].y,
                points[point_idx + 2].x, points[point_idx + 2].y);
            point_idx += 3;
            break;
        case VG_PATH_CLOSE:
            md_write_cstr(w, "Z ");
            break;
        default:
            break;
        }
    }
}

static void emit_svg_scene(md_writer_i* w, const vg_scene_t* scene) {
    md_write_cstr(w, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    md_writef(w,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"%.9g\" height=\"%.9g\" viewBox=\"0 0 %.9g %.9g\">\n",
        scene->size.x, scene->size.y, scene->size.x, scene->size.y);

    emit_defs(w, scene);

    for (size_t i = 0; i < md_array_size(scene->commands); ++i) {
        const vg_cmd_t* cmd = &scene->commands[i];
        switch (cmd->type) {
        case VG_CMD_RECT: {
            const md_float2_t size = { cmd->data.rect.max.x - cmd->data.rect.min.x, cmd->data.rect.max.y - cmd->data.rect.min.y };
            md_writef(w, "  <rect x=\"%.9g\" y=\"%.9g\" width=\"%.9g\" height=\"%.9g\"",
                cmd->data.rect.min.x, cmd->data.rect.min.y, size.x, size.y);
            if (cmd->data.rect.rounding > 0.0f) {
                md_writef(w, " rx=\"%.9g\" ry=\"%.9g\"", cmd->data.rect.rounding, cmd->data.rect.rounding);
            }
            emit_style(w, &cmd->style, i);
            md_write_cstr(w, " />\n");
            break;
        }
        case VG_CMD_BEZIER_CUBIC:
            md_writef(w, "  <path d=\"M %.9g %.9g C %.9g %.9g %.9g %.9g %.9g %.9g\"",
                cmd->data.bezier_cubic.p1.x, cmd->data.bezier_cubic.p1.y,
                cmd->data.bezier_cubic.p2.x, cmd->data.bezier_cubic.p2.y,
                cmd->data.bezier_cubic.p3.x, cmd->data.bezier_cubic.p3.y,
                cmd->data.bezier_cubic.p4.x, cmd->data.bezier_cubic.p4.y);
            emit_style(w, &cmd->style, i);
            md_write_cstr(w, " />\n");
            break;
        case VG_CMD_PATH:
            md_write_cstr(w, "  <path d=\"");
            emit_path_data(w, scene, cmd);
            md_write_char(w, '\"');
            emit_style(w, &cmd->style, i);
            md_write_cstr(w, " />\n");
            break;
        case VG_CMD_TEXT: {
            const md_str_t text = { scene->strings + cmd->data.text.text_offset, cmd->data.text.text_len };
            md_writef(w, "  <text x=\"%.9g\" y=\"%.9g\" font-size=\"%.9g\" font-family=\"sans-serif\" dominant-baseline=\"hanging\" xml:space=\"preserve\"",
                cmd->data.text.pos.x,
                cmd->data.text.pos.y,
                cmd->data.text.font_size);
            if (cmd->data.text.anchor == VG_TEXT_ANCHOR_MIDDLE)   md_write_cstr(w, " text-anchor=\"middle\"");
            else if (cmd->data.text.anchor == VG_TEXT_ANCHOR_END) md_write_cstr(w, " text-anchor=\"end\"");
            if (cmd->data.text.rotation != 0.0f) {
                md_writef(w, " transform=\"rotate(%.9g %.9g %.9g)\"", cmd->data.text.rotation, cmd->data.text.pos.x,
                    cmd->data.text.pos.y + cmd->data.text.font_size * 0.5f);
            }
            emit_style(w, &cmd->style, i);
            md_write_cstr(w, ">");
            emit_xml_escaped(w, text);
            md_write_cstr(w, "</text>\n");
            break;
        }
        default:
            break;
        }
    }

    md_write_cstr(w, "</svg>\n");
}

vg_color_t vg_color_from_vec4(md_float4_t color) {
    uint32_t r = (uint32_t)MD_CLAMP(color.x * 255.0f + 0.5f, 0.0f, 255.0f);
    uint32_t g = (uint32_t)MD_CLAMP(color.y * 255.0f + 0.5f, 0.0f, 255.0f);
    uint32_t b = (uint32_t)MD_CLAMP(color.z * 255.0f + 0.5f, 0.0f, 255.0f);
    uint32_t a = (uint32_t)MD_CLAMP(color.w * 255.0f + 0.5f, 0.0f, 255.0f);
    return VG_COLOR_RGBA(r, g, b, a);
}

void vg_scene_init(vg_scene_t* scene, md_float2_t size, md_allocator_i* alloc) {
    MD_ASSERT(scene);
    MD_ASSERT(alloc);
    *scene = vg_scene_t{};
    scene->size = size;
    scene->alloc = alloc;
}

bool vg_scene_write_svg_file(const vg_scene_t* scene, md_str_t path) {
    MD_ASSERT(scene);
    md_file_t file = {};
    if (!md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        return false;
    }
    // Streamed through a staging buffer, rather than built whole in memory first
    md_file_writer_t fw = md_file_writer_create(file, NULL, 0);
    emit_svg_scene(&fw.w, scene);
    const bool result = md_file_writer_close(&fw);
    md_file_close(&file);
    return result;
}

// A command of a type and a style, its data to be filled in
static vg_cmd_t make_cmd(vg_cmd_type_t type, const vg_style_t* style) {
    vg_cmd_t cmd = {};
    cmd.type  = type;
    cmd.style = safe_style(style);
    return cmd;
}

static vg_style_t style_fill(vg_color_t color) {
    vg_style_t style = vg_style_default();
    style.fill = vg_paint_color(color);
    return style;
}

void vg_add_rect(vg_scene_t* scene, md_float2_t min, md_float2_t max, vg_color_t color, float rounding, float thickness) {
    const vg_style_t style = vg_style_stroke(color, thickness);
    vg_add_rect_styled(scene, min, max, rounding, &style);
}

void vg_add_rect_filled(vg_scene_t* scene, md_float2_t min, md_float2_t max, vg_color_t color, float rounding) {
    const vg_style_t style = style_fill(color);
    vg_add_rect_styled(scene, min, max, rounding, &style);
}

void vg_add_rect_styled(vg_scene_t* scene, md_float2_t min, md_float2_t max, float rounding, const vg_style_t* style) {
    vg_cmd_t cmd = make_cmd(VG_CMD_RECT, style);
    cmd.data.rect.min = min;
    cmd.data.rect.max = max;
    cmd.data.rect.rounding = rounding;
    push_cmd(scene, &cmd);
}

void vg_add_bezier_cubic_styled(vg_scene_t* scene, md_float2_t p0, md_float2_t p1, md_float2_t p2, md_float2_t p3, const vg_style_t* style) {
    vg_cmd_t cmd = make_cmd(VG_CMD_BEZIER_CUBIC, style);
    cmd.data.bezier_cubic.p1 = p0;
    cmd.data.bezier_cubic.p2 = p1;
    cmd.data.bezier_cubic.p3 = p2;
    cmd.data.bezier_cubic.p4 = p3;
    cmd.style.fill = vg_paint_none();
    push_cmd(scene, &cmd);
}

void vg_add_text(vg_scene_t* scene, md_float2_t pos, float font_size, vg_color_t color, md_str_t text) {
    const vg_style_t style = style_fill(color);
    vg_add_text_styled(scene, pos, font_size, &style, text);
}

void vg_add_text_styled(vg_scene_t* scene, md_float2_t pos, float font_size, const vg_style_t* style, md_str_t text) {
    if (!text.ptr || text.len == 0) return;
    vg_cmd_t cmd = make_cmd(VG_CMD_TEXT, style);
    cmd.data.text.pos = pos;
    cmd.data.text.font_size = font_size;
    cmd.data.text.text_offset = push_string(scene, text);
    cmd.data.text.text_len = (uint32_t)text.len;
    push_cmd(scene, &cmd);
}

void vg_add_text_anchored(vg_scene_t* scene, md_float2_t pos, float font_size, vg_color_t color, vg_text_anchor_t anchor, md_str_t text) {
    const size_t before = md_array_size(scene->commands);
    vg_add_text(scene, pos, font_size, color, text);
    if (md_array_size(scene->commands) > before) scene->commands[before].data.text.anchor = (uint32_t)anchor;
}

void vg_add_text_rotated(vg_scene_t* scene, md_float2_t pos, float font_size, vg_color_t color, vg_text_anchor_t anchor, float rotation, md_str_t text) {
    const size_t before = md_array_size(scene->commands);
    vg_add_text_anchored(scene, pos, font_size, color, anchor, text);
    if (md_array_size(scene->commands) > before) scene->commands[before].data.text.rotation = rotation;
}

void vg_path_clear(vg_scene_t* scene) {
    MD_ASSERT(scene);
    md_array_shrink(scene->path_points, 0);
    md_array_shrink(scene->path_verbs, 0);
}

void vg_path_move_to(vg_scene_t* scene, md_float2_t pos) {
    MD_ASSERT(scene);
    md_array_push(scene->path_verbs, VG_PATH_MOVE_TO, scene->alloc);
    md_array_push(scene->path_points, pos, scene->alloc);
}

void vg_path_line_to(vg_scene_t* scene, md_float2_t pos) {
    MD_ASSERT(scene);
    if (md_array_size(scene->path_verbs) == 0) {
        vg_path_move_to(scene, pos);
        return;
    }
    md_array_push(scene->path_verbs, VG_PATH_LINE_TO, scene->alloc);
    md_array_push(scene->path_points, pos, scene->alloc);
}

void vg_path_bezier_cubic_curve_to(vg_scene_t* scene, md_float2_t p1, md_float2_t p2, md_float2_t p3) {
    MD_ASSERT(scene);
    md_array_push(scene->path_verbs, VG_PATH_CUBIC_TO, scene->alloc);
    md_array_push(scene->path_points, p1, scene->alloc);
    md_array_push(scene->path_points, p2, scene->alloc);
    md_array_push(scene->path_points, p3, scene->alloc);
}

void vg_path_close(vg_scene_t* scene) {
    MD_ASSERT(scene);
    if (md_array_size(scene->path_verbs) == 0) return;
    md_array_push(scene->path_verbs, VG_PATH_CLOSE, scene->alloc);
}

static void commit_path(vg_scene_t* scene, const vg_style_t* style) {
    if (!scene || md_array_size(scene->path_verbs) == 0) return;
    vg_cmd_t cmd = make_cmd(VG_CMD_PATH, style);
    cmd.data.path.verb_offset  = push_verbs(scene, scene->path_verbs, md_array_size(scene->path_verbs));
    cmd.data.path.verb_count   = (uint32_t)md_array_size(scene->path_verbs);
    cmd.data.path.point_offset = push_points(scene, scene->path_points, md_array_size(scene->path_points));
    cmd.data.path.point_count  = (uint32_t)md_array_size(scene->path_points);
    push_cmd(scene, &cmd);
    vg_path_clear(scene);
}

void vg_path_fill_styled(vg_scene_t* scene, const vg_style_t* style) {
    commit_path(scene, style);
}
