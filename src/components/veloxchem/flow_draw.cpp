#define IMGUI_DEFINE_MATH_OPERATORS

#include "flow_draw.h"

#include <md_vector_graphics.h>

#include <core/md_allocator.h>
#include <core/md_common.h>

#include <imgui_internal.h>

#include <math.h>
#include <string.h>

ImVec2 flow_view_to_screen(const FlowView& view, vec2_t graph_point) {
    const ImVec2 size = view.area_max - view.area_min;
    return ImVec2(
        view.area_min.x + (graph_point.x * view.zoom + view.pan.x) * size.x,
        view.area_min.y + (graph_point.y * view.zoom + view.pan.y) * size.y);
}

// A ribbon edge is a cubic with horizontal tangents at both ends: it leaves its source flat and
// arrives at its destination flat, which is what makes a Sankey read as flow rather than as
// wires. 'curvature' slides between a straight line and a full S.
static inline ImVec2 flow_edge_point(ImVec2 a, ImVec2 b, float curvature, float t) {
    const float dx = (b.x - a.x) * 0.5f * curvature;
    const ImVec2 c0 = { a.x + dx, a.y };
    const ImVec2 c1 = { b.x - dx, b.y };

    const float u = 1.0f - t;
    const float w0 = u * u * u;
    const float w1 = 3.0f * u * u * t;
    const float w2 = 3.0f * u * t * t;
    const float w3 = t * t * t;

    return ImVec2(
        a.x * w0 + c0.x * w1 + c1.x * w2 + b.x * w3,
        a.y * w0 + c0.y * w1 + c1.y * w2 + b.y * w3);
}

// Blend toward white. Used for node interiors: bright enough to read as a surface, tinted enough
// to still say which node it is.
static inline vec4_t flow_toward_white(vec4_t c, float keep) {
    const float k = ImClamp(keep, 0.0f, 1.0f);
    return vec4_t{ 1.0f - (1.0f - c.x) * k, 1.0f - (1.0f - c.y) * k, 1.0f - (1.0f - c.z) * k, c.w };
}

// Darken a colour until it stands out against a near-white interior, by scaling toward black until
// its luminance is at most 'target'. Scaling keeps the hue and spends only what is needed - a
// fixed multiplier would either leave pale colours invisible or crush dark ones to black.
//
// The two targets below were picked by measuring WCAG contrast against the tinted interior across
// the palettes this actually draws with (the role colours, the ImPlot Deep colormap a group can
// take, and the pale extremes a user could choose). 0.30 keeps every label at or above 5.4:1,
// comfortably past AA; 0.62 keeps every border at or above 2.2:1, which is what stops a pale group
// colour - or white hydrogen, back when these were CPK - from drawing an invisible box.
static inline vec4_t flow_cap_luminance(vec4_t c, float target) {
    const float lum = 0.2126f * c.x + 0.7152f * c.y + 0.0722f * c.z;
    if (lum <= target || lum <= 1.0e-4f) return c;
    const float k = target / lum;
    return vec4_t{ c.x * k, c.y * k, c.z * k, c.w };
}

static inline vec4_t flow_readable_text(vec4_t c)   { return flow_cap_luminance(c, 0.30f); }
static inline vec4_t flow_visible_border(vec4_t c)  { return flow_cap_luminance(c, 0.62f); }

static inline ImU32 flow_color_u32(vec4_t color, float alpha_scale) {
    ImVec4 c = { color.x, color.y, color.z, color.w * alpha_scale };
    return ImGui::ColorConvertFloat4ToU32(c);
}

static inline md_vg_color_t flow_color_vg(vec4_t color, float alpha_scale) {
    const float a = CLAMP(color.w * alpha_scale, 0.0f, 1.0f);
    return MD_VG_COLOR_RGBA(
        (uint32_t)(CLAMP(color.x, 0.0f, 1.0f) * 255.0f + 0.5f),
        (uint32_t)(CLAMP(color.y, 0.0f, 1.0f) * 255.0f + 0.5f),
        (uint32_t)(CLAMP(color.z, 0.0f, 1.0f) * 255.0f + 0.5f),
        (uint32_t)(a * 255.0f + 0.5f));
}

// Emits the ribbon as a gradient-shaded triangle strip. A stroked bezier (what the old code did)
// cannot express a band that changes thickness or colour along its length; a strip can, and it
// costs the same.
static void flow_ribbon_imgui(ImDrawList* draw_list, ImVec2 p0, ImVec2 p1, ImVec2 q0, ImVec2 q1,
                              ImU32 col_src, ImU32 col_dst, float curvature, int segments)
{
    segments = ImClamp(segments, 2, 128);

    const int num_vtx = (segments + 1) * 2;
    const int num_idx = segments * 6;
    draw_list->PrimReserve(num_idx, num_vtx);

    const ImVec2 uv = draw_list->_Data->TexUvWhitePixel;
    const unsigned int base = draw_list->_VtxCurrentIdx;

    for (int i = 0; i <= segments; ++i) {
        const float t = (float)i / (float)segments;
        const ImU32 col = ImGui::ColorConvertFloat4ToU32(
            ImLerp(ImGui::ColorConvertU32ToFloat4(col_src), ImGui::ColorConvertU32ToFloat4(col_dst), t));
        draw_list->PrimWriteVtx(flow_edge_point(p0, q0, curvature, t), uv, col);
        draw_list->PrimWriteVtx(flow_edge_point(p1, q1, curvature, t), uv, col);
    }

    for (int i = 0; i < segments; ++i) {
        const unsigned int a = base + (unsigned int)(i * 2);
        draw_list->PrimWriteIdx((ImDrawIdx)(a + 0));
        draw_list->PrimWriteIdx((ImDrawIdx)(a + 1));
        draw_list->PrimWriteIdx((ImDrawIdx)(a + 3));
        draw_list->PrimWriteIdx((ImDrawIdx)(a + 0));
        draw_list->PrimWriteIdx((ImDrawIdx)(a + 3));
        draw_list->PrimWriteIdx((ImDrawIdx)(a + 2));
    }
}

// x is monotonic along a ribbon, so the hit test is a walk to the segment containing the cursor
// and a vertical span test there. No sampling of the whole polygon, no winding rule.
static bool flow_ribbon_contains(ImVec2 p0, ImVec2 p1, ImVec2 q0, ImVec2 q1, float curvature,
                                 int segments, ImVec2 point)
{
    segments = ImClamp(segments, 2, 128);
    if (point.x < ImMin(p0.x, q0.x) || point.x > ImMax(p0.x, q0.x)) {
        return false;
    }

    ImVec2 prev_top = flow_edge_point(p0, q0, curvature, 0.0f);
    ImVec2 prev_bot = flow_edge_point(p1, q1, curvature, 0.0f);
    for (int i = 1; i <= segments; ++i) {
        const float t = (float)i / (float)segments;
        const ImVec2 top = flow_edge_point(p0, q0, curvature, t);
        const ImVec2 bot = flow_edge_point(p1, q1, curvature, t);

        const float x0 = ImMin(prev_top.x, top.x);
        const float x1 = ImMax(prev_top.x, top.x);
        if (point.x >= x0 && point.x <= x1) {
            const float span = (top.x - prev_top.x);
            const float u = (fabsf(span) > 1.0e-6f) ? (point.x - prev_top.x) / span : 0.0f;
            const float y_top = prev_top.y + (top.y - prev_top.y) * u;
            const float y_bot = prev_bot.y + (bot.y - prev_bot.y) * u;
            if (point.y >= ImMin(y_top, y_bot) && point.y <= ImMax(y_top, y_bot)) {
                return true;
            }
        }
        prev_top = top;
        prev_bot = bot;
    }
    return false;
}

static inline uint8_t flow_state_of(const FlowEmphasis& e, uint32_t cut_idx) {
    return (e.node_state && cut_idx < e.num_state_entries) ? e.node_state[cut_idx] : 0u;
}

// "In the foreground" for dimming purposes. Partial counts: a group holding a selected atom should
// not fade out just because its other atoms are unselected.
static inline bool flow_is_emphasized(const FlowEmphasis& e, uint32_t cut_idx) {
    return (flow_state_of(e, cut_idx) & (FlowNodeState_Selected | FlowNodeState_Partial | FlowNodeState_Highlight)) != 0;
}

// Text that fits inside a node reads as a label ON the thing; text beside it reads as an
// annotation NEXT to the thing, and at three columns those annotations collide with the ribbons
// leaving the next column over. So: inside when it fits, and nothing when it does not.
//
// Lines are dropped from the BOTTOM, so the order they arrive in is an order of importance: what
// the node IS, then what it is MADE OF, then how big it is.
static void flow_label_inside(ImDrawList* draw_list, ImVec2 min, ImVec2 max,
                              const char* const* lines, const ImU32* colors, int num_lines, float pad)
{
    const float avail_w = (max.x - min.x) - pad * 2.0f;
    const float avail_h = (max.y - min.y) - 2.0f;
    if (avail_w <= 1.0f || avail_h <= 1.0f || num_lines <= 0) return;

    int fit = 0;
    float total_h = 0.0f;
    for (int i = 0; i < num_lines; ++i) {
        if (!lines[i] || !lines[i][0]) continue;
        const ImVec2 sz = ImGui::CalcTextSize(lines[i]);
        if (sz.x > avail_w || total_h + sz.y > avail_h) break;
        total_h += sz.y;
        fit = i + 1;
    }
    if (fit == 0) return;

    const ImVec2 center = { (min.x + max.x) * 0.5f, (min.y + max.y) * 0.5f };
    float y = center.y - total_h * 0.5f;
    for (int i = 0; i < fit; ++i) {
        if (!lines[i] || !lines[i][0]) continue;
        const ImVec2 sz = ImGui::CalcTextSize(lines[i]);
        draw_list->AddText({ center.x - sz.x * 0.5f, y }, colors[i], lines[i]);
        y += sz.y;
    }
}

FlowHit flow_draw_imgui(ImDrawList* draw_list, const FlowView& view, const md_flow_graph_t* graph,
                        const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                        const FlowDrawStyle& style, const FlowEmphasis& emphasis)
{
    FlowHit hit = {};
    if (!draw_list || !graph || !cut || !layout) return hit;

    const ImVec2 mouse = ImGui::GetMousePos();
    const bool has_highlight = emphasis.any_selected || emphasis.hover_node >= 0 || emphasis.hover_link >= 0;

    // Ribbons first, so nodes read as sitting on top of the flow rather than behind it.
    const size_t num_ribbons = md_array_size(layout->ribbons);
    for (size_t i = 0; i < num_ribbons; ++i) {
        const md_flow_layout_ribbon_t* r = layout->ribbons + i;
        const md_flow_link_t* link = cut->links + r->link_idx;

        const md_flow_node_t* src = md_flow_cut_node(cut, graph, link->src);
        const md_flow_node_t* dst = md_flow_cut_node(cut, graph, link->dst);
        if (!src || !dst) continue;

        const ImVec2 p0 = flow_view_to_screen(view, r->p0);
        const ImVec2 p1 = flow_view_to_screen(view, r->p1);
        const ImVec2 q0 = flow_view_to_screen(view, r->q0);
        const ImVec2 q1 = flow_view_to_screen(view, r->q1);

        // Too thin to read as a quantity: skip it rather than lay another hair over the picture.
        // The hit test goes with it - a ribbon nobody can see should not be hoverable either.
        if (ImAbs(p1.y - p0.y) < style.ribbon_min_px && ImAbs(q1.y - q0.y) < style.ribbon_min_px) {
            continue;
        }

        const bool touched = !has_highlight
            || (int32_t)r->link_idx == emphasis.hover_link
            || (int32_t)link->src == emphasis.hover_node
            || (int32_t)link->dst == emphasis.hover_node
            || flow_is_emphasized(emphasis, link->src)
            || flow_is_emphasized(emphasis, link->dst);

        const float alpha = style.ribbon_alpha * (touched ? 1.0f : style.dim_factor);
        flow_ribbon_imgui(draw_list, p0, p1, q0, q1,
                          flow_color_u32(src->color, alpha), flow_color_u32(dst->color, alpha),
                          style.ribbon_curvature, style.ribbon_segments);

        if (hit.link < 0 && flow_ribbon_contains(p0, p1, q0, q1, style.ribbon_curvature, style.ribbon_segments, mouse)) {
            hit.link = (int32_t)r->link_idx;
        }
    }

    // Group bands. Outboard of the first and last columns, so they read as an outer layer wrapping
    // the atoms rather than as another thing in the same row. A middle column has no free side, so
    // its bands are skipped - nothing there is grouped in this design anyway.
    const size_t num_bands = md_array_size(layout->bands);
    for (size_t i = 0; i < num_bands; ++i) {
        const md_flow_layout_band_t* b = layout->bands + i;
        const md_flow_node_t* node = graph->nodes + b->node;

        const bool leftmost  = (b->column == 0);
        const bool rightmost = (b->column + 1 == graph->num_columns);
        if (!leftmost && !rightmost) continue;

        const ImVec2 min = flow_view_to_screen(view, b->min);
        const ImVec2 max = flow_view_to_screen(view, b->max);

        // Nested bands step further out, so a group of groups does not sit on top of its members.
        const float out = style.band_offset_px + (float)b->depth * (style.band_offset_px + style.band_width_px);
        const float x1 = leftmost ? (min.x - out) : (max.x + out);
        const float x0 = leftmost ? (x1 - style.band_width_px) : (x1 + style.band_width_px);

        const ImVec2 band_min = { ImMin(x0, x1), min.y };
        const ImVec2 band_max = { ImMax(x0, x1), max.y };

        const bool touched = !has_highlight || (int32_t)b->node == emphasis.hover_node
                             || flow_is_emphasized(emphasis, b->node);
        const float alpha = touched ? 1.0f : style.dim_factor;
        draw_list->AddRectFilled(band_min, band_max, flow_color_u32(node->color, alpha), style.band_width_px * 0.5f);

        if (mouse.x >= band_min.x && mouse.x <= band_max.x && mouse.y >= band_min.y && mouse.y <= band_max.y) {
            hit.node = (int32_t)b->node;
            hit.link = -1;
        }

        if (!style.show_band_labels || (band_max.y - band_min.y) < style.label_min_px) continue;

        char buf[128];
        snprintf(buf, sizeof(buf), "%.*s", (int)node->label.len, node->label.ptr);
        const ImVec2 text_size = ImGui::CalcTextSize(buf);
        // Vertical centring against the band, pushed to its outboard side.
        const ImVec2 text_pos = {
            leftmost ? (band_min.x - style.band_label_gap_px - text_size.x) : (band_max.x + style.band_label_gap_px),
            (band_min.y + band_max.y) * 0.5f - text_size.y * 0.5f,
        };
        draw_list->AddText(text_pos, style.text_color, buf);
    }

    const size_t num_nodes = md_array_size(layout->nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        const md_flow_layout_node_t* ln = layout->nodes + i;
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, ln->cut_idx);
        if (!node) continue;

        const ImVec2 min = flow_view_to_screen(view, ln->min);
        const ImVec2 max = flow_view_to_screen(view, ln->max);

        const uint8_t node_state = flow_state_of(emphasis, ln->cut_idx);
        const bool selected = (node_state & FlowNodeState_Selected) != 0;
        const bool partial  = (node_state & FlowNodeState_Partial) != 0;
        const bool hovered  = ((node_state & FlowNodeState_Highlight) != 0) || (int32_t)ln->cut_idx == emphasis.hover_node;
        const bool touched  = !has_highlight || selected || partial || hovered;
        const float alpha = touched ? 1.0f : style.dim_factor;

        // Dimmed nodes lose their tint as well as their alpha: the interior is nearly white
        // anyway, so fading only the border would leave a row of identical pale boxes.
        const vec4_t fill = flow_toward_white(node->color, style.node_fill_tint * (touched ? 1.0f : 0.5f));
        draw_list->AddRectFilled(min, max, flow_color_u32(fill, 1.0f), style.node_rounding);

        // Colour lives on the border now, and state is its weight: selection reads as a heavier
        // outline of the SAME colour rather than a different colour, so state and identity do not
        // compete for the one channel.
        // Fully selected reads heaviest, then partial, then merely hovered. One channel, ordered
        // by how much the user has committed to the thing.
        const float width = selected ? style.border_width_selected
                          : partial  ? style.border_width_partial
                          : hovered  ? style.border_width_hover
                                     : style.border_width;
        const vec4_t border = flow_visible_border(node->color.w > 0.0f ? node->color : vec4_t{0, 0, 0, 0.35f});
        draw_list->AddRect(min, max, flow_color_u32(border, alpha), style.node_rounding, 0, width);

        if (mouse.x >= min.x && mouse.x <= max.x && mouse.y >= min.y && mouse.y <= max.y) {
            hit.node = (int32_t)ln->cut_idx;
            hit.link = -1;  // a node under the cursor wins over the ribbons behind it
        }

        // Semantic LOD: a node too short for its label does not get one.
        const float height = max.y - min.y;
        if (height < style.label_min_px) continue;

        char name[96];
        char sub[96];
        char pct[32];
        snprintf(name, sizeof(name), "%.*s", (int)node->label.len, node->label.ptr);
        snprintf(sub,  sizeof(sub),  "%.*s", (int)node->sublabel.len, node->sublabel.ptr ? node->sublabel.ptr : "");
        snprintf(pct,  sizeof(pct),  "%.1f%%", node->weight * 100.0f);

        // The name carries the node's colour, darkened to stay legible. The sublabel is what the
        // node is MADE of, so it sits quieter; the percentage stays neutral, because it is a
        // magnitude and not an identity.
        const ImU32 name_col = style.color_labels
            ? flow_color_u32(flow_readable_text(node->color), alpha)
            : style.text_color;
        const ImU32 sub_col = (style.text_color & 0x00FFFFFFu) | 0xA0000000u;

        const char* lines[3]   = { name, sub, (height >= style.percent_min_px) ? pct : nullptr };
        const ImU32  colors[3] = { name_col, sub_col, style.text_color };
        flow_label_inside(draw_list, min, max, lines, colors, 3, style.label_pad_px);
    }

    return hit;
}

void flow_draw_vg(md_vg_scene_t* scene, const FlowView& view, const md_flow_graph_t* graph,
                  const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                  const FlowDrawStyle& style)
{
    if (!scene || !graph || !cut || !layout) return;

    auto to_vg = [&view](vec2_t p) -> vec2_t {
        const ImVec2 s = flow_view_to_screen(view, p);
        return vec2_t{ s.x, s.y };
    };

    const size_t num_ribbons = md_array_size(layout->ribbons);
    for (size_t i = 0; i < num_ribbons; ++i) {
        const md_flow_layout_ribbon_t* r = layout->ribbons + i;
        const md_flow_link_t* link = cut->links + r->link_idx;

        const md_flow_node_t* src = md_flow_cut_node(cut, graph, link->src);
        const md_flow_node_t* dst = md_flow_cut_node(cut, graph, link->dst);
        if (!src || !dst) continue;

        const vec2_t p0 = to_vg(r->p0);
        const vec2_t p1 = to_vg(r->p1);
        const vec2_t q0 = to_vg(r->q0);
        const vec2_t q1 = to_vg(r->q1);

        // Native cubics rather than the sampled strip the screen backend uses: an SVG viewer
        // resamples for its own resolution, so handing it curves beats handing it triangles.
        const float dx = (q0.x - p0.x) * 0.5f * style.ribbon_curvature;

        md_vg_path_clear(scene);
        md_vg_path_move_to(scene, p0);
        md_vg_path_bezier_cubic_curve_to(scene, vec2_t{p0.x + dx, p0.y}, vec2_t{q0.x - dx, q0.y}, q0);
        md_vg_path_line_to(scene, q1);
        md_vg_path_bezier_cubic_curve_to(scene, vec2_t{q1.x - dx, q1.y}, vec2_t{p1.x + dx, p1.y}, p1);
        md_vg_path_close(scene);

        md_vg_style_t vg_style = {};
        vg_style.fill.type = MD_VG_PAINT_LINEAR_GRADIENT;
        vg_style.fill.data.linear_gradient.p0   = p0;
        vg_style.fill.data.linear_gradient.p1   = q0;
        vg_style.fill.data.linear_gradient.col0 = flow_color_vg(src->color, style.ribbon_alpha);
        vg_style.fill.data.linear_gradient.col1 = flow_color_vg(dst->color, style.ribbon_alpha);
        md_vg_path_fill_styled(scene, &vg_style);
    }

    const size_t num_bands = md_array_size(layout->bands);
    for (size_t i = 0; i < num_bands; ++i) {
        const md_flow_layout_band_t* b = layout->bands + i;
        const md_flow_node_t* node = graph->nodes + b->node;

        const bool leftmost  = (b->column == 0);
        const bool rightmost = (b->column + 1 == graph->num_columns);
        if (!leftmost && !rightmost) continue;

        const vec2_t min = to_vg(b->min);
        const vec2_t max = to_vg(b->max);
        const float out = style.band_offset_px + (float)b->depth * (style.band_offset_px + style.band_width_px);
        const float x1 = leftmost ? (min.x - out) : (max.x + out);
        const float x0 = leftmost ? (x1 - style.band_width_px) : (x1 + style.band_width_px);

        md_vg_add_rect_filled(scene, vec2_t{MIN(x0, x1), min.y}, vec2_t{MAX(x0, x1), max.y},
                              flow_color_vg(node->color, 1.0f), style.band_width_px * 0.5f);

        if (!style.show_band_labels || (max.y - min.y) < style.label_min_px) continue;

        char buf[128];
        snprintf(buf, sizeof(buf), "%.*s", (int)node->label.len, node->label.ptr);
        const float font_size = 12.0f;
        const float text_w = font_size * 0.5f * (float)strlen(buf);
        const vec2_t pos = {
            leftmost ? (MIN(x0, x1) - style.band_label_gap_px - text_w) : (MAX(x0, x1) + style.band_label_gap_px),
            (min.y + max.y) * 0.5f + font_size * 0.35f,
        };
        md_vg_add_text(scene, pos, font_size, MD_VG_COLOR_RGB(0, 0, 0), str_from_cstr(buf));
    }

    const size_t num_nodes = md_array_size(layout->nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        const md_flow_layout_node_t* ln = layout->nodes + i;
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, ln->cut_idx);
        if (!node) continue;

        const vec2_t min = to_vg(ln->min);
        const vec2_t max = to_vg(ln->max);

        // Same treatment as the screen: near-white interior, colour on the border.
        md_vg_add_rect_filled(scene, min, max, flow_color_vg(flow_toward_white(node->color, style.node_fill_tint), 1.0f),
                              style.node_rounding);
        md_vg_add_rect(scene, min, max, flow_color_vg(flow_visible_border(node->color), 1.0f), style.node_rounding, style.border_width);

        if (max.y - min.y < style.label_min_px) continue;

        char buf[128];
        snprintf(buf, sizeof(buf), "%.*s  %.1f%%", (int)node->label.len, node->label.ptr, node->weight * 100.0f);
        const float font_size = 12.0f;
        // Centred inside the node, matching the screen backend. Width is estimated rather than
        // measured: md_vg has no font metrics, and the estimate only decides whether it fits.
        const float text_w = font_size * 0.5f * (float)strlen(buf);
        if (text_w > (max.x - min.x) - 4.0f) continue;
        const vec2_t pos = { (min.x + max.x) * 0.5f - text_w * 0.5f, (min.y + max.y) * 0.5f + font_size * 0.35f };
        md_vg_add_text(scene, pos, font_size,
                       style.color_labels ? flow_color_vg(flow_readable_text(node->color), 1.0f) : MD_VG_COLOR_RGB(0, 0, 0),
                       str_from_cstr(buf));
    }
}
