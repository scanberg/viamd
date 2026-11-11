#define IMGUI_DEFINE_MATH_OPERATORS

#include "flow_draw.h"

#include <md_vector_graphics.h>

#include <core/md_allocator.h>
#include <core/md_common.h>

#include <imgui_internal.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
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
    const float a = MD_CLAMP(color.w * alpha_scale, 0.0f, 1.0f);
    return MD_VG_COLOR_RGBA(
        (uint32_t)(MD_CLAMP(color.x, 0.0f, 1.0f) * 255.0f + 0.5f),
        (uint32_t)(MD_CLAMP(color.y, 0.0f, 1.0f) * 255.0f + 0.5f),
        (uint32_t)(MD_CLAMP(color.z, 0.0f, 1.0f) * 255.0f + 0.5f),
        (uint32_t)(a * 255.0f + 0.5f));
}

// Emits the ribbon as a gradient-shaded triangle strip. A stroked bezier (what the old code did)
// cannot express a band that changes thickness or colour along its length; a strip can, and it
// costs the same.
//
// Its long edges are anti-aliased the way ImGui's own fills are: a fringe of about a pixel along each, from the
// colour to nothing, the solid part pulled in by half of it so the ribbon keeps its width. 'fringe' says which edges
// get one (1 the top, 2 the bottom): the stripes of a ribbon meet edge to edge, and a fringe between two would
// show as a seam. The ends need none: they sit under the nodes.
enum { FLOW_FRINGE_TOP = 1, FLOW_FRINGE_BOTTOM = 2 };

static void flow_ribbon_imgui(ImDrawList* draw_list, ImVec2 p0, ImVec2 p1, ImVec2 q0, ImVec2 q1,
                              ImU32 col_src, ImU32 col_dst, float curvature, int segments, int fringe)
{
    segments = ImClamp(segments, 2, 128);
    const bool aa = (draw_list->Flags & ImDrawListFlags_AntiAliasedFill) != 0;
    if (!aa) fringe = 0;
    const float aa_size = draw_list->_FringeScale;
    const bool top_aa = (fringe & FLOW_FRINGE_TOP) != 0;
    const bool bot_aa = (fringe & FLOW_FRINGE_BOTTOM) != 0;

    // per sample: [top fringe] top bottom [bottom fringe]
    const int per = 2 + (top_aa ? 1 : 0) + (bot_aa ? 1 : 0);
    const int num_vtx = (segments + 1) * per;
    const int num_idx = segments * 6 * (per - 1);
    draw_list->PrimReserve(num_idx, num_vtx);

    const ImVec2 uv = draw_list->_Data->TexUvWhitePixel;
    const unsigned int base = draw_list->_VtxCurrentIdx;
    const ImVec4 c_src = ImGui::ColorConvertU32ToFloat4(col_src);
    const ImVec4 c_dst = ImGui::ColorConvertU32ToFloat4(col_dst);

    for (int i = 0; i <= segments; ++i) {
        const float t = (float)i / (float)segments;
        const ImVec4 cf = ImLerp(c_src, c_dst, t);
        const ImU32 col = ImGui::ColorConvertFloat4ToU32(cf);
        const ImU32 col_trans = col & ~IM_COL32_A_MASK;
        const ImVec2 top = flow_edge_point(p0, q0, curvature, t);
        const ImVec2 bot = flow_edge_point(p1, q1, curvature, t);
        // the normals of the edges, from the tangents, pointing away from the other edge
        const float dt = 1.0f / (float)segments * 0.5f;
        const float ta = ImMax(t - dt, 0.0f), tb = ImMin(t + dt, 1.0f);
        ImVec2 tan_top = flow_edge_point(p0, q0, curvature, tb) - flow_edge_point(p0, q0, curvature, ta);
        ImVec2 tan_bot = flow_edge_point(p1, q1, curvature, tb) - flow_edge_point(p1, q1, curvature, ta);
        ImVec2 n_top = ImVec2(tan_top.y, -tan_top.x);
        ImVec2 n_bot = ImVec2(tan_bot.y, -tan_bot.x);
        const float lt = ImSqrt(ImLengthSqr(n_top)), lb = ImSqrt(ImLengthSqr(n_bot));
        n_top = lt > 1.0e-6f ? n_top / lt : ImVec2(0.0f, -1.0f);
        n_bot = lb > 1.0e-6f ? n_bot / lb : ImVec2(0.0f, 1.0f);
        const ImVec2 across = top - bot;
        if (n_top.x * across.x + n_top.y * across.y < 0.0f) n_top = -n_top;
        if (n_bot.x * across.x + n_bot.y * across.y > 0.0f) n_bot = -n_bot;
        // pulled in by half the fringe, but never past the middle of a thin ribbon
        const float thick = ImSqrt(ImLengthSqr(across));
        const float in_top = top_aa ? ImMin(aa_size * 0.5f, thick * 0.5f) : 0.0f;
        const float in_bot = bot_aa ? ImMin(aa_size * 0.5f, thick * 0.5f) : 0.0f;
        if (top_aa) draw_list->PrimWriteVtx(top + n_top * (aa_size * 0.5f), uv, col_trans);
        draw_list->PrimWriteVtx(top - n_top * in_top, uv, col);
        draw_list->PrimWriteVtx(bot - n_bot * in_bot, uv, col);
        if (bot_aa) draw_list->PrimWriteVtx(bot + n_bot * (aa_size * 0.5f), uv, col_trans);
    }

    for (int i = 0; i < segments; ++i) {
        const unsigned int a = base + (unsigned int)(i * per);
        const unsigned int b = a + (unsigned int)per;
        for (int k = 0; k + 1 < per; ++k) {
            draw_list->PrimWriteIdx((ImDrawIdx)(a + k));
            draw_list->PrimWriteIdx((ImDrawIdx)(a + k + 1));
            draw_list->PrimWriteIdx((ImDrawIdx)(b + k + 1));
            draw_list->PrimWriteIdx((ImDrawIdx)(a + k));
            draw_list->PrimWriteIdx((ImDrawIdx)(b + k + 1));
            draw_list->PrimWriteIdx((ImDrawIdx)(b + k));
        }
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

// The stripes of link l to draw, or 0 when it is drawn whole
static uint32_t flow_stripe_count(const FlowStripes* st, uint32_t l) {
    if (!st || !st->offset || !st->frac || l >= st->num_links) return 0;
    const uint32_t n = st->offset[l + 1] - st->offset[l];
    if (n < 2 || st->frac[st->offset[l] + 1] < st->min_second) return 0;
    return n;
}

static inline ImVec2 flow_lerp(ImVec2 a, ImVec2 b, float t) { return ImVec2(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t); }


static inline vec4_t flow_with_alpha(vec4_t c, float a) { c.w *= a; return c; }
static inline vec4_t flow_scale_rgb(vec4_t c, float k)  { return vec4_t{ c.x * k, c.y * k, c.z * k, c.w }; }
static inline vec4_t flow_vec4(ImU32 c) { const ImVec4 v = ImGui::ColorConvertU32ToFloat4(c); return vec4_t{ v.x, v.y, v.z, v.w }; }

// Every other stripe a little lighter, so neighbours read apart without new hues
static inline vec4_t flow_stripe_tone(vec4_t c, uint32_t j) {
    if ((j & 1u) == 0) return c;
    return vec4_t{ c.x + (1.0f - c.x) * 0.45f, c.y + (1.0f - c.y) * 0.45f, c.z + (1.0f - c.z) * 0.45f, c.w };
}

// A share as the diagram prints it: whole percents from 1%, one decimal below, and a floor
static void flow_format_share(char* buf, size_t cap, float w) {
    const float pct = w * 100.0f;
    if (pct >= 0.995f)     snprintf(buf, cap, "%.0f%%", pct);
    else if (pct >= 0.05f) snprintf(buf, cap, "%.1f%%", pct);
    else                   snprintf(buf, cap, "<0.1%%");
}

// ---------------------------------------------------------------------------------------------------------
// The canvas: the few primitives the diagram is drawn with. One drawing routine (flow_draw_common) speaks
// to it, and each backend implements it, so the screen and the export cannot drift apart.
// ---------------------------------------------------------------------------------------------------------

struct FlowCanvas {
    void* user;
    float base_size;    // text sizes, in pixels
    float small_size;
    float (*measure)(void* user, const char* text, float size);
    void  (*fill_rect)(void* user, ImVec2 min, ImVec2 max, vec4_t col, float rounding);
    void  (*stroke_rect)(void* user, ImVec2 min, ImVec2 max, vec4_t col, float rounding, float width);
    // pos.y is the top of the line; pos.x its start (align -1), middle (0) or end (+1)
    void  (*text)(void* user, ImVec2 pos, float size, vec4_t col, int align, const char* text);
    // 'fringe' (FLOW_FRINGE_*): which long edges are outer edges of the ribbon, to anti-alias
    void  (*ribbon)(void* user, ImVec2 p0, ImVec2 p1, ImVec2 q0, ImVec2 q1, vec4_t col_src, vec4_t col_dst, int fringe);
    // a line turned to read upwards, centred on 'centre'
    void  (*text_up)(void* user, ImVec2 centre, float size, vec4_t col, const char* text);
    // null where nothing is clicked (the export): the controls are left out
    void  (*fill_tri)(void* user, ImVec2 a, ImVec2 b, ImVec2 c, vec4_t col);
};

static inline float flow_line_h(float size) { return size * 1.18f; }

// ---- ImGui ----

struct FlowImCanvas { ImDrawList* dl; float curvature; int segments; };

static float im_measure(void*, const char* t, float size) { return ImGui::GetFont()->CalcTextSizeA(size, FLT_MAX, 0.0f, t).x; }
static void im_fill(void* u, ImVec2 a, ImVec2 b, vec4_t c, float r) { ((FlowImCanvas*)u)->dl->AddRectFilled(a, b, flow_color_u32(c, 1.0f), r); }
static void im_stroke(void* u, ImVec2 a, ImVec2 b, vec4_t c, float r, float w) { ((FlowImCanvas*)u)->dl->AddRect(a, b, flow_color_u32(c, 1.0f), r, 0, w); }
static void im_text(void* u, ImVec2 p, float size, vec4_t c, int align, const char* t) {
    if (align >= 0) p.x -= im_measure(u, t, size) * (align == 0 ? 0.5f : 1.0f);
    ((FlowImCanvas*)u)->dl->AddText(ImGui::GetFont(), size, ImVec2(floorf(p.x + 0.5f), floorf(p.y + 0.5f)), flow_color_u32(c, 1.0f), t);
}
static void im_ribbon(void* u, ImVec2 p0, ImVec2 p1, ImVec2 q0, ImVec2 q1, vec4_t c0, vec4_t c1, int fringe) {
    const FlowImCanvas* cv = (const FlowImCanvas*)u;
    flow_ribbon_imgui(cv->dl, p0, p1, q0, q1, flow_color_u32(c0, 1.0f), flow_color_u32(c1, 1.0f), cv->curvature, cv->segments, fringe);
}
static void im_tri(void* u, ImVec2 a, ImVec2 b, ImVec2 c, vec4_t col) { ((FlowImCanvas*)u)->dl->AddTriangleFilled(a, b, c, flow_color_u32(col, 1.0f)); }
// ImGui draws text upright only: the glyphs' vertices, as written, turned a quarter about the centre
static void im_text_up(void* u, ImVec2 c, float size, vec4_t col, const char* t) {
    ImDrawList* dl = ((FlowImCanvas*)u)->dl;
    const float w = im_measure(u, t, size);
    const int first = dl->VtxBuffer.Size;
    dl->AddText(ImGui::GetFont(), size, ImVec2(c.x - w * 0.5f, c.y - size * 0.5f), flow_color_u32(col, 1.0f), t);
    for (int i = first; i < dl->VtxBuffer.Size; ++i) {
        ImDrawVert& v = dl->VtxBuffer[i];
        const ImVec2 d = v.pos - c;
        v.pos = ImVec2(c.x + d.y, c.y - d.x);
    }
}

// ---- md_vg (SVG) ----

#define FLOW_VG_TEXT_SIZE 13.0f

struct FlowVgCanvas { md_vg_scene_t* scene; float curvature; };

// An estimate, on the generous side: it only spaces things (the viewer aligns the text itself, by its
// anchor), and too wide leaves a little air where too narrow would collide.
static float vg_measure(void*, const char* t, float size) {
    size_t n = 0;     // code points, at a sans serif's average advance
    for (const unsigned char* c = (const unsigned char*)t; *c; ++c) n += ((*c & 0xC0) != 0x80);
    return 0.6f * size * (float)n;
}
static inline vec2_t flow_v2(ImVec2 p) { return vec2_t{ p.x, p.y }; }
static void vg_fill(void* u, ImVec2 a, ImVec2 b, vec4_t c, float r) { md_vg_add_rect_filled(((FlowVgCanvas*)u)->scene, flow_v2(a), flow_v2(b), flow_color_vg(c, 1.0f), r); }
static void vg_stroke(void* u, ImVec2 a, ImVec2 b, vec4_t c, float r, float w) { md_vg_add_rect(((FlowVgCanvas*)u)->scene, flow_v2(a), flow_v2(b), flow_color_vg(c, 1.0f), r, w); }
static void vg_text(void* u, ImVec2 p, float size, vec4_t c, int align, const char* t) {
    const md_vg_text_anchor_t anchor = align < 0 ? MD_VG_TEXT_ANCHOR_START : (align == 0 ? MD_VG_TEXT_ANCHOR_MIDDLE : MD_VG_TEXT_ANCHOR_END);
    md_vg_add_text_anchored(((FlowVgCanvas*)u)->scene, flow_v2(p), size, flow_color_vg(c, 1.0f), anchor, md_str_from_cstr(t));
}
static void vg_text_up(void* u, ImVec2 c, float size, vec4_t col, const char* t) {
    md_vg_add_text_rotated(((FlowVgCanvas*)u)->scene, vec2_t{ c.x, c.y - size * 0.5f }, size, flow_color_vg(col, 1.0f), MD_VG_TEXT_ANCHOR_MIDDLE, -90.0f, md_str_from_cstr(t));
}
static void vg_ribbon(void* u, ImVec2 p0, ImVec2 p1, ImVec2 q0, ImVec2 q1, vec4_t c0, vec4_t c1, int) {
    // Native cubics rather than the sampled strip of the screen: an SVG viewer resamples for its own
    // resolution, so handing it curves beats handing it triangles.
    const FlowVgCanvas* cv = (const FlowVgCanvas*)u;
    md_vg_scene_t* scene = cv->scene;
    const float dx = (q0.x - p0.x) * 0.5f * cv->curvature;
    md_vg_path_clear(scene);
    md_vg_path_move_to(scene, flow_v2(p0));
    md_vg_path_bezier_cubic_curve_to(scene, vec2_t{ p0.x + dx, p0.y }, vec2_t{ q0.x - dx, q0.y }, flow_v2(q0));
    md_vg_path_line_to(scene, flow_v2(q1));
    md_vg_path_bezier_cubic_curve_to(scene, vec2_t{ q1.x - dx, q1.y }, vec2_t{ p1.x + dx, p1.y }, flow_v2(p1));
    md_vg_path_close(scene);
    md_vg_style_t st = {};
    st.fill.type = MD_VG_PAINT_LINEAR_GRADIENT;
    st.fill.data.linear_gradient.p0   = flow_v2(p0);
    st.fill.data.linear_gradient.p1   = flow_v2(q0);
    st.fill.data.linear_gradient.col0 = flow_color_vg(c0, 1.0f);
    st.fill.data.linear_gradient.col1 = flow_color_vg(c1, 1.0f);
    md_vg_path_fill_styled(scene, &st);
}

// ---------------------------------------------------------------------------------------------------------
// Labels
// ---------------------------------------------------------------------------------------------------------

enum { FLOW_LABEL_LINES = 3, FLOW_LABEL_CAP = 96 };

struct FlowLabel {
    bool   visible;
    int    side;            // -1: left of the bar, +1: right of it, 0: on it (a middle column)
    float  anchor_x;        // right edge for -1, left edge for +1, centre for 0
    float  top, width, height;
    int    num_lines;
    char   line[FLOW_LABEL_LINES][FLOW_LABEL_CAP];
    float  size[FLOW_LABEL_LINES];
    bool   secondary[FLOW_LABEL_LINES];
};

static int flow_side_index(int side) { return side < 0 ? 0 : (side > 0 ? 2 : 1); }

static int flow_column_side(uint32_t column, uint32_t num_columns) {
    if (num_columns <= 1) return +1;
    if (column == 0) return -1;
    if (column + 1 == num_columns) return +1;
    return 0;
}

// Ancestors in the graph of a node (an "Others" counts from the parent it folds the children of)
static uint32_t flow_depth(const md_flow_graph_t* graph, const md_flow_node_t* node) {
    const size_t n = md_array_size(graph->nodes);
    uint32_t d = 0;
    for (uint32_t p = node->parent; p < n && d < 64; p = graph->nodes[p].parent) ++d;
    return d;
}

bool flow_node_is_container(const md_flow_graph_t* graph, const md_flow_cut_t* cut, uint32_t cut_idx) {
    const size_t n = md_array_size(graph->nodes);
    if (cut_idx >= n) return false;     // an "Others" opens by its own rules
    (void)cut;
    for (size_t j = 0; j < n; ++j) {
        const uint32_t p = graph->nodes[j].parent;
        if (p < n && graph->nodes[p].parent == cut_idx) return true;    // a grandchild
    }
    return false;
}

// The same for every graph node at once: out[i] for graph node i
static void flow_containers(bool* out, const md_flow_graph_t* graph) {
    const size_t n = md_array_size(graph->nodes);
    for (size_t i = 0; i < n; ++i) out[i] = false;
    for (size_t j = 0; j < n; ++j) {
        const uint32_t p = graph->nodes[j].parent;
        if (p < n && graph->nodes[p].parent < n) out[graph->nodes[p].parent] = true;
    }
}

void flow_layout_indent(md_flow_layout_t* layout, const md_flow_graph_t* graph, const md_flow_cut_t* cut, float gap) {
    if (!layout || !graph || !cut) return;
    const uint32_t C = graph->num_columns;
    const size_t num_nodes = md_array_size(layout->nodes);
    if (num_nodes == 0 || C < 2) return;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    // The deepest row of each side, and a level's pitch (a bar and the gap)
    uint32_t* depth = md_temp_alloc_array(temp, uint32_t, num_nodes);
    int32_t max_depth[3] = { -1, -1, -1 };
    float pitch[3] = { 0, 0, 0 };
    for (size_t i = 0; i < num_nodes; ++i) {
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, layout->nodes[i].cut_idx);
        depth[i] = node ? flow_depth(graph, node) : 0;
        const int s = node ? flow_side_index(flow_column_side(node->column, C)) : 1;
        max_depth[s] = MD_MAX(max_depth[s], (int32_t)depth[i]);
        pitch[s] = (layout->nodes[i].max.x - layout->nodes[i].min.x) * (1.0f + gap);
    }
    for (size_t i = 0; i < md_array_size(layout->bands); ++i) {
        const md_flow_layout_band_t* b = layout->bands + i;
        const int s = flow_side_index(flow_column_side(b->column, C));
        max_depth[s] = MD_MAX(max_depth[s], (int32_t)b->depth);
    }
    auto shift = [&](int side, uint32_t d) -> float {
        if (side == 0) return 0.0f;
        const int s = flow_side_index(side);
        const float levels = (float)MD_MAX(max_depth[s] - (int32_t)d, 0);
        return (float)side * levels * pitch[s];
    };

    // The rows, and where they are by cut index for the ribbons
    const size_t cut_space = md_array_size(graph->nodes) + md_array_size(cut->other);
    float* dx = md_temp_alloc_zero_array(temp, float, cut_space ? cut_space : 1);
    for (size_t i = 0; i < num_nodes; ++i) {
        md_flow_layout_node_t* ln = layout->nodes + i;
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, ln->cut_idx);
        if (!node) continue;
        const float d = shift(flow_column_side(node->column, C), depth[i]);
        ln->min.x += d;
        ln->max.x += d;
        if (ln->cut_idx < cut_space) dx[ln->cut_idx] = d;
    }
    for (size_t i = 0; i < md_array_size(layout->ribbons); ++i) {
        md_flow_layout_ribbon_t* r = layout->ribbons + i;
        const md_flow_link_t* l = cut->links + r->link_idx;
        const float ds = l->src < cut_space ? dx[l->src] : 0.0f;
        const float dd = l->dst < cut_space ? dx[l->dst] : 0.0f;
        r->p0.x += ds; r->p1.x += ds;
        r->q0.x += dd; r->q1.x += dd;
    }
    for (size_t i = 0; i < md_array_size(layout->bands); ++i) {
        md_flow_layout_band_t* b = layout->bands + i;
        const float d = shift(flow_column_side(b->column, C), b->depth);
        b->min.x += d;
        b->max.x += d;
    }
}


// The lines of a node's label for a node 'h' pixels tall: name over share (and what it is made of, when
// there is room), or one small line when there is not.
static void flow_label_lines(FlowLabel* lb, const FlowCanvas& cv, const md_flow_node_t* node, float h, const char* share_suffix = nullptr) {
    char name[FLOW_LABEL_CAP], share[FLOW_LABEL_CAP], sub[FLOW_LABEL_CAP];
    snprintf(name, sizeof(name), "%.*s", (int)node->label.len, node->label.ptr ? node->label.ptr : "");
    snprintf(sub, sizeof(sub), "%.*s", (int)node->sublabel.len, node->sublabel.ptr ? node->sublabel.ptr : "");
    flow_format_share(share, sizeof(share), node->weight);
    if (share_suffix && share_suffix[0]) {
        const size_t len = strlen(share);
        snprintf(share + len, sizeof(share) - len, " \xC2\xB7 %s", share_suffix);
    }

    lb->num_lines = 0;
    const float two = flow_line_h(cv.base_size) + flow_line_h(cv.small_size);
    if (h >= two * 0.85f) {
        snprintf(lb->line[0], FLOW_LABEL_CAP, "%s", name);   lb->size[0] = cv.base_size;  lb->secondary[0] = false;
        snprintf(lb->line[1], FLOW_LABEL_CAP, "%s", share);  lb->size[1] = cv.small_size; lb->secondary[1] = true;
        lb->num_lines = 2;
        if (sub[0] && h >= two + flow_line_h(cv.small_size)) {
            snprintf(lb->line[2], FLOW_LABEL_CAP, "%s", sub); lb->size[2] = cv.small_size; lb->secondary[2] = true;
            lb->num_lines = 3;
        }
    } else {
        snprintf(lb->line[0], FLOW_LABEL_CAP, "%s %s", name, share);
        lb->size[0] = cv.small_size;
        lb->secondary[0] = true;
        lb->num_lines = 1;
    }
    lb->width = 0.0f;
    lb->height = 0.0f;
    for (int i = 0; i < lb->num_lines; ++i) {
        lb->width = MD_MAX(lb->width, cv.measure(cv.user, lb->line[i], lb->size[i]));
        lb->height += flow_line_h(lb->size[i]);
    }
}

// The widest a node's label can be, whichever form it takes
static float flow_label_width_any(const FlowCanvas& cv, const md_flow_node_t* node, const char* share_suffix = nullptr) {
    FlowLabel a = {}, b = {};
    flow_label_lines(&a, cv, node, 1.0e6f, share_suffix);
    flow_label_lines(&b, cv, node, 0.0f, share_suffix);
    return MD_MAX(a.width, b.width);
}

// What a collapsed group holds, after its share: "4 fragments"
static void flow_member_count(char* buf, size_t cap, const FlowDrawStyle& style, uint32_t n) {
    const char* noun = style.member_noun[n == 1 ? 0 : 1];
    snprintf(buf, cap, "%u %s", n, noun ? noun : "");
}

// The children of every graph node, counted
static void flow_child_counts(uint32_t* out, const md_flow_graph_t* graph) {
    const size_t n = md_array_size(graph->nodes);
    for (size_t i = 0; i < n; ++i) out[i] = 0;
    for (size_t i = 0; i < n; ++i) {
        const uint32_t p = graph->nodes[i].parent;
        if (p < n) out[p] += 1;
    }
}

// ---------------------------------------------------------------------------------------------------------
// The diagram
// ---------------------------------------------------------------------------------------------------------

static vec4_t flow_node_base_color(const FlowDrawStyle& style, const md_flow_node_t* node, int side) {
    if (style.color_by_side) return style.side_color[flow_side_index(side)];
    return node->color.w > 0.0f ? node->color : flow_vec4(style.node_border_color);
}

// The size of the text along a band 'width' pixels wide, and whether a name fits along 'length' of it
static float flow_band_text_size(const FlowCanvas& cv, float width) { return MD_MIN(cv.small_size, width - 3.0f); }
static bool flow_band_label_fits(const FlowCanvas& cv, const FlowDrawStyle& style, const char* name, float width, float length) {
    if (!style.show_band_labels || !cv.text_up) return false;
    const float size = flow_band_text_size(cv, width);
    return size >= 6.0f && cv.measure(cv.user, name, size) + 8.0f <= length;
}

static void flow_node_name(char* buf, size_t cap, const md_flow_node_t* node) {
    snprintf(buf, cap, "%.*s", (int)node->label.len, node->label.ptr ? node->label.ptr : "");
}

// A collapsed group's arrow beside its row, between the bars and the label: its size, at a text size
static float flow_arrow_box(const FlowCanvas& cv) { return MD_MAX(8.0f, cv.small_size * 0.9f); }
static const float FLOW_ARROW_PAD = 4.0f;

// Whether any row of a side has an arrow beside it (container: flow_containers): that side's labels stand out by its room
static void flow_sides_with_arrows(bool out[3], const md_flow_graph_t* graph, const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                                   const bool* container) {
    out[0] = out[1] = out[2] = false;
    const size_t n = md_array_size(graph->nodes);
    for (size_t i = 0; i < md_array_size(layout->nodes); ++i) {
        const uint32_t c = layout->nodes[i].cut_idx;
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, c);
        if (!node || c >= n || !container[c]) continue;
        const int side = flow_column_side(node->column, graph->num_columns);
        if (side != 0) out[flow_side_index(side)] = true;
    }
}

static FlowHit flow_draw_common(const FlowCanvas& cv, const FlowView& view, const md_flow_graph_t* graph,
                                const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                                const FlowDrawStyle& style, const FlowEmphasis& emphasis,
                                const FlowStripes* stripes, const FlowAnnotations* notes, const ImVec2* mouse_ptr)
{
    FlowHit hit = {};
    if (!graph || !cut || !layout) return hit;
    const ImVec2 mouse = mouse_ptr ? *mouse_ptr : ImVec2(-1.0e9f, -1.0e9f);
    const bool has_highlight = emphasis.any_selected || emphasis.hover_node >= 0 || emphasis.hover_link >= 0;
    const vec4_t text_col = flow_vec4(style.text_color);
    const vec4_t text_sec = flow_vec4(style.text_secondary);

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

    const size_t num_nodes   = md_array_size(layout->nodes);
    const size_t num_ribbons = md_array_size(layout->ribbons);
    const uint32_t num_columns = graph->num_columns;

    // ---- the nodes on screen, their sides, their labels ----
    FlowLabel* labels = md_temp_alloc_zero_array(temp, FlowLabel, num_nodes ? num_nodes : 1);
    float label_w[3] = { 0, 0, 0 };    // the widest label placed, per side

    // How far out each outer side reaches (rows at their levels, bands): the labels stand beyond, in one line, after
    // the arrows of the collapsed groups where a side has any
    float outer[3] = { 1.0e30f, 0.0f, -1.0e30f };
    for (size_t i = 0; i < num_nodes; ++i) {
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, layout->nodes[i].cut_idx);
        if (!node) continue;
        const int side = flow_column_side(node->column, num_columns);
        if (side < 0) outer[0] = ImMin(outer[0], flow_view_to_screen(view, layout->nodes[i].min).x);
        if (side > 0) outer[2] = ImMax(outer[2], flow_view_to_screen(view, layout->nodes[i].max).x);
    }
    for (size_t i = 0; i < md_array_size(layout->bands); ++i) {
        const int side = flow_column_side(layout->bands[i].column, num_columns);
        if (side < 0) outer[0] = ImMin(outer[0], flow_view_to_screen(view, layout->bands[i].min).x);
        if (side > 0) outer[2] = ImMax(outer[2], flow_view_to_screen(view, layout->bands[i].max).x);
    }
    const size_t num_graph_nodes = md_array_size(graph->nodes);
    bool* container = md_temp_alloc_array(temp, bool, num_graph_nodes ? num_graph_nodes : 1);
    flow_containers(container, graph);
    bool arrows[3];
    flow_sides_with_arrows(arrows, graph, cut, layout, container);
    const float arrow_box = flow_arrow_box(cv);
    auto label_out = [&](int s) { return style.label_gap_px + ((arrows[s] && cv.fill_tri) ? arrow_box + FLOW_ARROW_PAD : 0.0f); };

    // A collapsed group's rows, counted for its label
    uint32_t* num_children = md_temp_alloc_array(temp, uint32_t, num_graph_nodes ? num_graph_nodes : 1);
    flow_child_counts(num_children, graph);
    {
        int32_t prev_col = -1;
        float prev_bottom = -1.0e30f;
        for (size_t i = 0; i < num_nodes; ++i) {
            const md_flow_layout_node_t* ln = layout->nodes + i;
            const md_flow_node_t* node = md_flow_cut_node(cut, graph, ln->cut_idx);
            if (!node) continue;
            const ImVec2 min = flow_view_to_screen(view, ln->min);
            const ImVec2 max = flow_view_to_screen(view, ln->max);
            const float h = max.y - min.y;
            FlowLabel* lb = labels + i;
            lb->side = flow_column_side(node->column, num_columns);
            char count[FLOW_LABEL_CAP] = "";
            if (lb->side != 0 && ln->cut_idx < num_graph_nodes && container[ln->cut_idx]) {
                flow_member_count(count, sizeof(count), style, num_children[ln->cut_idx]);
            }
            flow_label_lines(lb, cv, node, h, count);
            const float centre = (min.y + max.y) * 0.5f;
            if (lb->side == 0) {
                // on the bar: only when it fits its height
                if (h < flow_line_h(lb->size[0]) * 0.9f) continue;
                lb->anchor_x = (min.x + max.x) * 0.5f;
                lb->top = centre - lb->height * 0.5f;
                lb->visible = true;
                continue;
            }
            if ((int32_t)node->column != prev_col) { prev_col = (int32_t)node->column; prev_bottom = -1.0e30f; }
            float top = MD_MAX(centre - lb->height * 0.5f, prev_bottom + 1.0f);
            // pushed off its own node: it would name the wrong row
            if (top > max.y + 1.0f) continue;
            lb->top = top;
            const int s = flow_side_index(lb->side);
            lb->anchor_x = lb->side < 0 ? outer[0] - label_out(s) : outer[2] + label_out(s);
            lb->visible = true;
            prev_bottom = top + lb->height;
            label_w[flow_side_index(lb->side)] = MD_MAX(label_w[flow_side_index(lb->side)], lb->width);
        }
    }

    // ---- ribbons, under everything ----
    struct RibbonLabel { ImVec2 mid; float thick; uint32_t link; };
    RibbonLabel* rlabels = md_temp_alloc_array(temp, RibbonLabel, num_ribbons ? num_ribbons : 1);
    size_t num_rlabels = 0;

    for (size_t i = 0; i < num_ribbons; ++i) {
        const md_flow_layout_ribbon_t* r = layout->ribbons + i;
        if (r->flags & MD_FLOW_RIBBON_FLAG_SELF) continue;     // drawn with the nodes
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
        if (ImAbs(p1.y - p0.y) < style.ribbon_min_px && ImAbs(q1.y - q0.y) < style.ribbon_min_px) continue;

        const bool direct = (int32_t)r->link_idx == emphasis.hover_link
            || (int32_t)link->src == emphasis.hover_node || (int32_t)link->dst == emphasis.hover_node
            || flow_is_emphasized(emphasis, link->src) || flow_is_emphasized(emphasis, link->dst);
        // pale ribbons need a lift to read as the ones in focus, not only the others' fading
        const float alpha = style.ribbon_alpha * (!has_highlight ? 1.0f : (direct ? 1.6f : style.dim_factor));

        vec4_t c_src, c_dst;
        if (style.color_by_side) {
            c_src = c_dst = flow_with_alpha(style.ribbon_color, alpha);
        } else {
            c_src = flow_with_alpha(src->color, alpha);
            c_dst = flow_with_alpha(dst->color, alpha);
        }

        const uint32_t num_stripes = flow_stripe_count(stripes, r->link_idx);
        if (num_stripes == 0) {
            cv.ribbon(cv.user, p0, p1, q0, q1, c_src, c_dst, FLOW_FRINGE_TOP | FLOW_FRINGE_BOTTOM);
            if (mouse_ptr && hit.link < 0 && flow_ribbon_contains(p0, p1, q0, q1, style.ribbon_curvature, style.ribbon_segments, mouse)) {
                hit.link = (int32_t)r->link_idx;
            }
        } else {
            const float* frac = stripes->frac + stripes->offset[r->link_idx];
            float t0 = 0.0f;
            for (uint32_t j = 0; j < num_stripes; ++j) {
                const float t1 = (j + 1 == num_stripes) ? 1.0f : MD_MIN(1.0f, t0 + frac[j]);
                const ImVec2 a0 = flow_lerp(p0, p1, t0), a1 = flow_lerp(p0, p1, t1);
                const ImVec2 b0 = flow_lerp(q0, q1, t0), b1 = flow_lerp(q0, q1, t1);
                const int fringe = (j == 0 ? FLOW_FRINGE_TOP : 0) | (j + 1 == num_stripes ? FLOW_FRINGE_BOTTOM : 0);
                cv.ribbon(cv.user, a0, a1, b0, b1, flow_stripe_tone(c_src, j), flow_stripe_tone(c_dst, j), fringe);
                if (mouse_ptr && hit.link < 0 && flow_ribbon_contains(a0, a1, b0, b1, style.ribbon_curvature, style.ribbon_segments, mouse)) {
                    hit.link = (int32_t)r->link_idx;
                    hit.stripe = (int32_t)j;
                }
                t0 = t1;
            }
        }

        if (style.show_ribbon_labels) {
            const ImVec2 top = flow_edge_point(p0, q0, style.ribbon_curvature, 0.5f);
            const ImVec2 bot = flow_edge_point(p1, q1, style.ribbon_curvature, 0.5f);
            const float thick = ImAbs(bot.y - top.y);
            if (thick >= style.ribbon_label_min_px) {
                rlabels[num_rlabels++] = { ImVec2((top.x + bot.x) * 0.5f, (top.y + bot.y) * 0.5f), thick, r->link_idx };
            }
        }
    }

    // ---- group bands: the molecule (or set, or user group) a column's rows belong to, as a bar of its own against
    // theirs and spanning them, its name along it. Inside the labels, so the hierarchy reads as nesting: a molecule
    // holds its fragments the way a fragment holds its share. ----
    const size_t num_bands = md_array_size(layout->bands);
    for (size_t i = 0; i < num_bands; ++i) {
        const md_flow_layout_band_t* b = layout->bands + i;
        const md_flow_node_t* node = graph->nodes + b->node;
        const int side = flow_column_side(b->column, num_columns);
        if (side == 0) continue;     // a middle column has no free side

        // at its level (flow_layout_indent), as wide as a bar
        ImVec2 band_min = flow_view_to_screen(view, b->min);
        ImVec2 band_max = flow_view_to_screen(view, b->max);
        const float band_w = band_max.x - band_min.x;
        // never shorter than its arrow: a group of one small row is still a thing to close
        if (cv.fill_tri && band_max.y - band_min.y < band_w) {
            const float c = (band_min.y + band_max.y) * 0.5f;
            band_min.y = c - band_w * 0.5f;
            band_max.y = c + band_w * 0.5f;
        }

        const bool hovered = (int32_t)b->node == emphasis.hover_node;
        const bool touched = !has_highlight || hovered || flow_is_emphasized(emphasis, b->node);
        const float alpha = touched ? 1.0f : style.dim_factor;
        // the side's colour, as the bars, a little paler: a container, not a share of its own
        const vec4_t base = flow_node_base_color(style, node, side);
        cv.fill_rect(cv.user, band_min, band_max, flow_with_alpha(flow_toward_white(base, style.band_fill_tint), alpha), style.node_rounding);
        const vec4_t border = style.color_by_side ? flow_scale_rgb(base, 0.78f) : flow_visible_border(base);
        cv.stroke_rect(cv.user, band_min, band_max, flow_with_alpha(border, alpha), style.node_rounding,
                       hovered ? style.border_width_hover : style.border_width);

        if (mouse.x >= band_min.x && mouse.x <= band_max.x && mouse.y >= band_min.y && mouse.y <= band_max.y) {
            hit.node = (int32_t)b->node;
            hit.link = -1;
            hit.stripe = -1;
        }

        // Its arrow, at the top, pointing down along it: open. One click on it closes the group, on both sides.
        float label_top = band_min.y;
        if (cv.fill_tri) {
            const ImVec2 a0 = band_min, a1 = { band_max.x, band_min.y + band_w };
            const bool arrow_hovered = mouse.x >= a0.x && mouse.x <= a1.x && mouse.y >= a0.y && mouse.y <= a1.y;
            if (arrow_hovered) hit.toggle = (int32_t)b->node;
            const ImVec2 c = (a0 + a1) * 0.5f;
            const float r = MD_MIN(band_w, 16.0f) * (arrow_hovered ? 0.30f : 0.25f);
            const vec4_t col = flow_with_alpha(arrow_hovered ? text_col : flow_scale_rgb(border, 0.85f), alpha);
            cv.fill_tri(cv.user, ImVec2(c.x - r, c.y - r * 0.6f), ImVec2(c.x + r, c.y - r * 0.6f), ImVec2(c.x, c.y + r * 0.9f), col);
            label_top = a1.y;
        }

        // its name along it, where the band is long enough to hold it (else the tooltip has it)
        char buf[FLOW_LABEL_CAP];
        flow_node_name(buf, sizeof(buf), node);
        if (!flow_band_label_fits(cv, style, buf, band_w, band_max.y - label_top)) continue;
        const ImVec2 c = { (band_min.x + band_max.x) * 0.5f, (label_top + band_max.y) * 0.5f };
        cv.text_up(cv.user, c, flow_band_text_size(cv, band_w), flow_with_alpha(text_col, alpha), buf);
    }

    // ---- nodes ----
    for (size_t i = 0; i < num_nodes; ++i) {
        const md_flow_layout_node_t* ln = layout->nodes + i;
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, ln->cut_idx);
        if (!node) continue;
        const ImVec2 min = flow_view_to_screen(view, ln->min);
        const ImVec2 max = flow_view_to_screen(view, ln->max);
        const int side = labels[i].side;

        const uint8_t node_state = flow_state_of(emphasis, ln->cut_idx);
        const bool selected = (node_state & FlowNodeState_Selected) != 0;
        const bool partial  = (node_state & FlowNodeState_Partial) != 0;
        const bool hovered  = ((node_state & FlowNodeState_Highlight) != 0) || (int32_t)ln->cut_idx == emphasis.hover_node;
        const bool touched  = !has_highlight || selected || partial || hovered;
        const float alpha = touched ? 1.0f : style.dim_factor;

        const vec4_t base = flow_node_base_color(style, node, side);
        cv.fill_rect(cv.user, min, max, flow_toward_white(base, style.node_fill_tint * (touched ? 1.0f : 0.5f)), style.node_rounding);

        // The node's self flow, solid over its stub: the part that stays put
        for (size_t j = 0; j < num_ribbons; ++j) {
            const md_flow_layout_ribbon_t* r = layout->ribbons + j;
            if (!(r->flags & MD_FLOW_RIBBON_FLAG_SELF)) continue;
            const md_flow_link_t* link = cut->links + r->link_idx;
            vec2_t a, b;
            if      (link->src == ln->cut_idx) { a = r->p0; b = r->p1; }
            else if (link->dst == ln->cut_idx) { a = r->q0; b = r->q1; }
            else continue;
            const float y0 = flow_view_to_screen(view, a).y;
            const float y1 = flow_view_to_screen(view, b).y;
            cv.fill_rect(cv.user, { min.x, ImMin(y0, y1) }, { max.x, ImMax(y0, y1) }, flow_with_alpha(base, style.self_alpha * alpha), 0.0f);
        }

        // State is the border's weight, in the node's own colour: fully selected heaviest, then partial,
        // then hovered - one channel, ordered by how much the user has committed to the thing.
        const float width = selected ? style.border_width_selected
                          : partial  ? style.border_width_partial
                          : hovered  ? style.border_width_hover
                                     : style.border_width;
        const vec4_t border = style.color_by_side ? flow_scale_rgb(base, 0.78f) : flow_visible_border(base);
        cv.stroke_rect(cv.user, min, max, flow_with_alpha(border, alpha), style.node_rounding, width);

        if (mouse.x >= min.x && mouse.x <= max.x && mouse.y >= min.y && mouse.y <= max.y) {
            hit.node = (int32_t)ln->cut_idx;
            hit.link = -1;  // a node under the cursor wins over the ribbons behind it
            hit.stripe = -1;
        }

        // ---- a collapsed group's arrow: between the bars and its label, pointing the way it opens (in, towards the
        // ribbons). One click on it opens the group, on both sides. ----
        const FlowLabel* lb = labels + i;
        if (cv.fill_tri && side != 0 && ln->cut_idx < num_graph_nodes && container[ln->cut_idx]) {
            const int s = flow_side_index(side);
            const float cx = side < 0 ? outer[0] - style.label_gap_px * 0.5f - arrow_box * 0.5f
                                      : outer[2] + style.label_gap_px * 0.5f + arrow_box * 0.5f;
            const float cy = lb->visible ? lb->top + flow_line_h(lb->size[0]) * 0.5f : (min.y + max.y) * 0.5f;
            const ImVec2 a0 = { cx - arrow_box * 0.5f, cy - arrow_box * 0.5f }, a1 = { cx + arrow_box * 0.5f, cy + arrow_box * 0.5f };
            const bool arrow_hovered = mouse.x >= a0.x && mouse.x <= a1.x && mouse.y >= a0.y && mouse.y <= a1.y;
            if (arrow_hovered) {
                hit.toggle = (int32_t)ln->cut_idx;
                hit.node = (int32_t)ln->cut_idx;
                hit.link = -1;
                hit.stripe = -1;
            }
            const float r = arrow_box * (arrow_hovered ? 0.42f : 0.34f);
            const float d = (float)-side;   // the way in: right on the left side, left on the right
            const vec4_t col = flow_with_alpha(arrow_hovered ? text_col : flow_scale_rgb(border, 0.85f), alpha);
            cv.fill_tri(cv.user, ImVec2(cx - d * r * 0.6f, cy - r), ImVec2(cx + d * r * 0.9f, cy), ImVec2(cx - d * r * 0.6f, cy + r), col);
            (void)s;
        }

        // ---- its label ----
        if (!lb->visible) continue;
        const vec4_t name_col = (!style.color_by_side && style.color_labels) ? flow_readable_text(node->color) : text_col;
        if (lb->side == 0) {
            // on the bar, over a light plate so the ribbons behind do not bleed into it
            const float pad = 4.0f;
            const ImVec2 pmin = { lb->anchor_x - lb->width * 0.5f - pad, lb->top - 1.0f };
            const ImVec2 pmax = { lb->anchor_x + lb->width * 0.5f + pad, lb->top + lb->height + 1.0f };
            cv.fill_rect(cv.user, pmin, pmax, vec4_t{ 1.0f, 1.0f, 1.0f, 0.72f * alpha }, 3.0f);
        }
        const int align = lb->side < 0 ? +1 : (lb->side > 0 ? -1 : 0);
        float y = lb->top;
        for (int k = 0; k < lb->num_lines; ++k) {
            cv.text(cv.user, ImVec2(lb->anchor_x, y), lb->size[k], flow_with_alpha(lb->secondary[k] ? text_sec : name_col, alpha), align, lb->line[k]);
            y += flow_line_h(lb->size[k]);
        }
    }

    // ---- ribbon labels: the share of the big transfers, at their middle ----
    {
        ImVec2* placed_min = md_temp_alloc_array(temp, ImVec2, num_rlabels ? num_rlabels : 1);
        ImVec2* placed_max = md_temp_alloc_array(temp, ImVec2, num_rlabels ? num_rlabels : 1);
        size_t num_placed = 0;
        for (size_t i = 0; i < num_rlabels; ++i) {
            const RibbonLabel& rl = rlabels[i];
            const md_flow_link_t* link = cut->links + rl.link;
            char main[32], sub[64] = "";
            snprintf(main, sizeof(main), "%.1f%%", link->weight * 100.0f);
            if (rl.thick >= style.ribbon_sublabel_min_px && stripes && stripes->offset && rl.link < stripes->num_links) {
                const uint32_t n = stripes->offset[rl.link + 1] - stripes->offset[rl.link];
                if (n == 1)     snprintf(sub, sizeof(sub), "1 block NTO");
                else if (n > 1) snprintf(sub, sizeof(sub), "%u block NTOs, largest %.0f%%", n, 100.0f * stripes->frac[stripes->offset[rl.link]]);
            }
            const float wm = cv.measure(cv.user, main, cv.base_size);
            const float ws = sub[0] ? cv.measure(cv.user, sub, cv.small_size) : 0.0f;
            const float h = flow_line_h(cv.base_size) + (sub[0] ? flow_line_h(cv.small_size) : 0.0f);
            const float w = MD_MAX(wm, ws);
            const ImVec2 bmin = { rl.mid.x - w * 0.5f, rl.mid.y - h * 0.5f };
            const ImVec2 bmax = { rl.mid.x + w * 0.5f, rl.mid.y + h * 0.5f };
            bool clash = false;
            for (size_t k = 0; k < num_placed && !clash; ++k) {
                clash = bmin.x < placed_max[k].x && bmax.x > placed_min[k].x && bmin.y < placed_max[k].y && bmax.y > placed_min[k].y;
            }
            if (clash) continue;
            placed_min[num_placed] = bmin;
            placed_max[num_placed] = bmax;
            ++num_placed;
            cv.text(cv.user, ImVec2(rl.mid.x, bmin.y), cv.base_size, text_col, 0, main);
            if (sub[0]) cv.text(cv.user, ImVec2(rl.mid.x, bmin.y + flow_line_h(cv.base_size)), cv.small_size, text_sec, 0, sub);
        }
    }

    // ---- annotations: captions over the columns, the subtitle, the legend under it all ----
    if (notes && num_nodes > 0) {
        float col_x0[3] = { 0, 0, 0 }, col_x1[3] = { 0, 0, 0 };
        bool col_found[3] = { false, false, false };
        float top = 1.0e30f, bottom = -1.0e30f, left = 1.0e30f;
        for (size_t i = 0; i < num_nodes; ++i) {
            const md_flow_layout_node_t* ln = layout->nodes + i;
            const md_flow_node_t* node = md_flow_cut_node(cut, graph, ln->cut_idx);
            if (!node) continue;
            const ImVec2 min = flow_view_to_screen(view, ln->min);
            const ImVec2 max = flow_view_to_screen(view, ln->max);
            const int s = flow_side_index(labels[i].side);
            col_x0[s] = min.x; col_x1[s] = max.x; col_found[s] = true;
            top = MD_MIN(top, min.y);
            bottom = MD_MAX(bottom, max.y);
            left = MD_MIN(left, min.x);
        }
        const float cap_y = top - 10.0f - cv.base_size;
        for (int s = 0; s < 3; ++s) {
            if (!col_found[s] || !notes->caption[s]) continue;
            cv.text(cv.user, ImVec2((col_x0[s] + col_x1[s]) * 0.5f, cap_y), cv.base_size, text_col, 0, notes->caption[s]);
        }
        if (notes->subtitle && col_found[0] && col_found[2]) {
            const float w = cv.measure(cv.user, notes->subtitle, cv.small_size);
            const float cx = (col_x1[0] + col_x0[2]) * 0.5f;
            // a middle column's caption takes the centre of the captions' line: the subtitle goes above it
            const float y = col_found[1] ? cap_y - flow_line_h(cv.base_size) - 2.0f : cap_y + (cv.base_size - cv.small_size) * 0.5f;
            const ImVec2 pmin = { cx - w * 0.5f, y };
            cv.text(cv.user, ImVec2(cx, y), cv.small_size, text_sec, 0, notes->subtitle);
            if (mouse.x >= pmin.x && mouse.x <= pmin.x + w && mouse.y >= y && mouse.y <= y + flow_line_h(cv.small_size)) hit.subtitle = true;
        }
        if (notes->legend_ribbon || notes->legend_local) {
            const float y = bottom + 16.0f;
            const float sw = 12.0f, sh = 10.0f;
            const float ty = y + sh * 0.5f - cv.small_size * 0.5f;
            float x = left;
            if (notes->legend_ribbon) {
                cv.fill_rect(cv.user, { x, y }, { x + sw, y + sh }, flow_with_alpha(style.color_by_side ? style.ribbon_color : text_sec, style.ribbon_alpha), 0.0f);
                cv.text(cv.user, ImVec2(x + sw + 6.0f, ty), cv.small_size, text_sec, -1, notes->legend_ribbon);
                x += sw + 6.0f + cv.measure(cv.user, notes->legend_ribbon, cv.small_size) + 24.0f;
            }
            if (notes->legend_local) {
                const vec4_t solid = style.color_by_side ? style.side_color[2] : text_sec;
                cv.fill_rect(cv.user, { x, y }, { x + sw, y + sh }, flow_with_alpha(solid, style.self_alpha), 0.0f);
                cv.text(cv.user, ImVec2(x + sw + 6.0f, ty), cv.small_size, text_sec, -1, notes->legend_local);
            }
        }
    }

    return hit;
}

FlowHit flow_draw_imgui(ImDrawList* draw_list, const FlowView& view, const md_flow_graph_t* graph,
                        const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                        const FlowDrawStyle& style, const FlowEmphasis& emphasis,
                        const FlowStripes* stripes, const FlowAnnotations* notes)
{
    if (!draw_list) return FlowHit{};
    FlowImCanvas im = { draw_list, style.ribbon_curvature, style.ribbon_segments };
    const float base = ImGui::GetFontSize();
    const FlowCanvas cv = { &im, base, base * style.small_text_scale, im_measure, im_fill, im_stroke, im_text, im_ribbon, im_text_up, im_tri };
    const ImVec2 mouse = ImGui::GetMousePos();
    return flow_draw_common(cv, view, graph, cut, layout, style, emphasis, stripes, notes, &mouse);
}

void flow_draw_vg(md_vg_scene_t* scene, const FlowView& view, const md_flow_graph_t* graph,
                  const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                  const FlowDrawStyle& style, const FlowStripes* stripes, const FlowAnnotations* notes)
{
    if (!scene) return;
    FlowVgCanvas vg = { scene, style.ribbon_curvature };
    const float base = FLOW_VG_TEXT_SIZE;
    const FlowCanvas cv = { &vg, base, base * style.small_text_scale, vg_measure, vg_fill, vg_stroke, vg_text, vg_ribbon, vg_text_up, nullptr };
    flow_draw_common(cv, view, graph, cut, layout, style, FlowEmphasis{}, stripes, notes, nullptr);
}

static FlowMargins flow_measure_common(const FlowCanvas& cv, const md_flow_graph_t* graph, const md_flow_cut_t* cut,
                                       const md_flow_layout_t* layout, const FlowDrawStyle& style, const FlowAnnotations* notes)
{
    FlowMargins m = {};
    if (!graph || !cut || !layout) return m;

    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };
    const size_t num_graph_nodes = md_array_size(graph->nodes);
    bool*     container    = md_temp_alloc_array(temp, bool, num_graph_nodes ? num_graph_nodes : 1);
    uint32_t* num_children = md_temp_alloc_array(temp, uint32_t, num_graph_nodes ? num_graph_nodes : 1);
    flow_containers(container, graph);
    flow_child_counts(num_children, graph);

    float label_w[3] = { 0, 0, 0 };
    for (size_t i = 0; i < md_array_size(layout->nodes); ++i) {
        const uint32_t c = layout->nodes[i].cut_idx;
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, c);
        if (!node) continue;
        const int side = flow_column_side(node->column, graph->num_columns);
        char count[FLOW_LABEL_CAP] = "";
        if (side != 0 && c < num_graph_nodes && container[c]) flow_member_count(count, sizeof(count), style, num_children[c]);
        const int s = flow_side_index(side);
        label_w[s] = MD_MAX(label_w[s], flow_label_width_any(cv, node, count));
    }
    // the levels are in the layout (flow_layout_indent); outside them the arrows and the labels
    bool arrows[3];
    flow_sides_with_arrows(arrows, graph, cut, layout, container);
    for (int s = 0; s < 3; s += 2) {
        float w = style.label_gap_px + label_w[s];
        if (arrows[s] && cv.fill_tri) w += flow_arrow_box(cv) + FLOW_ARROW_PAD;
        (s == 0 ? m.left : m.right) = w;
    }
    if (notes) {
        bool middle = false;
        for (size_t i = 0; i < md_array_size(layout->nodes) && !middle; ++i) {
            const md_flow_node_t* node = md_flow_cut_node(cut, graph, layout->nodes[i].cut_idx);
            middle = node && flow_column_side(node->column, graph->num_columns) == 0;
        }
        const bool any_caption = notes->caption[0] || notes->caption[1] || notes->caption[2];
        m.top = any_caption ? 10.0f + flow_line_h(cv.base_size) : 0.0f;
        if (notes->subtitle && middle) m.top += flow_line_h(cv.base_size) + 2.0f;
        else if (notes->subtitle && !any_caption) m.top = 10.0f + flow_line_h(cv.base_size);
        if (notes->legend_ribbon || notes->legend_local) m.bottom = 16.0f + flow_line_h(cv.small_size);
    }
    return m;
}

FlowMargins flow_measure_margins(const md_flow_graph_t* graph, const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                                 const FlowDrawStyle& style, const FlowAnnotations* notes)
{
    const float base = ImGui::GetFontSize();
    // fill_tri only to say the screen has arrows: nothing is drawn
    const FlowCanvas cv = { nullptr, base, base * style.small_text_scale, im_measure, nullptr, nullptr, nullptr, nullptr, nullptr, im_tri };
    return flow_measure_common(cv, graph, cut, layout, style, notes);
}

FlowMargins flow_measure_margins_vg(const md_flow_graph_t* graph, const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                                    const FlowDrawStyle& style, const FlowAnnotations* notes)
{
    const float base = FLOW_VG_TEXT_SIZE;
    const FlowCanvas cv = { nullptr, base, base * style.small_text_scale, vg_measure, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };
    return flow_measure_common(cv, graph, cut, layout, style, notes);
}
