#pragma once

// Drawing backends for md_flow layouts.
//
// The layout (md_flow_layout_compute) is a pure function producing geometry in graph space.
// Everything here is a thin consumer of that geometry: one drawing routine, behind two backends - the
// screen (ImGui draw list) and export (md_vg_scene_t) - that only supply rectangles, ribbons and text.
// Neither computes a position of its own, which is why they cannot drift apart the way the hand-rolled
// im_sankey_diagram / vg_sankey_diagram pair did.

#include <md_flow.h>

#include <imgui.h>

struct md_vg_scene_t;

// Graph space -> screen. Pan and zoom live here rather than in the layout, so changing the view
// never triggers a relayout.
struct FlowView {
    ImVec2 area_min = {0, 0};
    ImVec2 area_max = {0, 0};
    ImVec2 pan      = {0, 0};   // in graph units
    float  zoom     = 1.0f;
};

struct FlowDrawStyle {
    // Colour by side (the default): the first column in the hole colour, the last in the particle colour,
    // a middle column in the orbital colour, and every ribbon in one neutral colour. What a node IS is said
    // by its label and by the mirrored order of the columns; colour says only which side it is on, so the
    // picture has two hues and the ribbons read as one material. Off: each node in its own colour
    // (md_flow_node_t.color), ribbons blending from source to destination.
    bool   color_by_side   = true;
    vec4_t side_color[3]   = { { 0.847f, 0.353f, 0.188f, 1.0f },     // first column: coral
                               { 0.498f, 0.467f, 0.867f, 1.0f },     // middle: the ribbons' violet
                               { 0.114f, 0.620f, 0.459f, 1.0f } };   // last column: teal
    vec4_t ribbon_color    = { 0.498f, 0.467f, 0.867f, 1.0f };

    float ribbon_alpha     = 0.32f;
    // A self flow (MD_FLOW_RIBBON_FLAG_SELF) is not a ribbon: it is the part of a node that stays where
    // it is - local excitation in the charge transfer diagram - drawn solid over its stub inside the node.
    float self_alpha       = 0.9f;
    int   ribbon_segments  = 24;
    float ribbon_curvature = 1.0f;   // 0 = straight, 1 = fully S-shaped
    float node_rounding    = 2.0f;

    // Below this on-screen thickness a ribbon carries no readable quantity and only adds hair to
    // the picture. Purely a draw-time filter - the model still holds the flow, and the tooltip
    // still reports it, so nothing silently disappears from the arithmetic.
    float ribbon_min_px    = 1.5f;

    // Labels sit OUTSIDE the outer columns - left of the first, right of the last - name over share, so
    // the bars can be thin and the space between the columns belongs to the ribbons. A node too short
    // for two lines gets one small line; a label that would collide with the one above is dropped
    // (the tooltip still has it). A middle column's labels sit on its bars.
    float label_gap_px     = 8.0f;
    float small_text_scale = 0.86f;
    bool  show_ribbon_labels      = true;
    float ribbon_label_min_px     = 22.0f;  // a ribbon this thick at its middle carries its share
    float ribbon_sublabel_min_px  = 44.0f;  // and this thick, what it is made of

    // Group bands: an expanded molecule (set, user group) as a bar of its own against the bars of its rows, spanning
    // them, its name along it; nested ones outward, the labels outside them all.
    float band_offset_px    = 3.0f;   // gap between the bars and the band, and between nested bands
    float band_width_px     = 15.0f;
    float band_fill_tint    = 0.12f;  // of the side's colour over white: paler than the bars
    bool  show_band_labels  = true;

    // A node's interior is this much of its colour over white, its border its colour darkened.
    float node_fill_tint    = 0.22f;
    float border_width      = 1.0f;
    float border_width_hover    = 1.75f;
    float border_width_selected = 2.5f;
    float border_width_partial  = 1.75f;
    bool  color_labels      = false;  // own colours only: the name in the node's colour

    ImU32 text_color        = IM_COL32(28, 28, 30, 255);
    ImU32 text_secondary    = IM_COL32(28, 28, 30, 150);
    ImU32 node_border_color = IM_COL32(0, 0, 0, 90);   // fallback for nodes with no identity

    // Applied to anything outside the emphasis. Deliberately close to 1: the point is to let the
    // selection read as foreground, not to hide the rest of the diagram.
    float dim_factor        = 0.55f;
};

// What the diagram's columns mean, which the backend cannot know: captions over them, a subtitle, and
// a legend under the diagram. Null entries are not drawn.
struct FlowAnnotations {
    const char* caption[3]    = {};       // over the first column, a middle one, the last
    const char* subtitle      = nullptr;  // centred between the outer columns, on the captions' line
    const char* legend_ribbon = nullptr;  // legend: what a ribbon is
    const char* legend_local  = nullptr;  // legend: what the solid part of a node is
};

// Room, in pixels outside the box of the nodes, that the labels, bands and annotations take: what a fit
// of the view leaves free.
struct FlowMargins {
    float left = 0, right = 0, top = 0, bottom = 0;
};

// Per-node state, as bit flags. A node can be several of these at once - selected AND hovered, or
// partially selected AND highlighted - so they are flags rather than an enum.
enum FlowNodeState : uint8_t {
    FlowNodeState_None      = 0,
    FlowNodeState_Selected  = 1 << 0,  // every atom this node stands for is selected
    FlowNodeState_Partial   = 1 << 1,  // some but not all of them are
    FlowNodeState_Highlight = 1 << 2,  // in the application's highlight (hover) mask
};

// What should read as foreground. An empty emphasis draws everything at full strength.
struct FlowEmphasis {
    // Per cut index. May be null.
    const uint8_t* node_state = nullptr;
    size_t         num_state_entries = 0;
    bool           any_selected = false;

    int32_t hover_node = -1;
    int32_t hover_link = -1;
};

struct FlowHit {
    int32_t node = -1;  // cut index space
    int32_t link = -1;  // index into cut->links
    int32_t stripe = -1;    // the stripe of that link under the cursor (FlowStripes), -1 for none
    bool subtitle = false;  // the cursor is over the annotations' subtitle
};

// What a ribbon is made of, drawn as stripes along it: for cut link l, the shares frac[offset[l] ..
// offset[l + 1]) of its weight, largest first, summing to 1. In the charge transfer diagram these are the
// orbital pairs (block NTOs) that carry one fragment-to-fragment transfer. A ribbon is striped only when its
// second share reaches 'min_second', so a transfer carried by one pair looks as it always did.
struct FlowStripes {
    const uint32_t* offset = nullptr;   // [num_links + 1]
    const float*    frac   = nullptr;
    size_t          num_links = 0;
    float           min_second = 0.05f;
};

ImVec2 flow_view_to_screen(const FlowView& view, vec2_t graph_point);

// Draws the diagram and returns what the mouse is over.
FlowHit flow_draw_imgui(ImDrawList* draw_list, const FlowView& view, const md_flow_graph_t* graph,
                        const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                        const FlowDrawStyle& style, const FlowEmphasis& emphasis,
                        const FlowStripes* stripes = nullptr, const FlowAnnotations* notes = nullptr);

// The same picture into a vector-graphics scene, for SVG export: one drawing routine behind both, so the
// export is the screen by construction. Text widths are estimated (md_vg has no font metrics).
void flow_draw_vg(md_vg_scene_t* scene, const FlowView& view, const md_flow_graph_t* graph,
                  const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                  const FlowDrawStyle& style, const FlowStripes* stripes = nullptr,
                  const FlowAnnotations* notes = nullptr);

// What flow_draw_imgui will need outside the nodes at the current font, for any zoom
FlowMargins flow_measure_margins(const md_flow_graph_t* graph, const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                                 const FlowDrawStyle& style, const FlowAnnotations* notes);
// The same for flow_draw_vg's text
FlowMargins flow_measure_margins_vg(const md_flow_graph_t* graph, const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                                    const FlowDrawStyle& style, const FlowAnnotations* notes);
