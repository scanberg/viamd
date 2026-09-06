#pragma once

// Drawing backends for md_flow layouts.
//
// The layout (md_flow_layout_compute) is a pure function producing geometry in graph space.
// Everything here is a thin consumer of that geometry: one backend for the screen (ImGui draw
// list) and one for export (md_vg_scene_t). Neither computes a position of its own - that is the
// whole point, and it is why these two cannot drift apart the way the hand-rolled
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
    float ribbon_alpha     = 0.55f;
    int   ribbon_segments   = 24;
    float ribbon_curvature  = 0.5f;   // 0 = straight, 1 = fully S-shaped
    float node_rounding     = 3.0f;

    // Below this on-screen thickness a ribbon carries no readable quantity and only adds hair to
    // the picture. Purely a draw-time filter - the model still holds the flow, and the tooltip
    // still reports it, so nothing silently disappears from the arithmetic.
    float ribbon_min_px     = 1.5f;

    // Semantic level of detail, in on-screen pixels of node height. A node too short to carry its
    // own label does not get one; a diagram of clipped text is worse than a diagram of bars.
    float label_min_px      = 11.0f;
    float percent_min_px    = 26.0f;
    float label_pad_px      = 5.0f;

    // Group bands, drawn outboard of the column they bracket. A band is the group as an
    // organizing layer AROUND its atoms, rather than an alternative to showing them.
    float band_offset_px    = 10.0f;  // gap between the column and the band
    float band_width_px     = 7.0f;
    float band_label_gap_px = 4.0f;
    bool  show_band_labels  = true;

    // Nodes are outlined, not filled. Colour identifies a node through its BORDER and its label;
    // the interior stays near-white so the ribbons crossing behind it are the only saturated thing
    // on screen. A tint rather than pure white keeps the box tied to its colour at a glance.
    float node_fill_tint    = 0.12f;   // 0 = white, 1 = the node's own colour
    float border_width      = 1.5f;
    float border_width_hover    = 2.5f;
    float border_width_selected = 3.0f;
    float border_width_partial  = 2.0f;
    bool  color_labels      = true;

    ImU32 text_color        = IM_COL32(0, 0, 0, 255);
    ImU32 node_border_color = IM_COL32(0, 0, 0, 90);   // fallback for nodes with no identity

    // Applied to anything outside the emphasis. Deliberately close to 1: the point is to let the
    // selection read as foreground, not to hide the rest of the diagram.
    float dim_factor        = 0.55f;
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
};

ImVec2 flow_view_to_screen(const FlowView& view, vec2_t graph_point);

// Draws the diagram and returns what the mouse is over.
FlowHit flow_draw_imgui(ImDrawList* draw_list, const FlowView& view, const md_flow_graph_t* graph,
                        const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                        const FlowDrawStyle& style, const FlowEmphasis& emphasis);

// Same geometry into a vector-graphics scene, for SVG export.
void flow_draw_vg(md_vg_scene_t* scene, const FlowView& view, const md_flow_graph_t* graph,
                  const md_flow_cut_t* cut, const md_flow_layout_t* layout,
                  const FlowDrawStyle& style);
