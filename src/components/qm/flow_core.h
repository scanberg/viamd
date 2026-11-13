#pragma once

// Layered flow graphs: the model behind VIAMD's charge-transfer / transition diagram.
//
// A flow graph is a layered DAG. Nodes live in columns; links only ever connect column c to
// column c+1. Within a column, nodes form a hierarchy (group -> atom -> orbital), and a node's
// children partition its weight exactly.
//
// The point of the split between flow_graph_t and flow_cut_t: the GRAPH is built once,
// at the finest level available, from whatever produced it. The CUT is what the user is looking
// at right now - which nodes they expanded, what they thresholded away. Every coarser view is
// derived from the graph by aggregation, never by asking the source for different numbers. That
// is what makes it impossible for the collapsed view and the expanded view to disagree.
//
// Layout is computed in GRAPH SPACE (x = column, y in [0,1] by cumulative weight). Pan and zoom
// belong to whoever draws it, so changing the zoom never triggers a relayout.

#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_str.h>
#include <core/md_vec.h>

#include <stdbool.h>
#include <stdint.h>


struct md_allocator_i;

#define FLOW_INVALID_INDEX 0xFFFFFFFFu
#define FLOW_MAX_COLUMNS   8

typedef enum flow_node_flag_t {
    FLOW_NODE_FLAG_NONE  = 0,
    // Synthesized by the cut to absorb every sibling below the threshold. Never present in a
    // graph, only ever in a cut, and it is expandable like any other node.
    FLOW_NODE_FLAG_OTHER = 1,
} flow_node_flag_t;

typedef struct flow_node_t {
    uint32_t column;    // 0 .. num_columns-1
    uint32_t parent;    // coarser node in the SAME column, or FLOW_INVALID_INDEX
    uint32_t level;     // 0 = coarsest. Purely descriptive; the hierarchy is 'parent'.
    float    weight;    // share of this node's column total, [0,1]
    md_float4_t   color;
    md_str_t    label;     // NOT owned - must outlive the graph
    md_str_t    sublabel;  // optional second line, same ownership rule; {0} for none
    uint64_t key;       // stable identity across rebuilds; interaction state keys off this
    uint32_t flags;     // flow_node_flag_t

    // What the node stands for REGARDLESS of its column; 0 for nothing in particular. Two nodes in
    // adjacent columns with the same non-zero entity are one thing seen from two sides - a fragment
    // as the donor and as the acceptor - and a link between them is a SELF flow: it does not go
    // anywhere. The layout marks those ribbons (FLOW_RIBBON_FLAG_SELF) so that a drawing backend
    // can show them as part of the nodes rather than as a band crossing the diagram.
    uint64_t entity;
} flow_node_t;

typedef struct flow_link_t {
    uint32_t src;       // node in column c
    uint32_t dst;       // node in column c+1
    float    weight;    // share of the grand total
} flow_link_t;

typedef struct flow_graph_t {
    md_array(flow_node_t) nodes;
    md_array(flow_link_t) links;
    uint32_t num_columns;
    struct md_allocator_i* alloc;
} flow_graph_t;

// ---------------------------------------------------------------------------
// Building
// ---------------------------------------------------------------------------

void     flow_graph_init (flow_graph_t* graph, uint32_t num_columns, struct md_allocator_i* alloc);
void     flow_graph_free (flow_graph_t* graph);
void     flow_graph_clear(flow_graph_t* graph);   // keeps the allocation

// Returns the new node's index, or FLOW_INVALID_INDEX if the node is malformed (column out of
// range, or a parent that is not an earlier node in the same column).
uint32_t flow_graph_add_node(flow_graph_t* graph, const flow_node_t* node);

// Links are only meaningful between LEAF nodes: coarser links are derived by the cut. A link with
// non-positive weight is dropped.
bool flow_graph_add_link(flow_graph_t* graph, uint32_t src, uint32_t dst, float weight);

bool flow_node_is_leaf(const flow_graph_t* graph, uint32_t node_idx);

// ---------------------------------------------------------------------------
// Validation - the four invariants from the design doc
// ---------------------------------------------------------------------------

typedef enum flow_result_t {
    FLOW_OK = 0,
    FLOW_ERR_HIERARCHY,      // parent in another column, forward reference, or a cycle
    FLOW_ERR_LAYERING,       // (1) a link that does not connect two leaves in adjacent columns
    FLOW_ERR_COLUMN_SUM,     // (2) leaf weights within a column do not sum to 1
    FLOW_ERR_BAND_SUM,       // (2) link weights across a band do not sum to 1
    FLOW_ERR_NODE_FLOW,      // (2) a node's links do not sum to its own weight
    FLOW_ERR_PARTITION,      // (3) a node's children do not partition its weight
} flow_result_t;

const char* flow_result_str(flow_result_t result);

// eps is an absolute tolerance on the weight sums. 1.0e-4f is a sane default.
flow_result_t flow_graph_validate(const flow_graph_t* graph, float eps);

// ---------------------------------------------------------------------------
// The cut: what is actually on screen
// ---------------------------------------------------------------------------
//
// A cut is one antichain per column - every leaf of the graph has exactly one ancestor-or-self
// in the cut. Expanding one node while its siblings stay collapsed is the normal case.
//
// Visible nodes are addressed in CUT INDEX SPACE: an index below md_array_size(graph->nodes) is
// a graph node, anything at or above it is a synthesized "Others". Use flow_cut_node() rather
// than indexing either array directly.

typedef struct flow_cut_t {
    // What the USER opened and closed, by graph node index. A pinned node stays as pinned whatever the
    // budget below decides. flow_cut_set_expanded sets the first directly; flow_cut_pin both.
    md_bitmask_t expanded;                  // pinned open
    md_bitmask_t collapsed;                 // pinned closed

    // Row budget per column, 0 for none. Resolving then opens nodes on its own, the heaviest first, as long
    // as no column would show more rows than this (before "Others" folding, which only removes rows): the
    // detail goes where the weight is, and the diagram stays readable however large the graph. Pins win.
    uint32_t max_rows;

    // Treat nodes of the same entity (flow_node_t.entity) in different columns as one: they open, close
    // and fold together, so mirrored columns stay mirrored.
    bool mirror;
    float threshold;                        // fold nodes below this share of their COLUMN, [0,1)

    // An "Others" never grows past this share of its column: folding proceeds from the smallest
    // node upward and stops once taking the next would cross the cap.
    //
    // DEFAULT 1 (uncapped), and it should usually stay there. A cap makes the threshold stop
    // meaning what it says - folding halts early and the LARGEST sub-threshold nodes stay on
    // screen, so a reader who asked for "above 1%" sees 0.5% rows and has no way to tell why.
    // Since an "Others" can now be opened (below), the cap is a preference rather than the fix it
    // was introduced as.
    float other_max;

    // Bit per column: that column's "Others" is open, so nothing in it folds. This is what lets
    // the threshold be honest - hide everything below it, and let the reader look anyway.
    uint32_t other_expanded;

    // Folding is per PARENT: the small children of a node fold into an "Others" inside it, which sits in
    // its parent's slot and band, and the small roots into one of their own. This opens the "Others" of
    // the children of one node (graph index), leaving the rest folded. See flow_cut_open_other.
    md_bitmask_t other_open;

    // Output of flow_cut_resolve, graph node index: open as resolved, by a pin or by the budget.
    md_bitmask_t open;

    // Outputs of flow_cut_resolve()
    md_array(uint32_t)       visible;       // cut index space, grouped by column
    md_array(flow_link_t) links;         // aggregated, endpoints in cut index space
    md_array(flow_node_t) other;         // synthesized "Others" nodes
    md_array(uint32_t)       other_count;   // how many nodes each "Others" absorbed

    struct md_allocator_i* alloc;
} flow_cut_t;

void flow_cut_init(flow_cut_t* cut, struct md_allocator_i* alloc);
void flow_cut_free(flow_cut_t* cut);

void flow_cut_set_expanded(flow_cut_t* cut, uint32_t node_idx, bool expanded);

// What a user's open or close of a node means: pin it (and, with 'mirror', the nodes of its entity) open
// or closed.
void flow_cut_pin(flow_cut_t* cut, const flow_graph_t* graph, uint32_t node_idx, bool open);
void flow_cut_clear_pins(flow_cut_t* cut);

// Open after flow_cut_resolve, whichever way it came to be
bool flow_cut_is_open(const flow_cut_t* cut, uint32_t node_idx);

// Open or close an "Others" (a node of cut->other): its siblings are then drawn individually. With
// 'mirror' the "Others" of the same parent entity in the other columns go with it.
void flow_cut_open_other(flow_cut_t* cut, const flow_graph_t* graph, const flow_node_t* other, bool open);
// Whether a node's siblings are shown individually because their "Others" was opened (a column's, or its
// parent's)
bool flow_cut_is_other_open_for(const flow_cut_t* cut, const flow_node_t* node);

// Open or close a column's "Others". An open one folds nothing, so its members are drawn
// individually however small they are.
void flow_cut_set_other_expanded(flow_cut_t* cut, uint32_t column, bool expanded);
bool flow_cut_is_expanded (const flow_cut_t* cut, uint32_t node_idx);
void flow_cut_collapse_all(flow_cut_t* cut);
void flow_cut_expand_all  (flow_cut_t* cut, const flow_graph_t* graph);

// Recompute 'visible', 'links' and 'other'. Call after any change to expansion, threshold or the
// graph itself. Cheap - it is an aggregation over the graph, no physics involved.
void flow_cut_resolve(flow_cut_t* cut, const flow_graph_t* graph);

// Resolve a cut index to a node. Never returns NULL for an index that appears in 'visible'.
const flow_node_t* flow_cut_node(const flow_cut_t* cut, const flow_graph_t* graph, uint32_t cut_idx);

// Invariant (4): weight is conserved on every valid cut - for any mix of expanded and collapsed
// nodes. This is the check that says an expansion cannot make the diagram lie.
flow_result_t flow_cut_validate(const flow_cut_t* cut, const flow_graph_t* graph, float eps);

// ---------------------------------------------------------------------------
// Layout - pure, in graph space, no drawing
// ---------------------------------------------------------------------------

// How the rows of a column are ordered. Both policies respect the hierarchy first - siblings stay
// contiguous within their ancestor's slot - and differ only in how they rank within it.
typedef enum flow_order_t {
    // Largest contribution at the top. Makes the column readable as a ranking, at the cost of some
    // ribbon crossings.
    FLOW_ORDER_WEIGHT = 0,
    // Barycentre sweeps: order by the mean position of the neighbours a node connects to. Fewer
    // crossings, but the rows are in no order the reader can name.
    FLOW_ORDER_BARYCENTRE,
    // The order the nodes were added to the graph in. A caller that adds the same entities to two
    // columns in the same order gets two MIRRORED columns, where a flow from a row to itself runs
    // straight across and every slanted band is a genuine transfer. "Others" go last.
    FLOW_ORDER_GRAPH,
} flow_order_t;

typedef struct flow_layout_params_t {
    float node_thickness;       // node bar thickness along the flow axis, graph units
    float node_gap;             // gap between adjacent nodes in a column, as a share of the column
    flow_order_t order;
    uint32_t barycentre_sweeps; // only for FLOW_ORDER_BARYCENTRE; 2 is plenty

    // Per column, the exponent applied to a node's weight before it is turned into a height:
    // h ∝ w^e, renormalized within the column. 1 = strictly proportional, the Sankey default.
    // Below 1 compresses the range, so a column whose weight is nearly all in one node still gives
    // the small ones a body you can read, label and click. 0 would make every row identical.
    //
    // This makes a column's scale LOCAL, so a ribbon spanning two columns with different exponents
    // tapers. That is not a lie: each end is proportional within its own column's scale, and the
    // stubs still exactly fill the node at both ends. It IS a reason not to compress a column whose
    // absolute proportions are the story - the atom columns - while compressing one that is really
    // a channel picker, like the NTO pairs.
    float weight_exponent[FLOW_MAX_COLUMNS];
} flow_layout_params_t;

flow_layout_params_t flow_layout_params_default(void);

typedef struct flow_layout_node_t {
    uint32_t cut_idx;
    md_float2_t   min;               // graph space; x is the flow axis, y is the weight axis
    md_float2_t   max;
} flow_layout_node_t;

// A ribbon is a quad with two curved edges. p0/p1 are the near and far edge of the band where it
// leaves the source; q0/q1 where it arrives at the destination. The curve between them is the
// drawing backend's business.
typedef enum flow_ribbon_flag_t {
    FLOW_RIBBON_FLAG_NONE = 0,
    // Source and destination are the same entity (flow_node_t.entity): the stubs say how much of
    // each node stays where it is, and there is nothing to connect.
    FLOW_RIBBON_FLAG_SELF = 1,
} flow_ribbon_flag_t;

typedef struct flow_layout_ribbon_t {
    uint32_t link_idx;          // index into cut->links
    uint32_t flags;             // flow_ribbon_flag_t
    md_float2_t   p0, p1;
    md_float2_t   q0, q1;
} flow_layout_ribbon_t;

// An ancestor that is currently EXPANDED, and therefore not drawn as a node itself. It still has
// an extent - the span its visible descendants occupy - and a diagram that shows that extent is
// one where the user can see which atoms belong to which group while looking at the atoms.
//
// Bands are what makes a group an organizing layer OUTSIDE the column rather than an alternative
// to it. The layout reports the extent and the column's own x range; which side of the column to
// put it on, and how far out, is the drawing backend's decision.
typedef struct flow_layout_band_t {
    uint32_t node;      // graph node index of the ancestor
    uint32_t column;
    uint32_t depth;     // 0 = outermost ancestor, 1 = its child, ...
    md_float2_t   min;       // y spans the visible descendants; x matches the column's node rects
    md_float2_t   max;
} flow_layout_band_t;

typedef struct flow_layout_t {
    md_array(flow_layout_node_t)   nodes;
    md_array(flow_layout_ribbon_t) ribbons;
    md_array(flow_layout_band_t)   bands;
    // Row order per column, cut index space. Kept so a caller can animate between two layouts.
    md_array(uint32_t) order;
    struct md_allocator_i* alloc;
} flow_layout_t;

void flow_layout_init(flow_layout_t* layout, struct md_allocator_i* alloc);
void flow_layout_free(flow_layout_t* layout);

// Lays out into the unit square [0,1] x [0,1] in graph space. Pan, zoom and the mapping to pixels
// belong to the caller, so zooming never needs a relayout.
void flow_layout_compute(flow_layout_t* layout, const flow_graph_t* graph, const flow_cut_t* cut,
                         const flow_layout_params_t* params);


