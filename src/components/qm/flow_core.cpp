#include "flow_core.h"

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

static bool flow_cut_is_other_expanded(const flow_cut_t* cut, uint32_t column);

// Key namespace for the synthesized "Others" nodes, kept clear of any caller's key scheme so the
// cut can be keyed by node identity without colliding with the graph's own keys.
#define FLOW_OTHER_KEY_TAG 0xF10A0000ULL

// ---------------------------------------------------------------------------
// Graph
// ---------------------------------------------------------------------------

void flow_graph_init(flow_graph_t* graph, uint32_t num_columns, struct md_allocator_i* alloc) {
    MD_ASSERT(graph);
    MD_ASSERT(alloc);
    md_memset(graph, 0, sizeof(flow_graph_t));
    graph->num_columns = num_columns;
    graph->alloc = alloc;
}

void flow_graph_free(flow_graph_t* graph) {
    MD_ASSERT(graph);
    if (graph->alloc) {
        md_array_free(graph->nodes, graph->alloc);
        md_array_free(graph->links, graph->alloc);
    }
    md_memset(graph, 0, sizeof(flow_graph_t));
}

void flow_graph_clear(flow_graph_t* graph) {
    MD_ASSERT(graph);
    md_array_shrink(graph->nodes, 0);
    md_array_shrink(graph->links, 0);
}

uint32_t flow_graph_add_node(flow_graph_t* graph, const flow_node_t* node) {
    MD_ASSERT(graph);
    MD_ASSERT(node);

    if (node->column >= graph->num_columns) {
        MD_LOG_ERROR("flow: node column %u is out of range (graph has %u columns)", node->column, graph->num_columns);
        return FLOW_INVALID_INDEX;
    }

    // A parent must be an EARLIER node in the SAME column. Requiring 'earlier' is what makes the
    // hierarchy acyclic by construction rather than by a check.
    if (node->parent != FLOW_INVALID_INDEX) {
        const uint32_t count = (uint32_t)md_array_size(graph->nodes);
        if (node->parent >= count) {
            MD_LOG_ERROR("flow: node parent %u is not an earlier node", node->parent);
            return FLOW_INVALID_INDEX;
        }
        if (graph->nodes[node->parent].column != node->column) {
            MD_LOG_ERROR("flow: node parent %u lives in column %u, child in column %u", node->parent,
                graph->nodes[node->parent].column, node->column);
            return FLOW_INVALID_INDEX;
        }
    }

    const uint32_t idx = (uint32_t)md_array_size(graph->nodes);
    md_array_push(graph->nodes, *node, graph->alloc);
    return idx;
}

bool flow_graph_add_link(flow_graph_t* graph, uint32_t src, uint32_t dst, float weight) {
    MD_ASSERT(graph);

    if (!(weight > 0.0f)) {
        // Zero and NaN links carry no information and would only add degenerate ribbons.
        return false;
    }

    const uint32_t count = (uint32_t)md_array_size(graph->nodes);
    if (src >= count || dst >= count) {
        MD_LOG_ERROR("flow: link (%u -> %u) references a node that does not exist", src, dst);
        return false;
    }

    const flow_link_t link = { .src = src, .dst = dst, .weight = weight };
    md_array_push(graph->links, link, graph->alloc);
    return true;
}

bool flow_node_is_leaf(const flow_graph_t* graph, uint32_t node_idx) {
    MD_ASSERT(graph);
    const size_t count = md_array_size(graph->nodes);
    for (size_t i = node_idx + 1; i < count; ++i) {
        // Children always come after their parent, so the scan can start just past the node.
        if (graph->nodes[i].parent == node_idx) {
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

const char* flow_result_str(flow_result_t result) {
    switch (result) {
    case FLOW_OK:                return "ok";
    case FLOW_ERR_HIERARCHY:     return "a node's parent is not an earlier node in the same column";
    case FLOW_ERR_LAYERING:      return "a link does not connect two leaves in adjacent columns";
    case FLOW_ERR_COLUMN_SUM:    return "the leaf weights of a column do not sum to one";
    case FLOW_ERR_BAND_SUM:      return "the link weights across a band do not sum to one";
    case FLOW_ERR_NODE_FLOW:     return "a node's links do not sum to its own weight";
    case FLOW_ERR_PARTITION:     return "a node's children do not partition its weight";
    default:                        return "unknown";
    }
}

// Marks, for every node, whether anything claims it as a parent. One pass instead of the
// quadratic flow_node_is_leaf() when the answer is needed for the whole graph.
static void flow_mark_internal(bool* out_has_children, const flow_graph_t* graph) {
    const size_t count = md_array_size(graph->nodes);
    md_memset(out_has_children, 0, sizeof(bool) * count);
    for (size_t i = 0; i < count; ++i) {
        const uint32_t parent = graph->nodes[i].parent;
        if (parent != FLOW_INVALID_INDEX) {
            out_has_children[parent] = true;
        }
    }
}

flow_result_t flow_graph_validate(const flow_graph_t* graph, float eps) {
    MD_ASSERT(graph);

    const size_t num_nodes = md_array_size(graph->nodes);
    const size_t num_links = md_array_size(graph->links);

    if (num_nodes == 0) {
        return FLOW_OK;
    }

    md_temp_scope_t temp = md_temp_begin();

    bool*  has_children  = md_temp_alloc_array(temp, bool,  num_nodes);
    double* child_sum    = md_temp_alloc_zero_array(temp, double, num_nodes);
    double* out_sum      = md_temp_alloc_zero_array(temp, double, num_nodes);
    double* in_sum       = md_temp_alloc_zero_array(temp, double, num_nodes);
    double* column_sum   = md_temp_alloc_zero_array(temp, double, graph->num_columns);
    double* band_sum     = md_temp_alloc_zero_array(temp, double, graph->num_columns);

    flow_mark_internal(has_children, graph);

    flow_result_t result = FLOW_OK;

    // (hierarchy) + accumulate the child sums and the per-column leaf sums
    for (size_t i = 0; i < num_nodes && result == FLOW_OK; ++i) {
        const flow_node_t* node = graph->nodes + i;

        if (node->column >= graph->num_columns) {
            result = FLOW_ERR_HIERARCHY;
            break;
        }
        if (node->parent != FLOW_INVALID_INDEX) {
            if (node->parent >= i || graph->nodes[node->parent].column != node->column) {
                result = FLOW_ERR_HIERARCHY;
                break;
            }
            child_sum[node->parent] += node->weight;
        }
        if (!has_children[i]) {
            column_sum[node->column] += node->weight;
        }
    }

    // (1) layering - links go c -> c+1, and only between leaves, because everything coarser is
    // derived by aggregation and a link on an internal node would be counted twice.
    for (size_t i = 0; i < num_links && result == FLOW_OK; ++i) {
        const flow_link_t* link = graph->links + i;
        const uint32_t src_col = graph->nodes[link->src].column;
        const uint32_t dst_col = graph->nodes[link->dst].column;
        if (src_col + 1 != dst_col || has_children[link->src] || has_children[link->dst]) {
            result = FLOW_ERR_LAYERING;
            break;
        }
        band_sum[src_col] += link->weight;
        out_sum[link->src] += link->weight;
        in_sum[link->dst]  += link->weight;
    }

    // (3) hierarchy partitions
    for (size_t i = 0; i < num_nodes && result == FLOW_OK; ++i) {
        if (has_children[i] && fabs(child_sum[i] - (double)graph->nodes[i].weight) > (double)eps) {
            result = FLOW_ERR_PARTITION;
        }
    }

    // A leaf's ribbons must exactly fill it. This is not decoration: the layout sub-allocates
    // ribbon stubs inside a node's own span, so a node whose links overrun its weight would draw
    // ribbons spilling out of the bar they belong to.
    for (size_t i = 0; i < num_nodes && result == FLOW_OK; ++i) {
        if (has_children[i]) continue;
        const flow_node_t* node = graph->nodes + i;
        if (node->column + 1 < graph->num_columns && fabs(out_sum[i] - (double)node->weight) > (double)eps) {
            result = FLOW_ERR_NODE_FLOW;
        }
        if (node->column > 0 && fabs(in_sum[i] - (double)node->weight) > (double)eps) {
            result = FLOW_ERR_NODE_FLOW;
        }
    }

    // (2) conservation
    for (uint32_t c = 0; c < graph->num_columns && result == FLOW_OK; ++c) {
        if (fabs(column_sum[c] - 1.0) > (double)eps) {
            result = FLOW_ERR_COLUMN_SUM;
        }
    }
    for (uint32_t c = 0; c + 1 < graph->num_columns && result == FLOW_OK; ++c) {
        if (fabs(band_sum[c] - 1.0) > (double)eps) {
            result = FLOW_ERR_BAND_SUM;
        }
    }

    md_temp_end(temp);
    return result;
}

// ---------------------------------------------------------------------------
// Cut
// ---------------------------------------------------------------------------

void flow_cut_init(flow_cut_t* cut, struct md_allocator_i* alloc) {
    MD_ASSERT(cut);
    MD_ASSERT(alloc);
    md_memset(cut, 0, sizeof(flow_cut_t));
    md_bitmask_init(&cut->expanded, 0, MD_BITSET_MAX_INDEX, alloc);  // words are allocated as nodes are expanded
    md_bitmask_init(&cut->collapsed, 0, MD_BITSET_MAX_INDEX, alloc);
    md_bitmask_init(&cut->other_open, 0, MD_BITSET_MAX_INDEX, alloc);
    md_bitmask_init(&cut->open, 0, MD_BITSET_MAX_INDEX, alloc);
    cut->other_max = 1.0f;      // uncapped unless the caller says otherwise
    cut->alloc = alloc;
}

void flow_cut_free(flow_cut_t* cut) {
    MD_ASSERT(cut);
    if (cut->alloc) {
        md_bitmask_free(&cut->expanded);
        md_bitmask_free(&cut->collapsed);
        md_bitmask_free(&cut->other_open);
        md_bitmask_free(&cut->open);
        md_array_free(cut->visible, cut->alloc);
        md_array_free(cut->links,   cut->alloc);
        md_array_free(cut->other,   cut->alloc);
        md_array_free(cut->other_count, cut->alloc);
    }
    md_memset(cut, 0, sizeof(flow_cut_t));
}

void flow_cut_set_expanded(flow_cut_t* cut, uint32_t node_idx, bool expanded) {
    MD_ASSERT(cut);
    if (expanded) {
        md_bitmask_set_bit(&cut->expanded, node_idx);
    } else {
        md_bitmask_clear_bit(&cut->expanded, node_idx);
    }
}

void flow_cut_pin(flow_cut_t* cut, const flow_graph_t* graph, uint32_t node_idx, bool open) {
    MD_ASSERT(cut);
    MD_ASSERT(graph);
    const size_t count = md_array_size(graph->nodes);
    if (node_idx >= count) return;
    const uint64_t entity = graph->nodes[node_idx].entity;
    for (size_t i = 0; i < count; ++i) {
        if (i != node_idx && !(cut->mirror && entity != 0 && graph->nodes[i].entity == entity)) continue;
        if (open) {
            md_bitmask_set_bit(&cut->expanded, (uint32_t)i);
            md_bitmask_clear_bit(&cut->collapsed, (uint32_t)i);
        } else {
            md_bitmask_clear_bit(&cut->expanded, (uint32_t)i);
            md_bitmask_set_bit(&cut->collapsed, (uint32_t)i);
        }
    }
}

void flow_cut_clear_pins(flow_cut_t* cut) {
    MD_ASSERT(cut);
    md_bitmask_clear(&cut->expanded);
    md_bitmask_clear(&cut->collapsed);
}

bool flow_cut_is_open(const flow_cut_t* cut, uint32_t node_idx) {
    MD_ASSERT(cut);
    return md_bitmask_test_bit(&cut->open, node_idx);
}

void flow_cut_open_other(flow_cut_t* cut, const flow_graph_t* graph, const flow_node_t* other, bool open) {
    MD_ASSERT(cut);
    MD_ASSERT(graph);
    if (!other) return;
    if (other->parent == FLOW_INVALID_INDEX) {
        // the roots' "Others": the column's, and with 'mirror' every column's
        for (uint32_t c = 0; c < graph->num_columns; ++c) {
            if (c == other->column || cut->mirror) flow_cut_set_other_expanded(cut, c, open);
        }
        return;
    }
    const size_t count = md_array_size(graph->nodes);
    if (other->parent >= count) return;
    const uint64_t entity = graph->nodes[other->parent].entity;
    for (size_t i = 0; i < count; ++i) {
        if (i != other->parent && !(cut->mirror && entity != 0 && graph->nodes[i].entity == entity)) continue;
        if (open) md_bitmask_set_bit(&cut->other_open, (uint32_t)i);
        else      md_bitmask_clear_bit(&cut->other_open, (uint32_t)i);
    }
}

bool flow_cut_is_other_open_for(const flow_cut_t* cut, const flow_node_t* node) {
    MD_ASSERT(cut);
    if (!node) return false;
    if (flow_cut_is_other_expanded(cut, node->column)) return true;
    return node->parent != FLOW_INVALID_INDEX && md_bitmask_test_bit(&cut->other_open, node->parent);
}

void flow_cut_set_other_expanded(flow_cut_t* cut, uint32_t column, bool expanded) {
    MD_ASSERT(cut);
    if (column >= FLOW_MAX_COLUMNS) return;
    if (expanded) cut->other_expanded |=  (1u << column);
    else          cut->other_expanded &= ~(1u << column);
}

static bool flow_cut_is_other_expanded(const flow_cut_t* cut, uint32_t column) {
    MD_ASSERT(cut);
    return column < FLOW_MAX_COLUMNS && (cut->other_expanded & (1u << column)) != 0;
}

bool flow_cut_is_expanded(const flow_cut_t* cut, uint32_t node_idx) {
    MD_ASSERT(cut);
    return md_bitmask_test_bit(&cut->expanded, node_idx);
}

void flow_cut_collapse_all(flow_cut_t* cut) {
    MD_ASSERT(cut);
    md_bitmask_clear(&cut->expanded);
    md_bitmask_clear(&cut->collapsed);
}

void flow_cut_expand_all(flow_cut_t* cut, const flow_graph_t* graph) {
    MD_ASSERT(cut);
    MD_ASSERT(graph);
    const size_t count = md_array_size(graph->nodes);
    if (count) md_bitmask_set_range(&cut->expanded, 0, (uint32_t)count);
}

const flow_node_t* flow_cut_node(const flow_cut_t* cut, const flow_graph_t* graph, uint32_t cut_idx) {
    MD_ASSERT(cut);
    MD_ASSERT(graph);
    const uint32_t num_graph_nodes = (uint32_t)md_array_size(graph->nodes);
    if (cut_idx < num_graph_nodes) {
        return graph->nodes + cut_idx;
    }
    const uint32_t other_idx = cut_idx - num_graph_nodes;
    if (other_idx < (uint32_t)md_array_size(cut->other)) {
        return cut->other + other_idx;
    }
    return NULL;
}

static int flow_link_cmp(const void* a, const void* b) {
    const flow_link_t* la = (const flow_link_t*)a;
    const flow_link_t* lb = (const flow_link_t*)b;
    if (la->src != lb->src) return la->src < lb->src ? -1 : 1;
    if (la->dst != lb->dst) return la->dst < lb->dst ? -1 : 1;
    return 0;
}

// The weight a node is ranked and thresholded by: its own, or with 'mirror' the largest of its entity's
static float flow_entity_weight(const flow_cut_t* cut, const flow_graph_t* graph, const uint32_t* peer_beg,
                                const uint32_t* peer_list, const uint32_t* peer_group, uint32_t i) {
    if (!cut->mirror || peer_group[i] == FLOW_INVALID_INDEX) return graph->nodes[i].weight;
    float w = 0.0f;
    const uint32_t g = peer_group[i];
    for (uint32_t k = peer_beg[g]; k < peer_beg[g + 1]; ++k) w = MD_MAX(w, graph->nodes[peer_list[k]].weight);
    return w;
}

void flow_cut_resolve(flow_cut_t* cut, const flow_graph_t* graph) {
    MD_ASSERT(cut);
    MD_ASSERT(graph);

    md_array_shrink(cut->visible, 0);
    md_array_shrink(cut->links,   0);
    md_array_shrink(cut->other,   0);
    md_array_shrink(cut->other_count, 0);
    md_bitmask_clear(&cut->open);

    const size_t num_nodes = md_array_size(graph->nodes);
    if (num_nodes == 0) {
        return;
    }

    md_temp_scope_t temp = md_temp_begin();

    bool* has_children = md_temp_alloc_array(temp, bool, num_nodes);
    flow_mark_internal(has_children, graph);

    uint32_t* num_children = md_temp_alloc_zero_array(temp, uint32_t, num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        if (graph->nodes[i].parent != FLOW_INVALID_INDEX) num_children[graph->nodes[i].parent]++;
    }

    // Peers: the nodes of each non-zero entity, grouped. peer_group[i] is the group of node i, INVALID for
    // a node without an entity (or when not mirroring).
    uint32_t* peer_group = md_temp_alloc_array(temp, uint32_t, num_nodes);
    uint32_t* peer_beg   = md_temp_alloc_zero_array(temp, uint32_t, num_nodes + 1);
    uint32_t* peer_list  = md_temp_alloc_array(temp, uint32_t, num_nodes);
    {
        for (size_t i = 0; i < num_nodes; ++i) peer_group[i] = FLOW_INVALID_INDEX;
        uint32_t num_groups = 0;
        if (cut->mirror) {
            // order by entity; a group is a run of equal entities
            uint32_t* order = md_temp_alloc_array(temp, uint32_t, num_nodes);
            uint32_t num = 0;
            for (size_t i = 0; i < num_nodes; ++i) if (graph->nodes[i].entity != 0) order[num++] = (uint32_t)i;
            // insertion sort is quadratic; nodes can be thousands of atoms, so a shell sort
            for (uint32_t gap = num / 2; gap > 0; gap /= 2) {
                for (uint32_t a = gap; a < num; ++a) {
                    const uint32_t v = order[a];
                    uint32_t b = a;
                    while (b >= gap && graph->nodes[order[b - gap]].entity > graph->nodes[v].entity) {
                        order[b] = order[b - gap];
                        b -= gap;
                    }
                    order[b] = v;
                }
            }
            for (uint32_t a = 0; a < num;) {
                uint32_t b = a + 1;
                while (b < num && graph->nodes[order[b]].entity == graph->nodes[order[a]].entity) b++;
                for (uint32_t k = a; k < b; ++k) {
                    peer_group[order[k]] = num_groups;
                    peer_list[k] = order[k];
                }
                peer_beg[num_groups + 1] = b;
                num_groups++;
                a = b;
            }
        }
        (void)num_groups;
    }

    // Pass 1: which nodes are open. Pins first: open where pinned open and not pinned closed.
    bool* open = md_temp_alloc_zero_array(temp, bool, num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        open[i] = has_children[i] && md_bitmask_test_bit(&cut->expanded, (uint32_t)i) &&
                  !md_bitmask_test_bit(&cut->collapsed, (uint32_t)i);
    }

    // The antichain for an open set: a node is in it when every ancestor is open and it is not
    bool* in_cut = md_temp_alloc_zero_array(temp, bool, num_nodes);
    uint32_t* rows = md_temp_alloc_zero_array(temp, uint32_t, graph->num_columns ? graph->num_columns : 1);
    // Reachability needs every ancestor, not only the parent: computed in pre-order
    bool* reach = md_temp_alloc_zero_array(temp, bool, num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        const uint32_t p = graph->nodes[i].parent;
        reach[i] = (p == FLOW_INVALID_INDEX) || (reach[p] && open[p]);
        if (reach[i] && !open[i]) {
            in_cut[i] = true;
            rows[graph->nodes[i].column]++;
        }
    }

    // The budget: open the heaviest node (with its peers) whose opening keeps every column within it
    if (cut->max_rows > 0) {
        bool* tried = md_temp_alloc_zero_array(temp, bool, num_nodes);
        uint32_t* group = md_temp_alloc_array(temp, uint32_t, num_nodes);
        uint32_t* delta = md_temp_alloc_zero_array(temp, uint32_t, graph->num_columns ? graph->num_columns : 1);
        for (;;) {
            uint32_t best = FLOW_INVALID_INDEX;
            float best_w = -1.0f;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (!in_cut[i] || !has_children[i] || tried[i]) continue;
                if (md_bitmask_test_bit(&cut->collapsed, (uint32_t)i)) continue;
                const float w = flow_entity_weight(cut, graph, peer_beg, peer_list, peer_group, (uint32_t)i);
                if (w > best_w) { best_w = w; best = (uint32_t)i; }
            }
            if (best == FLOW_INVALID_INDEX) break;

            uint32_t num_group = 0;
            if (cut->mirror && peer_group[best] != FLOW_INVALID_INDEX) {
                const uint32_t g = peer_group[best];
                for (uint32_t k = peer_beg[g]; k < peer_beg[g + 1]; ++k) {
                    const uint32_t j = peer_list[k];
                    if (in_cut[j] && has_children[j] && !md_bitmask_test_bit(&cut->collapsed, j)) group[num_group++] = j;
                }
            } else {
                group[num_group++] = best;
            }

            for (uint32_t c = 0; c < graph->num_columns; ++c) delta[c] = 0;
            for (uint32_t k = 0; k < num_group; ++k) delta[graph->nodes[group[k]].column] += num_children[group[k]] - 1;
            bool fits = true;
            for (uint32_t c = 0; c < graph->num_columns; ++c) fits &= rows[c] + delta[c] <= cut->max_rows;

            for (uint32_t k = 0; k < num_group; ++k) tried[group[k]] = true;
            if (!fits) continue;

            for (uint32_t k = 0; k < num_group; ++k) {
                const uint32_t j = group[k];
                open[j] = true;
                in_cut[j] = false;
                rows[graph->nodes[j].column] += num_children[j] - 1;
            }
            // the children of what was opened join the antichain (children follow their parent)
            for (uint32_t k = 0; k < num_group; ++k) {
                const uint32_t j = group[k];
                for (size_t i = j + 1; i < num_nodes; ++i) {
                    if (graph->nodes[i].parent == j) in_cut[i] = true;
                }
            }
        }
    }

    // The final antichain, from the final open set
    for (size_t i = 0; i < num_nodes; ++i) {
        const uint32_t p = graph->nodes[i].parent;
        reach[i] = (p == FLOW_INVALID_INDEX) || (reach[p] && open[p]);
        in_cut[i] = reach[i] && !open[i];
        if (reach[i] && open[i]) md_bitmask_set_bit(&cut->open, (uint32_t)i);
    }

    // Pass 2: threshold, per column and parent. The small children of a node fold into an "Others" of
    // their own, inside the parent's slot; the small roots into one per column. Smallest first, up to the
    // cap, and never a single node.
    const uint32_t num_graph_nodes = (uint32_t)num_nodes;

    uint32_t* rep = md_temp_alloc_array(temp, uint32_t, num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) rep[i] = FLOW_INVALID_INDEX;
    bool* folded = md_temp_alloc_zero_array(temp, bool, num_nodes);

    typedef struct fold_cand_t { uint32_t column, parent, node; float w; } fold_cand_t;
    fold_cand_t* cand = md_temp_alloc_array(temp, fold_cand_t, num_nodes);
    uint32_t num_cand = 0;
    if (cut->threshold > 0.0f) {
        for (size_t i = 0; i < num_nodes; ++i) {
            if (!in_cut[i]) continue;
            const flow_node_t* node = graph->nodes + i;
            if (flow_cut_is_other_open_for(cut, node)) continue;     // opened: fold nothing here
            const float w = flow_entity_weight(cut, graph, peer_beg, peer_list, peer_group, (uint32_t)i);
            if (w < cut->threshold) {
                cand[num_cand++] = fold_cand_t{ node->column, node->parent, (uint32_t)i, node->weight };
            }
        }
        // by (column, parent, weight): a shell sort, as above
        for (uint32_t gap = num_cand / 2; gap > 0; gap /= 2) {
            for (uint32_t a = gap; a < num_cand; ++a) {
                const fold_cand_t v = cand[a];
                uint32_t b = a;
                while (b >= gap) {
                    const fold_cand_t* u = cand + b - gap;
                    const bool greater = u->column != v.column ? u->column > v.column
                                       : u->parent != v.parent ? u->parent > v.parent
                                       : u->w > v.w;
                    if (!greater) break;
                    cand[b] = *u;
                    b -= gap;
                }
                cand[b] = v;
            }
        }
    }

    typedef struct fold_group_t { uint32_t column, parent, count; double weight; } fold_group_t;
    md_array(fold_group_t) groups = 0;
    {
        const float cap = (cut->other_max > 0.0f) ? cut->other_max : 1.0f;
        for (uint32_t a = 0; a < num_cand;) {
            uint32_t b = a;
            while (b < num_cand && cand[b].column == cand[a].column && cand[b].parent == cand[a].parent) b++;
            double sum = 0.0;
            uint32_t count = 0;
            for (uint32_t k = a; k < b; ++k) {
                const double w = (double)graph->nodes[cand[k].node].weight;
                if (sum + w > (double)cap) break;
                sum += w;
                folded[cand[k].node] = true;
                count++;
            }
            // One survivor in an Others is just that node wearing a worse label.
            if (count == 1) {
                folded[cand[a].node] = false;
            } else if (count > 1) {
                const fold_group_t g = { cand[a].column, cand[a].parent, count, sum };
                md_array_push(groups, g, md_temp_allocator(temp));
            }
            a = b;
        }
    }

    // The Others nodes, after the graph's in cut index space
    uint32_t* other_of = md_temp_alloc_array(temp, uint32_t, md_array_size(groups) ? md_array_size(groups) : 1);
    for (size_t g = 0; g < md_array_size(groups); ++g) {
        const fold_group_t* fg = groups + g;
        const flow_node_t other = {
            .column = fg->column,
            .parent = fg->parent,
            .level  = fg->parent == FLOW_INVALID_INDEX ? 0 : graph->nodes[fg->parent].level + 1,
            .weight = (float)fg->weight,
            .color  = {0.6f, 0.6f, 0.6f, 1.0f},
            .label  = MD_STR_INIT("Others"),
            .key    = (FLOW_OTHER_KEY_TAG << 32) | ((uint64_t)((fg->parent + 1) & 0xFFFFFFu) << 8) | (uint64_t)(fg->column & 0xFFu),
            .flags  = FLOW_NODE_FLAG_OTHER,
        };
        other_of[g] = num_graph_nodes + (uint32_t)md_array_size(cut->other);
        md_array_push(cut->other, other, cut->alloc);
        md_array_push(cut->other_count, fg->count, cut->alloc);
    }

    // Visible: the graph's in order, each Others right after the last of its siblings
    for (size_t i = 0; i < num_nodes; ++i) {
        if (!in_cut[i]) continue;
        if (folded[i]) {
            for (size_t g = 0; g < md_array_size(groups); ++g) {
                if (groups[g].column == graph->nodes[i].column && groups[g].parent == graph->nodes[i].parent) {
                    rep[i] = other_of[g];
                    break;
                }
            }
            continue;
        }
        rep[i] = (uint32_t)i;
        md_array_push(cut->visible, (uint32_t)i, cut->alloc);
    }
    for (size_t g = 0; g < md_array_size(groups); ++g) {
        // after the last visible node of the same column and parent, else after the last of the column's
        // nodes under that parent's ancestors, else at the end
        size_t pos = md_array_size(cut->visible);
        bool found = false;
        for (size_t v = md_array_size(cut->visible); v-- > 0;) {
            const uint32_t idx = cut->visible[v];
            const flow_node_t* n = flow_cut_node(cut, graph, idx);
            if (n && n->column == groups[g].column && n->parent == groups[g].parent) { pos = v + 1; found = true; break; }
        }
        if (!found) {
            // every sibling folded: where the first of them was, before the next visible node of the column
            uint32_t first = FLOW_INVALID_INDEX;
            for (size_t i = 0; i < num_nodes && first == FLOW_INVALID_INDEX; ++i) {
                if (folded[i] && rep[i] == other_of[g]) first = (uint32_t)i;
            }
            for (size_t v = 0; v < md_array_size(cut->visible); ++v) {
                const uint32_t idx = cut->visible[v];
                if (idx < num_graph_nodes && idx > first) { pos = v; break; }
            }
        }
        md_array_push(cut->visible, 0, cut->alloc);
        for (size_t v = md_array_size(cut->visible) - 1; v > pos; --v) cut->visible[v] = cut->visible[v - 1];
        cut->visible[pos] = other_of[g];
    }

    // Pass 3: every leaf attributes to its visible ancestor-or-self, so the links can be
    // aggregated. This is the step that guarantees a collapsed view and an expanded view tell
    // the same story - both are sums over the same leaf links.
    uint32_t* leaf_rep = md_temp_alloc_array(temp, uint32_t, num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        uint32_t n = (uint32_t)i;
        while (n != FLOW_INVALID_INDEX && !in_cut[n]) {
            n = graph->nodes[n].parent;
        }
        leaf_rep[i] = (n == FLOW_INVALID_INDEX) ? FLOW_INVALID_INDEX : rep[n];
    }

    const size_t num_links = md_array_size(graph->links);
    if (num_links > 0) {
        flow_link_t* scratch = md_temp_alloc_array(temp, flow_link_t, num_links);
        size_t count = 0;
        for (size_t i = 0; i < num_links; ++i) {
            const flow_link_t* link = graph->links + i;
            const uint32_t s = leaf_rep[link->src];
            const uint32_t d = leaf_rep[link->dst];
            if (s == FLOW_INVALID_INDEX || d == FLOW_INVALID_INDEX) continue;
            scratch[count].src    = s;
            scratch[count].dst    = d;
            scratch[count].weight = link->weight;
            count++;
        }

        qsort(scratch, count, sizeof(flow_link_t), flow_link_cmp);

        for (size_t i = 0; i < count;) {
            size_t j = i + 1;
            double sum = (double)scratch[i].weight;
            while (j < count && scratch[j].src == scratch[i].src && scratch[j].dst == scratch[i].dst) {
                sum += (double)scratch[j].weight;
                j++;
            }
            const flow_link_t merged = { .src = scratch[i].src, .dst = scratch[i].dst, .weight = (float)sum };
            md_array_push(cut->links, merged, cut->alloc);
            i = j;
        }
    }

    md_temp_end(temp);
}

flow_result_t flow_cut_validate(const flow_cut_t* cut, const flow_graph_t* graph, float eps) {
    MD_ASSERT(cut);
    MD_ASSERT(graph);

    if (md_array_size(cut->visible) == 0) {
        return FLOW_OK;
    }

    md_temp_scope_t temp = md_temp_begin();
    const uint32_t cut_space = (uint32_t)(md_array_size(graph->nodes) + md_array_size(cut->other));
    double* column_sum = md_temp_alloc_zero_array(temp, double, graph->num_columns);
    double* band_sum   = md_temp_alloc_zero_array(temp, double, graph->num_columns);
    double* out_sum    = md_temp_alloc_zero_array(temp, double, cut_space);
    double* in_sum     = md_temp_alloc_zero_array(temp, double, cut_space);

    flow_result_t result = FLOW_OK;

    const size_t num_visible = md_array_size(cut->visible);
    for (size_t i = 0; i < num_visible; ++i) {
        const flow_node_t* node = flow_cut_node(cut, graph, cut->visible[i]);
        column_sum[node->column] += (double)node->weight;
    }
    const size_t num_links = md_array_size(cut->links);
    for (size_t i = 0; i < num_links; ++i) {
        const flow_node_t* src = flow_cut_node(cut, graph, cut->links[i].src);
        band_sum[src->column] += (double)cut->links[i].weight;
        out_sum[cut->links[i].src] += (double)cut->links[i].weight;
        in_sum[cut->links[i].dst]  += (double)cut->links[i].weight;
    }

    for (size_t i = 0; i < num_visible && result == FLOW_OK; ++i) {
        const uint32_t n = cut->visible[i];
        const flow_node_t* node = flow_cut_node(cut, graph, n);
        if (node->column + 1 < graph->num_columns && fabs(out_sum[n] - (double)node->weight) > (double)eps) {
            result = FLOW_ERR_NODE_FLOW;
        }
        if (node->column > 0 && fabs(in_sum[n] - (double)node->weight) > (double)eps) {
            result = FLOW_ERR_NODE_FLOW;
        }
    }

    for (uint32_t c = 0; c < graph->num_columns && result == FLOW_OK; ++c) {
        if (fabs(column_sum[c] - 1.0) > (double)eps) {
            result = FLOW_ERR_COLUMN_SUM;
        }
    }
    for (uint32_t c = 0; c + 1 < graph->num_columns && result == FLOW_OK; ++c) {
        if (fabs(band_sum[c] - 1.0) > (double)eps) {
            result = FLOW_ERR_BAND_SUM;
        }
    }

    md_temp_end(temp);
    return result;
}

// ---------------------------------------------------------------------------
// Layout
// ---------------------------------------------------------------------------
//
// Everything below works in graph space: the flow axis (x) runs 0..1 across the columns and the
// weight axis (y) runs 0..1 down each column. Pan, zoom and pixels are the caller's business,
// which is what lets a zoom happen without a relayout.

flow_layout_params_t flow_layout_params_default(void) {
    flow_layout_params_t params = {
        .node_thickness    = 0.12f,
        .node_gap          = 0.02f,
        .order             = FLOW_ORDER_WEIGHT,
        .barycentre_sweeps = 2,
    };
    for (uint32_t c = 0; c < FLOW_MAX_COLUMNS; ++c) params.weight_exponent[c] = 1.0f;
    return params;
}

void flow_layout_init(flow_layout_t* layout, struct md_allocator_i* alloc) {
    MD_ASSERT(layout);
    MD_ASSERT(alloc);
    md_memset(layout, 0, sizeof(flow_layout_t));
    layout->alloc = alloc;
}

void flow_layout_free(flow_layout_t* layout) {
    MD_ASSERT(layout);
    if (layout->alloc) {
        md_array_free(layout->nodes,   layout->alloc);
        md_array_free(layout->ribbons, layout->alloc);
        md_array_free(layout->bands,   layout->alloc);
        md_array_free(layout->order,   layout->alloc);
    }
    md_memset(layout, 0, sizeof(flow_layout_t));
}

// Stable insertion sort over an index array, keyed by a parallel float array. Columns and bands
// hold tens of entries at most, so this beats reaching for qsort with no context parameter.
static void flow_sort_by_key(uint32_t* indices, const float* key, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        const uint32_t v = indices[i];
        const float    k = key[v];
        size_t j = i;
        while (j > 0 && key[indices[j - 1]] > k) {
            indices[j] = indices[j - 1];
            j--;
        }
        indices[j] = v;
    }
}


// Sort keys for one column: (which family, then where inside it).
//
// The families are ranked by an INTEGER rather than by their raw score. That matters: two families
// with identical scores - two groups of equal weight, say - would interleave under a float key,
// because their members' scores then decide the global order and siblings end up separated. An
// integer rank cannot tie, so a family always occupies one contiguous run.
static void flow_column_keys(float* out_key, const uint32_t* nodes, uint32_t count, const uint32_t* root_of,
                             const float* family_score, const float* member_score, float* scratch_score,
                             uint32_t* scratch_root)
{
    uint32_t num_roots = 0;
    for (uint32_t i = 0; i < count; ++i) {
        const uint32_t r = root_of[nodes[i]];
        bool seen = false;
        for (uint32_t j = 0; j < num_roots; ++j) {
            if (scratch_root[j] == r) { seen = true; break; }
        }
        if (!seen) {
            scratch_root[num_roots] = r;
            scratch_score[num_roots] = family_score[nodes[i]];
            num_roots++;
        }
    }

    for (uint32_t i = 0; i < count; ++i) {
        const uint32_t r = root_of[nodes[i]];
        uint32_t rank = 0;
        float score = 0.0f;
        for (uint32_t j = 0; j < num_roots; ++j) {
            if (scratch_root[j] == r) { score = scratch_score[j]; break; }
        }
        // Rank counts the families that score HIGHER, so the highest-scoring family gets rank 0
        // and lands at the top. Callers that want "small value at the top" negate their score.
        for (uint32_t j = 0; j < num_roots; ++j) {
            if (scratch_score[j] > score || (scratch_score[j] == score && scratch_root[j] < r)) {
                rank++;
            }
        }
        // Member scores are in [-1, 1] for weights and [0, count) for barycentres, so the family
        // stride has to clear the larger of the two.
        out_key[nodes[i]] = (float)rank * (float)(count + 2) + member_score[nodes[i]];
    }
}

void flow_layout_compute(flow_layout_t* layout, const flow_graph_t* graph, const flow_cut_t* cut,
                         const flow_layout_params_t* params)
{
    MD_ASSERT(layout);
    MD_ASSERT(graph);
    MD_ASSERT(cut);

    md_array_shrink(layout->nodes,   0);
    md_array_shrink(layout->ribbons, 0);
    md_array_shrink(layout->bands,   0);
    md_array_shrink(layout->order,   0);

    const size_t num_visible = md_array_size(cut->visible);
    const uint32_t num_columns = graph->num_columns;
    if (num_visible == 0 || num_columns == 0) {
        return;
    }

    const flow_layout_params_t p = params ? *params : flow_layout_params_default();

    md_temp_scope_t temp = md_temp_begin();

    // Cut index space is dense: graph nodes first, then the synthesized "Others".
    const uint32_t cut_space = (uint32_t)(md_array_size(graph->nodes) + md_array_size(cut->other));

    float*    row      = md_temp_alloc_array(temp, float,    cut_space);
    float*    bary     = md_temp_alloc_array(temp, float,    cut_space);
    uint32_t* root_of  = md_temp_alloc_array(temp, uint32_t, cut_space);
    float*    y_min    = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    y_max    = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    key      = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    fam_key  = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    mem_key  = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    out_cur  = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    in_cur   = md_temp_alloc_zero_array(temp, float, cut_space);
    uint32_t* col_of   = md_temp_alloc_array(temp, uint32_t, cut_space);

    for (uint32_t i = 0; i < cut_space; ++i) {
        row[i] = 0.0f;
        col_of[i] = FLOW_INVALID_INDEX;
        root_of[i] = i;
    }

    // A visible node's outermost ancestor. Crossing reduction orders by this first, so that
    // expanding a group shows its members side by side in the group's own slot instead of
    // scattering them down the column - an expansion that shuffles the whole column reads as a
    // different diagram rather than as a closer look at the same one.
    {
        const uint32_t num_graph_nodes = (uint32_t)md_array_size(graph->nodes);
        for (uint32_t i = 0; i < num_graph_nodes; ++i) {
            uint32_t r = i;
            while (graph->nodes[r].parent != FLOW_INVALID_INDEX) {
                r = graph->nodes[r].parent;
            }
            root_of[i] = r;
        }
        // An "Others" belongs to the family of its parent, so it sits in the parent's slot
        for (size_t k = 0; k < md_array_size(cut->other); ++k) {
            const uint32_t p = cut->other[k].parent;
            if (p != FLOW_INVALID_INDEX && p < num_graph_nodes) root_of[num_graph_nodes + k] = root_of[p];
        }
    }

    // Bucket the visible nodes by column, preserving insertion order as the initial row order.
    uint32_t* col_count = md_temp_alloc_zero_array(temp, uint32_t, num_columns);
    for (size_t i = 0; i < num_visible; ++i) {
        const uint32_t cut_idx = cut->visible[i];
        const flow_node_t* node = flow_cut_node(cut, graph, cut_idx);
        if (!node) continue;
        col_of[cut_idx] = node->column;
        row[cut_idx] = (float)col_count[node->column]++;
    }

    uint32_t max_col_count = 0;
    for (uint32_t c = 0; c < num_columns; ++c) {
        max_col_count = MD_MAX(max_col_count, col_count[c]);
    }
    if (max_col_count == 0) {
        md_temp_end(temp);
        return;
    }

    // Column-major index lists, so a sweep can reorder one column at a time.
    uint32_t** col_nodes = md_temp_alloc_array(temp, uint32_t*, num_columns);
    for (uint32_t c = 0; c < num_columns; ++c) {
        col_nodes[c] = col_count[c] ? md_temp_alloc_array(temp, uint32_t, col_count[c]) : NULL;
    }
    {
        uint32_t* fill = md_temp_alloc_zero_array(temp, uint32_t, num_columns);
        for (size_t i = 0; i < num_visible; ++i) {
            const uint32_t cut_idx = cut->visible[i];
            const uint32_t c = col_of[cut_idx];
            if (c == FLOW_INVALID_INDEX) continue;
            col_nodes[c][fill[c]++] = cut_idx;
        }
    }

    // Crossing reduction: order each column by the weighted mean row of its neighbours in the
    // column the sweep came from. Two sweeps is plenty at these sizes and it is the single
    // largest readability win over a fixed order.
    const size_t num_cut_links = md_array_size(cut->links);

    // Rank by contribution: heaviest first, hierarchy respected. The ancestor's own weight ranks
    // the families and a node's weight ranks it inside its family, so a column reads as a ranking
    // at every level of expansion rather than only at the top one.
    if (p.order == FLOW_ORDER_WEIGHT) {
        float*    scratch_score = md_temp_alloc_array(temp, float,    cut_space);
        uint32_t* scratch_root  = md_temp_alloc_array(temp, uint32_t, cut_space);

        for (uint32_t c = 0; c < num_columns; ++c) {
            if (col_count[c] < 2) continue;
            for (uint32_t i = 0; i < col_count[c]; ++i) {
                const uint32_t n = col_nodes[c][i];
                const uint32_t r = root_of[n];
                // A family's score is the ancestor's own weight when it has one, and the node's
                // own weight when the node is its own root.
                fam_key[n] = (r != n && r < (uint32_t)md_array_size(graph->nodes))
                           ? graph->nodes[r].weight
                           : flow_cut_node(cut, graph, n)->weight;
                // Negated, because the sort is ascending and the heaviest belongs at the top.
                mem_key[n] = -flow_cut_node(cut, graph, n)->weight;
            }
            flow_column_keys(key, col_nodes[c], col_count[c], root_of, fam_key, mem_key,
                             scratch_score, scratch_root);
            flow_sort_by_key(col_nodes[c], key, col_count[c]);
            for (uint32_t i = 0; i < col_count[c]; ++i) {
                row[col_nodes[c][i]] = (float)i;
            }
        }
    }

    float*    sweep_score = md_temp_alloc_array(temp, float,    cut_space);
    uint32_t* sweep_root  = md_temp_alloc_array(temp, uint32_t, cut_space);
    const uint32_t sweeps = (p.order == FLOW_ORDER_BARYCENTRE) ? p.barycentre_sweeps : 0u;
    for (uint32_t sweep = 0; sweep < sweeps; ++sweep) {
        for (int dir = 0; dir < 2; ++dir) {
            const bool forward = (dir == 0);
            for (uint32_t step = 0; step + 1 < num_columns; ++step) {
                const uint32_t c = forward ? (step + 1) : (num_columns - 2 - step);
                if (col_count[c] < 2) continue;

                for (uint32_t i = 0; i < col_count[c]; ++i) {
                    bary[col_nodes[c][i]] = row[col_nodes[c][i]];
                }

                for (uint32_t i = 0; i < col_count[c]; ++i) {
                    const uint32_t n = col_nodes[c][i];
                    double num = 0.0;
                    double den = 0.0;
                    for (size_t l = 0; l < num_cut_links; ++l) {
                        const flow_link_t* link = cut->links + l;
                        const uint32_t self  = forward ? link->dst : link->src;
                        const uint32_t other = forward ? link->src : link->dst;
                        if (self != n) continue;
                        num += (double)row[other] * (double)link->weight;
                        den += (double)link->weight;
                    }
                    if (den > 0.0) {
                        bary[n] = (float)(num / den);
                    }
                }

                // Order by (which family, own barycentre), the family ranked by the weighted mean
                // barycentre of its members. Negated below because a family's rank counts how many
                // score LOWER, and here a smaller barycentre belongs nearer the top.
                for (uint32_t i = 0; i < col_count[c]; ++i) {
                    const uint32_t n = col_nodes[c][i];
                    const uint32_t r = root_of[n];
                    double num = 0.0;
                    double den = 0.0;
                    for (uint32_t j = 0; j < col_count[c]; ++j) {
                        const uint32_t m = col_nodes[c][j];
                        if (root_of[m] != r) continue;
                        const float w = flow_cut_node(cut, graph, m)->weight;
                        num += (double)bary[m] * (double)w;
                        den += (double)w;
                    }
                    fam_key[n] = -((den > 0.0) ? (float)(num / den) : bary[n]);  // small bary -> top
                    mem_key[n] = bary[n];
                }
                flow_column_keys(key, col_nodes[c], col_count[c], root_of, fam_key, mem_key,
                                 sweep_score, sweep_root);

                flow_sort_by_key(col_nodes[c], key, col_count[c]);
                for (uint32_t i = 0; i < col_count[c]; ++i) {
                    row[col_nodes[c][i]] = (float)i;
                }
            }
        }
    }

    // One weight scale shared by every column, so a ribbon has the same thickness at both ends.
    // The gap budget is set by the busiest column; shorter columns are centred in the area.
    float gap = p.node_gap;
    float usable = 1.0f - gap * (float)(max_col_count > 0 ? max_col_count - 1 : 0);
    if (usable < 0.1f) {
        // Too many nodes for the requested gap - shrink the gap rather than invert the layout.
        usable = 0.1f;
        gap = (max_col_count > 1) ? (1.0f - usable) / (float)(max_col_count - 1) : 0.0f;
    }

    const float thickness = MD_CLAMP(p.node_thickness, 0.0f, 1.0f);
    const float span = (num_columns > 1) ? (1.0f - thickness) / (float)(num_columns - 1) : 0.0f;

    // A node's drawn share of its column, after the column's own exponent. Kept separately from
    // the node's weight, which stays the truth the tooltips and the invariants speak in.
    float* draw_share = md_temp_alloc_zero_array(temp, float, cut_space);

    for (uint32_t c = 0; c < num_columns; ++c) {
        if (col_count[c] == 0) continue;

        const float e = (c < FLOW_MAX_COLUMNS) ? p.weight_exponent[c] : 1.0f;

        double share_sum = 0.0;
        for (uint32_t i = 0; i < col_count[c]; ++i) {
            const uint32_t n = col_nodes[c][i];
            const float w = flow_cut_node(cut, graph, n)->weight;
            const double sh = (e == 1.0f) ? (double)w : pow((double)MD_MAX(w, 0.0f), (double)e);
            draw_share[n] = (float)sh;
            share_sum += sh;
        }
        if (share_sum <= 0.0) continue;
        for (uint32_t i = 0; i < col_count[c]; ++i) {
            draw_share[col_nodes[c][i]] /= (float)share_sum;
        }

        // Compressed or not, a column fills the same band - only the split between its rows moves.
        const float height = usable + gap * (float)(col_count[c] - 1);
        float y = (1.0f - height) * 0.5f;

        const float x0 = (num_columns > 1) ? (float)c * span : (1.0f - thickness) * 0.5f;
        const float x1 = x0 + thickness;

        for (uint32_t i = 0; i < col_count[c]; ++i) {
            const uint32_t n = col_nodes[c][i];
            const float h = draw_share[n] * usable;

            y_min[n] = y;
            y_max[n] = y + h;
            out_cur[n] = y;
            in_cur[n] = y;

            const flow_layout_node_t ln = {
                .cut_idx = n,
                .min = { x0, y },
                .max = { x1, y + h },
            };
            md_array_push(layout->nodes, ln, layout->alloc);
            md_array_push(layout->order, n, layout->alloc);

            y += h + gap;
        }
    }

    // Bands: the extent of every expanded ancestor, taken from where its descendants ended up.
    // Derived rather than laid out, so a band can never disagree with the nodes it brackets.
    {
        const uint32_t num_graph_nodes = (uint32_t)md_array_size(graph->nodes);
        float* band_min = md_temp_alloc_array(temp, float, num_graph_nodes ? num_graph_nodes : 1);
        float* band_max = md_temp_alloc_array(temp, float, num_graph_nodes ? num_graph_nodes : 1);
        bool*  band_any = md_temp_alloc_zero_array(temp, bool, num_graph_nodes ? num_graph_nodes : 1);

        for (uint32_t i = 0; i < num_graph_nodes; ++i) {
            band_min[i] =  1.0e30f;
            band_max[i] = -1.0e30f;
        }

        for (size_t i = 0; i < num_visible; ++i) {
            const uint32_t n = cut->visible[i];
            if (col_of[n] == FLOW_INVALID_INDEX) continue;
            const flow_node_t* vn = flow_cut_node(cut, graph, n);   // an "Others" has its parent's ancestry
            if (!vn) continue;

            for (uint32_t p = vn->parent; p != FLOW_INVALID_INDEX; p = graph->nodes[p].parent) {
                band_any[p] = true;
                band_min[p] = MD_MIN(band_min[p], y_min[n]);
                band_max[p] = MD_MAX(band_max[p], y_max[n]);
            }
        }

        for (uint32_t i = 0; i < num_graph_nodes; ++i) {
            if (!band_any[i]) continue;

            const uint32_t c = graph->nodes[i].column;
            const float x0 = (num_columns > 1) ? (float)c * span : (1.0f - thickness) * 0.5f;

            uint32_t depth = 0;
            for (uint32_t p = graph->nodes[i].parent; p != FLOW_INVALID_INDEX; p = graph->nodes[p].parent) {
                depth++;
            }

            const flow_layout_band_t band = {
                .node   = i,
                .column = c,
                .depth  = depth,
                .min    = { x0, band_min[i] },
                .max    = { x0 + thickness, band_max[i] },
            };
            md_array_push(layout->bands, band, layout->alloc);
        }
    }

    // Ribbon stubs. Outgoing are allocated in destination-row order and incoming in source-row
    // order; that alone removes most self-crossings and costs nothing.
    if (num_cut_links > 0) {
        uint32_t* by_out  = md_temp_alloc_array(temp, uint32_t, num_cut_links);
        uint32_t* by_in   = md_temp_alloc_array(temp, uint32_t, num_cut_links);
        float*    link_key = md_temp_alloc_array(temp, float,   num_cut_links);

        // Sort keys pack the two row indices into one float. Rows are small integers, so a
        // scale of 4096 keeps the pairs ordered without needing a comparator with context.
        for (size_t l = 0; l < num_cut_links; ++l) {
            by_out[l] = (uint32_t)l;
            by_in[l]  = (uint32_t)l;
        }

        for (size_t l = 0; l < num_cut_links; ++l) {
            link_key[l] = row[cut->links[l].src] * 4096.0f + row[cut->links[l].dst];
        }
        flow_sort_by_key(by_out, link_key, num_cut_links);

        // A stub takes the link's share OF ITS OWN NODE, scaled by that node's drawn height. With
        // every exponent at 1 this is identical to weight * usable; with a compressed column it is
        // what keeps the stubs exactly filling the node instead of overflowing it.
        float* out_span_min = md_temp_alloc_array(temp, float, num_cut_links);
        float* out_span_max = md_temp_alloc_array(temp, float, num_cut_links);
        for (size_t i = 0; i < num_cut_links; ++i) {
            const uint32_t l = by_out[i];
            const uint32_t s = cut->links[l].src;
            const flow_node_t* src = flow_cut_node(cut, graph, s);
            const float denom = (src && src->weight > 0.0f) ? src->weight : 1.0f;
            const float h = (cut->links[l].weight / denom) * (y_max[s] - y_min[s]);
            out_span_min[l] = out_cur[s];
            out_span_max[l] = out_cur[s] + h;
            out_cur[s] += h;
        }

        for (size_t l = 0; l < num_cut_links; ++l) {
            link_key[l] = row[cut->links[l].dst] * 4096.0f + row[cut->links[l].src];
        }
        flow_sort_by_key(by_in, link_key, num_cut_links);

        for (size_t i = 0; i < num_cut_links; ++i) {
            const uint32_t l = by_in[i];
            const uint32_t s = cut->links[l].src;
            const uint32_t d = cut->links[l].dst;
            const flow_node_t* dst = flow_cut_node(cut, graph, d);
            const float denom = (dst && dst->weight > 0.0f) ? dst->weight : 1.0f;
            const float h = (cut->links[l].weight / denom) * (y_max[d] - y_min[d]);

            const uint32_t src_col = col_of[s];
            const uint32_t dst_col = col_of[d];
            if (src_col == FLOW_INVALID_INDEX || dst_col == FLOW_INVALID_INDEX) continue;

            const float src_x = ((num_columns > 1) ? (float)src_col * span : (1.0f - thickness) * 0.5f) + thickness;
            const float dst_x = (num_columns > 1) ? (float)dst_col * span : (1.0f - thickness) * 0.5f;

            const flow_node_t* src = flow_cut_node(cut, graph, s);
            const bool self = src && dst && src->entity != 0 && src->entity == dst->entity;

            const flow_layout_ribbon_t ribbon = {
                .link_idx = l,
                .flags = self ? FLOW_RIBBON_FLAG_SELF : FLOW_RIBBON_FLAG_NONE,
                .p0 = { src_x, out_span_min[l] },
                .p1 = { src_x, out_span_max[l] },
                .q0 = { dst_x, in_cur[d] },
                .q1 = { dst_x, in_cur[d] + h },
            };
            md_array_push(layout->ribbons, ribbon, layout->alloc);
            in_cur[d] += h;
        }
    }

    md_temp_end(temp);
}

static const flow_layout_node_t* flow_find_node(const flow_layout_t* layout, uint32_t cut_idx) {
    const size_t count = md_array_size(layout->nodes);
    for (size_t i = 0; i < count; ++i) {
        if (layout->nodes[i].cut_idx == cut_idx) return layout->nodes + i;
    }
    return NULL;
}

static const flow_layout_ribbon_t* flow_find_ribbon(const flow_layout_t* layout, uint32_t link_idx) {
    const size_t count = md_array_size(layout->ribbons);
    for (size_t i = 0; i < count; ++i) {
        if (layout->ribbons[i].link_idx == link_idx) return layout->ribbons + i;
    }
    return NULL;
}

static md_float2_t flow_lerp2(md_float2_t a, md_float2_t b, float t) {
    md_float2_t r = { a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t };
    return r;
}
