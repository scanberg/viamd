#include "channels.h"

#include <core/md_allocator.h>
#include <core/md_common.h>

#include <float.h>
#include <math.h>
#include <string.h>

#define CHANNEL_FLAG_TOP    0x1
#define CHANNEL_FLAG_BOTTOM 0x2

typedef struct sweep_ctx_t {
    // Union find over the qualifying voxels, with ids allocated in sweep order so that an older
    // component always carries the smaller id.
    md_array(uint32_t) parent;
    md_array(uint32_t) size;
    md_array(uint8_t)  flags;     // Which faces the component touches
    md_array(uint32_t) node_at;   // Branch currently carrying this component; valid at the root only

    md_array(int32_t)  node_slab; // Scratch: last slab a branch was seen in

    channel_tree_t* tree;
    struct md_allocator_i* alloc;
    float z_now;
    bool  want_tree;
} sweep_ctx_t;

static uint32_t uf_make(sweep_ctx_t* c, uint8_t flags) {
    const uint32_t id = (uint32_t)md_array_size(c->parent);
    md_array_push(c->parent, id, c->alloc);
    md_array_push(c->size,   1u, c->alloc);
    md_array_push(c->flags,  flags, c->alloc);
    if (c->want_tree) {
        md_array_push(c->node_at, CHANNEL_INVALID_INDEX, c->alloc);
    }
    return id;
}

static uint32_t uf_find(sweep_ctx_t* c, uint32_t x) {
    while (c->parent[x] != x) {
        c->parent[x] = c->parent[c->parent[x]];   // path halving
        x = c->parent[x];
    }
    return x;
}

static uint32_t node_new(sweep_ctx_t* c, float z, bool at_top) {
    channel_node_t n;
    MEMSET(&n, 0, sizeof(n));
    n.parent        = CHANNEL_INVALID_INDEX;
    n.first_child   = CHANNEL_INVALID_INDEX;
    n.next_sibling  = CHANNEL_INVALID_INDEX;
    n.z_top         = z;
    n.z_bot         = z;
    n.max_clearance = 0.0f;
    n.min_clearance = FLT_MAX;
    n.reaches_top   = at_top;

    md_array_push(c->tree->nodes, n, c->alloc);
    md_array_push(c->node_slab, INT32_MIN, c->alloc);
    return (uint32_t)md_array_size(c->tree->nodes) - 1;
}

static void node_attach(sweep_ctx_t* c, uint32_t child, uint32_t parent, float z) {
    channel_node_t* nodes = c->tree->nodes;
    nodes[child].parent  = parent;
    nodes[child].z_bot   = MIN(nodes[child].z_bot, z);
    nodes[child].next_sibling = nodes[parent].first_child;
    nodes[parent].first_child = child;
}

static void uf_union(sweep_ctx_t* c, uint32_t a, uint32_t b) {
    uint32_t ra = uf_find(c, a);
    uint32_t rb = uf_find(c, b);
    if (ra == rb) return;

    const uint32_t na = c->want_tree ? c->node_at[ra] : CHANNEL_INVALID_INDEX;
    const uint32_t nb = c->want_tree ? c->node_at[rb] : CHANNEL_INVALID_INDEX;
    const uint8_t  fl = (uint8_t)(c->flags[ra] | c->flags[rb]);

    if (c->size[ra] < c->size[rb]) {
        const uint32_t t = ra; ra = rb; rb = t;
    }
    c->parent[rb] = ra;
    c->size[ra]  += c->size[rb];
    c->flags[ra]  = fl;

    if (!c->want_tree) return;

    uint32_t node;
    if (na == CHANNEL_INVALID_INDEX) {
        node = nb;
    } else if (nb == CHANNEL_INVALID_INDEX || na == nb) {
        node = na;
    } else {
        // Two branches which already have an identity meet here, so this is where they join.
        //
        // Three or more branches meeting in the same slab is one junction, not a chain of them. The
        // unions arrive pairwise, so a merge node already made at this z is folded into rather than
        // nested under, which keeps the junction n-ary. Nesting instead leaves the lower node with no
        // z extent and - because the per slab accumulation below only ever reaches the node the
        // component currently carries - with no voxels, no path and no clearance either, which is
        // what used to fill the diagram with zero rows.
        //
        // A node made at this z is necessarily one of those junctions: branch nodes for this slab are
        // not created until every union below has been applied.
        const channel_node_t* nodes = c->tree->nodes;
        const bool na_fresh = (nodes[na].z_top == c->z_now) && (nodes[na].num_voxels == 0);
        const bool nb_fresh = (nodes[nb].z_top == c->z_now) && (nodes[nb].num_voxels == 0);

        if (na_fresh && nb_fresh) {
            // Two junctions made at this z, reached from different sides: one adopts the other's
            // branches. Reading each next_sibling before reattaching keeps the list well formed. The
            // emptied node is left unreferenced and is never visited again, since every traversal
            // starts from a root.
            node = na;
            uint32_t child = c->tree->nodes[nb].first_child;
            c->tree->nodes[nb].first_child = CHANNEL_INVALID_INDEX;
            while (child != CHANNEL_INVALID_INDEX) {
                const uint32_t next = c->tree->nodes[child].next_sibling;
                node_attach(c, child, node, c->z_now);
                child = next;
            }
        } else if (na_fresh) {
            node = na;
            node_attach(c, nb, node, c->z_now);
        } else if (nb_fresh) {
            node = nb;
            node_attach(c, na, node, c->z_now);
        } else {
            node = node_new(c, c->z_now, false);
            node_attach(c, na, node, c->z_now);
            node_attach(c, nb, node, c->z_now);
        }
    }
    c->node_at[ra] = node;
}

void channel_tree_free(channel_tree_t* tree) {
    if (!tree || !tree->alloc) return;
    for (size_t i = 0; i < md_array_size(tree->nodes); ++i) {
        md_array_free(tree->nodes[i].path, tree->alloc);
    }
    md_array_free(tree->nodes, tree->alloc);
    md_array_free(tree->roots, tree->alloc);
    MEMSET(tree, 0, sizeof(channel_tree_t));
}

bool channel_sweep(channel_tree_t* out, const channel_field_t* field, double probe_radius, bool want_tree, struct md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(field);
    ASSERT(alloc);

    MEMSET(out, 0, sizeof(channel_tree_t));
    out->alloc = alloc;
    out->probe_radius = probe_radius;

    if (!field->data) return false;
    if (field->pbc[2]) return false;                       // The sweep axis has to have two faces
    const int nx = field->dim[0], ny = field->dim[1], nz = field->dim[2];
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;

    const float r = (float)probe_radius;
    const size_t plane = (size_t)nx * (size_t)ny;

    sweep_ctx_t ctx;
    MEMSET(&ctx, 0, sizeof(ctx));
    ctx.tree      = out;
    ctx.alloc     = alloc;
    ctx.want_tree = want_tree;

    md_array(uint32_t) id_prev = 0;
    md_array(uint32_t) id_cur  = 0;
    md_array_resize(id_prev, plane, alloc);
    md_array_resize(id_cur,  plane, alloc);
    for (size_t i = 0; i < plane; ++i) id_prev[i] = CHANNEL_INVALID_INDEX;

    // Sweep from the far face towards the near one. Only the two slabs in flight are ever indexed by
    // position, so the per position storage is one plane rather than one volume.
    for (int k = nz - 1; k >= 0; --k) {
        const float z = field->origin[2] + ((float)k + 0.5f) * field->spacing[2];
        ctx.z_now = z;

        const float* slab = field->data + (size_t)k * plane;

        uint8_t face = 0;
        if (k == nz - 1) face |= CHANNEL_FLAG_TOP;
        if (k == 0)      face |= CHANNEL_FLAG_BOTTOM;

        for (size_t i = 0; i < plane; ++i) {
            id_cur[i] = (slab[i] >= r) ? uf_make(&ctx, face) : CHANNEL_INVALID_INDEX;
        }

        // In plane neighbours
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const size_t p = (size_t)y * nx + x;
                if (id_cur[p] == CHANNEL_INVALID_INDEX) continue;
                if (x > 0 && id_cur[p - 1] != CHANNEL_INVALID_INDEX)  uf_union(&ctx, id_cur[p], id_cur[p - 1]);
                if (y > 0 && id_cur[p - nx] != CHANNEL_INVALID_INDEX) uf_union(&ctx, id_cur[p], id_cur[p - nx]);
            }
        }
        // Periodic seams, joined after the interior so the raster scan above stays a plain two neighbour pass
        if (field->pbc[0] && nx > 1) {
            for (int y = 0; y < ny; ++y) {
                const size_t a = (size_t)y * nx + (nx - 1), b = (size_t)y * nx;
                if (id_cur[a] != CHANNEL_INVALID_INDEX && id_cur[b] != CHANNEL_INVALID_INDEX) uf_union(&ctx, id_cur[a], id_cur[b]);
            }
        }
        if (field->pbc[1] && ny > 1) {
            for (int x = 0; x < nx; ++x) {
                const size_t a = (size_t)(ny - 1) * nx + x, b = (size_t)x;
                if (id_cur[a] != CHANNEL_INVALID_INDEX && id_cur[b] != CHANNEL_INVALID_INDEX) uf_union(&ctx, id_cur[a], id_cur[b]);
            }
        }
        // Onto the slab already swept
        for (size_t i = 0; i < plane; ++i) {
            if (id_cur[i] != CHANNEL_INVALID_INDEX && id_prev[i] != CHANNEL_INVALID_INDEX) {
                uf_union(&ctx, id_cur[i], id_prev[i]);
            }
        }

        if (want_tree) {
            // A component with no branch yet is one which appears at this slab
            for (size_t i = 0; i < plane; ++i) {
                if (id_cur[i] == CHANNEL_INVALID_INDEX) continue;
                const uint32_t root = uf_find(&ctx, id_cur[i]);
                if (ctx.node_at[root] == CHANNEL_INVALID_INDEX) {
                    ctx.node_at[root] = node_new(&ctx, z, k == nz - 1);
                }
            }

            // Representative point per branch per slab: its widest voxel in this slab
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    const size_t i = (size_t)y * nx + x;
                    if (id_cur[i] == CHANNEL_INVALID_INDEX) continue;

                    const uint32_t node = ctx.node_at[uf_find(&ctx, id_cur[i])];
                    channel_node_t* n = ctx.tree->nodes + node;

                    n->num_voxels    += 1;
                    n->max_clearance  = MAX(n->max_clearance, slab[i]);
                    n->z_bot          = MIN(n->z_bot, z);

                    if (ctx.node_slab[node] != k) {
                        ctx.node_slab[node] = k;
                        vec4_t pt = {0, 0, 0, -FLT_MAX};
                        md_array_push(n->path, pt, alloc);
                        n = ctx.tree->nodes + node;   // push may have moved the array
                    }
                    vec4_t* pt = md_array_last(n->path);
                    if (slab[i] > pt->w) {
                        pt->x = field->origin[0] + ((float)x + 0.5f) * field->spacing[0];
                        pt->y = field->origin[1] + ((float)y + 0.5f) * field->spacing[1];
                        pt->z = z;
                        pt->w = slab[i];
                    }
                }
            }
        }

        md_array(uint32_t) tmp = id_prev;
        id_prev = id_cur;
        id_cur  = tmp;
    }

    // Components, and which of them get all the way through
    for (uint32_t i = 0; i < (uint32_t)md_array_size(ctx.parent); ++i) {
        if (uf_find(&ctx, i) != i) continue;
        out->num_components += 1;
        const bool spans = (ctx.flags[i] & CHANNEL_FLAG_TOP) && (ctx.flags[i] & CHANNEL_FLAG_BOTTOM);
        if (spans) out->num_spanning += 1;
        if (want_tree) {
            const uint32_t node = ctx.node_at[i];
            if (node != CHANNEL_INVALID_INDEX) {
                out->nodes[node].reaches_bottom = (ctx.flags[i] & CHANNEL_FLAG_BOTTOM) != 0;
                md_array_push(out->roots, node, alloc);
            }
        }
    }

    if (want_tree) {
        // A merge node is always created after the branches it joins, so one ascending pass carries
        // "something above me reaches the far face" all the way to the roots.
        for (size_t i = 0; i < md_array_size(out->nodes); ++i) {
            channel_node_t* n = out->nodes + i;
            for (size_t p = 0; p < md_array_size(n->path); ++p) {
                n->min_clearance = MIN(n->min_clearance, n->path[p].w);
            }
            if (n->min_clearance == FLT_MAX) n->min_clearance = 0.0f;
            if (n->reaches_top && n->parent != CHANNEL_INVALID_INDEX) {
                out->nodes[n->parent].reaches_top = true;
            }
        }
    }

    md_array_free(id_prev, alloc);
    md_array_free(id_cur,  alloc);
    md_array_free(ctx.parent, alloc);
    md_array_free(ctx.size, alloc);
    md_array_free(ctx.flags, alloc);
    md_array_free(ctx.node_at, alloc);
    md_array_free(ctx.node_slab, alloc);

    return true;
}

float channel_tree_layout(float* out_slot, const channel_tree_t* tree, struct md_allocator_i* alloc) {
    ASSERT(out_slot);
    ASSERT(tree);

    const size_t num_nodes = md_array_size(tree->nodes);
    for (size_t i = 0; i < num_nodes; ++i) out_slot[i] = -1.0f;

    // An explicit traversal stack. Post order: a node's column is only known once every child has
    // one, so a frame stays on the stack accumulating its children's columns until they are done.
    typedef struct {
        uint32_t node;
        uint32_t next_child;
        float    sum;
        uint32_t count;
    } frame_t;

    md_array(frame_t) stack = 0;
    float slot = 0.0f;

    for (size_t i = 0; i < md_array_size(tree->roots); ++i) {
        const uint32_t root = tree->roots[i];
        if (!(tree->nodes[root].reaches_top && tree->nodes[root].reaches_bottom)) continue;

        md_array_shrink(stack, 0);
        frame_t root_frame = { root, tree->nodes[root].first_child, 0.0f, 0 };
        md_array_push(stack, root_frame, alloc);

        while (md_array_size(stack) > 0) {
            frame_t* top = md_array_last(stack);

            if (top->next_child != CHANNEL_INVALID_INDEX) {
                const uint32_t child = top->next_child;
                top->next_child = tree->nodes[child].next_sibling;
                frame_t child_frame = { child, tree->nodes[child].first_child, 0.0f, 0 };
                // The push may move the array, so nothing may be read through top after this
                md_array_push(stack, child_frame, alloc);
                continue;
            }

            const float column = top->count ? (top->sum / (float)top->count) : slot++;
            out_slot[top->node] = column;
            md_array_pop(stack);

            if (md_array_size(stack) > 0) {
                frame_t* parent = md_array_last(stack);
                parent->sum   += column;
                parent->count += 1;
            }
        }
    }

    md_array_free(stack, alloc);
    return slot;
}

void channel_spanning_counts(uint32_t* out_counts, const double* radii, size_t num_radii, const channel_field_t* field, struct md_allocator_i* alloc) {
    for (size_t i = 0; i < num_radii; ++i) {
        channel_tree_t t;
        out_counts[i] = channel_sweep(&t, field, radii[i], false, alloc) ? t.num_spanning : 0;
        channel_tree_free(&t);
    }
}

double channel_critical_radius(const channel_field_t* field, double r_lo, double r_hi, double tol, struct md_allocator_i* alloc) {
    channel_tree_t t;

    // Nothing gets through even at the loosest radius asked for
    if (!channel_sweep(&t, field, r_lo, false, alloc)) return r_lo - tol;
    const bool lo_spans = t.num_spanning > 0;
    channel_tree_free(&t);
    if (!lo_spans) return r_lo - tol;

    if (channel_sweep(&t, field, r_hi, false, alloc) && t.num_spanning > 0) {
        channel_tree_free(&t);
        return r_hi;
    }
    channel_tree_free(&t);

    // Spanning is monotone in the radius: the superlevel set only shrinks as the probe grows, so a
    // radius which gets through implies every smaller one does.
    while (r_hi - r_lo > tol) {
        const double mid = 0.5 * (r_lo + r_hi);
        const bool spans = channel_sweep(&t, field, mid, false, alloc) && t.num_spanning > 0;
        channel_tree_free(&t);
        if (spans) r_lo = mid; else r_hi = mid;
    }
    return r_lo;
}
