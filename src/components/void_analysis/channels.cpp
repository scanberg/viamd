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

// ---------------------------------------------------------------------------------------------
// Percolation: one pass in order of decreasing clearance
// ---------------------------------------------------------------------------------------------

#define PERC_INVALID 0xFFFFFFFFu
#define CHANNEL_FLAG_BOTH (CHANNEL_FLAG_TOP | CHANNEL_FLAG_BOTTOM)

typedef struct perc_ctx_t {
    // Union find over the active voxels, indexed by position in the clearance-ordered list. That
    // ordering is what does the work: an id smaller than the one being inserted is, by construction,
    // a voxel that is already in - so "have I seen this neighbour" is a comparison, not a lookup.
    uint32_t* parent;
    uint32_t* size;
    uint8_t*  flags;
    uint16_t* min_k;
    uint16_t* max_k;

    // Running totals, maintained through the unions rather than recomputed per sample. A component
    // leaves the totals before it is merged and the merged one re-enters, so every sample is a read
    // of four integers.
    uint64_t open_voxels;
    uint64_t span_voxels;
    uint32_t num_roots;
    uint32_t num_spanning;

    // Monotone: a component's reach only grows and components only merge, so these never need to be
    // taken back out.
    int32_t deepest_top;
    int32_t highest_bot;
} perc_ctx_t;

static uint32_t perc_find(perc_ctx_t* c, uint32_t x) {
    while (c->parent[x] != x) {
        c->parent[x] = c->parent[c->parent[x]];   // path halving
        x = c->parent[x];
    }
    return x;
}

static void perc_enter(perc_ctx_t* c, uint32_t r) {
    const uint8_t f = c->flags[r];
    if (f) c->open_voxels += c->size[r];
    if (f == CHANNEL_FLAG_BOTH) {
        c->span_voxels += c->size[r];
        c->num_spanning += 1;
    }
    if (f & CHANNEL_FLAG_TOP)    c->deepest_top = MIN(c->deepest_top, (int32_t)c->min_k[r]);
    if (f & CHANNEL_FLAG_BOTTOM) c->highest_bot = MAX(c->highest_bot, (int32_t)c->max_k[r]);
}

static void perc_leave(perc_ctx_t* c, uint32_t r) {
    const uint8_t f = c->flags[r];
    if (f) c->open_voxels -= c->size[r];
    if (f == CHANNEL_FLAG_BOTH) {
        c->span_voxels -= c->size[r];
        c->num_spanning -= 1;
    }
}

// Returns true when this union is what made a component touch both faces for the first time.
static bool perc_union(perc_ctx_t* c, uint32_t a, uint32_t b) {
    uint32_t ra = perc_find(c, a);
    uint32_t rb = perc_find(c, b);
    if (ra == rb) return false;

    const bool was_spanning = (c->flags[ra] == CHANNEL_FLAG_BOTH) || (c->flags[rb] == CHANNEL_FLAG_BOTH);

    perc_leave(c, ra);
    perc_leave(c, rb);

    if (c->size[ra] < c->size[rb]) { const uint32_t t = ra; ra = rb; rb = t; }
    c->parent[rb] = ra;
    c->size[ra]  += c->size[rb];
    c->flags[ra]  = (uint8_t)(c->flags[ra] | c->flags[rb]);
    c->min_k[ra]  = MIN(c->min_k[ra], c->min_k[rb]);
    c->max_k[ra]  = MAX(c->max_k[ra], c->max_k[rb]);
    c->num_roots -= 1;

    perc_enter(c, ra);
    return !was_spanning && c->flags[ra] == CHANNEL_FLAG_BOTH;
}

void channel_percolation_free(channel_percolation_t* perc) {
    if (!perc || !perc->alloc) return;
    md_array_free(perc->radius, perc->alloc);
    md_array_free(perc->frac_void, perc->alloc);
    md_array_free(perc->frac_open, perc->alloc);
    md_array_free(perc->frac_spanning, perc->alloc);
    md_array_free(perc->z_from_top, perc->alloc);
    md_array_free(perc->z_from_bottom, perc->alloc);
    md_array_free(perc->num_components, perc->alloc);
    md_array_free(perc->num_spanning, perc->alloc);
    MEMSET(perc, 0, sizeof(channel_percolation_t));
}

bool channel_percolate(channel_percolation_t* out, const channel_field_t* field, double r_min, uint32_t num_samples, struct md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(field);
    ASSERT(alloc);

    MEMSET(out, 0, sizeof(channel_percolation_t));
    out->alloc = alloc;

    if (!field->data) return false;
    if (field->pbc[2]) return false;                       // The sweep axis has to have two faces
    const int nx = field->dim[0], ny = field->dim[1], nz = field->dim[2];
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;
    if (nz > 65535) return false;                          // min_k / max_k are uint16

    const size_t plane = (size_t)nx * (size_t)ny;
    const size_t N     = plane * (size_t)nz;
    if (N > 0xFFFFFFF0u) return false;                     // Voxel indices are uint32 throughout

    const float rmin = (float)r_min;

    // What the run is sized by, and the widest clearance the curves have to reach.
    float  d_max = rmin;
    size_t M     = 0;
    for (size_t i = 0; i < N; ++i) {
        const float d = field->data[i];
        if (d >= rmin) {
            M += 1;
            if (d > d_max) d_max = d;
        }
    }
    out->num_active = M;
    if (M == 0) return false;

    const uint32_t ns  = (uint32_t)CLAMP((int)num_samples, 2, 4096);
    // Internal buckets are far finer than the reported samples. The bucket width is what bounds the
    // error on r_c, and at these counts it lands orders of magnitude below the grid's own limit -
    // half a voxel of unresolved gap - so the grid stays the thing that limits the answer.
    const uint32_t per = MAX(1u, (4096u + ns - 1u) / ns);
    const uint32_t NB  = ns * per;

    const double span  = MAX(1.0e-6, (double)d_max - (double)rmin);
    const double width = span / (double)NB;

    const size_t bytes = N * sizeof(uint32_t)
                       + M * (sizeof(uint32_t) * 3 + sizeof(uint8_t) + sizeof(uint16_t) * 2)
                       + (size_t)NB * sizeof(uint32_t) * 2;
    out->bytes = bytes;

    uint32_t* cid = (uint32_t*)md_alloc(alloc, N * sizeof(uint32_t));
    uint32_t* vox = (uint32_t*)md_alloc(alloc, M * sizeof(uint32_t));
    uint32_t* off = (uint32_t*)md_alloc(alloc, (size_t)NB * sizeof(uint32_t));
    uint32_t* cur = (uint32_t*)md_alloc(alloc, (size_t)NB * sizeof(uint32_t));

    perc_ctx_t c;
    MEMSET(&c, 0, sizeof(c));
    c.parent = (uint32_t*)md_alloc(alloc, M * sizeof(uint32_t));
    c.size   = (uint32_t*)md_alloc(alloc, M * sizeof(uint32_t));
    c.flags  = (uint8_t*) md_alloc(alloc, M * sizeof(uint8_t));
    c.min_k  = (uint16_t*)md_alloc(alloc, M * sizeof(uint16_t));
    c.max_k  = (uint16_t*)md_alloc(alloc, M * sizeof(uint16_t));
    c.deepest_top = INT32_MAX;
    c.highest_bot = INT32_MIN;

    if (!cid || !vox || !off || !cur || !c.parent || !c.size || !c.flags || !c.min_k || !c.max_k) {
        md_free(alloc, cid, N * sizeof(uint32_t));
        md_free(alloc, vox, M * sizeof(uint32_t));
        md_free(alloc, off, (size_t)NB * sizeof(uint32_t));
        md_free(alloc, cur, (size_t)NB * sizeof(uint32_t));
        md_free(alloc, c.parent, M * sizeof(uint32_t));
        md_free(alloc, c.size,   M * sizeof(uint32_t));
        md_free(alloc, c.flags,  M * sizeof(uint8_t));
        md_free(alloc, c.min_k,  M * sizeof(uint16_t));
        md_free(alloc, c.max_k,  M * sizeof(uint16_t));
        return false;
    }

    MEMSET(off, 0, (size_t)NB * sizeof(uint32_t));

    // Counting sort by clearance, descending. A comparison sort of a hundred million voxels would
    // cost more than the union find it feeds.
    for (size_t i = 0; i < N; ++i) {
        cid[i] = PERC_INVALID;
        const float d = field->data[i];
        if (d < rmin) continue;
        int b = (int)(((double)d - (double)rmin) / width);
        b = CLAMP(b, 0, (int)NB - 1);
        off[b] += 1;
    }
    {   // Bucket NB-1 is widest and goes first, so the offsets accumulate downward
        uint32_t acc = 0;
        for (int b = (int)NB - 1; b >= 0; --b) {
            const uint32_t n = off[b];
            off[b] = acc;
            cur[b] = acc;
            acc += n;
        }
    }
    for (size_t i = 0; i < N; ++i) {
        const float d = field->data[i];
        if (d < rmin) continue;
        int b = (int)(((double)d - (double)rmin) / width);
        b = CLAMP(b, 0, (int)NB - 1);
        const uint32_t p = cur[b]++;
        vox[p] = (uint32_t)i;
        cid[i] = p;
    }

    md_array_resize(out->radius,         ns, alloc);
    md_array_resize(out->frac_void,      ns, alloc);
    md_array_resize(out->frac_open,      ns, alloc);
    md_array_resize(out->frac_spanning,  ns, alloc);
    md_array_resize(out->z_from_top,     ns, alloc);
    md_array_resize(out->z_from_bottom,  ns, alloc);
    md_array_resize(out->num_components, ns, alloc);
    md_array_resize(out->num_spanning,   ns, alloc);

    const double z_lo   = (double)field->origin[2];
    const double z_hi   = (double)field->origin[2] + (double)field->spacing[2] * (double)nz;
    const double inv_N  = 1.0 / (double)N;

    uint32_t p = 0;
    for (int b = (int)NB - 1; b >= 0; --b) {
        const uint32_t end = (b > 0) ? off[b - 1] : (uint32_t)M;

        for (; p < end; ++p) {
            const uint32_t v = vox[p];
            const int k = (int)((size_t)v / plane);
            const int y = (int)(((size_t)v % plane) / (size_t)nx);
            const int x = (int)((size_t)v % (size_t)nx);

            c.parent[p] = p;
            c.size[p]   = 1;
            c.min_k[p]  = (uint16_t)k;
            c.max_k[p]  = (uint16_t)k;
            c.flags[p]  = (uint8_t)(((k == nz - 1) ? CHANNEL_FLAG_TOP : 0) | ((k == 0) ? CHANNEL_FLAG_BOTTOM : 0));
            c.num_roots += 1;
            perc_enter(&c, p);

            bool connected = (!out->has_r_c && c.flags[p] == CHANNEL_FLAG_BOTH);

            // Six neighbours. An id below p is a voxel already inserted, since the list is in
            // insertion order - no separate "visited" mark is needed or kept.
            //
            // Both wrap directions are probed even though either alone would do: a seam pair is the
            // same pair seen from both sides, and whichever voxel is inserted second makes the
            // union. Keeping the loop symmetric costs two comparisons and means the periodic case
            // reads the same as the interior one.
            uint32_t nb[6];
            int      nn = 0;
            if (x > 0)               nb[nn++] = cid[(size_t)v - 1];
            else if (field->pbc[0] && nx > 1) nb[nn++] = cid[(size_t)v + (nx - 1)];
            if (x < nx - 1)          nb[nn++] = cid[(size_t)v + 1];
            else if (field->pbc[0] && nx > 1) nb[nn++] = cid[(size_t)v - (nx - 1)];
            if (y > 0)               nb[nn++] = cid[(size_t)v - nx];
            else if (field->pbc[1] && ny > 1) nb[nn++] = cid[(size_t)v + (size_t)(ny - 1) * nx];
            if (y < ny - 1)          nb[nn++] = cid[(size_t)v + nx];
            else if (field->pbc[1] && ny > 1) nb[nn++] = cid[(size_t)v - (size_t)(ny - 1) * nx];
            if (k > 0)               nb[nn++] = cid[(size_t)v - plane];
            if (k < nz - 1)          nb[nn++] = cid[(size_t)v + plane];

            for (int q = 0; q < nn; ++q) {
                if (nb[q] == PERC_INVALID || nb[q] >= p) continue;
                if (perc_union(&c, p, nb[q])) connected = true;
            }

            if (connected && !out->has_r_c) {
                // The voxel whose insertion joined the two faces is the tightest point of the widest
                // route, and its clearance is the critical radius. Bisecting on "does anything span"
                // brackets this number and never tells you which voxel it was.
                out->has_r_c   = true;
                out->r_c       = (double)field->data[v];
                out->throat[0] = field->origin[0] + ((float)x + 0.5f) * field->spacing[0];
                out->throat[1] = field->origin[1] + ((float)y + 0.5f) * field->spacing[1];
                out->throat[2] = field->origin[2] + ((float)k + 0.5f) * field->spacing[2];
            }
        }

        if ((uint32_t)b % per == 0) {
            const uint32_t s = (uint32_t)b / per;
            out->radius[s]         = (double)rmin + (double)b * width;
            out->frac_void[s]      = (double)p * inv_N;
            out->frac_open[s]      = (double)c.open_voxels * inv_N;
            out->frac_spanning[s]  = (double)c.span_voxels * inv_N;
            out->num_components[s] = c.num_roots;
            out->num_spanning[s]   = c.num_spanning;
            out->z_from_top[s]     = (c.deepest_top == INT32_MAX) ? z_hi
                                   : (double)field->origin[2] + ((double)c.deepest_top + 0.5) * (double)field->spacing[2];
            out->z_from_bottom[s]  = (c.highest_bot == INT32_MIN) ? z_lo
                                   : (double)field->origin[2] + ((double)c.highest_bot + 0.5) * (double)field->spacing[2];
        }
    }

    md_free(alloc, cid, N * sizeof(uint32_t));
    md_free(alloc, vox, M * sizeof(uint32_t));
    md_free(alloc, off, (size_t)NB * sizeof(uint32_t));
    md_free(alloc, cur, (size_t)NB * sizeof(uint32_t));
    md_free(alloc, c.parent, M * sizeof(uint32_t));
    md_free(alloc, c.size,   M * sizeof(uint32_t));
    md_free(alloc, c.flags,  M * sizeof(uint8_t));
    md_free(alloc, c.min_k,  M * sizeof(uint16_t));
    md_free(alloc, c.max_k,  M * sizeof(uint16_t));

    return true;
}

// ---------------------------------------------------------------------------------------------
// Tracing a route
// ---------------------------------------------------------------------------------------------

namespace {

const int PATH_DIR[6][3] = { {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1} };
const uint8_t PATH_SEED  = 7;

// Voxel index of an integer lattice point, wrapping x and y where the field is periodic. Returns
// false for a point off a non periodic face.
inline bool path_voxel(const channel_field_t* f, long ix, long iy, long ik, size_t* out) {
    const long nx = f->dim[0], ny = f->dim[1], nz = f->dim[2];
    if (f->pbc[0]) { ix = ((ix % nx) + nx) % nx; } else if (ix < 0 || ix >= nx) return false;
    if (f->pbc[1]) { iy = ((iy % ny) + ny) % ny; } else if (iy < 0 || iy >= ny) return false;
    if (ik < 0 || ik >= nz) return false;
    *out = ((size_t)ik * (size_t)ny + (size_t)iy) * (size_t)nx + (size_t)ix;
    return true;
}

// Does the straight segment between two lattice points stay inside {d >= r}? Sampled at half a
// voxel, which is the finest statement the field itself supports.
bool path_segment_clear(const channel_field_t* f, float r, const double a[3], const double b[3]) {
    const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
    const double len = sqrt(dx * dx + dy * dy + dz * dz);
    const int steps = (int)(2.0 * len) + 1;
    for (int i = 0; i <= steps; ++i) {
        const double t = (double)i / (double)steps;
        size_t v;
        if (!path_voxel(f, lround(a[0] + t * dx), lround(a[1] + t * dy), lround(a[2] + t * dz), &v)) return false;
        if (f->data[v] < r) return false;
    }
    return true;
}

// One sample of the route: the wrapped world position, with the clearance the field actually
// reports there rather than an interpolation between two distant anchors.
void path_emit(md_array(vec4_t)* out, const channel_field_t* f, double sx, double sy, double sz, struct md_allocator_i* alloc) {
    size_t v = 0;
    if (!path_voxel(f, lround(sx), lround(sy), lround(sz), &v)) return;
    double wx = sx, wy = sy;
    if (f->pbc[0]) { wx = fmod(sx, (double)f->dim[0]); if (wx < 0.0) wx += (double)f->dim[0]; }
    if (f->pbc[1]) { wy = fmod(sy, (double)f->dim[1]); if (wy < 0.0) wy += (double)f->dim[1]; }
    vec4_t pt;
    pt.x = f->origin[0] + (float)(wx + 0.5) * f->spacing[0];
    pt.y = f->origin[1] + (float)(wy + 0.5) * f->spacing[1];
    pt.z = f->origin[2] + (float)(sz + 0.5) * f->spacing[2];
    pt.w = f->data[v];
    md_array_push(*out, pt, alloc);
}

}  // namespace

bool channel_trace_path(md_array(vec4_t)* out_path, double* out_length, const channel_field_t* field, double r, struct md_allocator_i* alloc) {
    ASSERT(out_path);
    ASSERT(field);
    ASSERT(alloc);

    if (!field->data) return false;
    if (field->pbc[2]) return false;
    const int nx = field->dim[0], ny = field->dim[1], nz = field->dim[2];
    if (nx <= 0 || ny <= 0 || nz <= 0) return false;

    const size_t plane = (size_t)nx * (size_t)ny;
    const size_t N     = plane * (size_t)nz;
    const float  rr    = (float)r;

    // One byte per voxel: which way we arrived. A parent index would be four, and the direction is
    // all a backtrace needs.
    uint8_t* from = (uint8_t*)md_alloc(alloc, N);
    if (!from) return false;
    MEMSET(from, 0, N);

    md_array(uint32_t) queue = 0;
    size_t head = 0;
    size_t found = N;

    for (size_t i = 0; i < plane; ++i) {
        const size_t v = (size_t)(nz - 1) * plane + i;
        if (field->data[v] < rr) continue;
        from[v] = PATH_SEED;
        if (nz == 1) { found = v; break; }
        md_array_push(queue, (uint32_t)v, alloc);
    }

    while (found == N && head < md_array_size(queue)) {
        const size_t v = (size_t)queue[head++];
        const int k = (int)(v / plane);
        const int y = (int)((v % plane) / (size_t)nx);
        const int x = (int)(v % (size_t)nx);

        for (int d = 0; d < 6; ++d) {
            size_t w;
            if (!path_voxel(field, (long)x + PATH_DIR[d][0], (long)y + PATH_DIR[d][1], (long)k + PATH_DIR[d][2], &w)) continue;
            if (from[w] || field->data[w] < rr) continue;
            from[w] = (uint8_t)(d + 1);
            if (w < plane) { found = w; break; }        // Reached the bottom face
            md_array_push(queue, (uint32_t)w, alloc);
        }
    }

    md_array_free(queue, alloc);

    if (found == N) {
        md_free(alloc, from, N);
        return false;
    }

    // Backtrace, then reverse, carrying unwrapped lattice coordinates so a route that leaves through
    // a periodic face keeps going in a straight line instead of jumping the box.
    md_array(int32_t) chain = 0;                        // Directions, bottom to top
    {
        size_t v = found;
        while (from[v] != PATH_SEED) {
            const int d = (int)from[v] - 1;
            md_array_push(chain, (int32_t)d, alloc);
            size_t w;
            const int k = (int)(v / plane);
            const int y = (int)((v % plane) / (size_t)nx);
            const int x = (int)(v % (size_t)nx);
            if (!path_voxel(field, (long)x - PATH_DIR[d][0], (long)y - PATH_DIR[d][1], (long)k - PATH_DIR[d][2], &w)) break;
            v = w;
        }
        // v is now the seed at the top face; rebuild the route downward from it
        md_array(double) pts = 0;
        long ix = (long)(v % (size_t)nx);
        long iy = (long)((v % plane) / (size_t)nx);
        long ik = (long)(v / plane);
        md_array_push(pts, (double)ix, alloc);
        md_array_push(pts, (double)iy, alloc);
        md_array_push(pts, (double)ik, alloc);
        for (size_t i = md_array_size(chain); i > 0; --i) {
            const int d = (int)chain[i - 1];
            ix += PATH_DIR[d][0];
            iy += PATH_DIR[d][1];
            ik += PATH_DIR[d][2];
            md_array_push(pts, (double)ix, alloc);
            md_array_push(pts, (double)iy, alloc);
            md_array_push(pts, (double)ik, alloc);
        }
        md_array_free(chain, alloc);
        md_free(alloc, from, N);

        const size_t n = md_array_size(pts) / 3;

        // Straighten it. Breadth first through a six connected grid returns a staircase: it is the
        // shortest route in voxel steps, which overstates the length of anything not axis aligned by
        // up to sqrt(3). Replacing a run by the straight segment between its ends, where that segment
        // stays inside the set, recovers a length worth quoting and a route worth looking at.
        const size_t LOOKAHEAD = 64;
        md_array(double) anchor = 0;
        for (size_t i = 0; ; ) {
            md_array_push(anchor, pts[i*3+0], alloc);
            md_array_push(anchor, pts[i*3+1], alloc);
            md_array_push(anchor, pts[i*3+2], alloc);
            if (i + 1 >= n) break;

            size_t best = i + 1;
            const size_t far = MIN(n - 1, i + LOOKAHEAD);
            for (size_t j = far; j > i + 1; --j) {
                const double a[3] = { pts[i*3+0], pts[i*3+1], pts[i*3+2] };
                const double b[3] = { pts[j*3+0], pts[j*3+1], pts[j*3+2] };
                if (path_segment_clear(field, rr, a, b)) { best = j; break; }
            }
            i = best;
        }
        md_array_free(pts, alloc);

        const size_t na = md_array_size(anchor) / 3;

        double length = 0.0;
        for (size_t i = 1; i < na; ++i) {
            const double dx = (anchor[i*3+0] - anchor[(i-1)*3+0]) * (double)field->spacing[0];
            const double dy = (anchor[i*3+1] - anchor[(i-1)*3+1]) * (double)field->spacing[1];
            const double dz = (anchor[i*3+2] - anchor[(i-1)*3+2]) * (double)field->spacing[2];
            length += sqrt(dx * dx + dy * dy + dz * dz);
        }
        if (out_length) *out_length = length;

        // Resample the straightened route at half a voxel, reading the clearance at each sample.
        //
        // The anchors alone are not something to draw with. Straightening is free to leave them
        // sixty voxels apart, and the clearance between two of them is whatever the field says, not
        // the interpolation of its endpoints - so anything drawn from the anchors is at its widest
        // exactly where its width was never measured. Sampling also keeps every step short, which is
        // what lets a caller drop the one segment that crosses a periodic face without leaving a
        // visible gap in the rest.
        md_array_shrink(*out_path, 0);
        for (size_t i = 0; i + 1 < na; ++i) {
            const double ax = anchor[i*3+0],     ay = anchor[i*3+1],     az = anchor[i*3+2];
            const double dx = anchor[(i+1)*3+0] - ax;
            const double dy = anchor[(i+1)*3+1] - ay;
            const double dz = anchor[(i+1)*3+2] - az;
            const int steps = MAX(1, (int)(2.0 * sqrt(dx * dx + dy * dy + dz * dz)));
            for (int t = 0; t < steps; ++t) {
                const double u = (double)t / (double)steps;
                path_emit(out_path, field, ax + u * dx, ay + u * dy, az + u * dz, alloc);
            }
        }
        path_emit(out_path, field, anchor[(na-1)*3+0], anchor[(na-1)*3+1], anchor[(na-1)*3+2], alloc);
        md_array_free(anchor, alloc);
    }

    return md_array_size(*out_path) > 1;
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
