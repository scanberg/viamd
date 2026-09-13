#pragma once

// Channels through a structure along the sweep axis.
//
// Given a clearance field - for every voxel, the radius of the largest probe sphere centred there -
// a probe of radius r can occupy exactly the superlevel set {d >= r}. A channel through the
// structure is a connected component of that set which touches both faces of the sweep axis.
//
// Sweeping a plane from the far face towards the near one, the components of {d >= r} restricted to
// the part already swept can only appear or merge, never split, because the swept region only ever
// grows. That makes the sweep history a forest: a genuine tree rooted at the near face, whose
// branches are stretches of channel and whose leaves are where a branch first appears. It is the
// structure behind both the count and the diagram.
//
// The sweep axis is z (axis 2) and must be non periodic. x and y may be periodic.

#include <core/md_array.h>
#include <core/md_vec_math.h>

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

struct md_allocator_i;

#define CHANNEL_INVALID_INDEX 0xFFFFFFFFu

// A clearance field on a regular grid, indexed x fastest.
typedef struct channel_field_t {
    const float* data;
    int   dim[3];
    float spacing[3];
    float origin[3];
    bool  pbc[3];       // Periodicity per axis. pbc[2] must be false.
} channel_field_t;

// One branch of the sweep forest: a component over the z range in which it had its own identity,
// from where it appeared down to where it merged into its parent.
typedef struct channel_node_t {
    uint32_t parent;            // Merge node below this branch, CHANNEL_INVALID_INDEX for a root
    uint32_t first_child;
    uint32_t next_sibling;

    float    z_top;             // Where the branch appears
    float    z_bot;             // Where it merges into its parent, or the bottom of the grid
    float    max_clearance;     // Largest clearance anywhere in the branch
    float    min_clearance;     // Smallest clearance along its own centreline, i.e. its tightest point
    uint32_t num_voxels;

    bool     reaches_top;       // This branch or one below it touches the far face
    bool     reaches_bottom;    // Only a root can, and only if its component touches the near face

    // Representative centreline: one point per slab the branch spans, from its top down. xyz is the
    // world position of the widest voxel of the branch in that slab, w is the clearance there. It is
    // a representative route, not the exact widest geodesic - see channel_critical_radius for the
    // number which is exact.
    md_array(vec4_t) path;
} channel_node_t;

typedef struct channel_tree_t {
    md_array(channel_node_t) nodes;
    md_array(uint32_t) roots;

    uint32_t num_components;    // Connected components of {d >= r} in the whole grid
    uint32_t num_spanning;      // Of those, the ones touching both faces: the channels
    double   probe_radius;

    struct md_allocator_i* alloc;
} channel_tree_t;

// Percolation and reachability, from one pass in order of decreasing clearance.
//
// Insert voxels into a union-find from the widest outward. At the moment some component first
// touches both z faces, the clearance of the voxel just inserted IS the critical radius - it is the
// tightest point of the widest route, which is what a maximum-capacity path gives and what bisecting
// on "does anything span" only brackets. That voxel is also *where* the route is limited, which no
// amount of bisection recovers: a scalar r_c says a probe of that size gets through and nothing
// about what stops a larger one.
//
// Everything else falls out of the same order, at no extra cost:
//
// - A component touching a z face is reachable from outside, so open and closed porosity separate
//   without the second flood fill the design document still lists as owed. A cavity that fits the
//   probe but connects to nothing is counted by the accessible volume and not by this.
// - The lowest z any top-touching component reaches is how far a probe of that radius penetrates
//   from above, and the same from below. They approach each other as the probe shrinks and meet at
//   r_c, which is the percolation threshold arrived at from the other direction - and unlike a
//   single number, the pair of curves says how much of the film a probe too large to cross can
//   still get into.
//
// Counting connected components, which is what channel_sweep reports, is nearly useless on a real
// network: below r_c the void space is one component and the count is 1, above it 0. The curves
// here are what that step function was standing in for.
typedef struct channel_percolation_t {
    md_array(double)   radius;          // Ascending, from r_min to the widest clearance in the field
    md_array(double)   frac_void;       // Vol{d >= r} / V, the whole box
    md_array(double)   frac_open;       // Of the box: in a component which touches a z face
    md_array(double)   frac_spanning;   // Of the box: in a component which touches both
    md_array(double)   z_from_top;      // Lowest z reached from the top face; the top of the grid if none
    md_array(double)   z_from_bottom;   // Highest z reached from the bottom face
    md_array(uint32_t) num_components;
    md_array(uint32_t) num_spanning;

    bool   has_r_c;
    double r_c;                         // Clearance at which the two faces first connect
    float  throat[3];                   // Where: the voxel whose insertion connected them

    size_t num_active;                  // Voxels at or above r_min, which is what the run was sized by
    size_t bytes;                       // What it allocated, so the cost of a larger field is legible

    struct md_allocator_i* alloc;
} channel_percolation_t;

#ifdef __cplusplus
extern "C" {
#endif

void channel_tree_free(channel_tree_t* tree);
void channel_percolation_free(channel_percolation_t* perc);

// One pass. r_min bounds the memory - voxels below it never enter the structure - and is therefore
// also the lowest radius the curves reach. num_samples is the number of points on them.
//
// Costs 4 bytes per voxel of the field plus 17 per voxel at or above r_min, against 9 per active
// voxel for a single channel_sweep. It replaces a bisection of about a dozen of those plus one sweep
// per reported radius - measured at 6.5x faster than that pair on a 7.7 Mvoxel disordered network,
// and exact where they were bracketed.
bool channel_percolate(channel_percolation_t* out, const channel_field_t* field, double r_min, uint32_t num_samples, struct md_allocator_i* alloc);

// A route a probe of radius r can actually follow from the top face to the bottom: breadth first
// through {d >= r}, so it is the shortest such route in voxel steps, then straightened wherever the
// segment between two of its points stays inside the set.
//
// Unlike the representative centreline on a tree branch - the widest voxel per slab, which need not
// be connected to the widest voxel of the slab below - every point of this is reachable from the
// previous one at radius r. Trace it at r_c and it is the route the critical radius belongs to.
//
// out_path is xyz in world space with w the clearance there, wrapped into the box on a periodic
// axis. out_length is the length of the unwrapped route, so length / |z span| is the tortuosity.
// Returns false when nothing gets through at r.
bool channel_trace_path(md_array(vec4_t)* out_path, double* out_length, const channel_field_t* field, double r, struct md_allocator_i* alloc);

// Sweep the field at one probe radius. With want_tree false only the counts are produced, which is
// what the radius sweeps below use.
bool channel_sweep(channel_tree_t* out, const channel_field_t* field, double probe_radius, bool want_tree, struct md_allocator_i* alloc);

// Number of channels at each supplied radius.
void channel_spanning_counts(uint32_t* out_counts, const double* radii, size_t num_radii, const channel_field_t* field, struct md_allocator_i* alloc);

// Assign a dendrogram column to every branch belonging to a channel which gets all the way through.
// Leaves take consecutive columns in traversal order and a merge node sits at the mean of its
// children, so a branch is drawn above the span of what feeds it. Branches of any other component
// are left at -1: a blind pore would only crowd the diagram.
//
// out_slot must have room for md_array_size(tree->nodes) entries. Returns the number of columns used.
//
// The traversal carries its own stack. A merge tree has one node per pair of branches that join, so
// a real network produces a chain thousands of levels deep and recursing over it overflows.
float channel_tree_layout(float* out_slot, const channel_tree_t* tree, struct md_allocator_i* alloc);

// The largest probe radius which still gets through: bisected until the bracket is below tol.
//
// channel_percolate answers this exactly and in one pass, and is what the component uses. This stays
// because it asks a completely different question of the field - a threshold test repeated, rather
// than an ordering - so the two agreeing on r_c is evidence about the answer rather than about one
// implementation. The tests cross check them.
// This is a property of the tightest throat along the best route, not of any cavity along it, and it
// is the number to quote - the branch clearances above are read off a representative centreline.
// Returns a value below r_lo when nothing gets through even at r_lo.
double channel_critical_radius(const channel_field_t* field, double r_lo, double r_hi, double tol, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
