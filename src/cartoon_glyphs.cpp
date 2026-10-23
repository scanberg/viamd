#include <cartoon_glyphs.h>

#include <md_system.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_log.h>

#include <math.h>
#include <float.h>

// Radii of the spheres and capsules, Å
static constexpr float HEAD_RADIUS_SCALE        = 1.6f;   // A head of n heavy atoms is 1.6 cbrt(n): 4.3 for the 20 of a phosphatidyl choline
static constexpr float HEAD_RADIUS_MIN          = 1.8f;
static constexpr float HEAD_RADIUS_MAX          = 7.0f;
static constexpr float HEAVY_ATOMS_PER_BEAD     = 4.0f;   // What a Martini bead stands for
static constexpr float TAIL_RADIUS              = 1.1f;
// The ring system of a sterol, a slab: as thick as a chain, about twice as wide in its profile, its sections along it
// about a ring apart. Without its rings to lay it out (no reference positions, fewer than three ring atoms) a round tube.
static constexpr float STEROL_HALF_THICKNESS    = 1.0f;
static constexpr float STEROL_HALF_WIDTH        = 2.0f;
static constexpr float STEROL_SECTION           = 3.2f;
static constexpr float SKELETON_RADIUS          = 1.5f;
// The backbone of a nucleic acid: a flat band 2 (r + l) wide and 2 r thick, its faces towards the bases. The r is that of
// a representation (CARTOON_NUCLEIC_BAND_RADIUS is the default, which it scales), l stays.
static constexpr float NUCLEIC_BACKBONE_RADIUS  = CARTOON_NUCLEIC_BAND_RADIUS;
static constexpr float NUCLEIC_BACKBONE_FLAT    = 1.6f;

static inline bool atom_is_bead(const md_system_t& sys, size_t i) {
    return md_atom_particle_kind(&sys.atom, i) != MD_PARTICLE_ATOM && md_atom_atomic_number(&sys.atom, i) == 0;
}

// The residue of an atom
static int32_t residue_of_atom(const md_system_t& sys, int32_t atom) {
    return (int32_t)md_system_residue_find_by_atom_idx(&sys, (size_t)atom);
}

// The phosphorus of a phosphate, or a Martini phosphate bead: PO and digits (PO4, PO1, PO41), not the P1 .. P5 of the
// phosphates on an inositol, which belong to the head group
static bool is_phosphate_center(const md_system_t& sys, int32_t atom) {
    const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, atom);
    if (z != 0) return z == 15;
    const str_t name = md_atom_name(&sys.atom, atom);
    if (name.len < 3 || name.ptr[0] != 'P' || name.ptr[1] != 'O') return false;
    for (size_t i = 2; i < name.len; ++i) {
        if (name.ptr[i] < '0' || name.ptr[i] > '9') return false;
    }
    return true;
}

// The atom a head takes its colour from: its heaviest element (the phosphorus of a phospholipid, the oxygen of a
// sterol), or for beads the phosphate bead if there is one
static int32_t head_color_atom(const md_system_t& sys, const int32_t* atoms, size_t n) {
    int32_t best = atoms[0];
    md_atomic_number_t best_z = md_atom_atomic_number(&sys.atom, atoms[0]);
    for (size_t k = 0; k < n; ++k) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, atoms[k]);
        if (z > best_z) {
            best = atoms[k];
            best_z = z;
        }
        if (best_z == 0 && is_phosphate_center(sys, atoms[k])) best = atoms[k];
    }
    return best;
}

// The radius of the sphere of a head (or a part of one) of n heavy atoms or beads
static float head_radius(const md_system_t& sys, const int32_t* atoms, size_t n) {
    const float heavy = (float)n * (atom_is_bead(sys, atoms[0]) ? HEAVY_ATOMS_PER_BEAD : 1.0f);
    return CLAMP(HEAD_RADIUS_SCALE * cbrtf(heavy), HEAD_RADIUS_MIN, HEAD_RADIUS_MAX);
}

// A head split at its phosphates into the head group beyond them, the phosphates with it (a phosphocholine, a
// phosphoserine, the bridging phosphoglycerol of a cardiolipin), and the backbone the tails hang off (a glycerol and its
// esters, the polar part of a sphingosine): a group per connected part of either. The Martini beads of a phospholipid
// split the same way, NC3 PO4 | GL1 GL2.
struct HeadSplit {
    md_array(int32_t)  atoms = nullptr;     // The atoms of the groups one after the other, ascending within each
    md_array(uint32_t) offset = nullptr;    // [count + 1]
    md_array(uint8_t)  outer = nullptr;     // [count] Of the head group, else of the backbone
};

static inline uint32_t head_split_count(const HeadSplit& split) {
    return split.offset ? (uint32_t)md_array_size(split.offset) - 1 : 0;
}

// The local index of an atom of the head (its atoms ascend), -1 for none
static int32_t head_local(const int32_t* atoms, size_t n, int32_t atom) {
    size_t lo = 0, hi = n;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (atoms[mid] < atom) lo = mid + 1;
        else hi = mid;
    }
    return (lo < n && atoms[lo] == atom) ? (int32_t)lo : -1;
}

// Labels the connected parts of the atoms of the head that keep (keep[i]), through bonds between such atoms. Returns
// their count; label[i] is -1 for those that do not keep.
static int32_t head_label_parts(int32_t* label, const uint8_t* keep, const md_system_t& sys, const int32_t* atoms, size_t n, int32_t* stack) {
    for (size_t i = 0; i < n; ++i) label[i] = -1;
    int32_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!keep[i] || label[i] != -1) continue;
        size_t top = 0;
        stack[top++] = (int32_t)i;
        label[i] = count;
        while (top > 0) {
            const int32_t a = stack[--top];
            for (md_bond_iter_t it = md_bond_iter(&sys.bond, atoms[a]); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
                const int32_t b = head_local(atoms, n, md_bond_iter_atom_index(&it));
                if (b >= 0 && keep[b] && label[b] == -1) {
                    label[b] = count;
                    stack[top++] = b;
                }
            }
        }
        ++count;
    }
    return count;
}

// Splits the head (its atoms ascend) whose tails hang off the anchors. No split for a head without a phosphate, with
// nothing beyond it or nothing before it, or with a tail that hangs off a phosphate itself. The split and the scratch
// space come from alloc, a temporary allocator (a scope of its own here would take back what the split grows into).
static bool split_head(HeadSplit* out, const md_system_t& sys, const int32_t* atoms, size_t n, const int32_t* anchors, size_t num_anchors, md_allocator_i* alloc) {
    if (n < 2 || num_anchors == 0) return false;
    uint8_t* center = (uint8_t*)md_alloc(alloc, n);
    uint8_t* keep   = (uint8_t*)md_alloc(alloc, n);
    int32_t* label  = (int32_t*)md_alloc(alloc, n * sizeof(int32_t));
    int32_t* outer_label = (int32_t*)md_alloc(alloc, n * sizeof(int32_t));
    int32_t* stack  = (int32_t*)md_alloc(alloc, n * sizeof(int32_t));
    uint8_t* anchored = (uint8_t*)md_alloc(alloc, n);
    bool ok = false;
    bool any_center = false;
    for (size_t i = 0; i < n; ++i) {
        center[i] = is_phosphate_center(sys, atoms[i]) ? 1 : 0;
        keep[i] = !center[i];
        any_center |= center[i] != 0;
    }
    int32_t num_inner = 0, num_outer = 0;
    if (any_center) {
        // The parts without the phosphates: those the tails hang off are the backbone, the rest and the phosphates the head group
        const int32_t num_parts = head_label_parts(label, keep, sys, atoms, n, stack);
        MEMSET(anchored, 0, n);
        bool anchors_ok = true;
        for (size_t k = 0; k < num_anchors; ++k) {
            const int32_t a = head_local(atoms, n, anchors[k]);
            if (a < 0 || center[a]) {
                anchors_ok = false;
                break;
            }
            anchored[label[a]] = 1;
        }
        if (anchors_ok) {
            for (size_t i = 0; i < n; ++i) keep[i] = center[i] || !anchored[label[i]];
            num_outer = head_label_parts(outer_label, keep, sys, atoms, n, stack);
            for (int32_t p = 0; p < num_parts; ++p) num_inner += anchored[p] ? 1 : 0;
            ok = num_outer > 0 && num_inner > 0;

            if (ok) {
                md_array_push(out->offset, 0, alloc);
                for (int32_t p = 0; p < num_parts; ++p) {
                    if (!anchored[p]) continue;
                    for (size_t i = 0; i < n; ++i) {
                        if (!center[i] && label[i] == p) md_array_push(out->atoms, atoms[i], alloc);
                    }
                    md_array_push(out->offset, (uint32_t)md_array_size(out->atoms), alloc);
                    md_array_push(out->outer, 0, alloc);
                }
                for (int32_t p = 0; p < num_outer; ++p) {
                    for (size_t i = 0; i < n; ++i) {
                        if (outer_label[i] == p) md_array_push(out->atoms, atoms[i], alloc);
                    }
                    md_array_push(out->offset, (uint32_t)md_array_size(out->atoms), alloc);
                    md_array_push(out->outer, 1, alloc);
                }
            }
        }
    }
    return ok;
}


struct GlyphBuilder {
    CartoonGlyphs* g;
    md_allocator_i* alloc;

    // The mean of the atoms, or their weighted sum with the weights (n of them)
    uint32_t point(const int32_t* atoms, size_t n, int32_t ref, int32_t color_atom, int32_t res, float radius, CartoonGlyphKind kind, CartoonGlyphVariant variant = CartoonGlyphVariant_Always, const float* weights = nullptr) {
        const uint32_t idx = (uint32_t)md_array_size(g->point_ref);
        md_array_push(g->point_variant, (uint8_t)variant, alloc);
        md_array_push_array(g->point_atom, atoms, n, alloc);
        for (size_t k = 0; k < n; ++k) md_array_push(g->point_weight, weights ? weights[k] : 1.0f / (float)n, alloc);
        md_array_push(g->point_offset, (uint32_t)md_array_size(g->point_atom), alloc);
        md_array_push(g->point_ref, ref, alloc);
        md_array_push(g->point_color_atom, color_atom, alloc);
        md_array_push(g->point_res, res, alloc);
        md_array_push(g->point_radius, radius, alloc);
        md_array_push(g->point_kind, (uint8_t)kind, alloc);
        return idx;
    }

    // The next control point of the current tube, which tube_end closes
    void tube_point(uint32_t p, float radius) {
        md_array_push(g->tube_point, p, alloc);
        md_array_push(g->tube_radius, radius, alloc);
        md_array_push(g->tube_flat, 0.0f, alloc);
        md_array_push(g->tube_face, UINT32_MAX, alloc);
        md_array_push(g->tube_side, UINT32_MAX, alloc);
    }

    // flat: half the length of the flat sides of the profile, face: the point they face, side: the point the long axis
    // of the profile points at (UINT32_MAX for none)
    void tube_point(uint32_t p, float radius, float flat, uint32_t face, uint32_t side = UINT32_MAX) {
        tube_point(p, radius);
        g->tube_flat[md_array_size(g->tube_flat) - 1] = flat;
        g->tube_face[md_array_size(g->tube_face) - 1] = face;
        g->tube_side[md_array_size(g->tube_side) - 1] = side;
    }

    // through: the tube passes through its control points (a nucleic backbone through its C4', which its sugars hang
    // off), rather than near them (a lipid tail, smoothed)
    void tube_end(bool through = false) {
        md_array_push(g->tube_offset, (uint32_t)md_array_size(g->tube_point), alloc);
        md_array_push(g->tube_through, (uint8_t)(through ? 1 : 0), alloc);
    }

    void slab(const uint32_t* corners, size_t n) {
        md_array_push_array(g->slab_point, corners, n, alloc);
        md_array_push(g->slab_offset, (uint32_t)md_array_size(g->slab_point), alloc);
    }
};

// Where a tube starts at a head: the point of the head, or of the part of the backbone of a split head the tail hangs
// off. It is a point of the tube, never hidden by the number of head spheres.
struct TubeRoot {
    uint32_t head;
    int32_t  anchor;
    uint32_t point;
    uint32_t size;      // Atoms of the head or the part
};

// The ring system of a sterol as a slab: a tube through the middle of its rings in sections along it, flat in their plane
// (STEROL_HALF_THICKNESS, STEROL_HALF_WIDTH), from the centre of its head and on along its side chain as a chain is. The
// rings are the ring atoms of the skeleton, or its ring beads (R1 .. R6) for a coarse grained sterol, whose bonds need not
// make rings. They are laid out once, from the reference positions of the system: the long axis and the width of the
// ring system are the first two principal axes of their positions; a section is a share of the ring atoms on either
// side of the long axis, about STEROL_SECTION along it, and its profile lies towards those on the one side. The ring
// system is rigid, so what its reference positions lay out holds in every frame. False where the
// rings cannot lay it out (fewer than three ring atoms, no reference positions): a round tube then.
static bool sterol_slab(GlyphBuilder& b, const md_system_t& sys, const md_lipid_data_t& L, uint32_t t, uint32_t head_point, int32_t ref,
                        const uint8_t* in_ring, md_allocator_i* temp_alloc) {
    const int32_t* atoms = md_index_range_beg(&L.tail.atoms, t);
    const size_t num_atoms = md_index_range_size(&L.tail.atoms, t);
    const int32_t* spine = md_index_range_beg(&L.tail.spine, t);
    const size_t len = md_index_range_size(&L.tail.spine, t);
    const md_system_state_t& refs = sys.reference;
    if (!refs.xyz || refs.num_atoms != sys.atom.count || len == 0) return false;

    auto is_ring = [&](int32_t i) -> bool {
        if (atom_is_bead(sys, i)) {
            const str_t name = md_atom_name(&sys.atom, i);
            return name.len >= 2 && name.ptr[0] == 'R' && name.ptr[1] >= '0' && name.ptr[1] <= '9';
        }
        return in_ring && in_ring[i] && md_atom_atomic_number(&sys.atom, i) != MD_Z_H;
    };
    md_array(int32_t) ring = nullptr;
    for (size_t k = 0; k < num_atoms; ++k) {
        if (is_ring(atoms[k])) md_array_push(ring, atoms[k], temp_alloc);
    }
    const size_t n = md_array_size(ring);
    if (n < 3) return false;

    // The positions about the anchor, each in its image nearest to it
    const int32_t anchor = L.tail.anchor_beg[t];
    vec3_t* d = (vec3_t*)md_alloc(temp_alloc, n * sizeof(vec3_t));
    for (size_t k = 0; k < n; ++k) d[k] = vec3_sub(refs.xyz[ring[k]], refs.xyz[anchor]);
    if (md_unitcell_flags(&refs.unitcell)) md_util_min_image_vec3(d, n, &refs.unitcell);
    vec3_t c = vec3_zero();
    for (size_t k = 0; k < n; ++k) c = vec3_add(c, d[k]);
    c = vec3_mul1(c, 1.0f / (float)n);
    mat3_t M = {};
    for (size_t k = 0; k < n; ++k) {
        const vec3_t x = vec3_sub(d[k], c);
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M.elem[i][j] += x.elem[i] * x.elem[j];
    }
    const mat3_eigen_t eig = mat3_eigen(M);
    vec3_t axis = eig.vectors.col[0];
    const vec3_t width = eig.vectors.col[1];
    // From the head: the anchor (the origin) at the low end
    if (vec3_dot(c, axis) < 0.0f) axis = vec3_mul1(axis, -1.0f);

    float* along = (float*)md_alloc(temp_alloc, n * sizeof(float));
    float* across = (float*)md_alloc(temp_alloc, n * sizeof(float));
    float lo = FLT_MAX, hi = -FLT_MAX;
    for (size_t k = 0; k < n; ++k) {
        const vec3_t x = vec3_sub(d[k], c);
        along[k] = vec3_dot(x, axis);
        across[k] = vec3_dot(x, width);
        lo = MIN(lo, along[k]);
        hi = MAX(hi, along[k]);
    }

    // The ring atoms on either side of the axis, each side in order along it. A section takes a share of each side in
    // turn, so that its atoms on the one side lie across from those on the other: the middle of the section is on the
    // axis, and the direction from it to those on the one side is across the ring system, in its plane.
    md_array(int32_t) pos = nullptr;
    md_array(int32_t) neg = nullptr;
    for (size_t k = 0; k < n; ++k) {
        if (across[k] > 0.0f) md_array_push(pos, (int32_t)k, temp_alloc);
        else md_array_push(neg, (int32_t)k, temp_alloc);
    }
    auto sort_along = [&](int32_t* a, size_t m) {
        for (size_t i = 1; i < m; ++i) {
            const int32_t x = a[i];
            size_t j = i;
            while (j > 0 && along[a[j - 1]] > along[x]) {
                a[j] = a[j - 1];
                --j;
            }
            a[j] = x;
        }
    };
    const size_t num_pos = md_array_size(pos), num_neg = md_array_size(neg);
    if (num_pos == 0 || num_neg == 0) return false;
    sort_along(pos, num_pos);
    sort_along(neg, num_neg);
    const size_t num_sections = (size_t)CLAMP((int)((hi - lo) / STEROL_SECTION + 0.5f), 1, (int)MIN(num_pos, num_neg));

    const float flat = STEROL_HALF_WIDTH - STEROL_HALF_THICKNESS;
    md_array(int32_t) group = nullptr;
    md_array(uint32_t) points = nullptr;
    md_array(uint32_t) sides = nullptr;
    for (size_t sec = 0; sec < num_sections; ++sec) {
        const size_t pb = sec * num_pos / num_sections, pe = (sec + 1) * num_pos / num_sections;
        const size_t nb = sec * num_neg / num_sections, ne = (sec + 1) * num_neg / num_sections;
        md_array_shrink(group, 0);
        for (size_t i = pb; i < pe; ++i) md_array_push(group, ring[pos[i]], temp_alloc);
        const size_t m = md_array_size(group);
        for (size_t i = nb; i < ne; ++i) md_array_push(group, ring[neg[i]], temp_alloc);
        const size_t g = md_array_size(group);
        const uint32_t p = b.point(group, g, ref, group[0], residue_of_atom(sys, group[0]), 0.0f, CartoonGlyph_SterolSkeleton);
        const uint32_t side = b.point(group, m, ref, group[0], residue_of_atom(sys, group[0]), 0.0f, CartoonGlyph_SterolSkeleton);
        md_array_push(points, p, temp_alloc);
        md_array_push(sides, side, temp_alloc);
    }
    const uint32_t first_side = sides[0];
    const uint32_t last_side = sides[num_sections - 1];

    // The dummies of a coarse grained sterol, which carry the mass of its rings and do not interact (VOH V1 V2 V3 of
    // the optimal LINCS cholesterol of Martini 2, its ring beads virtual sites of them): of the ring system, which its
    // spine runs through to reach the side chain, though the slab is laid out from the ring beads alone
    auto is_dummy = [&](int32_t i) -> bool {
        return atom_is_bead(sys, i) && md_atom_name(&sys.atom, i).len >= 1 && md_atom_name(&sys.atom, i).ptr[0] == 'V';
    };

    // The side chain: the spine beyond its last ring atom (or dummy)
    size_t chain_beg = 0;
    for (size_t k = 0; k < len; ++k) {
        if (is_ring(spine[k]) || is_dummy(spine[k])) chain_beg = k + 1;
    }

    b.tube_point(head_point, STEROL_HALF_THICKNESS, 0.0f, UINT32_MAX, first_side);
    for (size_t i = 0; i < md_array_size(points); ++i) {
        b.tube_point(points[i], STEROL_HALF_THICKNESS, flat, UINT32_MAX, sides[i]);
    }
    const int32_t chain_res = residue_of_atom(sys, spine[len - 1]);
    size_t k = chain_beg;
    if (k < len && !atom_is_bead(sys, spine[k])) {
        for (; k + 1 < len; k += 2) {
            b.tube_point(b.point(spine + k, 2, ref, spine[k], chain_res, 0.0f, CartoonGlyph_SterolSkeleton), TAIL_RADIUS, 0.0f, UINT32_MAX, last_side);
        }
        k = len - 1;
    }
    for (; k < len; ++k) {
        b.tube_point(b.point(spine + k, 1, ref, spine[k], chain_res, 0.0f, CartoonGlyph_SterolSkeleton), TAIL_RADIUS, 0.0f, UINT32_MAX, last_side);
    }
    b.tube_end();
    return true;
}

static void build_lipids(GlyphBuilder& b, const md_system_t& sys) {
    const md_lipid_data_t& L = sys.lipid;
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    // The atoms of rings, which lay out the ring system of a sterol
    uint8_t* in_ring = nullptr;
    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    if (num_rings > 0) {
        in_ring = (uint8_t*)md_alloc(temp_alloc, sys.atom.count);
        MEMSET(in_ring, 0, sys.atom.count);
        for (size_t r = 0; r < num_rings; ++r) {
            const int32_t* ra = md_index_range_beg(&sys.ring, r);
            for (size_t k = 0; k < md_index_range_size(&sys.ring, r); ++k) in_ring[ra[k]] = 1;
        }
    }
    md_array(TubeRoot) roots = nullptr;     // Of the heads of the lipid
    md_array(int32_t) anchors = nullptr;
    md_array(int32_t) whole = nullptr;
    for (size_t l = 0; l < L.count; ++l) {
        const int32_t* residues = md_index_range_beg(&L.residues, l);
        const size_t num_residues = md_index_range_size(&L.residues, l);
        const uint32_t h0 = L.head_offset[l];
        const uint32_t h1 = L.head_offset[l + 1];
        const uint32_t t0 = L.tail_offset[l];
        const uint32_t t1 = L.tail_offset[l + 1];
        if (h1 == h0) {
            // A lipid whose heads and tails are not known (no bonds to find them by, a molecule of the lipidomes whose
            // beads tell no head from a tail: a hopanoid, a bile acid, a part of a lipid A on its own): one sphere for all
            // of it, at the centre of its heavy atoms, so that it is not missing from the cartoon
            md_array_shrink(whole, 0);
            for (size_t c = 0; c < num_residues; ++c) {
                const md_urange_t r = md_system_residue_atom_range(&sys, (size_t)residues[c]);
                for (uint32_t i = r.beg; i < r.end; ++i) {
                    if (md_atom_atomic_number(&sys.atom, i) != MD_Z_H) md_array_push(whole, (int32_t)i, temp_alloc);
                }
            }
            const size_t n = md_array_size(whole);
            if (n == 0) continue;
            const int32_t color = head_color_atom(sys, whole, n);
            b.point(whole, n, whole[0], color, residue_of_atom(sys, color), head_radius(sys, whole, n), CartoonGlyph_LipidHead);
            continue;
        }
        md_array_shrink(roots, 0);
        // Every point of the lipid in the image of its first head atom
        const int32_t ref = md_index_range_beg(&L.head.atoms, h0)[0];

        for (uint32_t h = h0; h < h1; ++h) {
            const int32_t* atoms = md_index_range_beg(&L.head.atoms, h);
            const size_t n = md_index_range_size(&L.head.atoms, h);
            const int32_t color = head_color_atom(sys, atoms, n);
            const int32_t res = residue_of_atom(sys, color);

            // The atoms of the head its tails hang off
            md_array_shrink(anchors, 0);
            for (uint32_t t = t0; t < t1; ++t) {
                if (L.tail.head_beg[t] == h) md_array_push(anchors, L.tail.anchor_beg[t], temp_alloc);
                if (L.tail.head_end[t] == h && L.tail.anchor_end[t] != -1) md_array_push(anchors, L.tail.anchor_end[t], temp_alloc);
            }

            HeadSplit split = {};
            if (!split_head(&split, sys, atoms, n, anchors, md_array_size(anchors), temp_alloc)) {
                const uint32_t p = b.point(atoms, n, ref, color, res, head_radius(sys, atoms, n), CartoonGlyph_LipidHead);
                md_array_push(roots, (TubeRoot{h, -1, p, (uint32_t)n}), temp_alloc);
                continue;
            }

            // One sphere for the head, or two (or more): a sphere for each group of the split, and a point without a
            // sphere at each part of the backbone where its tails start
            b.g->num_split_heads += 1;
            b.point(atoms, n, ref, color, res, head_radius(sys, atoms, n), CartoonGlyph_LipidHead, CartoonGlyphVariant_OneHeadSphere);
            for (uint32_t k = 0; k < head_split_count(split); ++k) {
                const int32_t* g_atoms = split.atoms + split.offset[k];
                const size_t   g_n     = split.offset[k + 1] - split.offset[k];
                const int32_t  g_color = head_color_atom(sys, g_atoms, g_n);
                const int32_t  g_res   = residue_of_atom(sys, g_color);
                const CartoonGlyphKind kind = split.outer[k] ? CartoonGlyph_LipidHead : CartoonGlyph_LipidBackbone;
                b.point(g_atoms, g_n, ref, g_color, g_res, head_radius(sys, g_atoms, g_n), kind, CartoonGlyphVariant_TwoHeadSpheres);
                if (split.outer[k]) continue;
                const uint32_t root = b.point(g_atoms, g_n, ref, g_color, g_res, 0.0f, kind);
                for (size_t a = 0; a < md_array_size(anchors); ++a) {
                    if (head_local(g_atoms, g_n, anchors[a]) >= 0) {
                        md_array_push(roots, (TubeRoot{h, anchors[a], root, (uint32_t)g_n}), temp_alloc);
                    }
                }
            }
        }

        // The root of a tube at a head: the part its anchor is in, else the head
        auto root_of = [&](uint32_t h, int32_t anchor) -> TubeRoot {
            TubeRoot head_root = {h, -1, UINT32_MAX, 0};
            for (size_t i = 0; i < md_array_size(roots); ++i) {
                if (roots[i].head != h) continue;
                if (roots[i].anchor == anchor) return roots[i];
                if (roots[i].anchor == -1) head_root = roots[i];
            }
            return head_root;
        };

        for (uint32_t t = t0; t < t1; ++t) {
            const int32_t* spine = md_index_range_beg(&L.tail.spine, t);
            const size_t len = md_index_range_size(&L.tail.spine, t);
            if (len == 0) continue;
            const bool skeleton = L.tail.kind[t] == MD_LIPID_TAIL_SKELETON;
            const CartoonGlyphKind kind = skeleton ? CartoonGlyph_SterolSkeleton : CartoonGlyph_LipidTail;
            const float radius = skeleton ? SKELETON_RADIUS : TAIL_RADIUS;
            const int32_t res = residue_of_atom(sys, spine[0]);
            const bool beads = atom_is_bead(sys, spine[0]);

            // The control points of the tube are points without a sphere of their own (a radius of 0)

            // From the centre of the head (or the part of it the tail hangs off) over the atom it hangs off, unless that
            // atom is all of it (the hydroxyl bead of a coarse grained sterol)
            const int32_t anchor_beg = L.tail.anchor_beg[t];
            const TubeRoot beg = root_of(L.tail.head_beg[t], anchor_beg);
            if (beg.point == UINT32_MAX) continue;
            if (skeleton && sterol_slab(b, sys, L, t, beg.point, ref, in_ring, temp_alloc)) continue;
            b.tube_point(beg.point, radius);
            if (beg.size > 1) {
                b.tube_point(b.point(&anchor_beg, 1, ref, spine[0], res, 0.0f, kind), radius);
            }

            // Along the spine: a control point per bead, or per two atoms, and the tip
            size_t k = 0;
            if (!beads) {
                for (; k + 1 < len; k += 2) {
                    b.tube_point(b.point(spine + k, 2, ref, spine[k], res, 0.0f, kind), radius);
                }
                k = len - 1;
            }
            for (; k < len; ++k) {
                b.tube_point(b.point(spine + k, 1, ref, spine[k], res, 0.0f, kind), radius);
            }

            // A spanning chain goes on to its other head
            if (L.tail.anchor_end[t] != -1 && L.tail.head_end[t] != UINT32_MAX) {
                const int32_t anchor_end = L.tail.anchor_end[t];
                const TubeRoot end = root_of(L.tail.head_end[t], anchor_end);
                if (end.point != UINT32_MAX) {
                    if (end.size > 1) {
                        b.tube_point(b.point(&anchor_end, 1, ref, spine[len - 1], res, 0.0f, kind), radius);
                    }
                    b.tube_point(end.point, radius);
                }
            }
            b.tube_end();
        }
    }
    md_temp_end(temp);
}

// The rings of each component which lie within it, as (component, ring) sorted by component
struct ResRing {
    int32_t  res;
    uint32_t ring;
};

static md_array(ResRing) rings_by_residue(const md_system_t& sys, md_allocator_i* alloc) {
    md_array(ResRing) out = nullptr;
    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* atoms = md_index_range_beg(&sys.ring, r);
        const size_t n = md_index_range_size(&sys.ring, r);
        if (n < 3) continue;
        const md_residue_idx_t res = md_system_residue_find_by_atom_idx(&sys, (size_t)atoms[0]);
        if (res < 0) continue;
        const md_urange_t range = md_system_residue_atom_range(&sys, (size_t)res);
        bool within = true;
        for (size_t k = 0; k < n; ++k) {
            within &= range.beg <= (uint32_t)atoms[k] && (uint32_t)atoms[k] < range.end;
        }
        if (within) md_array_push(out, (ResRing{res, (uint32_t)r}), alloc);
    }
    // Rings come in the order of their atoms, which is nearly that of the residues: sort to be sure
    const size_t n = md_array_size(out);
    for (size_t i = 1; i < n; ++i) {
        const ResRing x = out[i];
        size_t j = i;
        while (j > 0 && out[j - 1].res > x.res) {
            out[j] = out[j - 1];
            --j;
        }
        out[j] = x;
    }
    return out;
}

static bool ring_has_atom(const md_system_t& sys, uint32_t ring, int32_t atom) {
    const int32_t* atoms = md_index_range_beg(&sys.ring, ring);
    const size_t n = md_index_range_size(&sys.ring, ring);
    for (size_t k = 0; k < n; ++k) {
        if (atoms[k] == atom) return true;
    }
    return false;
}

// The rings of 3 to 7 atoms of each atom: atom a is in ring[offset[a] .. offset[a + 1]). For the nucleotides without a
// residue (their backbones found from the bonds alone, as for an XYZ file), whose rings are found from their atoms.
struct AtomRings {
    md_array(uint32_t) offset = nullptr;    // [num_atoms + 1]
    md_array(uint32_t) ring = nullptr;
};

static AtomRings rings_by_atom(const md_system_t& sys, md_allocator_i* alloc) {
    AtomRings out;
    const size_t num_atoms = sys.atom.count;
    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    md_array_resize(out.offset, num_atoms + 1, alloc);
    MEMSET(out.offset, 0, sizeof(uint32_t) * (num_atoms + 1));
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* atoms = md_index_range_beg(&sys.ring, r);
        const size_t n = md_index_range_size(&sys.ring, r);
        if (n < 3 || n > 7) continue;
        for (size_t k = 0; k < n; ++k) out.offset[atoms[k] + 1] += 1;
    }
    for (size_t i = 0; i < num_atoms; ++i) out.offset[i + 1] += out.offset[i];
    md_array_resize(out.ring, out.offset[num_atoms], alloc);
    uint32_t* cursor = (uint32_t*)md_alloc(alloc, sizeof(uint32_t) * (num_atoms + 1));
    MEMCPY(cursor, out.offset, sizeof(uint32_t) * num_atoms);
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* atoms = md_index_range_beg(&sys.ring, r);
        const size_t n = md_index_range_size(&sys.ring, r);
        if (n < 3 || n > 7) continue;
        for (size_t k = 0; k < n; ++k) out.ring[cursor[atoms[k]]++] = (uint32_t)r;
    }
    return out;
}

static bool array_has(const uint32_t* arr, size_t n, uint32_t x) {
    for (size_t i = 0; i < n; ++i) {
        if (arr[i] == x) return true;
    }
    return false;
}

// The rings of a nucleotide from its atoms rather than from its residue: the sugar, the ring with its C1'; then those
// of its base, the ring with the atom bonded to the C1' outside the sugar (the glycosidic N) and the rings fused to it,
// which share a bond with it (the second ring of a purine).
static void nucleotide_rings(md_array(uint32_t)* out, const AtomRings& ar, const md_system_t& sys, int32_t c1, md_allocator_i* alloc) {
    md_array_shrink(*out, 0);
    if (c1 < 0) return;
    uint32_t sugar = UINT32_MAX;
    if (ar.offset[c1] < ar.offset[c1 + 1]) {
        sugar = ar.ring[ar.offset[c1]];
        md_array_push(*out, sugar, alloc);
    }

    int32_t glycosidic = -1;
    md_bond_iter_t it = md_bond_iter(&sys.bond, (size_t)c1);
    for (; md_bond_iter_has_next(&it) && glycosidic < 0; md_bond_iter_next(&it)) {
        const int32_t a = md_bond_iter_atom_index(&it);
        if (md_atom_atomic_number(&sys.atom, (size_t)a) == MD_Z_H) continue;
        if (sugar != UINT32_MAX && ring_has_atom(sys, sugar, a)) continue;
        if (ar.offset[a] < ar.offset[a + 1]) glycosidic = a;
    }
    if (glycosidic < 0) return;

    const size_t first = md_array_size(*out);
    for (uint32_t k = ar.offset[glycosidic]; k < ar.offset[glycosidic + 1]; ++k) {
        if (!array_has(*out, md_array_size(*out), ar.ring[k])) md_array_push(*out, ar.ring[k], alloc);
    }
    for (size_t i = first; i < md_array_size(*out); ++i) {
        const uint32_t ring = (*out)[i];
        const int32_t* atoms = md_index_range_beg(&sys.ring, ring);
        const size_t n = md_index_range_size(&sys.ring, ring);
        for (size_t k = 0; k < n; ++k) {
            for (uint32_t j = ar.offset[atoms[k]]; j < ar.offset[atoms[k] + 1]; ++j) {
                const uint32_t other = ar.ring[j];
                if (other == sugar || array_has(*out, md_array_size(*out), other)) continue;
                size_t shared = 0;
                for (size_t l = 0; l < n; ++l) shared += ring_has_atom(sys, other, atoms[l]) ? 1 : 0;
                if (shared >= 2) md_array_push(*out, other, alloc);
            }
        }
    }
}

// The rings of the bases of the Martini nucleic acids, each ring atom a weighted sum of the real beads of its base
// (cartoon_nucleobases.inl, which tools/gen_cartoon_nucleobases.py generates and checks): the atoms N1 C2 N3 C4 C5 C6,
// then N7 C8 N9 for a purine. The bases of the models are told apart by the number of their beads, virtual sites
// included (3 and 4 in Martini 2, 4 to 6 in Martini 3), which name the same.
static constexpr uint32_t NUCLEOBASE_MAX_REAL = 4;
static constexpr uint32_t NUCLEOBASE_MAX_CORNERS = 9;

struct CartoonNucleobase {
    char        base;           // A, C, G, T or U
    uint8_t     num_beads;      // Of the base, real and virtual
    uint8_t     num_real;       // The beads the ring atoms are weighted sums of
    const char* real[NUCLEOBASE_MAX_REAL];
    uint8_t     num_corners;    // 6 for a pyrimidine, 9 for a purine
    float       weight[NUCLEOBASE_MAX_CORNERS][NUCLEOBASE_MAX_REAL];   // [corner][real bead]
};

#include "cartoon_nucleobases.inl"

static constexpr uint32_t NUCLEOBASE_RING6[] = { 0, 1, 2, 3, 4, 5 };    // N1 C2 N3 C4 C5 C6
static constexpr uint32_t NUCLEOBASE_RING5[] = { 3, 4, 6, 7, 8 };       // C4 C5 N7 C8 N9

// The base of a nucleotide by its name: DA, A, RA, and those of a terminal one (DA5, A3). 0 for none.
static char nucleotide_base(str_t name) {
    while (name.len > 1 && (name.ptr[name.len - 1] == '3' || name.ptr[name.len - 1] == '5')) name.len -= 1;
    char c = 0;
    if (name.len == 1) c = name.ptr[0];
    else if (name.len == 2 && (name.ptr[0] == 'D' || name.ptr[0] == 'R')) c = name.ptr[1];
    return (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U') ? c : 0;
}

// The base of a coarse grained nucleotide with the rings of an ideal one: a point for each ring atom, a weighted sum of
// the real beads of the base, taken whole about ref, and a slab for each ring. The point of the glycosidic N, or
// UINT32_MAX for a base of none of the models (another name, number of beads, a bead missing), which is not drawn.
static uint32_t coarse_grained_base(GlyphBuilder& b, const md_system_t& sys, int32_t res, int32_t ref) {
    const char base = nucleotide_base(md_system_residue_name(&sys, res));
    if (!base) return UINT32_MAX;
    const md_urange_t range = md_system_residue_atom_range(&sys, res);
    uint32_t num_beads = 0;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        num_beads += (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_NUCLEOBASE) ? 1 : 0;
    }
    const CartoonNucleobase* t = nullptr;
    for (const CartoonNucleobase& e : cartoon_nucleobases) {
        if (e.base == base && e.num_beads == num_beads) {
            t = &e;
            break;
        }
    }
    if (!t) return UINT32_MAX;

    int32_t real[NUCLEOBASE_MAX_REAL];
    for (uint32_t k = 0; k < t->num_real; ++k) {
        real[k] = -1;
        for (uint32_t i = range.beg; i < range.end && real[k] < 0; ++i) {
            if (str_eq_cstr(md_atom_name(&sys.atom, i), t->real[k])) real[k] = (int32_t)i;
        }
        if (real[k] < 0) return UINT32_MAX;
    }

    uint32_t corner[NUCLEOBASE_MAX_CORNERS];
    for (uint32_t k = 0; k < t->num_corners; ++k) {
        corner[k] = b.point(real, t->num_real, ref, real[0], res, 0.0f, CartoonGlyph_NucleicBase, CartoonGlyphVariant_Always, t->weight[k]);
    }
    const bool purine = t->num_corners == 9;
    uint32_t ring[6];
    for (uint32_t k = 0; k < ARRAY_SIZE(NUCLEOBASE_RING6); ++k) ring[k] = corner[NUCLEOBASE_RING6[k]];
    b.slab(ring, ARRAY_SIZE(NUCLEOBASE_RING6));
    if (purine) {
        for (uint32_t k = 0; k < ARRAY_SIZE(NUCLEOBASE_RING5); ++k) ring[k] = corner[NUCLEOBASE_RING5[k]];
        b.slab(ring, ARRAY_SIZE(NUCLEOBASE_RING5));
    }
    return purine ? corner[8] : corner[0];
}

// The cartoon of the nucleic acids, after Lindow et al. (Interactive Visualization of RNA and DNA Structures, 2019): a
// flat band along each backbone through the C4' of its nucleotides, its faces towards the bases; the sugar ring as a
// slab; a stick from its C1' to the base; and each ring of the base as a slab of its own (a purine is two, side by side).
// The rings are those of the system (its ring perception) which lie within the nucleotide: the sugar is the one with
// its C1', and the glycosidic nitrogen the atom of a base ring bonded to it. Each corner of a slab is a point of its own
// atom, without a sphere; the points of a nucleotide are taken whole about its C4'.
// A coarse grained nucleotide has no rings in the system: its band goes through the sugar bead with the C4', its base
// has the rings of an ideal one (coarse_grained_base), the stick goes from the band to it, and its sugar is not drawn.
static void build_nucleic_acids(GlyphBuilder& b, const md_system_t& sys) {
    const md_nucleic_backbone_data_t& bb = sys.nucleic_backbone;
    if (bb.range.count == 0) return;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_array(ResRing) res_rings = rings_by_residue(sys, temp_alloc);
    const size_t num_res_rings = md_array_size(res_rings);

    // Nucleotides without a residue have their rings found from their atoms
    bool without_res = false;
    for (size_t s = 0; s < bb.segment.count && !without_res; ++s) without_res = bb.segment.res_idx[s] < 0;
    AtomRings atom_rings = {};
    if (without_res) atom_rings = rings_by_atom(sys, temp_alloc);

    // The point of each atom of the current nucleotide, its rings, and the corners of a slab
    md_array(int32_t)  corner_atom = nullptr;
    md_array(uint32_t) corner_point = nullptr;
    md_array(uint32_t) seg_rings = nullptr;
    md_array(uint32_t) corners = nullptr;

    for (size_t r = 0; r < bb.range.count; ++r) {
        const bool beads = md_nucleic_backbone_range_is_coarse_grained(&bb, r);
        bool open = false;
        for (uint32_t s = bb.range.offset[r]; s < bb.range.offset[r + 1]; ++s) {
            const md_nucleic_acid_atoms_t& atoms = bb.segment.atoms[s];
            const int32_t res = bb.segment.res_idx[s];
            if (atoms.c4 < 0) {
                if (open) b.tube_end(true);
                open = false;
                continue;
            }
            const int32_t ref = atoms.c4;
            md_array_shrink(corner_atom, 0);
            md_array_shrink(corner_point, 0);
            auto point_of = [&](int32_t atom, CartoonGlyphKind kind) -> uint32_t {
                for (size_t i = 0; i < md_array_size(corner_atom); ++i) {
                    if (corner_atom[i] == atom) return corner_point[i];
                }
                const uint32_t p = b.point(&atom, 1, ref, atom, res, 0.0f, kind);
                md_array_push(corner_atom, atom, temp_alloc);
                md_array_push(corner_point, p, temp_alloc);
                return p;
            };

            const uint32_t c4 = point_of(atoms.c4, CartoonGlyph_NucleicBackbone);
            uint32_t base_n = UINT32_MAX;

            if (beads) {
                base_n = coarse_grained_base(b, sys, res, ref);
                if (base_n != UINT32_MAX) {
                    // The stick from the band, whose colour it takes
                    const uint32_t stick[2] = { c4, base_n };
                    b.slab(stick, 2);
                }
                b.tube_point(c4, NUCLEIC_BACKBONE_RADIUS, NUCLEIC_BACKBONE_FLAT, base_n != UINT32_MAX ? base_n : point_of(atoms.c1, CartoonGlyph_NucleicSugar));
                open = true;
                continue;
            }

            // The rings of the nucleotide: those within its residue (binary search for its first), or without one
            // those found from its C1'
            md_array_shrink(seg_rings, 0);
            if (res >= 0) {
                size_t lo = 0, hi = num_res_rings;
                while (lo < hi) {
                    const size_t mid = (lo + hi) / 2;
                    if (res_rings[mid].res < res) lo = mid + 1;
                    else hi = mid;
                }
                for (size_t i = lo; i < num_res_rings && res_rings[i].res == res; ++i) {
                    md_array_push(seg_rings, res_rings[i].ring, temp_alloc);
                }
            } else {
                nucleotide_rings(&seg_rings, atom_rings, sys, atoms.c1, temp_alloc);
            }
            uint32_t sugar_c1 = UINT32_MAX;
            for (size_t i = 0; i < md_array_size(seg_rings); ++i) {
                const uint32_t ring = seg_rings[i];
                const int32_t* ring_atoms = md_index_range_beg(&sys.ring, ring);
                const size_t n = md_index_range_size(&sys.ring, ring);
                if (n > 7) continue;
                const bool sugar = atoms.c1 >= 0 && ring_has_atom(sys, ring, atoms.c1);
                const CartoonGlyphKind kind = sugar ? CartoonGlyph_NucleicSugar : CartoonGlyph_NucleicBase;

                // The corners from a carbon (the C1' of the sugar), whose colour the slab takes
                size_t first = 0;
                for (size_t k = 0; k < n; ++k) {
                    const bool pick = sugar ? ring_atoms[k] == atoms.c1 : md_atom_atomic_number(&sys.atom, ring_atoms[k]) == 6;
                    if (pick) {
                        first = k;
                        break;
                    }
                }
                md_array_shrink(corners, 0);
                for (size_t k = 0; k < n; ++k) {
                    md_array_push(corners, point_of(ring_atoms[(first + k) % n], kind), temp_alloc);
                }
                b.slab(corners, n);

                if (sugar) {
                    sugar_c1 = point_of(atoms.c1, kind);
                } else if (atoms.c1 >= 0 && base_n == UINT32_MAX) {
                    // The glycosidic nitrogen: the atom of the base bonded to the C1'
                    for (size_t k = 0; k < n; ++k) {
                        if (md_bond_find(&sys.bond, ring_atoms[k], atoms.c1) != -1) {
                            base_n = point_of(ring_atoms[k], kind);
                        }
                    }
                }
            }
            if (sugar_c1 != UINT32_MAX && base_n != UINT32_MAX) {
                // The stick: a slab of two corners, from the sugar whose colour it takes
                const uint32_t stick[2] = { sugar_c1, base_n };
                b.slab(stick, 2);
            }

            // The band faces the base, or failing one the C1'
            uint32_t face = base_n;
            if (face == UINT32_MAX && atoms.c1 >= 0) face = point_of(atoms.c1, CartoonGlyph_NucleicSugar);
            b.tube_point(c4, NUCLEIC_BACKBONE_RADIUS, NUCLEIC_BACKBONE_FLAT, face);
            open = true;
        }
        if (open) b.tube_end(true);
    }
    md_temp_end(temp);
}

void cartoon_glyphs_init(CartoonGlyphs* g, const md_system_t& sys, md_gl_mol_t gl_mol, md_allocator_i* alloc) {
    ASSERT(g);
    ASSERT(alloc);
    cartoon_glyphs_free(g);
    g->alloc = alloc;
    md_array_push(g->point_offset, 0, alloc);
    md_array_push(g->tube_offset, 0, alloc);
    md_array_push(g->slab_offset, 0, alloc);

    GlyphBuilder b = { .g = g, .alloc = alloc };
    build_lipids(b, sys);
    build_nucleic_acids(b, sys);

    g->num_points      = (uint32_t)md_array_size(g->point_ref);
    g->num_tubes       = (uint32_t)md_array_size(g->tube_offset) - 1;
    g->num_tube_points = (uint32_t)md_array_size(g->tube_point);
    g->num_slabs       = (uint32_t)md_array_size(g->slab_offset) - 1;
    g->num_spheres     = 0;
    for (uint32_t i = 0; i < g->num_points; ++i) g->num_spheres += g->point_radius[i] > 0.0f ? 1 : 0;
    g->num_source_atoms = (uint32_t)sys.atom.count;
    if (g->num_points == 0 || gl_mol.id == 0) return;

    g->gl_mol = md_gl_mol_create_particles(g->num_points);
    md_gl_mol_set_atom_radius(g->gl_mol, 0, g->num_points, g->point_radius, 0);

    static_assert(sizeof(int32_t) == sizeof(uint32_t), "The atom indices are passed as unsigned");
    md_gl_mol_set_derived(g->gl_mol, gl_mol, g->point_offset, (const uint32_t*)g->point_atom, g->point_weight, (const uint32_t*)g->point_ref);

    if (g->num_tubes > 0) {
        const md_gl_tubes_t tubes = {
            .count    = g->num_tubes,
            .offset   = g->tube_offset,
            .particle = g->tube_point,
            .radius   = g->tube_radius,
            .flat     = g->tube_flat,
            .face     = g->tube_face,
            .side     = g->tube_side,
            .through  = g->tube_through,
        };
        md_gl_mol_set_tubes(g->gl_mol, &tubes);
    }
    if (g->num_slabs > 0) {
        md_gl_mol_set_slabs(g->gl_mol, g->num_slabs, g->slab_offset, g->slab_point);
    }
    g->positions_valid = false;
}

void cartoon_glyphs_free(CartoonGlyphs* g) {
    ASSERT(g);
    if (g->gl_mol.id) md_gl_mol_destroy(g->gl_mol);
    if (g->alloc) {
        md_array_free(g->point_offset, g->alloc);
        md_array_free(g->point_atom, g->alloc);
        md_array_free(g->point_weight, g->alloc);
        md_array_free(g->point_ref, g->alloc);
        md_array_free(g->point_color_atom, g->alloc);
        md_array_free(g->point_res, g->alloc);
        md_array_free(g->point_radius, g->alloc);
        md_array_free(g->point_kind, g->alloc);
        md_array_free(g->point_variant, g->alloc);
        md_array_free(g->tube_offset, g->alloc);
        md_array_free(g->tube_point, g->alloc);
        md_array_free(g->tube_radius, g->alloc);
        md_array_free(g->tube_flat, g->alloc);
        md_array_free(g->tube_face, g->alloc);
        md_array_free(g->tube_side, g->alloc);
        md_array_free(g->tube_through, g->alloc);
        md_array_free(g->slab_offset, g->alloc);
        md_array_free(g->slab_point, g->alloc);
    }
    *g = CartoonGlyphs{};
}

void cartoon_glyphs_update_positions(CartoonGlyphs* g, const md_system_state_t& state) {
    ASSERT(g);
    if (!cartoon_glyphs_valid(*g)) return;
    const md_unitcell_t* cell = md_unitcell_flags(&state.unitcell) ? &state.unitcell : nullptr;
    // The positions on the CPU only place the tubes for their culling, and must be those of every atom
    const vec3_t* xyz = (state.xyz && state.num_atoms >= g->num_source_atoms) ? state.xyz : nullptr;
    md_gl_mol_compute_derived(g->gl_mol, cell, xyz);
    g->positions_valid = true;
}

void cartoon_glyphs_update_flags(const CartoonGlyphs& g, const uint8_t* atom_flags, size_t num_atoms) {
    if (!cartoon_glyphs_valid(g) || !atom_flags) return;
    md_temp_scope_t temp = md_temp_begin();
    uint8_t* flags = (uint8_t*)md_temp_alloc(temp, g.num_points);
    for (uint32_t i = 0; i < g.num_points; ++i) {
        uint8_t f = 0;
        for (uint32_t j = g.point_offset[i]; j < g.point_offset[i + 1]; ++j) {
            const int32_t a = g.point_atom[j];
            if ((size_t)a < num_atoms) f |= atom_flags[a];
        }
        flags[i] = f;
    }
    md_gl_mol_set_atom_flags(g.gl_mol, 0, g.num_points, flags, 0);
    md_temp_end(temp);
}

md_gl_rep_t cartoon_glyphs_rep_create(const CartoonGlyphs& g) {
    if (!cartoon_glyphs_valid(g)) return {};
    return md_gl_rep_create(g.gl_mol);
}

bool cartoon_glyphs_rep_set_colors(const CartoonGlyphs& g, md_gl_rep_t rep, const uint32_t* atom_colors, size_t num_atoms, bool two_head_spheres) {
    if (!cartoon_glyphs_valid(g) || !rep.id || !atom_colors) return false;
    md_temp_scope_t temp = md_temp_begin();
    uint32_t* colors = (uint32_t*)md_temp_alloc(temp, g.num_points * sizeof(uint32_t));
    const uint8_t hidden_variant = two_head_spheres ? CartoonGlyphVariant_OneHeadSphere : CartoonGlyphVariant_TwoHeadSpheres;
    bool any = false;
    for (uint32_t i = 0; i < g.num_points; ++i) {
        bool shown = false;
        if (g.point_variant[i] == hidden_variant) {
            colors[i] = 0;
            continue;
        }
        for (uint32_t j = g.point_offset[i]; j < g.point_offset[i + 1] && !shown; ++j) {
            const int32_t a = g.point_atom[j];
            shown = (size_t)a < num_atoms && (atom_colors[a] & 0xFF000000U);
        }
        const int32_t c = g.point_color_atom[i];
        const uint32_t rgb = (size_t)c < num_atoms ? (atom_colors[c] & 0x00FFFFFFU) : 0x00FFFFFFU;
        colors[i] = rgb | (shown ? 0xFF000000U : 0U);
        any |= shown;
    }
    md_gl_rep_set_atom_colors(rep, 0, g.num_points, colors, 0);
    md_temp_end(temp);
    return any;
}
