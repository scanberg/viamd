#include <cartoon_glyphs.h>

#include <md_system.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_log.h>

#include <math.h>

// Radii of the spheres and capsules, Å
static constexpr float HEAD_RADIUS_SCALE        = 1.6f;   // A head of n heavy atoms is 1.6 cbrt(n): 4.3 for the 20 of a phosphatidyl choline
static constexpr float HEAD_RADIUS_MIN          = 1.8f;
static constexpr float HEAD_RADIUS_MAX          = 7.0f;
static constexpr float HEAVY_ATOMS_PER_BEAD     = 4.0f;   // What a Martini bead stands for
static constexpr float TAIL_RADIUS              = 1.1f;
static constexpr float SKELETON_RADIUS          = 2.0f;
// The backbone of a nucleic acid: a flat band 2 (r + l) wide and 2 r thick, its faces towards the bases. The r is that of
// a representation (CARTOON_NUCLEIC_BAND_RADIUS is the default, which it scales), l stays.
static constexpr float NUCLEIC_BACKBONE_RADIUS  = CARTOON_NUCLEIC_BAND_RADIUS;
static constexpr float NUCLEIC_BACKBONE_FLAT    = 1.6f;

static inline bool atom_is_bead(const md_system_t& sys, size_t i) {
    return md_atom_particle_kind(&sys.atom, i) != MD_PARTICLE_ATOM && md_atom_atomic_number(&sys.atom, i) == 0;
}

// The component of an atom among those of [beg, end)
static int32_t component_of_atom(const md_system_t& sys, uint32_t comp_beg, uint32_t comp_end, int32_t atom) {
    for (uint32_t c = comp_beg; c < comp_end; ++c) {
        const md_urange_t r = md_component_atom_range(&sys.component, c);
        if (r.beg <= (uint32_t)atom && (uint32_t)atom < r.end) return (int32_t)c;
    }
    return (int32_t)comp_beg;
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

    uint32_t point(const int32_t* atoms, size_t n, int32_t ref, int32_t color_atom, int32_t comp, float radius, CartoonGlyphKind kind, CartoonGlyphVariant variant = CartoonGlyphVariant_Always) {
        const uint32_t idx = (uint32_t)md_array_size(g->point_ref);
        md_array_push(g->point_variant, (uint8_t)variant, alloc);
        md_array_push_array(g->point_atom, atoms, n, alloc);
        md_array_push(g->point_offset, (uint32_t)md_array_size(g->point_atom), alloc);
        md_array_push(g->point_ref, ref, alloc);
        md_array_push(g->point_color_atom, color_atom, alloc);
        md_array_push(g->point_comp, comp, alloc);
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
    }

    // flat: half the length of the flat sides of the profile, face: the point they face (UINT32_MAX for none)
    void tube_point(uint32_t p, float radius, float flat, uint32_t face) {
        tube_point(p, radius);
        g->tube_flat[md_array_size(g->tube_flat) - 1] = flat;
        g->tube_face[md_array_size(g->tube_face) - 1] = face;
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

static void build_lipids(GlyphBuilder& b, const md_system_t& sys) {
    const md_lipid_data_t& L = sys.lipid;
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_array(TubeRoot) roots = nullptr;     // Of the heads of the lipid
    md_array(int32_t) anchors = nullptr;
    for (size_t l = 0; l < L.count; ++l) {
        const md_urange_t cr = L.comp_range[l];
        const uint32_t h0 = L.head_offset[l];
        const uint32_t h1 = L.head_offset[l + 1];
        const uint32_t t0 = L.tail_offset[l];
        const uint32_t t1 = L.tail_offset[l + 1];
        if (h1 == h0) continue;
        md_array_shrink(roots, 0);
        // Every point of the lipid in the image of its first head atom
        const int32_t ref = md_index_range_beg(&L.head.atoms, h0)[0];

        for (uint32_t h = h0; h < h1; ++h) {
            const int32_t* atoms = md_index_range_beg(&L.head.atoms, h);
            const size_t n = md_index_range_size(&L.head.atoms, h);
            const int32_t color = head_color_atom(sys, atoms, n);
            const int32_t comp = component_of_atom(sys, cr.beg, cr.end, color);

            // The atoms of the head its tails hang off
            md_array_shrink(anchors, 0);
            for (uint32_t t = t0; t < t1; ++t) {
                if (L.tail.head_beg[t] == h) md_array_push(anchors, L.tail.anchor_beg[t], temp_alloc);
                if (L.tail.head_end[t] == h && L.tail.anchor_end[t] != -1) md_array_push(anchors, L.tail.anchor_end[t], temp_alloc);
            }

            HeadSplit split = {};
            if (!split_head(&split, sys, atoms, n, anchors, md_array_size(anchors), temp_alloc)) {
                const uint32_t p = b.point(atoms, n, ref, color, comp, head_radius(sys, atoms, n), CartoonGlyph_LipidHead);
                md_array_push(roots, (TubeRoot{h, -1, p, (uint32_t)n}), temp_alloc);
                continue;
            }

            // One sphere for the head, or two (or more): a sphere for each group of the split, and a point without a
            // sphere at each part of the backbone where its tails start
            b.g->num_split_heads += 1;
            b.point(atoms, n, ref, color, comp, head_radius(sys, atoms, n), CartoonGlyph_LipidHead, CartoonGlyphVariant_OneHeadSphere);
            for (uint32_t k = 0; k < head_split_count(split); ++k) {
                const int32_t* g_atoms = split.atoms + split.offset[k];
                const size_t   g_n     = split.offset[k + 1] - split.offset[k];
                const int32_t  g_color = head_color_atom(sys, g_atoms, g_n);
                const int32_t  g_comp  = component_of_atom(sys, cr.beg, cr.end, g_color);
                const CartoonGlyphKind kind = split.outer[k] ? CartoonGlyph_LipidHead : CartoonGlyph_LipidBackbone;
                b.point(g_atoms, g_n, ref, g_color, g_comp, head_radius(sys, g_atoms, g_n), kind, CartoonGlyphVariant_TwoHeadSpheres);
                if (split.outer[k]) continue;
                const uint32_t root = b.point(g_atoms, g_n, ref, g_color, g_comp, 0.0f, kind);
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
            const int32_t comp = component_of_atom(sys, cr.beg, cr.end, spine[0]);
            const bool beads = atom_is_bead(sys, spine[0]);

            // The control points of the tube are points without a sphere of their own (a radius of 0)

            // From the centre of the head (or the part of it the tail hangs off) over the atom it hangs off, unless that
            // atom is all of it (the hydroxyl bead of a coarse grained sterol)
            const int32_t anchor_beg = L.tail.anchor_beg[t];
            const TubeRoot beg = root_of(L.tail.head_beg[t], anchor_beg);
            if (beg.point == UINT32_MAX) continue;
            b.tube_point(beg.point, radius);
            if (beg.size > 1) {
                b.tube_point(b.point(&anchor_beg, 1, ref, spine[0], comp, 0.0f, kind), radius);
            }

            // Along the spine: a control point per bead, or per two atoms, and the tip
            size_t k = 0;
            if (!beads) {
                for (; k + 1 < len; k += 2) {
                    b.tube_point(b.point(spine + k, 2, ref, spine[k], comp, 0.0f, kind), radius);
                }
                k = len - 1;
            }
            for (; k < len; ++k) {
                b.tube_point(b.point(spine + k, 1, ref, spine[k], comp, 0.0f, kind), radius);
            }

            // A spanning chain goes on to its other head
            if (L.tail.anchor_end[t] != -1 && L.tail.head_end[t] != UINT32_MAX) {
                const int32_t anchor_end = L.tail.anchor_end[t];
                const TubeRoot end = root_of(L.tail.head_end[t], anchor_end);
                if (end.point != UINT32_MAX) {
                    if (end.size > 1) {
                        b.tube_point(b.point(&anchor_end, 1, ref, spine[len - 1], comp, 0.0f, kind), radius);
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
struct CompRing {
    int32_t  comp;
    uint32_t ring;
};

static md_array(CompRing) rings_by_component(const md_system_t& sys, md_allocator_i* alloc) {
    md_array(CompRing) out = nullptr;
    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* atoms = md_index_range_beg(&sys.ring, r);
        const size_t n = md_index_range_size(&sys.ring, r);
        if (n < 3) continue;
        const md_component_idx_t comp = md_system_component_find_by_atom_idx(&sys, (size_t)atoms[0]);
        if (comp < 0) continue;
        const md_urange_t range = md_component_atom_range(&sys.component, (size_t)comp);
        bool within = true;
        for (size_t k = 0; k < n; ++k) {
            within &= range.beg <= (uint32_t)atoms[k] && (uint32_t)atoms[k] < range.end;
        }
        if (within) md_array_push(out, (CompRing{comp, (uint32_t)r}), alloc);
    }
    // Rings come in the order of their atoms, which is nearly that of the components: sort to be sure
    const size_t n = md_array_size(out);
    for (size_t i = 1; i < n; ++i) {
        const CompRing x = out[i];
        size_t j = i;
        while (j > 0 && out[j - 1].comp > x.comp) {
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

// The cartoon of the nucleic acids, after Lindow et al. (Interactive Visualization of RNA and DNA Structures, 2019): a
// flat band along each backbone through the C4' of its nucleotides, its faces towards the bases; the sugar ring as a
// slab; a stick from its C1' to the base; and each ring of the base as a slab of its own (a purine is two, side by side).
// The rings are those of the system (its ring perception) which lie within the nucleotide: the sugar is the one with
// its C1', and the glycosidic nitrogen the atom of a base ring bonded to it. Each corner of a slab is a point of its own
// atom, without a sphere; the points of a nucleotide are taken whole about its C4'.
static void build_nucleic_acids(GlyphBuilder& b, const md_system_t& sys) {
    const md_nucleic_backbone_data_t& bb = sys.nucleic_backbone;
    if (bb.range.count == 0) return;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_array(CompRing) comp_rings = rings_by_component(sys, temp_alloc);
    const size_t num_comp_rings = md_array_size(comp_rings);

    // The point of each atom of the current nucleotide, and the corners of a slab
    md_array(int32_t)  corner_atom = nullptr;
    md_array(uint32_t) corner_point = nullptr;
    md_array(uint32_t) corners = nullptr;

    for (size_t r = 0; r < bb.range.count; ++r) {
        bool open = false;
        for (uint32_t s = bb.range.offset[r]; s < bb.range.offset[r + 1]; ++s) {
            const md_nucleic_acid_atoms_t& atoms = bb.segment.atoms[s];
            const int32_t comp = bb.segment.comp_idx[s];
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
                const uint32_t p = b.point(&atom, 1, ref, atom, comp, 0.0f, kind);
                md_array_push(corner_atom, atom, temp_alloc);
                md_array_push(corner_point, p, temp_alloc);
                return p;
            };

            const uint32_t c4 = point_of(atoms.c4, CartoonGlyph_NucleicBackbone);

            // The rings of the nucleotide: binary search for its first
            size_t lo = 0, hi = num_comp_rings;
            while (lo < hi) {
                const size_t mid = (lo + hi) / 2;
                if (comp_rings[mid].comp < comp) lo = mid + 1;
                else hi = mid;
            }
            uint32_t sugar_c1 = UINT32_MAX;
            uint32_t base_n   = UINT32_MAX;
            for (size_t i = lo; i < num_comp_rings && comp_rings[i].comp == comp; ++i) {
                const uint32_t ring = comp_rings[i].ring;
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

    // Each point the mean of its atoms
    md_temp_scope_t temp = md_temp_begin();
    const size_t num_atoms = md_array_size(g->point_atom);
    float* weight = (float*)md_temp_alloc(temp, MAX(num_atoms, 1) * sizeof(float));
    for (uint32_t i = 0; i < g->num_points; ++i) {
        const uint32_t beg = g->point_offset[i];
        const uint32_t end = g->point_offset[i + 1];
        for (uint32_t j = beg; j < end; ++j) weight[j] = 1.0f / (float)(end - beg);
    }
    static_assert(sizeof(int32_t) == sizeof(uint32_t), "The atom indices are passed as unsigned");
    md_gl_mol_set_derived(g->gl_mol, gl_mol, g->point_offset, (const uint32_t*)g->point_atom, weight, (const uint32_t*)g->point_ref);
    md_temp_end(temp);

    if (g->num_tubes > 0) {
        const md_gl_tubes_t tubes = {
            .count    = g->num_tubes,
            .offset   = g->tube_offset,
            .particle = g->tube_point,
            .radius   = g->tube_radius,
            .flat     = g->tube_flat,
            .face     = g->tube_face,
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
        md_array_free(g->point_ref, g->alloc);
        md_array_free(g->point_color_atom, g->alloc);
        md_array_free(g->point_comp, g->alloc);
        md_array_free(g->point_radius, g->alloc);
        md_array_free(g->point_kind, g->alloc);
        md_array_free(g->point_variant, g->alloc);
        md_array_free(g->tube_offset, g->alloc);
        md_array_free(g->tube_point, g->alloc);
        md_array_free(g->tube_radius, g->alloc);
        md_array_free(g->tube_flat, g->alloc);
        md_array_free(g->tube_face, g->alloc);
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
