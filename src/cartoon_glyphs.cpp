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
static constexpr float NUCLEIC_BACKBONE_RADIUS  = 1.8f;
static constexpr float NUCLEIC_BASE_RADIUS      = 1.0f;

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
        if (best_z == 0 && str_eq_cstr(md_atom_name(&sys.atom, atoms[k]), "PO4")) best = atoms[k];
    }
    return best;
}

// The atom of the base of a nucleotide at its edge which pairs: N1 of a purine (it has an N9), N3 of a pyrimidine
static int32_t nucleic_base_tip(const md_system_t& sys, size_t comp_idx) {
    const md_urange_t r = md_component_atom_range(&sys.component, comp_idx);
    bool purine = false;
    for (uint32_t i = r.beg; i < r.end; ++i) {
        if (str_eq_cstr(md_atom_name(&sys.atom, i), "N9")) purine = true;
    }
    const char* tip = purine ? "N1" : "N3";
    for (uint32_t i = r.beg; i < r.end; ++i) {
        if (str_eq_cstr(md_atom_name(&sys.atom, i), tip)) return (int32_t)i;
    }
    return -1;
}

struct GlyphBuilder {
    CartoonGlyphs* g;
    md_allocator_i* alloc;
    md_array(uint32_t) seg_a = nullptr;
    md_array(uint32_t) seg_b = nullptr;

    uint32_t point(const int32_t* atoms, size_t n, int32_t ref, int32_t color_atom, int32_t comp, float radius, CartoonGlyphKind kind) {
        const uint32_t idx = (uint32_t)md_array_size(g->point_ref);
        md_array_push_array(g->point_atom, atoms, n, alloc);
        md_array_push(g->point_offset, (uint32_t)md_array_size(g->point_atom), alloc);
        md_array_push(g->point_ref, ref, alloc);
        md_array_push(g->point_color_atom, color_atom, alloc);
        md_array_push(g->point_comp, comp, alloc);
        md_array_push(g->point_radius, radius, alloc);
        md_array_push(g->point_kind, (uint8_t)kind, alloc);
        return idx;
    }

    void segment(uint32_t a, uint32_t b, float radius) {
        md_array_push(g->segment_point, a, alloc);
        md_array_push(g->segment_point, b, alloc);
        md_array_push(g->segment_radius, radius, alloc);
    }
};

static void build_lipids(GlyphBuilder& b, const md_system_t& sys) {
    const md_lipid_data_t& L = sys.lipid;
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_array(uint32_t) head_point = nullptr;    // The point of each head of the lipid
    for (size_t l = 0; l < L.count; ++l) {
        const md_urange_t cr = L.comp_range[l];
        const uint32_t h0 = L.head_offset[l];
        const uint32_t h1 = L.head_offset[l + 1];
        if (h1 == h0) continue;
        md_array_resize(head_point, h1 - h0, temp_alloc);
        // Every point of the lipid in the image of its first head atom
        const int32_t ref = md_index_range_beg(&L.head.atoms, h0)[0];

        for (uint32_t h = h0; h < h1; ++h) {
            const int32_t* atoms = md_index_range_beg(&L.head.atoms, h);
            const size_t n = md_index_range_size(&L.head.atoms, h);
            const float heavy = (float)n * (atom_is_bead(sys, atoms[0]) ? HEAVY_ATOMS_PER_BEAD : 1.0f);
            const float radius = CLAMP(HEAD_RADIUS_SCALE * cbrtf(heavy), HEAD_RADIUS_MIN, HEAD_RADIUS_MAX);
            const int32_t color = head_color_atom(sys, atoms, n);
            head_point[h - h0] = b.point(atoms, n, ref, color, component_of_atom(sys, cr.beg, cr.end, color), radius, CartoonGlyph_LipidHead);
        }

        for (uint32_t t = L.tail_offset[l]; t < L.tail_offset[l + 1]; ++t) {
            const int32_t* spine = md_index_range_beg(&L.tail.spine, t);
            const size_t len = md_index_range_size(&L.tail.spine, t);
            if (len == 0) continue;
            const bool skeleton = L.tail.kind[t] == MD_LIPID_TAIL_SKELETON;
            const CartoonGlyphKind kind = skeleton ? CartoonGlyph_SterolSkeleton : CartoonGlyph_LipidTail;
            const float radius = skeleton ? SKELETON_RADIUS : TAIL_RADIUS;
            const int32_t comp = component_of_atom(sys, cr.beg, cr.end, spine[0]);
            const bool beads = atom_is_bead(sys, spine[0]);

            // From the centre of the head over the atom it hangs off, unless that atom is all of the head (the hydroxyl
            // bead of a coarse grained sterol)
            uint32_t prev = head_point[L.tail.head_beg[t] - h0];
            const int32_t anchor_beg = L.tail.anchor_beg[t];
            if (md_index_range_size(&L.head.atoms, L.tail.head_beg[t]) > 1) {
                const uint32_t p = b.point(&anchor_beg, 1, ref, spine[0], comp, radius, kind);
                b.segment(prev, p, radius);
                prev = p;
            }

            // Along the spine: a point per bead, or per two atoms, and the tip
            size_t k = 0;
            if (!beads) {
                for (; k + 1 < len; k += 2) {
                    const uint32_t p = b.point(spine + k, 2, ref, spine[k], comp, radius, kind);
                    b.segment(prev, p, radius);
                    prev = p;
                }
                k = len - 1;
            }
            for (; k < len; ++k) {
                const uint32_t p = b.point(spine + k, 1, ref, spine[k], comp, radius, kind);
                b.segment(prev, p, radius);
                prev = p;
            }

            // A spanning chain goes on to its other head
            if (L.tail.anchor_end[t] != -1 && L.tail.head_end[t] != UINT32_MAX) {
                const int32_t anchor_end = L.tail.anchor_end[t];
                if (md_index_range_size(&L.head.atoms, L.tail.head_end[t]) > 1) {
                    const uint32_t p = b.point(&anchor_end, 1, ref, spine[len - 1], comp, radius, kind);
                    b.segment(prev, p, radius);
                    prev = p;
                }
                b.segment(prev, head_point[L.tail.head_end[t] - h0], radius);
            }
        }
    }
    md_temp_end(temp);
}

static void build_nucleic_acids(GlyphBuilder& b, const md_system_t& sys) {
    const md_nucleic_backbone_data_t& bb = sys.nucleic_backbone;
    for (size_t r = 0; r < bb.range.count; ++r) {
        uint32_t prev = UINT32_MAX;
        for (uint32_t s = bb.range.offset[r]; s < bb.range.offset[r + 1]; ++s) {
            const md_nucleic_acid_atoms_t& atoms = bb.segment.atoms[s];
            const int32_t comp = bb.segment.comp_idx[s];
            if (atoms.c4 < 0) {
                prev = UINT32_MAX;
                continue;
            }
            const uint32_t p = b.point(&atoms.c4, 1, atoms.c4, atoms.c4, comp, NUCLEIC_BACKBONE_RADIUS, CartoonGlyph_NucleicBackbone);
            if (prev != UINT32_MAX) b.segment(prev, p, NUCLEIC_BACKBONE_RADIUS);
            prev = p;

            // The base, over the C1' to its edge which pairs
            const int32_t tip = nucleic_base_tip(sys, (size_t)comp);
            if (atoms.c1 >= 0 && tip >= 0) {
                const uint32_t c1 = b.point(&atoms.c1, 1, atoms.c1, tip, comp, NUCLEIC_BASE_RADIUS, CartoonGlyph_NucleicBase);
                const uint32_t e  = b.point(&tip, 1, tip, tip, comp, NUCLEIC_BASE_RADIUS, CartoonGlyph_NucleicBase);
                b.segment(p, c1, NUCLEIC_BASE_RADIUS);
                b.segment(c1, e, NUCLEIC_BASE_RADIUS);
            }
        }
    }
}

void cartoon_glyphs_init(CartoonGlyphs* g, const md_system_t& sys, md_allocator_i* alloc) {
    ASSERT(g);
    ASSERT(alloc);
    cartoon_glyphs_free(g);
    g->alloc = alloc;
    md_array_push(g->point_offset, 0, alloc);

    GlyphBuilder b = { .g = g, .alloc = alloc };
    build_lipids(b, sys);
    build_nucleic_acids(b, sys);

    g->num_points   = (uint32_t)md_array_size(g->point_ref);
    g->num_segments = (uint32_t)md_array_size(g->segment_radius);
    if (g->num_points == 0) return;

    md_array_resize(g->xyz, g->num_points, alloc);
    MEMSET(g->xyz, 0, md_array_bytes(g->xyz));
    g->gl_mol = md_gl_mol_create_particles(g->num_points);
    md_gl_mol_set_atom_radius(g->gl_mol, 0, g->num_points, g->point_radius, 0);
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
        md_array_free(g->segment_point, g->alloc);
        md_array_free(g->segment_radius, g->alloc);
        md_array_free(g->xyz, g->alloc);
    }
    *g = CartoonGlyphs{};
}

void cartoon_glyphs_update_positions(CartoonGlyphs* g, const md_system_state_t& state, bool velocity) {
    ASSERT(g);
    if (!cartoon_glyphs_valid(*g) || state.num_atoms == 0) return;
    const md_unitcell_t* cell = md_unitcell_flags(&state.unitcell) ? &state.unitcell : nullptr;

    vec3_t dx[64];
    for (uint32_t i = 0; i < g->num_points; ++i) {
        const uint32_t beg = g->point_offset[i];
        const uint32_t end = g->point_offset[i + 1];
        const int32_t ref = g->point_ref[i];
        if (end - beg == 1 && g->point_atom[beg] == ref) {
            g->xyz[i] = state.xyz[ref];
            continue;
        }
        // The mean of the atoms at their images nearest to the reference, in batches
        const vec3_t r = state.xyz[ref];
        vec3_t sum = {0, 0, 0};
        for (uint32_t j = beg; j < end; j += (uint32_t)ARRAY_SIZE(dx)) {
            const uint32_t n = MIN(end - j, (uint32_t)ARRAY_SIZE(dx));
            for (uint32_t k = 0; k < n; ++k) dx[k] = vec3_sub(state.xyz[g->point_atom[j + k]], r);
            md_util_min_image_vec3(dx, n, cell);
            for (uint32_t k = 0; k < n; ++k) sum = vec3_add(sum, dx[k]);
        }
        g->xyz[i] = vec3_add(r, vec3_mul1(sum, 1.0f / (float)(end - beg)));
    }

    md_gl_mol_set_atom_position(g->gl_mol, 0, g->num_points, g->xyz);
    if (velocity) {
        vec3_t pbc_ext = {0, 0, 0};
        md_unitcell_diag_extract_float(pbc_ext.elem, &state.unitcell);
        md_gl_mol_compute_velocity(g->gl_mol, pbc_ext.elem);
    } else {
        md_gl_mol_zero_velocity(g->gl_mol);
    }
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
    md_gl_rep_t rep = md_gl_rep_create(g.gl_mol);
    if (g.num_segments > 0) {
        md_temp_scope_t temp = md_temp_begin();
        uint32_t* a = (uint32_t*)md_temp_alloc(temp, g.num_segments * sizeof(uint32_t));
        uint32_t* b = (uint32_t*)md_temp_alloc(temp, g.num_segments * sizeof(uint32_t));
        for (uint32_t i = 0; i < g.num_segments; ++i) {
            a[i] = g.segment_point[2 * i + 0];
            b[i] = g.segment_point[2 * i + 1];
        }
        // The radius of a capsule is the radius of the draw (the scale of the representation) times its weight
        md_gl_rep_set_weak_bonds(rep, g.num_segments, a, b, g.segment_radius);
        md_temp_end(temp);
    }
    return rep;
}

bool cartoon_glyphs_rep_set_colors(const CartoonGlyphs& g, md_gl_rep_t rep, const uint32_t* atom_colors, size_t num_atoms) {
    if (!cartoon_glyphs_valid(g) || !rep.id || !atom_colors) return false;
    md_temp_scope_t temp = md_temp_begin();
    uint32_t* colors = (uint32_t*)md_temp_alloc(temp, g.num_points * sizeof(uint32_t));
    bool any = false;
    for (uint32_t i = 0; i < g.num_points; ++i) {
        bool shown = false;
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
