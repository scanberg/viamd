#pragma once

#include <stdint.h>
#include <stddef.h>

#include <md_gl.h>
#include <core/md_array.h>
#include <core/md_vec_math.h>

struct md_system_t;
struct md_system_state_t;
struct md_allocator_i;

// THE CARTOON OF LIPIDS AND NUCLEIC ACIDS
//
// The cartoon of a protein is a spline through its backbone, extruded with a profile its secondary structure drives
// (md_gl). That of a lipid is made the same way: a sphere for each head and a spline for each tail, a tube of round
// profile. That of a nucleic acid is a flat band along its backbone, with the rings of its sugars and bases as slabs.
//
//   lipid      A sphere for each head, at the centre of its heavy atoms and sized by how many there are (a phosphatidyl
//              choline is larger than the hydroxyl of a sterol, an inositol phosphate larger still). Or two, where a
//              representation asks for them (CartoonGlyphVariant): the head of a phospholipid splits at its phosphate
//              into the head group with the phosphate (a phosphocholine) and the backbone the tails hang off (a glycerol
//              and its esters), a sphere for each, sized the same way (a cardiolipin has two backbones and one head
//              group between them). A tube for each tail
//              (md_lipid_data_t) from the centre of its head over the atom it hangs off and along its spine to its tip,
//              or on to the head a spanning chain joins (a bolalipid). An atomistic spine gets a control point per two
//              of its atoms, which is on the axis of a chain in its all trans zigzag and so smooths it away; a coarse
//              grained one a control point per bead. The ring system of a sterol is a slab: a tube flat in the plane
//              of its rings, wider than a chain and as thick, through the middle of the rings in sections along it, and
//              its side chain a tube as a chain is. The tubes have round ends.
//              A lipid whose heads and tails are not known (md_lipid_data_t has none for it) is one sphere for all of it.
//   nucleic    After Lindow et al. (Interactive Visualization of RNA and DNA Structures, 2019): a flat band along each
//              backbone (md_nucleic_backbone_data_t) through the C4' of its nucleotides, its flat faces towards the bases;
//              the sugar ring a slab; a stick from its C1' to the base; each ring of the base a slab of its own (a purine
//              two side by side, rather than one plate over both). A slab is the ring flattened onto its plane and made
//              solid by a sphere: capsules along its edges and a plate between them (md_gl_mol_set_slabs). The stick is a
//              slab of two corners, the C1' and the glycosidic N: one capsule, as thick as the rims of the rings it joins
//              and as round where it meets them. A representation sets the radius of the band (half its thickness, the
//              width stays) and that of the rims (CARTOON_NUCLEIC_BAND_RADIUS, CARTOON_NUCLEIC_RING_RADIUS).
//              A coarse grained nucleic acid (Martini, md_nucleic_backbone_range_is_coarse_grained) holds no more of a
//              base than where it is: its beads are a rigid body in the plane of its rings. Its base is drawn with the
//              rings of an ideal base where its beads are, each ring atom a fixed weighted sum of them (the one which
//              gives it from the beads of the ideal base, cartoon_nucleobases.inl); the band goes through the sugar
//              bead with the C4' (BB2), and the stick from it to the glycosidic N. Its sugar is not drawn: two beads
//              tell neither the pucker of the ring nor how it lies.
//
// The points (the centres of the spheres, the control points of the tubes, the ends of the capsules) are the particles
// of a molecule of their own (md_gl_mol_create_particles). Each is a weighted sum of the atoms it stands for, most of
// them the mean, derived on the GPU from the positions of the atoms of the system (md_gl_mol_set_derived,
// md_gl_mol_compute_derived), each atom in its periodic image nearest to the first head atom of its lipid, so a lipid
// stays whole across the periodic boundary. The spheres are drawn as such (MD_GL_REP_SPACE_FILL, a control point of a
// tube or a corner of a slab has no sphere: a radius of 0), the tubes as the splines of md_gl through the control
// points with a stadium profile, round or flat (md_gl_mol_set_tubes, MD_GL_REP_TUBES) and the slabs raycast through
// their corners (md_gl_mol_set_slabs, MD_GL_REP_SLABS). A point stands for a component, which is what picking it
// resolves to, takes its colour from one of its atoms, and is shown when any of its atoms is shown. A part of a tube is
// that of its control point, a slab that of its first corner.

// The radii of the cartoon of a nucleic acid, Å, by default (Representation): of the band of its backbone, which its
// tubes are made with (a representation scales them to its own, md_gl_draw_op_t args.tubes.band_scale), and of the rims
// of its slabs and the stick from each sugar to its base
static constexpr float CARTOON_NUCLEIC_BAND_RADIUS = 0.25f;
static constexpr float CARTOON_NUCLEIC_RING_RADIUS = 0.2f;

enum CartoonGlyphKind : uint8_t {
    CartoonGlyph_LipidHead = 0,     // A head, or the head group of a split head
    CartoonGlyph_LipidBackbone,     // The backbone of a split head
    CartoonGlyph_LipidTail,
    CartoonGlyph_SterolSkeleton,
    CartoonGlyph_NucleicBackbone,
    CartoonGlyph_NucleicSugar,
    CartoonGlyph_NucleicBase,
    CartoonGlyph_Count
};

// Which representations show a point: lipid heads are one sphere, or two where the head splits
enum CartoonGlyphVariant : uint8_t {
    CartoonGlyphVariant_Always = 0,
    CartoonGlyphVariant_OneHeadSphere,      // The sphere of a head which splits
    CartoonGlyphVariant_TwoHeadSpheres,     // The spheres of its parts
};

struct CartoonGlyphs {
    md_gl_mol_t gl_mol = {};
    uint32_t num_points = 0;
    uint32_t num_spheres = 0;       // Points with a sphere (radius > 0): the heads of lipids. A nucleic acid has none
    uint32_t num_tubes = 0;
    uint32_t num_tube_points = 0;   // Control points of all tubes
    uint32_t num_slabs = 0;
    uint32_t num_source_atoms = 0;  // Of the system the points are derived from
    uint32_t num_split_heads = 0;   // Lipid heads which are two spheres (or more) in a representation which asks for that

    // A point is the weighted sum of its atoms, point_atom[point_offset[i] .. point_offset[i + 1]) with the weights
    // point_weight[...] of the same range (their mean, but for the ring atoms of a coarse grained base), each atom
    // taken at its image nearest to point_ref[i]. The points of one lipid share their reference, so a lipid stays whole
    // across the periodic boundary.
    md_array(uint32_t) point_offset = nullptr;       // [num_points + 1]
    md_array(int32_t)  point_atom = nullptr;
    md_array(float)    point_weight = nullptr;       // Of point_atom
    md_array(int32_t)  point_ref = nullptr;          // [num_points]
    md_array(int32_t)  point_color_atom = nullptr;   // [num_points] The atom whose colour it takes
    md_array(int32_t)  point_comp = nullptr;         // [num_points] The component it stands for
    md_array(float)    point_radius = nullptr;       // [num_points] Of its sphere, Å, before the scale of the representation. 0 for none
    md_array(uint8_t)  point_kind = nullptr;         // [num_points] CartoonGlyphKind
    md_array(uint8_t)  point_variant = nullptr;      // [num_points] CartoonGlyphVariant

    // Tube t goes through the points tube_point[tube_offset[t] .. tube_offset[t + 1])
    md_array(uint32_t) tube_offset = nullptr;        // [num_tubes + 1]
    md_array(uint32_t) tube_point = nullptr;         // [num_tube_points]
    md_array(float)    tube_radius = nullptr;        // [num_tube_points] Å, before the scale of the representation
    md_array(float)    tube_flat = nullptr;          // [num_tube_points] Half the length of the flat sides of the profile, Å: 0 is round
    md_array(uint32_t) tube_face = nullptr;          // [num_tube_points] The point the flat sides face, UINT32_MAX for none
    md_array(uint32_t) tube_side = nullptr;          // [num_tube_points] The point the long axis of the profile points at, UINT32_MAX for none
    md_array(uint8_t)  tube_through = nullptr;       // [num_tubes] It passes through its control points rather than near them

    // Slab s has the corners slab_point[slab_offset[s] .. slab_offset[s + 1]), points in order around the ring, or the
    // two ends of a stick
    md_array(uint32_t) slab_offset = nullptr;        // [num_slabs + 1]
    md_array(uint32_t) slab_point = nullptr;

    bool positions_valid = false;                    // The points are derived from the current positions and velocities of the atoms

    md_allocator_i* alloc = nullptr;
};

// The glyphs of the lipids and nucleic acid backbones of the system, none when it has neither (cartoon_glyphs_valid).
// Their points are derived from the atoms of gl_mol, the molecule of the system, which must outlive them.
void cartoon_glyphs_init(CartoonGlyphs* glyphs, const md_system_t& sys, md_gl_mol_t gl_mol, md_allocator_i* alloc);
void cartoon_glyphs_free(CartoonGlyphs* glyphs);

static inline bool cartoon_glyphs_valid(const CartoonGlyphs& glyphs) {
    return glyphs.num_points > 0 && glyphs.gl_mol.id != 0;
}

// Places the points from the positions and velocities of the atoms on the GPU (those of the molecule of the system, set
// from the state first), on the GPU. The positions of the state on the CPU only place the tubes for their culling.
void cartoon_glyphs_update_positions(CartoonGlyphs* glyphs, const md_system_state_t& state);

// The flags of the points (highlighted, selected, visible, see AtomBit_): those of any of their atoms
void cartoon_glyphs_update_flags(const CartoonGlyphs& glyphs, const uint8_t* atom_flags, size_t num_atoms);

// A representation of the glyphs for a representation of the system. Destroy with md_gl_rep_destroy.
md_gl_rep_t cartoon_glyphs_rep_create(const CartoonGlyphs& glyphs);

// The colours of the points from those of the atoms of a representation (an alpha of 0 hides an atom), and whether the
// heads which split are two spheres or one. Whether any point is shown.
bool cartoon_glyphs_rep_set_colors(const CartoonGlyphs& glyphs, md_gl_rep_t rep, const uint32_t* atom_colors, size_t num_atoms, bool two_head_spheres);

// Picking: the spheres write the index of a point, the tubes num_points + the index of their control
// point, the slabs num_points + num_tube_points + their index (md_gl_picking_offset_t: atom_base at the beginning of the
// range, backbone_base num_points after it, ring_base num_tube_points after that)
static inline uint32_t cartoon_glyphs_picking_count(const CartoonGlyphs& glyphs) {
    return glyphs.num_points + glyphs.num_tube_points + glyphs.num_slabs;
}

// The point a picking index stands for, UINT32_MAX for none
static inline uint32_t cartoon_glyphs_picked_point(const CartoonGlyphs& glyphs, uint32_t picking_idx) {
    if (picking_idx < glyphs.num_points) return picking_idx;
    picking_idx -= glyphs.num_points;
    if (picking_idx < glyphs.num_tube_points) return glyphs.tube_point[picking_idx];
    picking_idx -= glyphs.num_tube_points;
    if (picking_idx < glyphs.num_slabs) return glyphs.slab_point[glyphs.slab_offset[picking_idx]];
    return UINT32_MAX;
}

// The component a point stands for, -1 for none
static inline int32_t cartoon_glyphs_point_component(const CartoonGlyphs& glyphs, uint32_t point_idx) {
    return point_idx < glyphs.num_points ? glyphs.point_comp[point_idx] : -1;
}
