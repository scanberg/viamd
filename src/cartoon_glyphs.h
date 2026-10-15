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
// profile. Nucleic acids have a cheaper one for now, of spheres and capsules.
//
//   lipid      A sphere for each head, at the centre of its heavy atoms and sized by how many there are (a phosphatidyl
//              choline is larger than the hydroxyl of a sterol, an inositol phosphate larger still). A tube for each tail
//              (md_lipid_data_t) from the centre of its head over the atom it hangs off and along its spine to its tip,
//              or on to the head a spanning chain joins (a bolalipid). An atomistic spine gets a control point per two
//              of its atoms, which is on the axis of a chain in its all trans zigzag and so smooths it away; a coarse
//              grained one a control point per bead. The skeleton of a sterol is a thicker tube.
//   nucleic    A capsule from the C4' of each nucleotide of a backbone (md_nucleic_backbone_data_t) to the next, and a
//              thinner stick from the C4' over the C1' to the edge of the base which pairs: N1 of a purine, N3 of a
//              pyrimidine. A sphere at each of these points.
//
// The points (the centres of the spheres, the control points of the tubes, the ends of the capsules) are the particles of
// a molecule of their own (md_gl_mol_create_particles). Each is the mean of the atoms it stands for, derived on the GPU
// from the positions of the atoms of the system (md_gl_mol_set_derived, md_gl_mol_compute_derived), each atom in its
// periodic image nearest to the first head atom of its lipid, so a lipid stays whole across the periodic boundary. The
// spheres are drawn as such (MD_GL_REP_SPACE_FILL, a control point of a tube has no sphere: a radius of 0), the tubes as
// the splines of md_gl through the control points (md_gl_mol_set_tubes, MD_GL_REP_TUBES) and the capsules as the weak
// bonds of each representation (md_gl_rep_set_weak_bonds, MD_GL_REP_WEAK_BONDS, one dash which fills it). A point stands
// for a component, which is what picking it resolves to, takes its colour from one of its atoms, and is shown when any
// of its atoms is shown. A part of a tube is that of its control point.
//
// What comes after: a spline with an oval profile for the backbones of nucleic acids.

enum CartoonGlyphKind : uint8_t {
    CartoonGlyph_LipidHead = 0,
    CartoonGlyph_LipidTail,
    CartoonGlyph_SterolSkeleton,
    CartoonGlyph_NucleicBackbone,
    CartoonGlyph_NucleicBase,
    CartoonGlyph_Count
};

struct CartoonGlyphs {
    md_gl_mol_t gl_mol = {};
    uint32_t num_points = 0;
    uint32_t num_tubes = 0;
    uint32_t num_tube_points = 0;   // Control points of all tubes
    uint32_t num_segments = 0;
    uint32_t num_source_atoms = 0;  // Of the system the points are derived from

    // A point is the mean of its atoms, point_atom[point_offset[i] .. point_offset[i + 1]), each atom taken at its
    // image nearest to point_ref[i]. The points of one lipid share their reference, so a lipid stays whole across the
    // periodic boundary.
    md_array(uint32_t) point_offset = nullptr;       // [num_points + 1]
    md_array(int32_t)  point_atom = nullptr;
    md_array(int32_t)  point_ref = nullptr;          // [num_points]
    md_array(int32_t)  point_color_atom = nullptr;   // [num_points] The atom whose colour it takes
    md_array(int32_t)  point_comp = nullptr;         // [num_points] The component it stands for
    md_array(float)    point_radius = nullptr;       // [num_points] Of its sphere, Å, before the scale of the representation. 0 for none
    md_array(uint8_t)  point_kind = nullptr;         // [num_points] CartoonGlyphKind

    // Tube t goes through the points tube_point[tube_offset[t] .. tube_offset[t + 1])
    md_array(uint32_t) tube_offset = nullptr;        // [num_tubes + 1]
    md_array(uint32_t) tube_point = nullptr;         // [num_tube_points]
    md_array(float)    tube_radius = nullptr;        // [num_tube_points] Å, before the scale of the representation

    md_array(uint32_t) segment_point = nullptr;      // [2 * num_segments] The points at the ends of each capsule
    md_array(float)    segment_radius = nullptr;     // [num_segments] Å, before the scale of the representation

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

// A representation of the glyphs, with their capsules, for a representation of the system. Destroy with md_gl_rep_destroy.
md_gl_rep_t cartoon_glyphs_rep_create(const CartoonGlyphs& glyphs);

// The colours of the points from those of the atoms of a representation (an alpha of 0 hides an atom). Whether any point
// is shown.
bool cartoon_glyphs_rep_set_colors(const CartoonGlyphs& glyphs, md_gl_rep_t rep, const uint32_t* atom_colors, size_t num_atoms);

// Picking: the spheres and capsules write the index of a point, the tubes num_points + the index of their control point
// (md_gl_picking_offset_t: atom_base at the beginning of the range, backbone_base num_points after it)
static inline uint32_t cartoon_glyphs_picking_count(const CartoonGlyphs& glyphs) {
    return glyphs.num_points + glyphs.num_tube_points;
}

// The point a picking index stands for, UINT32_MAX for none
static inline uint32_t cartoon_glyphs_picked_point(const CartoonGlyphs& glyphs, uint32_t picking_idx) {
    if (picking_idx < glyphs.num_points) return picking_idx;
    if (picking_idx - glyphs.num_points < glyphs.num_tube_points) return glyphs.tube_point[picking_idx - glyphs.num_points];
    return UINT32_MAX;
}

// The component a point stands for, -1 for none
static inline int32_t cartoon_glyphs_point_component(const CartoonGlyphs& glyphs, uint32_t point_idx) {
    return point_idx < glyphs.num_points ? glyphs.point_comp[point_idx] : -1;
}
