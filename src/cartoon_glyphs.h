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
// (md_gl). Lipids and nucleic acids have a cheaper one, of spheres and capsules: a sphere at every point and a capsule
// between two.
//
//   lipid      A sphere for each head, at the centre of its heavy atoms and sized by how many there are (a phosphatidyl
//              choline is larger than the hydroxyl of a sterol, an inositol phosphate larger still). A tube for each tail
//              (md_lipid_data_t) from the centre of its head over the atom it hangs off and along its spine to its tip,
//              or on to the head a spanning chain joins (a bolalipid). An atomistic spine gets a point per two of its
//              atoms, which is on the axis of a chain in its all trans zigzag and so smooths it away; a coarse grained
//              one a point per bead. The skeleton of a sterol is a thicker tube.
//   nucleic    A tube through the C4' of each nucleotide of a backbone (md_nucleic_backbone_data_t), and a thinner stick
//              from the C4' over the C1' to the edge of the base which pairs: N1 of a purine, N3 of a pyrimidine.
//
// The points are the particles of a molecule of their own (md_gl_mol_create_particles), placed every frame from the
// atoms they stand for, and drawn as spheres (MD_GL_REP_SPACE_FILL); the capsules are the weak bonds of each
// representation (md_gl_rep_set_weak_bonds, MD_GL_REP_WEAK_BONDS, one dash which fills it). A point stands for a
// component, which is what picking it resolves to, takes its colour from one of its atoms, and is shown when any of its
// atoms is shown.
//
// What comes after: splines with extruded profiles for the backbones of nucleic acids (an oval) and for the tails of
// lipids (a circle), as for the proteins. The points here, in order along each path, are what such a spline goes
// through.

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
    uint32_t num_segments = 0;

    // A point is the mean of its atoms, point_atom[point_offset[i] .. point_offset[i + 1]), each atom taken at its
    // image nearest to point_ref[i]. The points of one lipid share their reference, so a lipid stays whole across the
    // periodic boundary.
    md_array(uint32_t) point_offset = nullptr;       // [num_points + 1]
    md_array(int32_t)  point_atom = nullptr;
    md_array(int32_t)  point_ref = nullptr;          // [num_points]
    md_array(int32_t)  point_color_atom = nullptr;   // [num_points] The atom whose colour it takes
    md_array(int32_t)  point_comp = nullptr;         // [num_points] The component it stands for
    md_array(float)    point_radius = nullptr;       // [num_points] Of its sphere, Å, before the scale of the representation
    md_array(uint8_t)  point_kind = nullptr;         // [num_points] CartoonGlyphKind

    md_array(uint32_t) segment_point = nullptr;      // [2 * num_segments] The points at the ends of each capsule
    md_array(float)    segment_radius = nullptr;     // [num_segments] Å, before the scale of the representation

    md_array(vec3_t)   xyz = nullptr;                // [num_points] Where the points are
    bool positions_valid = false;                    // xyz and the positions of gl_mol are those of the current state
    bool velocity_valid = false;                     // They were valid before the state moved, so the move is a velocity

    md_allocator_i* alloc = nullptr;
};

// The glyphs of the lipids and nucleic acid backbones of the system, none when it has neither (cartoon_glyphs_valid)
void cartoon_glyphs_init(CartoonGlyphs* glyphs, const md_system_t& sys, md_allocator_i* alloc);
void cartoon_glyphs_free(CartoonGlyphs* glyphs);

static inline bool cartoon_glyphs_valid(const CartoonGlyphs& glyphs) {
    return glyphs.num_points > 0 && glyphs.gl_mol.id != 0;
}

// Places the points from the coordinates of the state and uploads them. With velocity, the motion from where they were
// is their velocity (for the temporal reprojection); without it they stand still.
void cartoon_glyphs_update_positions(CartoonGlyphs* glyphs, const md_system_state_t& state, bool velocity);

// The flags of the points (highlighted, selected, visible, see AtomBit_): those of any of their atoms
void cartoon_glyphs_update_flags(const CartoonGlyphs& glyphs, const uint8_t* atom_flags, size_t num_atoms);

// A representation of the glyphs, with their capsules, for a representation of the system. Destroy with md_gl_rep_destroy.
md_gl_rep_t cartoon_glyphs_rep_create(const CartoonGlyphs& glyphs);

// The colours of the points from those of the atoms of a representation (an alpha of 0 hides an atom). Whether any point
// is shown.
bool cartoon_glyphs_rep_set_colors(const CartoonGlyphs& glyphs, md_gl_rep_t rep, const uint32_t* atom_colors, size_t num_atoms);

// The component a point stands for, -1 for none
static inline int32_t cartoon_glyphs_point_component(const CartoonGlyphs& glyphs, uint32_t point_idx) {
    return point_idx < glyphs.num_points ? glyphs.point_comp[point_idx] : -1;
}
