#pragma once

#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_lru_cache.inl>
#include <core/md_hash.h>
#include <core/md_array.h>

#include <md_system.h>
#include <md_element.h>
#include <color_utils.h>
#include <md_run.h>
#include <plane.h>
#include <clock.h>
#include <transfer_function_widget.h>
#include <md_gto.h>
#include <core/md_grid.h>
#include <script_engine.h>
#include <string>
#include <md_gl.h>
#include <md_infer.h>
#include <md_chem.h>
#include <md_edr.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#include <core/md_vec_math.h>
#endif

#include <app/IconsFontAwesome6.h>
#include <app/application.h>
#include <gfx/camera.h>
#include <gfx/camera_utils.h>
#include <gfx/gbuffer.h>
#include <gfx/view_param.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <task_system.h>
#include <surface_field.h>
#include <glyphs.h>
#include <cartoon_glyphs.h>
#include <loader.h>
#include <file_list.h>
#include <drop_files.h>
#include <event.h>
#include <viamd_event.h>
#include <plot_series.h>
#include <timeline_periods.h>
#include <filter_input.h>
#include <representation_tree.h>

#define IMGUI_DEFINE_MATH_OPERATORS

#include <script_editor.h>
#include <implot.h>
#include <imgui_notify.h>

#include <bitset>
#include <stdint.h>
#include <stddef.h>

#define JITTER_SEQUENCE_MAX 128 // The jitter's positions at most: temporal_aa.jitter_samples, or more when upscaling (jitter_count)
#define DEFAULT_COLORMAP 5      // This corresponds to plasma colormap (Do not want to include implot.h just for this)
#define FRAME_CACHE_SIZE 4      // Number of frames used in the frame cache
#define PERIODIC_IMAGES_MAX 3   // Periodic images each way along a cell vector: (2 n + 1)^3 - 1 copies of the scene at most

#if FRAME_CACHE_SIZE == 4
#define FRAME_CACHE_LRU_TYPE md_lru_cache4_t
#elif FRAME_CACHE_SIZE == 8
#define FRAME_CACHE_LRU_TYPE md_lru_cache8_t
#else
#error "Unsupported frame cache size"
#endif

// Profiling: PUSH_GPU_SECTION / POP_GPU_SECTION / PROFILE_SCOPE (gfx/profiler.h). A section is a debug group and is
// timed on the GPU and the CPU; one that issues no GL work is timed on the CPU only, in effect.
#include <gfx/profiler.h>
#define PUSH_CPU_SECTION(lbl) PUSH_GPU_SECTION(lbl)
#define POP_CPU_SECTION()     POP_GPU_SECTION()

// For logging
#define VIAMD_LOG_INFO  MD_LOG_INFO
#define VIAMD_LOG_DEBUG MD_LOG_DEBUG
#define VIAMD_LOG_ERROR MD_LOG_ERROR
#define VIAMD_LOG_SUCCESS(...) ImGui::InsertNotification(ImGuiToast(ImGuiToastType_Success, 6000, __VA_ARGS__))

// One colour per script property (by its place in the script), shared by every view of it
inline constexpr uint32_t PROPERTY_COLORS[] = {4293119554, 4290017311, 4287291314, 4281114675, 4288256763, 4280031971, 4285513725, 4278222847, 4292260554, 4288298346, 4288282623, 4280834481};

#define INVALID_PICKING_IDX (~0U)

// The id of no system (SystemContext::id): of a picking range or a hit that is of none
#define SYSTEM_ID_NONE (~0U)

constexpr ImGuiKey KEY_PLAY_PAUSE               = ImGuiKey_Space;
constexpr ImGuiKey KEY_SKIP_TO_PREV_FRAME       = ImGuiKey_LeftArrow;
constexpr ImGuiKey KEY_SKIP_TO_NEXT_FRAME       = ImGuiKey_RightArrow;
constexpr ImGuiKey KEY_RECOMPILE_SHADERS        = ImGuiKey_F5;
constexpr ImGuiKey KEY_SHOW_DEBUG_WINDOW        = ImGuiKey_F11;
constexpr ImGuiKey KEY_RECENTER_ON_HIGHLIGHT    = ImGuiKey_F2;

constexpr ImGuiKeyChord KEY_SCRIPT_EVALUATE          = ImGuiMod_Shift | ImGuiKey_Enter;
constexpr ImGuiKeyChord KEY_SHOW_ATTRIBUTE_WINDOW = ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_A;

constexpr md_str_t WORKSPACE_FILE_EXTENSION = MD_STR_INIT("via");

typedef uint64_t PickingDomainID;
typedef uint64_t PickingSourceID;
typedef uint64_t InteractionSurfaceID;

constexpr PickingDomainID PickingDomain_Atom   = HASH_STR_LIT64("picking domain atom");
// A bond: what the middle half of a bond picks as, the quarters at its ends picking as their atoms, however it is drawn.
// Of the system (key 0): a bond of the graph, or one formed in the frame shown after them (picked_bond). Or a weak bond
// of a representation (key: its id), in a range of its own: one of a Bonds representation, an interaction over the
// sticks of a licorice.
constexpr PickingDomainID PickingDomain_Bond   = HASH_STR_LIT64("picking domain bond");
constexpr PickingDomainID PickingDomain_Dipole = HASH_STR_LIT64("picking domain dipole");
// A protein backbone segment (md_system_t::protein_backbone.segment): what the cartoon and ribbons pick as. They are
// drawn per segment, not per atom, so a hit names the segment and resolves to the residue it stems from.
constexpr PickingDomainID PickingDomain_BackboneSegment = HASH_STR_LIT64("picking domain backbone segment");
// A ring of the system (md_system_t::ring): what the circle of an aromatic ring picks as (licorice, ball and stick with
// bond orders). A hit names the ring, ring(local_idx + 1) in the script, and resolves to its atoms.
constexpr PickingDomainID PickingDomain_Ring = HASH_STR_LIT64("picking domain ring");
// A point of the cartoon of the lipids and nucleic acids (CartoonGlyphs): a head, a point along a tail or a backbone, or
// the part of a tube around a control point. A hit names the point (cartoon_glyphs_picked_point) and resolves to the
// residue it stands for.
constexpr PickingDomainID PickingDomain_CartoonGlyph = HASH_STR_LIT64("picking domain cartoon glyph");
// The gap of a chain over a break (md_gl_mol_set_chain_gaps): the dashed line where residues are missing. A hit names
// the residue before the break (MD_RESIDUE_FLAG_BREAK_AFTER) and resolves to it and the one after it.
constexpr PickingDomainID PickingDomain_ChainGap = HASH_STR_LIT64("picking domain chain gap");
// A point of a point set drawn as glyphs (glyphs.h), keyed by the representation: what its glyph picks as. The index is
// the point.
constexpr PickingDomainID PickingDomain_Glyph = HASH_STR_LIT64("picking domain point glyph");

// The domains whose hits stand for atoms of the system: hovering one highlights them, and a click commits the highlight
// to the selection
static inline bool picking_domain_resolves_to_atoms(PickingDomainID domain) {
    return domain == PickingDomain_Atom || domain == PickingDomain_Bond || domain == PickingDomain_BackboneSegment || domain == PickingDomain_Ring ||
           domain == PickingDomain_CartoonGlyph || domain == PickingDomain_ChainGap;
}

constexpr uint64_t interaction_surface_main = HASH_STR_LIT64("interaction surface main"); // This is the main interaction surface which corresponds to the main interaction window, but we want to keep it separate from the picking source and domain ids as we may want to have different picking sources/domains for different interaction surfaces in the future

enum class PlaybackMode { Stopped, Playing, Count };
// Where the transport is over the view: at its top or its bottom, in the middle of the part of the view uncovered by
// windows, or not shown. The values are those a workspace stores ([Animation] Transport), so a new one goes last.
enum TransportAnchor : int { TransportAnchor_Top = 0, TransportAnchor_Bottom, TransportAnchor_Hidden, TransportAnchor_Count };
// Part: the part of the atom (md_part_data_t), a covalent piece of its molecule: a strand of a polymer with its caps, a
// glycan, a ligand. Chain: the chain of the atom (md_chain_data_t), or its part for one in no chain (the data gives no
// chain ids: the strand it is of). Molecule: the molecule of the atom (md_molecule_data_t), atoms joined by bonds with
// their residues and polymer strands whole. Fragment: the fragment of the atom (md_fragment.h), a ring system, a
// functional group, a backbone, an alkyl chain. The values are those a workspace stores, so a new one goes last;
// selection_granularity_order is the order to offer them in, finest first.
enum class SelectionGranularity { Atom, Residue, Part, Molecule, Fragment, Count };
enum class SelectionOperator { Or, And, AndNot, Set, Clear, Count };
enum class SelectionGrowthMode { CovalentBond, Radial, Count };
enum class CameraMode { Perspective, Orthographic, Count };
enum class InterpolationMode { Nearest, Linear, CubicSpline, BSpline, Count };
// TAA: each frame blended with the ones before; FXAA: each frame on its own (nothing carried over time)
enum class AntiAliasingMode { TAA, FXAA, Count };

// Stored in workspaces by value: append, never insert or reorder. The order shown is color_mapping_order.
enum class ColorMapping {
    Uniform,
    Type,               // the colour of the particle type
    Serial,             // one per atom index
    ResName,            // one per residue name, of no meaning
    ResSeqId,
    ResIndex,
    ChainId,            // one per chain id the data gives (md_chain_data_t)
    ChainIndex,         // one per chain, by its position
    SecondaryStructure,
    Attribute,          // a per atom attribute of the system, through a colour scale
    ResType,            // the colour of the residue type, from its standard monomer unless overridden (md_residue_type_data_t)
    PartIndex,          // one per part (md_part_data_t), by its position: each strand, ligand and water
    Count
};

enum class ElectronicStructureSource {
    MolecularOrbital,
    NaturalTransitionOrbital,
    TransitionDensity,
    ElectronDensity,
    DensityProperty,
    ElectronLocalization,       // the electron localization function, ELF (md_qm_elf)
    LocalizedOrbitalLocator,    // the localized orbital locator, LOL (md_qm_lol)
    // Fields of the system as a whole, evaluated as the fields a surface is coloured by are (surface_field.h, the kind
    // of electronic_structure_source_field_kind): the volume all of the grid, on the GPU where mdlib evaluates them
    // there, else on the CPU
    ElectrostaticPotential,
    EmbeddingPotential,
    AverageLocalIonizationEnergy,
    LocalElectronAttachmentEnergy,
    DensityScaleLength,
    GostshypWallPotential,      // of the GOSTSHYP wall of a run under pressure (gostshyp.h)
    GostshypEffectivePotential,
    GostshypEnergyDensity,
    Count
};

enum class ElectronicStructureField {
    Amplitude,
    Density,
    Count
};

// What colours an isosurface: its own colour, the colours of the atoms blended over it, or a field
// evaluated on it and shown through a colour map (surface_field.h)
enum class SurfaceColoring {
    Uniform,
    AtomColors,
    Field,
    Count
};

inline const char* surface_coloring_str[(int)SurfaceColoring::Count] = {
    "Uniform",
    "Atom Colors",
    "Field",
};

enum class ElectronicStructureSpin {
    None,
    Alpha,
    Beta,
    Total,
    Difference,
    Count
};

// These bits are a compressed form of flags which are passed onto rendering as the rendering only supports 8-bits
enum AtomBit_ {
    AtomBit_None        = 0,
    AtomBit_Highlighted = 1u << 0,
    AtomBit_Selected    = 1u << 1,
    AtomBit_Visible     = 1u << 2,
    AtomBit_Marked      = 1u << 3,  // Highlighted or selected: the atoms the mark pass draws (md_gl's atom mask takes all its bits)
    AtomBit_Dimmed      = 1u << 4,  // Not selected while something is: drawn with less of its colour (selection.color.saturation)
    // The group of the highlight an atom is in, of HIGHLIGHT_GROUP_COLORS (selection.highlight_groups): which colour it
    // is marked in. 0 when the highlight has no groups.
    AtomBit_GroupShift  = 5,
    AtomBit_GroupMask   = 7u << AtomBit_GroupShift,
};

// The colours the groups of the highlight are marked in, in turn: two groups fewer than this apart differ in colour
#define HIGHLIGHT_GROUP_COLORS 8

enum MolBit_ {
    MolBit_None                     = 0,
    MolBit_DirtyPosition            = 1u << 0,
    MolBit_DirtyRadius              = 1u << 1,
    MolBit_DirtyFlags               = 1u << 3,
    MolBit_DirtyBonds               = 1u << 4,
    MolBit_ClearVelocity            = 1u << 5,
    MolBit_ResetBackboneHistory     = 1u << 6,  // The positions jumped: the backbone orientation starts over (md_gl_mol_reset_backbone_history)
};

enum class RepresentationType {
    SpaceFill = MD_GL_REP_SPACE_FILL,
    Licorice = MD_GL_REP_LICORICE,
    BallAndStick = MD_GL_REP_BALL_AND_STICK,
    Ribbons = MD_GL_REP_RIBBONS,
    Cartoon = MD_GL_REP_CARTOON,
    ElectronicStructure,
    DipoleMoment,
    SimulationBox,
    Bonds,              // The bonds its filter picks (script::filter_bonds), as sticks or dashed, without the atoms
    Glyphs,             // A glyph per point of a point set of the attributes (glyphs.h): a cavity's tesserae, ...
    Count
};

inline const char* selection_granularity_str[(int)SelectionGranularity::Count] = {
    "Atom",
    "Residue",
    "Part",
    "Molecule",
    "Fragment",
};

inline const SelectionGranularity selection_granularity_order[(int)SelectionGranularity::Count] = {
    SelectionGranularity::Atom,
    SelectionGranularity::Fragment,
    SelectionGranularity::Residue,
    SelectionGranularity::Part,
    SelectionGranularity::Molecule,
};

inline const char* interpolation_mode_str[(int)InterpolationMode::Count] = {
    "Nearest",
    "Linear",
    "Cubic Spline",
    "B-Spline",
};

inline const char* representation_type_str[(int)RepresentationType::Count] = {
    "Spacefill",
    "Licorice",
    "Ball And Stick",
    "Ribbons",
    "Cartoon",
    "Field",
    "Dipole Moment",
    "Simulation Box",
    "Bonds",
    "Glyphs",
};

// In the words of the files: residue, and chain for the chains the data gives (md_chain_data_t); a part is a strand, a
// ligand or a water (md_part_data_t), which tells the strands of a system without chain ids apart
inline const char* color_mapping_str[(int)ColorMapping::Count] = {
    "Uniform",
    "Atom Type",
    "Atom Index",
    "Res Name",
    "Res Id",
    "Res Index",
    "Chain Id",
    "Chain Index",
    "Secondary Structure",
    "Attribute",
    "Res Type",
    "Part Index",
};

// The order the colour mappings are listed in
inline const ColorMapping color_mapping_order[(int)ColorMapping::Count] = {
    ColorMapping::Uniform,
    ColorMapping::Type,
    ColorMapping::Serial,
    ColorMapping::ResName,
    ColorMapping::ResType,
    ColorMapping::ResSeqId,
    ColorMapping::ResIndex,
    ColorMapping::ChainId,
    ColorMapping::ChainIndex,
    ColorMapping::PartIndex,
    ColorMapping::SecondaryStructure,
    ColorMapping::Attribute,
};


inline const char* electronic_structure_source_str[(int)ElectronicStructureSource::Count] = {
    "Molecular Orbital",
    "Natural Transition Orbital",
    "Transition Density",
    "Electron Density",
    "Density Property",
    "Electron Localization (ELF)",
    "Localized Orbital Locator (LOL)",
    "Electrostatic Potential",
    "Embedding Potential",
    "Local Ionization Energy (ALIE)",
    "Local Electron Attachment Energy",
    "Density Scale Length",
    "GOSTSHYP Wall Potential",
    "GOSTSHYP Effective Potential",
    "GOSTSHYP Energy Density",
};

// The kind of surface field a source is evaluated as (surface_field.h): Count for those evaluated of orbitals and
// densities on the GPU (viamd.cpp, es_evaluate_into). ELF and LOL are both: by their own kernel on the GPU when
// there is one, else as a surface field on the CPU.
static inline SurfaceFieldKind electronic_structure_source_field_kind(ElectronicStructureSource s) {
    switch (s) {
    case ElectronicStructureSource::ElectronLocalization:          return SurfaceFieldKind::ElectronLocalization;
    case ElectronicStructureSource::LocalizedOrbitalLocator:       return SurfaceFieldKind::LocalizedOrbitalLocator;
    case ElectronicStructureSource::ElectrostaticPotential:        return SurfaceFieldKind::ElectrostaticPotential;
    case ElectronicStructureSource::EmbeddingPotential:            return SurfaceFieldKind::EmbeddingPotential;
    case ElectronicStructureSource::AverageLocalIonizationEnergy:  return SurfaceFieldKind::AverageLocalIonizationEnergy;
    case ElectronicStructureSource::LocalElectronAttachmentEnergy: return SurfaceFieldKind::LocalElectronAttachmentEnergy;
    case ElectronicStructureSource::DensityScaleLength:            return SurfaceFieldKind::DensityScaleLength;
    case ElectronicStructureSource::GostshypWallPotential:         return SurfaceFieldKind::GostshypWallPotential;
    case ElectronicStructureSource::GostshypEffectivePotential:    return SurfaceFieldKind::GostshypEffectivePotential;
    case ElectronicStructureSource::GostshypEnergyDensity:         return SurfaceFieldKind::GostshypEnergyDensity;
    default:                                                       return SurfaceFieldKind::Count;
    }
}

// Whether the source is evaluated as a surface field (surface_field.h) whatever the GPU: all but the orbitals,
// the densities and ELF and LOL
static inline bool electronic_structure_source_is_surface_field(ElectronicStructureSource s) {
    return electronic_structure_source_field_kind(s) != SurfaceFieldKind::Count &&
           s != ElectronicStructureSource::ElectronLocalization && s != ElectronicStructureSource::LocalizedOrbitalLocator;
}

// Whether isosurface meshes can be extracted from the source (iso_meshes.h): the orbitals and densities, which
// md_gto_iso evaluates, ELF and LOL, which it bounds through the density and the kinetic energy density, and the
// GOSTSHYP potentials, orbitals over the tesserae (gostshyp_wall_orbital). The others' isosurfaces are raycast
// through their volume, whatever the setting.
static inline bool electronic_structure_source_meshable(ElectronicStructureSource s) {
    return s == ElectronicStructureSource::MolecularOrbital || s == ElectronicStructureSource::NaturalTransitionOrbital ||
           s == ElectronicStructureSource::TransitionDensity || s == ElectronicStructureSource::ElectronDensity ||
           s == ElectronicStructureSource::DensityProperty || s == ElectronicStructureSource::ElectronLocalization ||
           s == ElectronicStructureSource::LocalizedOrbitalLocator || s == ElectronicStructureSource::GostshypWallPotential ||
           s == ElectronicStructureSource::GostshypEffectivePotential;
}

enum ElectronicStructureSourceFlag_ : uint32_t {
    ElectronicStructureSourceFlag_MolecularOrbital          = 1u << (int)ElectronicStructureSource::MolecularOrbital,
    ElectronicStructureSourceFlag_NaturalTransitionOrbital  = 1u << (int)ElectronicStructureSource::NaturalTransitionOrbital,
    ElectronicStructureSourceFlag_TransitionDensity         = 1u << (int)ElectronicStructureSource::TransitionDensity,
    ElectronicStructureSourceFlag_ElectronDensity           = 1u << (int)ElectronicStructureSource::ElectronDensity,
    ElectronicStructureSourceFlag_DensityProperty           = 1u << (int)ElectronicStructureSource::DensityProperty,
    ElectronicStructureSourceFlag_ElectronLocalization      = 1u << (int)ElectronicStructureSource::ElectronLocalization,
    ElectronicStructureSourceFlag_LocalizedOrbitalLocator   = 1u << (int)ElectronicStructureSource::LocalizedOrbitalLocator,
    ElectronicStructureSourceFlag_ElectrostaticPotential    = 1u << (int)ElectronicStructureSource::ElectrostaticPotential,
    ElectronicStructureSourceFlag_EmbeddingPotential        = 1u << (int)ElectronicStructureSource::EmbeddingPotential,
    ElectronicStructureSourceFlag_AverageLocalIonizationEnergy  = 1u << (int)ElectronicStructureSource::AverageLocalIonizationEnergy,
    ElectronicStructureSourceFlag_LocalElectronAttachmentEnergy = 1u << (int)ElectronicStructureSource::LocalElectronAttachmentEnergy,
    ElectronicStructureSourceFlag_DensityScaleLength        = 1u << (int)ElectronicStructureSource::DensityScaleLength,
    ElectronicStructureSourceFlag_GostshypWallPotential     = 1u << (int)ElectronicStructureSource::GostshypWallPotential,
    ElectronicStructureSourceFlag_GostshypEffectivePotential= 1u << (int)ElectronicStructureSource::GostshypEffectivePotential,
    ElectronicStructureSourceFlag_GostshypEnergyDensity     = 1u << (int)ElectronicStructureSource::GostshypEnergyDensity,
};

typedef uint32_t ElectronicStructureSourceFlags;

inline const char* electronic_structure_spin_str[(int)ElectronicStructureSpin::Count] = {
    "None",
    "Alpha",
    "Beta",
    "Total",
    "Difference",
};

enum class ElectronicStructureNtoComponent {
    Particle,
    Hole,
    Count
};

inline const char* electronic_structure_nto_component_str[(int)ElectronicStructureNtoComponent::Count] = {
    "Particle",
    "Hole",
};

enum class ElectronicStructureTransitionDensityComponent {
    Attachment,
    Detachment,
    Difference,
    Count
};

inline const char* electronic_structure_transition_density_component_str[(int)ElectronicStructureTransitionDensityComponent::Count] = {
    "Attachment",
    "Detachment",
    "Difference"
};

enum class ElectronicStructureLegacyType {
    MolecularOrbital,
    MolecularOrbitalDensity,
    NaturalTransitionOrbitalParticle,
    NaturalTransitionOrbitalHole,
    NaturalTransitionOrbitalDensityParticle,
    NaturalTransitionOrbitalDensityHole,
    AttachmentDensity,
    DetachmentDensity,
    ElectronDensity,
    Count,
};

struct LoadDatasetWindowState {
    char path_buf[1024] = "";
    char atom_format_buf[128] = "";
    char err_buf[128] = "";
    bool path_is_valid = false;
    bool path_changed = false;
    bool load_topology = false;
    bool load_trajectory = false;
    bool show_window = false;
    bool show_file_dialog = false;
    bool atom_format_valid = false;
    int  loader_idx = 0;
    int  atom_format_idx = 0;
    // Where it goes: into the system with this id, or into a new one (asked when several are loaded)
    uint32_t system = SYSTEM_ID_NONE;
    bool new_system = false;
};

// Hint flags for operations to be performed for the files
enum {
    FileFlags_None = 0,
    FileFlags_ShowDialogue = 1,
    FileFlags_DisableCacheWrite = 4,
    FileFlags_NewSystem = 8,        // Into a system of its own, added for it
    FileFlags_Dropped = 16,         // Dropped on the window: where it was dropped is told the frame after (FileQueue)
};

typedef uint32_t FileFlags;

// Bohr is what md_gto evaluates in, and what a GTO basis' exponents are expressed against; Angstrom
// is what the world, the camera and a system's coordinates are in. Every conversion between the two
// goes through these, so there is one value and not one per file that happened to need it.
inline constexpr double ANGSTROM_TO_BOHR = 1.8897261246257702;
inline constexpr double BOHR_TO_ANGSTROM = 0.5291772109029999;

enum class VolumeResolution {
    Low,
    Mid,
    High,
    Count,
};

inline const char* volume_resolution_str[(int)VolumeResolution::Count] = {
    "Low",
    "Mid",
    "High",
};


// How an electronic structure representation draws its field
enum class ElectronicStructureDisplay {
    Isosurface,     // surfaces where the field takes given values
    Volume,         // direct volume rendering through a transfer function, cut by planes
    Planes,         // the field on a set of planes, each of its own opacity
    Count
};

inline const char* electronic_structure_display_str[(int)ElectronicStructureDisplay::Count] = {
    "Isosurface",
    "Volume",
    "Planes",
};

#define ES_MAX_PLANES 8

// The spacing of the planes display's images (each plane evaluated on itself)
enum class PlaneDetail {
    Coarse,
    Fine,
    Finest,
    Count,
};

// Angstrom per pixel of a plane's image for each PlaneDetail; an image is at most GPU_PLANE_DIM pixels a side, the
// spacing made coarser where the plane is larger than that
inline const float plane_detail_spacing_angstrom[(int)PlaneDetail::Count] = {
    0.05f,
    0.02f,
    0.01f,
};

// A plane of an electronic structure representation: a cut of its volume display, or a slice of its planes display
struct ElectronicStructurePlane {
    Plane plane = {};
    bool  enabled = true;       // a cut cuts, a slice is drawn
    float opacity = 1.0f;       // a slice's
    // A slice's rectangle on its plane (a cut has none: it is the plane through the box): its centre, its first side's
    // direction (u, unit, along the plane; v = normal x u) and its half sizes along u and v (Angstrom). Not limited
    // to the box: outside it the volume's values at the box's faces carry on. An extent of 0: not set yet (the
    // rectangle around where the plane cuts the box, once the box is there: es_slices_conform).
    md_float3_t center = {};
    md_float3_t axis = {};
    md_float2_t extent = {};
};

struct ElectronicStructurePlaneSet {
    ElectronicStructurePlane items[ES_MAX_PLANES] = {};
    int count = 0;
};

// Samples per Angstrom for each VolumeResolution. Beside the enum, because a value per enumerator
// is part of what the enum MEANS and a second table elsewhere is how the two drift apart.
inline const double volume_resolution_samples_per_angstrom[(int)VolumeResolution::Count] = {
    4.0,
    8.0,
    16.0,
};

// Primitives whose contribution falls below this are pruned when a basis is expanded. One value for
// the whole application: the uploaded basis is built once per dataset with it, so an evaluation
// asking for a different one would silently get this one anyway.
#define DEFAULT_GTO_CUTOFF_VALUE 1.0e-6

enum class ScreenshotResolution {
    View,       // the main view's, as it is shown
    FHD,
    QHD,
    UHD_4K,
    UHD_8K,
    Custom,
    Count,
};

inline const char* screenshot_resolution_str[(int)ScreenshotResolution::Count] = {
    "View",
    "Full HD (1920x1080)",
    "Quad HD (2560x1440)",
    "Ultra HD 4K (3840x2160)",
    "Ultra HD 8K (7680x4320)",
    "Custom",
};

enum class BondColorMode {
    NearestAtom,
	SmoothAtom,
    Uniform,
    Count,
};

// Which bonds a Bonds representation draws dashed, each in the style of its kind (md_gl_dashes_t, BOND DASHES in md_gl.h)
enum class BondDashes {
    ByKind,     // Those which are not covalent: coordinate bonds, hydrogen bonds, interactions
    None,
    All,
    Count,
};

inline const char* bond_color_mode_str[(int)BondColorMode::Count] = {
    "Nearest Atom",
    "Smooth Atom",
    "Uniform Color",
};

struct FileQueue {
    struct Entry {
        md_str_t path;
        FileFlags flags;
        int prio;
        uint32_t system;    // The id of the system it is loaded into, SYSTEM_ID_NONE for where it goes (file_queue_process)
    };
    Entry arr[64] = {};     // A power of two: the ring's indices wrap with the unsigned arithmetic
    uint32_t head = 0;
    uint32_t tail = 0;

    // We use a ring alloc here to do a deep copy of the paths within the added entries
    // Because its a ring alloc, there is no need to free entries
    md_allocator_i* ring = NULL;

    // Files dropped on the window go where they were dropped: told while the windows of the frame they were dropped in
    // are drawn (the Scene window: onto a system's rows, into it; below them, a new system), and taken by
    // file_queue_process the frame after. Dropped anywhere else, they are placed by the rules.
    int      drop_frame = -1;               // The frame told in (ImGui's frame count)
    uint32_t drop_system = SYSTEM_ID_NONE;  // Onto this system
    bool     drop_new = false;              // Into a new system
};

// What load_data made of a file
enum LoadResult {
    LoadResult_Failed = 0,
    LoadResult_System,      // A new system (built onto one without atoms), with its run when the file holds several frames
    LoadResult_Joined,      // Of the system's atoms, joined to it (md_system_merge): its coordinates, names, bonds, ...
    LoadResult_Trajectory,  // A new run for the system
    LoadResult_Supplement,  // Data of the system's atoms, or along its run: the system and the run are as they were
};

struct SingleSelectionSequence {
    int32_t idx[4] = {-1, -1, -1, -1};
    uint32_t system = SYSTEM_ID_NONE;   // The id of the system the atoms are of: a pick in another starts it over
};

// One dipole group the producer published, as two attributes sharing an index space:
//     dipole/<group>/vector    the moment, in `unit`
//     dipole/<group>/origin    where to draw it from, Angstrom
// A group holds ONE moment (a ground state) or one per excited state, so (key, element index)
// addresses a dipole and key alone does not - which is why a representation stores both.
//
// This describes the group and carries no values: reading one moment is dipole_moment_read, and it
// takes the same two numbers a representation and a picking hit already hold. Nothing owns a
// DipoleGroup and 'label' points into the attribute table, so it must not outlive the system.
struct DipoleGroup {
    md_attribute_id_t key = MD_ATTRIBUTE_INVALID;   // id of the vector attribute, stable across a reload
    md_str_t  label  = { 0 };                          // group name as the path spells it
    uint32_t count = 0;                             // elements in the group; 1 means the index is decorative
    md_unit_t unit = md_unit_none();
};

// One density property out of the system's attribute table: a leaf under es_path::density_property,
// named by its own label or, failing that, by the last segment of its path.
// Build these with density_properties_gather; nothing owns or caches them, and 'label' points into
// the attribute table rather than at a copy - so it is good for as long as the system is loaded and
// must not outlive it.
struct DensityProperty {
    // How many a caller gathers into a buffer on its stack: well past what any file carries
    static constexpr size_t MaxCount = 512;

    uint64_t key    = 0;    // attribute id, stable across a reload
    md_str_t label     = { 0 };
};

// ---------------------------------------------------------------------------
// Volume geometry
// ---------------------------------------------------------------------------
// Fitting an object aligned box around a set of points and turning it into a grid and a volume.
// None of it knows what is being evaluated in that volume or who asked, which is why it sits here
// rather than in whichever component first needed it.

// Object aligned bounding box
struct OABB {
    md_float3x3_t orientation = md_float3x3_ident();
    md_float3_t min_ext = { 0 };
    md_float3_t max_ext = { 0 };
};

struct AABB {
    md_float3_t min_ext = { 0 };
    md_float3_t max_ext = { 0 };
};

md_float3x3_t mat3_PCA(const md_float4_t* xyzw, size_t count);

void calculate_bounds(float out_min[3], float out_max[3], const md_float4_t* xyzw, size_t count, const md_float3x3_t& orientation = md_float3x3_ident());

// Construct texture to world transformation matrix for Volume
// extent is the extent of the volume (dim * voxel_size)
md_float4x4_t compute_texture_to_world_mat(const md_float3x3_t& orientation, const md_float3_t& origin, const md_float3_t& extent);

md_float4x4_t compute_world_to_model_mat(const md_float3x3_t& orientation, const md_float3_t& origin);

md_float4x4_t compute_index_to_world_mat(const md_float3x3_t& orientation, const md_float3_t& in_origin, const md_float3_t& stepsize);

// Attempts to compute fitting volume dimensions given an input extent and a suggested number of samples per length unit
void compute_dim(int out_dim[3], const md_float3_t& in_ext, double samples_per_unit_length);

// Grid units are BOHR, matching what md_gto evaluates in; a Volume's transforms are Angstrom,
// matching the world the camera lives in. init_volume is where the two meet, and it is the only
// place that conversion belongs.
void init_grid(md_grid_t* grid, const md_float3x3_t& orientation, const md_float3_t& min_ext, const md_float3_t& max_ext, double samples_per_unit_length);

struct Volume {
    md_float4x4_t world_to_model   = {};   // Roto-translation into volume local axes, no scaling applied (preserves world length units)
    md_float4x4_t texture_to_world = {};   // Texture space [0,1] to world coordinates
    md_float3_t voxel_size  = {1,1,1};   // Size of each voxel in world units
    int dim[3] = {128, 128, 128};
    uint32_t tex_id = 0;
};

// Sizes and allocates the GL 3D texture for a grid, and fills in the transforms that place it in
// the world. Defined in viamd.cpp because it is the half that touches GL.
void init_volume(Volume* vol, const md_grid_t& grid, GLenum format = GL_R16F);

// How isosurfaces are drawn: a Field representation's (ElectronicStructureRepresentation::rendering), and the
// panels of the Orbital Grid and Transition Analysis windows (Settings::isosurface_rendering)
enum IsosurfaceRendering {
    IsosurfaceRendering_Raycast = 0,    // through the evaluated volume (volume::render_isosurfaces)
    IsosurfaceRendering_Mesh    = 1,    // meshes extracted from the field itself (iso_meshes.h, isomesh::render)
};

// The meshes of an electronic structure representation (iso_meshes.cpp)
struct EsIsoMesh;

// Descriptor for handling iso surfaces
struct IsoDesc {
    size_t count;
    float values[8];
    md_float4_t colors[8];
    float optical_densities[8];
};

// The image of a slice evaluated on its plane (the planes display; viamd.cpp, electronic_structure_plane_images_update): two, the one drawn and the one being evaluated, swapped when it is there. Drawn only while it is of the
// plane as it is (current): until then, and while the plane is dragged, the plane is drawn from the volume (and so is
// anything of it outside the image's rectangle).
struct ElectronicStructurePlaneImage {
    Volume    img[2] = {};      // W x H x 1 R32F, the texture to world of its rectangle (texture_to_world)
    uint64_t  key[2] = {};      // of what each is of (the source and the slice's rectangle): 0 none
    float     spacing[2] = {};  // each one's (Angstrom; it covers the rectangle exactly, its pixels near this square)
    int       front = 0;        // the one drawn
    uint32_t  serial = 0;       // counted up as the front changes: what the scene's hash takes
    bool      pending = false;  // the other is being evaluated
    bool      current = false;  // the front one is of the plane, the source and the frame as they are
    double    issued = 0.0;     // when the pending one was issued (s), and how long the front one took (ms)
    float     ms = 0.0f;
};

// A field over R^3 of the electronic structure: what it is (its path, field.h, names the source and a density
// property's attribute) and the parameters that tell the fields of one family apart. What a representation draws is one.
struct FieldRef {
    ElectronicStructureSource source = ElectronicStructureSource::MolecularOrbital;
    bool use_magnitude = false;
    ElectronicStructureSpin spin = ElectronicStructureSpin::Alpha;
    ElectronicStructureNtoComponent nto_component = ElectronicStructureNtoComponent::Particle;
    ElectronicStructureTransitionDensityComponent transition_density_component = ElectronicStructureTransitionDensityComponent::Attachment;

    int orbital_idx = 0;
    int excited_state_idx = 0;
    int nto_lambda_idx = 0;
    // The density property, named by the id of its attribute rather than by a position in a list. An index silently
    // re-points at a different property whenever the set changes across a reload; an id is derived from the path, so
    // it names the same thing or nothing at all.
    md_attribute_id_t density_property_key = MD_ATTRIBUTE_INVALID;
};

struct ElectronicStructureRepresentation {
    Volume density_vol = {};
    Volume color_vol   = {};
    // The analytic normals of density_vol (RGB10_A2, its dims), where its evaluation made them (es_normals_of): what the
    // raycast isosurfaces are shaded with
    uint32_t normal_tex = 0;
    EsIsoMesh* iso_mesh = nullptr;     // owned: the mesh rendering mode's surfaces, made on first use
    VolumeResolution resolution = VolumeResolution::Mid;

    // The isovalue of the orbitals and densities (a density's squared: electronic_structure_uses_density_iso_squared)
    double iso_value = 0.05;
    // and of the other fields (electronic_structure_source_field_kind), in their unit: both signs of it for one about
    // zero (electronic_structure_is_signed), else this one, which may be negative
    double field_iso_value = 0.8;

    // These are default values for different volume types.
    md_float4_t col_psi_pos = {0.f/255.f,75.f/255.f,135.f/255.f,0.75f};
    md_float4_t col_psi_neg = {255.f/255.f,205.f/255.f,0.f/255.f,0.75f};
    md_float4_t col_den     = {255.f/255.f,255.f/255.f,255.f/255.f,0.75f};
    md_float4_t col_att     = {0, 162.0f/255.0f, 135.0f/255.0f, 0.75f};
    md_float4_t col_det     = {162.0f/255.0f, 35.0f/255.0f, 135.0f/255.0f, 0.75f};

    md_float4_t tint_psi_pos = { 1.0f, 1.0f, 1.0f, 0.75f };
    md_float4_t tint_psi_neg = { 1.0f, 1.0f, 1.0f, 0.75f };
    md_float4_t tint_den     = { 1.0f, 1.0f, 1.0f, 0.75f };
    md_float4_t tint_att     = { 1.0f, 1.0f, 1.0f, 0.75f };
    md_float4_t tint_det     = { 1.0f, 1.0f, 1.0f, 0.75f };

    struct {
        int num_isos = 0;
        double values[8];
        md_float4_t colors[8];
    } density_property;

    // Scaling factor of *power* in the gaussians to splat the color volume (when using atom colors for volumes)
    double gaussian_splatting_power = 10.0;

    // Optical scaling factor which controls attenuation of light within iso surfaces.
    double iso_optical_density = 0.005;

    SurfaceColoring coloring = SurfaceColoring::Uniform;

    // SurfaceColoring::Field: the field, how it maps to colour, and its values on the grid of
    // density_vol - evaluated only around the isosurfaces this representation draws
    SurfaceFieldKind    field_kind = SurfaceFieldKind::EmbeddingPotential;
    ColorScale          field_map  = surface_field_default_scale(SurfaceFieldKind::EmbeddingPotential);
    ColorScaleIsolines  field_isolines = {};    // lines on the surfaces where the field takes given values
    SurfaceFieldVolume  field_vol  = {};

    // The grid density_vol was evaluated on (bohr), kept for what is evaluated on the same grid
    md_grid_t grid = {};

    // The field (FieldRef): the data the representation draws, in every display - its isosurfaces, its volume, its
    // planes. The isosurface display can colour its surfaces by another (coloring, field_kind).
    FieldRef shown;

    // How its isosurfaces are drawn (IsosurfaceRendering): meshes only of what they can be made of (es_mesh_mode),
    // the others raycast whatever this says
    int rendering = IsosurfaceRendering_Raycast;

    // The source evaluated as a surface field (electronic_structure_source_is_surface_field, and ELF and LOL without
    // their kernel): all of the grid, taken into density_vol once complete; its values' spread (1st to 99th
    // percentile, the extremes) for the isovalue's range and the colour map's
    SurfaceFieldVolume source_vol = {};
    ColorScaleSpan     source_span = {};

	uint64_t col_hash = 0;
	uint64_t vol_hash = 0;

    ElectronicStructureDisplay display = ElectronicStructureDisplay::Isosurface;

    // The planes of each display, apart, so that one is switched to and from without changing the other's
    // (electronic_structure_planes picks the display's). Volume: the cuts (what is drawn is the box of the volume on
    // the kept side of every enabled one). Planes: the slices, each drawn with its opacity, the direction of its
    // normal of no matter. Edited in the view with a gizmo (es_planes.h).
    ElectronicStructurePlaneSet cuts;
    ElectronicStructurePlaneSet slices;

    // Planes: each slice evaluated on itself (plane_images, at plane_detail's spacing); drawn from the volume, at Mid
    // (electronic_structure_volume_resolution), where its image is not there, while it is dragged, and without the GPU
    PlaneDetail   plane_detail = PlaneDetail::Fine;
    ElectronicStructurePlaneImage plane_images[ES_MAX_PLANES] = {};

    // The box the volume is (or would be) evaluated in: its unit cube to the world (Angstrom), as the volume's
    // texture_to_world; set whenever the representation is updated, with or without a volume (es_volume_box)
    md_float4x4_t box_to_world = {};
    bool   has_box = false;

    // Volume and Planes: what a value is drawn as, a colour map's colour over a range and an opacity curve over it
    // (transfer_function.h, edited with transfer_function_widget.h); the planes take the colours alone
    struct {
        TransferFunction tf = {};
        bool  map_init = false;        // tf was set up for a field (es_volume_map_update) ...
        bool  map_auto = false;        // ... and its range is to be taken from the values when they come (a potential's)
        int   map_kind = 0;            // ... of this kind (ElectronicStructureValueKind): set up again when that changes
        int   map_source = -1;         // ... and of this field where it is a range of its own (a surface field's source;
                                       // -2 for the others, -1 as loaded: taken to be the field it is of)
        uint32_t tf_tex = 0;           // tf baked into a texture, and the hash of what it was baked from
        uint64_t tf_hash = 0;
        TransferFunctionView view = {};// the editor's view of the values, and the histogram it shows (of the volume)
    } volume;
};

// The resolution the representation's volume is evaluated at: the planes display's is Mid (what a plane is drawn from
// until its own image is there, while it is dragged, and without the GPU), not a setting; else the representation's
static inline VolumeResolution electronic_structure_volume_resolution(const ElectronicStructureRepresentation& es) {
    return es.display == ElectronicStructureDisplay::Planes ? VolumeResolution::Mid : es.resolution;
}

// The planes of the representation's display: its slices in the planes display, its cuts otherwise
static inline ElectronicStructurePlaneSet& electronic_structure_planes(ElectronicStructureRepresentation& es) {
    return es.display == ElectronicStructureDisplay::Planes ? es.slices : es.cuts;
}
static inline const ElectronicStructurePlaneSet& electronic_structure_planes(const ElectronicStructureRepresentation& es) {
    return es.display == ElectronicStructureDisplay::Planes ? es.slices : es.cuts;
}

static inline ElectronicStructureSourceFlags electronic_structure_source_flag(ElectronicStructureSource source) {
    return 1u << (uint32_t)source;
}

static inline bool electronic_structure_source_supported(ElectronicStructureSourceFlags mask, ElectronicStructureSource source) {
    return (mask & electronic_structure_source_flag(source)) != 0;
}

static inline bool electronic_structure_is_density_property(const ElectronicStructureRepresentation& rep) {
    return rep.shown.source == ElectronicStructureSource::DensityProperty;
}

// What the values of a representation's field are, which its colour map is set up for (es_volume_map_update)
enum class ElectronicStructureValueKind {
    Density,        // >= 0, largest at the nuclei
    Signed,         // an orbital, a difference: about 0
    Localization,   // in [0, 1]: ELF, LOL
    Ranged,         // of a range of its own, set from the values, no value special: a potential that is not about
                    // its sign, a length, an energy density
    SignedRanged,   // likewise about 0, its sign the physics: an electrostatic potential
};

// ELF and LOL: in [0, 1] everywhere, of both spins together
static inline bool electronic_structure_is_localization(const ElectronicStructureRepresentation& rep) {
    return rep.shown.source == ElectronicStructureSource::ElectronLocalization ||
           rep.shown.source == ElectronicStructureSource::LocalizedOrbitalLocator;
}

static inline bool electronic_structure_uses_orbital_idx(const ElectronicStructureRepresentation& rep) {
    return rep.shown.source == ElectronicStructureSource::MolecularOrbital;
}

static inline bool electronic_structure_uses_excited_state_idx(const ElectronicStructureRepresentation& rep) {
    return rep.shown.source == ElectronicStructureSource::NaturalTransitionOrbital ||
           rep.shown.source == ElectronicStructureSource::TransitionDensity;
}

static inline bool electronic_structure_uses_nto_lambda_idx(const ElectronicStructureRepresentation& rep) {
    return rep.shown.source == ElectronicStructureSource::NaturalTransitionOrbital;
}

static inline bool electronic_structure_uses_magnitude_toggle(const ElectronicStructureRepresentation& rep) {
    return rep.shown.source == ElectronicStructureSource::MolecularOrbital ||
           rep.shown.source == ElectronicStructureSource::NaturalTransitionOrbital ||
           (rep.shown.source == ElectronicStructureSource::ElectronDensity && rep.shown.spin == ElectronicStructureSpin::Difference);
}

static inline bool electronic_structure_uses_spin(const ElectronicStructureRepresentation& rep) {
    return rep.shown.source == ElectronicStructureSource::MolecularOrbital ||
           rep.shown.source == ElectronicStructureSource::ElectronDensity;
}

static inline bool electronic_structure_is_signed(const ElectronicStructureRepresentation& rep) {
    if (electronic_structure_is_density_property(rep)) {
        return false;
    }
    const SurfaceFieldKind fk = electronic_structure_source_field_kind(rep.shown.source);
    if (fk != SurfaceFieldKind::Count) {
        return surface_field_kind_info(fk).scale == SurfaceFieldScale::Diverging;
    }
    if (rep.shown.source == ElectronicStructureSource::TransitionDensity &&
        rep.shown.transition_density_component == ElectronicStructureTransitionDensityComponent::Difference) {
        return true;
    }
    return electronic_structure_uses_magnitude_toggle(rep) && !rep.shown.use_magnitude;
}

static inline ElectronicStructureValueKind electronic_structure_value_kind(const ElectronicStructureRepresentation& rep) {
    if (electronic_structure_is_localization(rep)) return ElectronicStructureValueKind::Localization;
    if (electronic_structure_source_is_surface_field(rep.shown.source)) {
        return electronic_structure_is_signed(rep) ? ElectronicStructureValueKind::SignedRanged : ElectronicStructureValueKind::Ranged;
    }
    return electronic_structure_is_signed(rep) ? ElectronicStructureValueKind::Signed : ElectronicStructureValueKind::Density;
}

static inline bool electronic_structure_uses_density_iso_squared(const ElectronicStructureRepresentation& rep) {
    return rep.shown.source == ElectronicStructureSource::TransitionDensity ||
           rep.shown.source == ElectronicStructureSource::ElectronDensity;
}

static inline md_float4_t electronic_structure_density_property_default_color(const ElectronicStructureRepresentation& rep, int idx) {
    switch (idx) {
    case 0: return rep.col_psi_pos;
    case 1: return rep.col_psi_neg;
    case 2: return rep.col_den;
    case 3: return rep.col_att;
    case 4: return rep.col_det;
    default: return rep.col_den;
    }
}

static inline void electronic_structure_density_property_init_defaults(ElectronicStructureRepresentation* rep) {
    MD_ASSERT(rep);
    rep->density_property.num_isos = 2;
    rep->density_property.values[0] =  0.05;
    rep->density_property.values[1] = -0.05;
    rep->density_property.colors[0] = rep->col_psi_pos;
    rep->density_property.colors[1] = rep->col_psi_neg;
    for (int i = 2; i < (int)MD_COUNTOF(rep->density_property.values); ++i) {
        rep->density_property.values[i] = 0.05 * (double)(i + 1);
        rep->density_property.colors[i] = electronic_structure_density_property_default_color(*rep, i);
    }
}

static inline const char* electronic_structure_iso_value_label() {
    return "isovalue (*)";
}

static inline const char* electronic_structure_iso_value_tooltip(const ElectronicStructureRepresentation& rep) {
    if (electronic_structure_is_density_property(rep)) {
        return "User-defined isosurfaces for arbitrary density-property volumes.";
    }
    if (electronic_structure_source_field_kind(rep.shown.source) != SurfaceFieldKind::Count) {
        return electronic_structure_is_signed(rep) ?
            "Where the surfaces are: the field at +isovalue and -isovalue, in its unit" :
            "Where the surface is: the field at the isovalue, in its unit";
    }
    if (electronic_structure_uses_density_iso_squared(rep)) {
        return electronic_structure_is_signed(rep) ?
            "Visual surface threshold. Signed density fields render both +isovalue^2 and -isovalue^2 surfaces." :
            "Visual surface threshold. Density fields render at isovalue^2 for perceptual consistency with amplitude fields.";
    }
    return electronic_structure_is_signed(rep) ?
        "Visual surface threshold. Signed fields render both +isovalue and -isovalue surfaces." :
        "Visual surface threshold for the rendered field.";
}

static inline void electronic_structure_iso_desc_init(IsoDesc* iso, const ElectronicStructureRepresentation& rep) {
    MD_ASSERT(iso);
    *iso = {};    

    if (electronic_structure_is_density_property(rep)) {
        const int num_isos = MD_CLAMP(rep.density_property.num_isos, 0, (int)MD_COUNTOF(rep.density_property.values));
        iso->count = (size_t)num_isos;
        for (int i = 0; i < num_isos; ++i) {
            iso->values[i] = (float)rep.density_property.values[i];
            iso->colors[i] = rep.density_property.colors[i];
            iso->optical_densities[i] = (float)rep.iso_optical_density;
        }
        return;
    }

    const bool field = electronic_structure_source_field_kind(rep.shown.source) != SurfaceFieldKind::Count;
    const float iso_value = field ? (float)rep.field_iso_value : (float)rep.iso_value;
    const float iso_threshold = electronic_structure_uses_density_iso_squared(rep) ? iso_value * iso_value :
                                electronic_structure_is_signed(rep) ? fabsf(iso_value) : iso_value;

    if (electronic_structure_is_signed(rep)) {
        iso->count = 2;
        iso->values[0] =  iso_threshold;
        iso->values[1] = -iso_threshold;
        if (rep.coloring != SurfaceColoring::Uniform) {
            iso->colors[0] = rep.tint_psi_pos;
            iso->colors[1] = rep.tint_psi_neg;
        } else {
            iso->colors[0] = rep.col_psi_pos;
            iso->colors[1] = rep.col_psi_neg;
        }
        iso->optical_densities[0] = (float)rep.iso_optical_density;
        iso->optical_densities[1] = (float)rep.iso_optical_density;
    } else {
        iso->count = 1;
        iso->values[0] = iso_threshold;
        if (rep.coloring != SurfaceColoring::Uniform) {
            if (rep.shown.source == ElectronicStructureSource::TransitionDensity && rep.shown.transition_density_component == ElectronicStructureTransitionDensityComponent::Attachment) {
                iso->colors[0] = rep.tint_att;
            } else if (rep.shown.source == ElectronicStructureSource::TransitionDensity && rep.shown.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment) {
                iso->colors[0] = rep.tint_det;
            } else {
                iso->colors[0] = rep.tint_den;
            }
        } else {
            if (rep.shown.source == ElectronicStructureSource::TransitionDensity && rep.shown.transition_density_component == ElectronicStructureTransitionDensityComponent::Attachment) {
                iso->colors[0] = rep.col_att;
            } else if (rep.shown.source == ElectronicStructureSource::TransitionDensity && rep.shown.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment) {
                iso->colors[0] = rep.col_det;
            } else {
                iso->colors[0] = rep.col_den;
            }
        }
        iso->optical_densities[0] = (float)rep.iso_optical_density;
    }
}

static inline ElectronicStructureField electronic_structure_legacy_field(const ElectronicStructureRepresentation& rep) {
    return rep.shown.use_magnitude ? ElectronicStructureField::Density : ElectronicStructureField::Amplitude;
}

static inline void electronic_structure_set_legacy_field(ElectronicStructureRepresentation* rep, int field) {
    MD_ASSERT(rep);
    rep->shown.use_magnitude = (rep->shown.source == ElectronicStructureSource::MolecularOrbital ||
                          rep->shown.source == ElectronicStructureSource::NaturalTransitionOrbital) &&
                         (ElectronicStructureField)field == ElectronicStructureField::Density;
}

// What a field of its source starts out with: the spin it is of, and its magnitude off
static inline void field_ref_set_defaults(FieldRef* f) {
    MD_ASSERT(f);
    f->use_magnitude = false;
    switch (f->source) {
    case ElectronicStructureSource::MolecularOrbital:
        if (f->spin != ElectronicStructureSpin::Alpha && f->spin != ElectronicStructureSpin::Beta) {
            f->spin = ElectronicStructureSpin::Alpha;
        }
        break;
    case ElectronicStructureSource::ElectronDensity:
        f->spin = ElectronicStructureSpin::Total;
        break;
    default:
        // of no spin, or of both together
        f->spin = ElectronicStructureSpin::None;
        break;
    }
}

// What the representation starts out with for the field it shows: the field's defaults, and its isovalue and surface
static inline void electronic_structure_set_source_defaults(ElectronicStructureRepresentation* rep) {
    MD_ASSERT(rep);
    field_ref_set_defaults(&rep->shown);
    const SurfaceFieldKind fk = electronic_structure_source_field_kind(rep->shown.source);
    if (fk != SurfaceFieldKind::Count) {
        // a field of the system as a whole: an isovalue in its unit
        rep->field_iso_value = surface_field_kind_info(fk).iso;
        rep->source_span = {};
    }
    if (rep->shown.source == ElectronicStructureSource::DensityProperty) {
        // The atoms' colours are not offered for a density property; a field is, as it is evaluated on the same grid
        if (rep->coloring == SurfaceColoring::AtomColors) {
            rep->coloring = SurfaceColoring::Uniform;
        }
        if (rep->density_property.num_isos <= 0) {
            electronic_structure_density_property_init_defaults(rep);
        }
    }
}

static inline void electronic_structure_set_legacy_type(ElectronicStructureRepresentation* rep, int type) {
    MD_ASSERT(rep);
    switch ((ElectronicStructureLegacyType)type) {
    case ElectronicStructureLegacyType::MolecularOrbital:
        rep->shown.source = ElectronicStructureSource::MolecularOrbital;
        rep->shown.use_magnitude = false;
        break;
    case ElectronicStructureLegacyType::MolecularOrbitalDensity:
        rep->shown.source = ElectronicStructureSource::MolecularOrbital;
        rep->shown.use_magnitude = true;
        break;
    case ElectronicStructureLegacyType::NaturalTransitionOrbitalParticle:
        rep->shown.source = ElectronicStructureSource::NaturalTransitionOrbital;
        rep->shown.use_magnitude = false;
        rep->shown.nto_component = ElectronicStructureNtoComponent::Particle;
        break;
    case ElectronicStructureLegacyType::NaturalTransitionOrbitalHole:
        rep->shown.source = ElectronicStructureSource::NaturalTransitionOrbital;
        rep->shown.use_magnitude = false;
        rep->shown.nto_component = ElectronicStructureNtoComponent::Hole;
        break;
    case ElectronicStructureLegacyType::NaturalTransitionOrbitalDensityParticle:
        rep->shown.source = ElectronicStructureSource::NaturalTransitionOrbital;
        rep->shown.use_magnitude = true;
        rep->shown.nto_component = ElectronicStructureNtoComponent::Particle;
        break;
    case ElectronicStructureLegacyType::NaturalTransitionOrbitalDensityHole:
        rep->shown.source = ElectronicStructureSource::NaturalTransitionOrbital;
        rep->shown.use_magnitude = true;
        rep->shown.nto_component = ElectronicStructureNtoComponent::Hole;
        break;
    case ElectronicStructureLegacyType::AttachmentDensity:
        rep->shown.source = ElectronicStructureSource::TransitionDensity;
        rep->shown.use_magnitude = false;
        rep->shown.transition_density_component = ElectronicStructureTransitionDensityComponent::Attachment;
        break;
    case ElectronicStructureLegacyType::DetachmentDensity:
        rep->shown.source = ElectronicStructureSource::TransitionDensity;
        rep->shown.use_magnitude = false;
        rep->shown.transition_density_component = ElectronicStructureTransitionDensityComponent::Detachment;
        break;
    case ElectronicStructureLegacyType::ElectronDensity:
        rep->shown.source = ElectronicStructureSource::ElectronDensity;
        rep->shown.use_magnitude = false;
        rep->shown.spin = ElectronicStructureSpin::Total;
        break;
    default:
        MD_ASSERT(false);
        break;
    }
}


static inline int electronic_structure_legacy_type(const ElectronicStructureRepresentation& rep) {
    switch (rep.shown.source) {
    case ElectronicStructureSource::MolecularOrbital:
        return (int)(!rep.shown.use_magnitude ?
            ElectronicStructureLegacyType::MolecularOrbital :
            ElectronicStructureLegacyType::MolecularOrbitalDensity);
    case ElectronicStructureSource::NaturalTransitionOrbital:
        if (!rep.shown.use_magnitude) {
            return (int)(rep.shown.nto_component == ElectronicStructureNtoComponent::Particle ?
                ElectronicStructureLegacyType::NaturalTransitionOrbitalParticle :
                ElectronicStructureLegacyType::NaturalTransitionOrbitalHole);
        }
        return (int)(rep.shown.nto_component == ElectronicStructureNtoComponent::Particle ?
            ElectronicStructureLegacyType::NaturalTransitionOrbitalDensityParticle :
            ElectronicStructureLegacyType::NaturalTransitionOrbitalDensityHole);
    case ElectronicStructureSource::TransitionDensity:
        return (int)(rep.shown.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment ?
            ElectronicStructureLegacyType::DetachmentDensity :
            ElectronicStructureLegacyType::AttachmentDensity);
    case ElectronicStructureSource::ElectronDensity:
        return (int)ElectronicStructureLegacyType::ElectronDensity;
    default:
        // nothing of the legacy types: the electron density, as the key has always said for what it did not know
        return (int)ElectronicStructureLegacyType::ElectronDensity;
    }
}

// Colouring atoms by one of their attributes: a per atom scalar field of the system's attribute
// table (atom/...). The attribute itself is stored as its id and nothing else about it - no copied
// label, no cached list position. The id is a hash of the path and stable across a reload, so a
// representation keeps pointing at the same quantity when the data is reloaded; an index into a
// gathered list would not.
//
// The values map to colours by a colour scale, the same one an isosurface coloured by a field uses
// (color_scale.h). Its automatic range follows the span of the atoms the representation SHOWS - the
// atoms of its filter, over every variant so the colours do not jump as the variant changes.
struct AtomAttributeColoring {
    md_attribute_id_t key = MD_ATTRIBUTE_INVALID;
    int variant_idx = 0;                // position along the leading axis, when the field has one
    ColorScale scale = {};
    ColorScaleSpan span = {};           // of the shown atoms, refreshed with the colours
    uint64_t span_hash = 0;             // what span was measured over: the attribute, its version, the atoms
};

struct DipoleRepresentation {
    md_attribute_id_t dipole_key = MD_ATTRIBUTE_INVALID;   // the vector attribute, so it survives a reload
    uint32_t dipole_index = 0;                             // which element of that group; see DipoleGroup
    md_float4_t color = { 0, 0, 0, 1 };
    md_float3_t offset = { 0, 0, 0 };
    double scale = 1.0;
	float radius = 0.05f;
};

// The edges of the unit cell, of the frame shown (it may change with it)
struct BoxRepresentation {
    md_float4_t color = { 0.25f, 0.25f, 0.25f, 1.0f };
    float thickness = 1.5f;     // pt, the width of the edges on screen (tubes, lit), never thinner than a pixel; 0 for
                                // lines of a pixel, unlit
    bool axes = false;          // The edges along a, b and c from the origin of the cell red, green and blue
};

// The representations, the groups they can be put in and everything done with them is in representation.h.
struct Representation {
    RepresentationNode node = {};   // where it is in the Representations window, see representation_tree.h
    // The id of the instance of its system it is of (SystemInstance::id), as what is around it in the tree is: drawn
    // from its coordinates, in its place. Of the first where it names none there is.
    uint32_t instance = 0;

    char name[64] = "rep";
    // Evaluated within the context of the groups it is in (see RepresentationGroup). Its result says whether what it
    // picks depends on the frame, through that context too. Empty, as a new one is, it picks everything, as 'all' does.
    AtomFilter filter = {};

    RepresentationType type = RepresentationType::BallAndStick;
    ColorMapping color_mapping = ColorMapping::Type;
    // What it shows: what its filter picked when it was last valid, kept while the filter or its context is not
    // (representation_filter_invalid); empty when it has not been valid since the system was loaded
    md_bitset_t atom_mask = {};
    md_gl_rep_t md_rep = {};
    // The representation of the cartoon of the lipids and nucleic acids (CartoonGlyphs), drawn by a Cartoon, and whether
    // it shows any of them (a cartoon of proteins alone draws none)
    md_gl_rep_t glyph_rep = {};
    bool glyphs_shown = false;

    bool enabled = true;            // its own switch: it is shown when every group around it is too (representation_shown)
	bool needs_update = true;
    bool type_is_valid = false;
    // Drawn as spheres for now: what its type draws (bonds, a backbone) is not there yet, while the system is inferred in
    // the background (SystemInference)
    bool sphere_fallback = false;
    bool filt_is_dirty = true;
    bool dynamic_evaluation = true; // follow the frame when what the filter picks depends on it ("auto-update")

    // User defined base color used in uniform mode
    md_float4_t base_color = {1.0f, 1.0f, 1.0f, 1.0f};

    // Global post processing parameters applied to the final colors.
    md_float4_t tint_color = {1.0f, 0.0f, 0.0f, 1.0f};
    float tint_scale = 0.0f;
    float saturation = 1.0f;

    // scaling parameter (radius, width, height, etc depending on type). Cartoon: coil, sheet, helix, and the thickness of
    // the lipids and nucleic acids (w)
    md_float4_t scale = {1.0f, 1.0f, 1.0f, 1.0f};

	BondColorMode bond_color = BondColorMode::NearestAtom;
    float  bond_sharpness = 0.5f; // 0 = sharper, 1 = smoother
	md_float4_t bond_base_color = { 1.0f, 1.0f, 1.0f, 1.0f };
    bool   bond_orders = true;    // Licorice, ball and stick: double and triple bonds, aromatic and delocalized bonds (md_gl_mol_set_bond_orders)
    bool   chain_gaps = true;     // Cartoon, ribbons, licorice, ball and stick: a dashed line over each break of a chain (md_gl_mol_set_chain_gaps)
    bool   elastic_bonds = false; // Licorice, ball and stick: the elastic bonds (MD_BOND_KIND_ELASTIC), springs of the force field such as an elastic network, which are not bonds of the molecules
    bool   interaction_bonds = true; // Licorice, ball and stick: the interactions of the system among its atoms, dashed and thinner than its sticks. Off where a Bonds representation is to draw them instead
    bool   coordinate_bonds_dashed = true; // Licorice, ball and stick: the coordinate bonds (MD_BOND_KIND_COORDINATE) dashed, rather than sticks
    // Bonds: what the filter picks, as of its last evaluation (in the persistent allocator), and which are dashed. In
    // VIAL also the flags (md_bond_flags_t) of each pair's bond, the system's or one found for the frame (a hydrogen, ionic
    // or halogen bond bond() finds, a bond formed()): empty in v1, whose pairs are told by the system's bonds
    md_array(md_atom_pair_t) bond_pairs = nullptr;
    md_array(uint32_t) bond_flags = nullptr;
    BondDashes bond_dashes = BondDashes::ByKind;
    bool   open_properties = false; // The window opens its properties the next time it lists it (one made from the context menu, to be styled)
    double invalid_since = 0.0;     // ImGui::GetTime() when its filter was found not valid, for the notice in the viewport; 0 while valid
    bool   two_head_spheres = false;  // Cartoon: the head of a phospholipid as its head group and its backbone (CartoonGlyphs)
    // Cartoon of nucleic acids, Å (before the scale of the lipids and nucleic acids): the radius of the band of the
    // backbone, half its thickness, and that of the rims of the rings of the sugars and bases and of the sticks between them
    float  nucleic_band_radius = CARTOON_NUCLEIC_BAND_RADIUS;
    float  nucleic_ring_radius = CARTOON_NUCLEIC_RING_RADIUS;

    ElectronicStructureRepresentation electronic_structure = {};
    AtomAttributeColoring atom_attribute = {};
	DipoleRepresentation dipole = {};
    BoxRepresentation box = {};
    GlyphRepresentation glyphs = {};
    bool periodic_images = true;    // Repeated in the periodic images of the view (system.periodic_images), as they are shown
};

// A group of representations. What it holds is shown only while the group is, and its filter, when it has one, is the
// context of everything within it: their filters are evaluated within the atoms it picks, as the left hand side of 'in'
// is within its right in a script - residue(1) is the first residue of the group, all is the group. A group within a
// group has its filter evaluated within the outer one's in turn. Without a filter a group only groups, and what is in
// it sees whatever context the group itself is in.
struct RepresentationGroup {
    RepresentationNode node = {0, 0, 0, true};
    uint32_t instance = 0;          // As a representation's (Representation::instance), and so is all within it

    char name[64] = "group";
    // Empty for no context of its own. Its result says whether the context is valid and whether it depends on the
    // frame, the context it is evaluated in included.
    AtomFilter filter = {};

    bool enabled = true;
    bool filt_is_dirty = true;
    md_bitset_t atom_mask = {};     // what the filter picks within its context, empty without a filter
};

struct FrameCache {
    FRAME_CACHE_LRU_TYPE lru = {};
    md_system_state_t states[FRAME_CACHE_SIZE] = {};
    int32_t frame_idx[FRAME_CACHE_SIZE] = {};
};

// A block of picking indices belonging to one domain. 'key' names WHAT the block indexes, for a
// domain whose objects are not the system's own arrays: a dipole range is keyed by the attribute id
// of one dipole group, so an index within it is the element index inside that group and the hit
// needs nothing else to resolve. Domains that index the system directly - atoms, bonds - leave it 0.
struct PickingRange {
    PickingDomainID domain = 0;
    uint32_t beg = 0;
    uint32_t end = 0;
    uint64_t key = 0;
    uint32_t system = SYSTEM_ID_NONE;   // The id of the system whose objects it indexes, none for what is of no system
};

// Sized for the domains of every system (atoms, bonds, ...: six each, for at most SYSTEM_MAX) and one range per keyed
// object. A VeloxChem file publishes four dipole groups today (ground state plus the electric, magnetic and velocity
// transition sets), and each representation which draws weak bonds has a range of its own; reserving past the end fails
// rather than aliasing, which costs those objects their picking, nothing else.
struct PickingSpace {
    size_t num_ranges = 0;
    PickingRange ranges[128] = {};
};

struct PickingReadbackSlot {
    uint32_t color_pbo = 0;
    uint32_t depth_pbo = 0;
    
    uint32_t submitted_frame_idx = 0;
    bool pending = false;
    
    uint32_t viewport_width = 0;
    uint32_t viewport_height = 0;
    
    md_float2_t surface_coord = {0};
    md_float2_t screen_coord  = {0};
    md_float4x4_t clip_to_world = MD_MAT4_IDENT_INIT;
};

struct PickingSurface {
    PickingSourceID source = 0;
    uint32_t slot_cursor = 0;
    PickingReadbackSlot slots[2] = {};
};

struct PickingHandler {
    uint32_t frame_idx = 0;

    struct {
        PickingSpace space = {};
        uint32_t submitted_frame_idx = 0;
    } history[2];
};

struct PickingHit {
    PickingSourceID source = 0;
    PickingDomainID domain = 0;

    // What was hit, in the domain's own terms: 'key' is the range's key (0 for a domain that does
    // not use one) and 'local_idx' is the offset within that range. For a dipole that is exactly
    // (group attribute, element) - the identity, carried through, with nothing to reconstruct.
    uint64_t key = 0;

    // The id of the system it is of (picking_hit_system), SYSTEM_ID_NONE for what is of no system or no hit
    uint32_t system = SYSTEM_ID_NONE;

    uint32_t frame_idx = 0;    uint32_t raw_idx = INVALID_PICKING_IDX;
    uint32_t local_idx = 0;

    md_float2_t surface_coord = {0};
    md_float2_t screen_coord = {0};
    md_float3_t world_pos = {0};
    float depth = 1.0f;
};

struct PickingReadbackRequest {
    uint32_t fbo = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    md_float2_t surface_coord = {0};
    md_float2_t screen_coord = {0};
    md_float4x4_t clip_to_world = MD_MAT4_IDENT_INIT;
};

// The key of the script evaluation's slot of a distribution or a volume over a timeline period, by its id
// (script::aggregate_request): asked for by what is shown over the period (series_resolve_histogram,
// series_resolve_volume)
#define SCRIPT_SLOT_PERIOD(id) (((uint64_t)1 << 32) | (uint64_t)(id))

// An evaluation of the script over a run of the first system other than the one script.eval is of
struct ScriptRunEval {
    char          run[128] = "";    // "run/<name>"
    script::Eval* eval = nullptr;
    bool          queued = false;   // to be evaluated, once what runs now has ended
};

// What the script visualization is to show this frame (script_visualize_ref, script_visualize_str)
struct ScriptVisTarget {
    enum Kind : uint8_t { None = 0, Ref, Str };
    Kind kind = None;
    script::Ref ref = {};               // Ref
    md_str_t str = {};                     // Str: an expression, copied into the frame allocator
    int subidx = -1;
    uint32_t flags = 0;                 // script::DRAW_*
};

// A visualization computed off the main thread (viamd.cpp)
struct ScriptVisJob;

// Edge length of the device evaluation scratch (ApplicationState::gpu_volume), in voxels: the largest
// grid any evaluation on the device can write, which is the largest compute_dim makes.
inline constexpr uint32_t GPU_VOLUME_DIM = 512;

// Edge length of the device scratch the images of planes are evaluated into (ApplicationState::gpu_plane): a
// GPU_PLANE_DIM^2 x 1 R32F texture, the largest plane image (es_plane_images.h)
inline constexpr uint32_t GPU_PLANE_DIM = 2048;

// The largest grid the analytic normals of an evaluation are made for (ApplicationState::gpu_grad): 20 bytes a voxel
// on the device while it is evaluated (the gradient in float, its normals packed), 4 in the GL texture. Past it the
// isosurfaces take their normals from the volume's differences.
inline constexpr size_t GPU_GRAD_MAX_VOXELS = size_t(1) << 25;

// What is derived from the core a loader gives and the user edits (the bonds, what the residues are, the structures,
// the chemistry: md_system_infer), inferred in the background and put in place of what the system had between frames.
// Until then the system is shown as it was: as loaded, its atoms, as spheres where a representation needs what is not
// there yet; after an edit of its bonds, as inferred before. One system at a time: a load cancels the inference of the
// system before it.
struct SystemInference {
    task_system::Job    job;
    md_system_infer_work_t work = {};             // Inferred on: the core of the system shared, the derived in an arena of its own
    md_allocator_i*     arena = nullptr;        // The rest of the work: its copies of the types of the atoms
    std::atomic<bool>   done{false};            // The work is inferred: not cancelled half way
    bool                pending = false;        // The system shows its core: inference is under way, or to be started again
    bool                restart = false;        // Cancelled (the pool was interrupted) or of a topology since changed
    // The default representations wait for it: a placeholder shows the atoms meanwhile, replaced unless it was changed
    bool                defaults_pending = false;
    char                placeholder_name[64] = "";
};

// THE DATA OF A SYSTEM: the system itself (md_system_t), its state as shown, what is derived from it to draw it, its
// run and the frames read of it. ApplicationState::mold of the one system loaded (SystemContext::data).
struct SystemData {
    // The arena backing this dataset is sys.alloc - there is no second handle to it. It is
    // created once at startup and survives every load: free_system_data rewinds it and puts the
    // handle back, because zeroing a system also zeroes the allocator it was using.
    // The cartoon's secondary structure (blend weights between coil, helix and sheet) lives on the GPU only: it is
    // computed where it is produced and uploaded straight away (upload_secondary_structure_weights,
    // interpolate_system_state), never read back and kept nowhere else.
    md_system_t         sys = {};
    SystemInference     inference;

    // The version of the system's groups last seen (md_group_data_t::version): a change of them, a file of groups loaded
    // or the user's edited, is broadcast once (EventType_ViamdGroupsChanged), in the frame loop
    uint64_t            group_version = 0;

    // The system's own coordinates: what its files gave (a structure, a restart; the first frame of its first run when
    // they gave none), which it is inferred of and its instances without a run show. What is drawn is each instance's
    // (SystemInstance::state), made from these or from a run's frames.
    md_system_state_t   state = {};
    // The excited states whose natural transition orbitals have been read on a worker thread (viamd.cpp,
    // NtoWarm; es_nto_ready). Of the system: freed with it, once its tasks have ended.
    void*               nto_warm = nullptr;
    // The SCF densities the CPU paths read (the ESP, the critical point search), each kept as its packed
    // triangle in double once read (density_cache.h; density_cache_of). Its entries are of the system and go
    // with it; the cache itself is the application's.
    struct DensityCache* density_cache = nullptr;

#if MD_ENABLE_GPU
    // GPU side data derived from sys, sitting beside gl_mol and for the same reason: it is
    // built FROM the system so that something can draw or evaluate it, and mdlib never reads it
    // back. The test is simple - no mdlib function takes an md_system_t* and touches these -
    // which is what keeps them out of md_system_t itself, unlike the trajectory it does read.
    //
    // This whole block is per DATASET, so it is what gets replicated when several systems can
    // be loaded at once. That is also why the evaluation scratch is not here but on the device.
    md_gto_gpu_basis_t  gpu_basis = nullptr;   // built from the basis/ attributes, not from a loader
    md_gpu_addr_t       gpu_atoms = 0;         // packed float4 positions, xyz in Bohr
    // Hash of the positions currently in gpu_atoms, 0 when nothing has been uploaded. The
    // positions come from the system STATE, so they move with the trajectory; comparing what is
    // uploaded against what is wanted is the only test that cannot go stale, and it costs a hash
    // over the QM atoms against an upload of the same array.
    uint64_t            gpu_atoms_hash = 0;
    // The SCF density the evaluations read, as a matrix in tiles (md_gto_gpu_dmat): only the AO pairs
    // an evaluation reads, built on the device from the orbitals, so no matrix is made on the host.
    // Its pattern follows the atoms while they stay within its margin; the matrix is rebuilt when the
    // density asked for, or the orbitals it is made of, change.
    md_gto_gpu_dmat_t   gpu_dmat = nullptr;
    uint64_t            gpu_dmat_key = 0;      // the density and the orbitals' versions it holds, 0: none
    // The next one, made on a worker thread for the representations (viamd.cpp, DmatBuild), and the key of
    // one that could not be (made where it is asked for instead)
    void*               gpu_dmat_build = nullptr;
    uint64_t            gpu_dmat_failed_key = 0;
    // What the ELF and LOL kernel reads (viamd.cpp, LocalizationGpu): a GPU basis at their AO cutoff and the
    // spins' density matrices in tiles, made by a pool task; and the key of one that could not be made
    void*               gpu_localization = nullptr;
    uint64_t            gpu_localization_failed_key = 0;
    // Any other density a representation evaluates (a transition density, made by its provider on every
    // read; a density property, read from its file) as the packed triangle the kernel reads, extracted on a
    // worker thread straight into a device buffer of its own (viamd.cpp, DensityPackBuild) and kept: the
    // last one, until another is asked for
    md_gpu_addr_t       gpu_dpack = 0;
    uint64_t            gpu_dpack_key = 0;     // the density, slice and version it holds, 0: none
    void*               gpu_dpack_build = nullptr;
    uint64_t            gpu_dpack_failed_key = 0;
#endif
};

// The most representations a system has. A hard limit, kept far above what is used, so that what a renderer keeps per
// atom and representation fits a fixed word: which representations an atom is in, as one bit each (uint32_t). Groups
// are not counted: they draw nothing.
constexpr size_t REPRESENTATION_MAX = 32;

// THE REPRESENTATIONS OF A SYSTEM, and the groups they are in (representation.h). ApplicationState::representation of
// the one system loaded (SystemContext::representation). At most REPRESENTATION_MAX of them.
struct SystemRepresentations {
    // Both in the order the window lists them. Where each one is in the tree of groups is on the
    // representation or group itself (representation_tree.h).
    md_array(Representation) reps = 0;
    md_array(RepresentationGroup) groups = 0;
    uint32_t next_rep_id = 1;
    uint32_t next_group_id = 1;
    md_bitmask_t visibility_mask = {};     // the union of the masks of the shown representations
    uint64_t visibility_mask_hash = 0;
    bool atom_visibility_mask_dirty = false;
    uint64_t gl_version = 0;    // bumped whenever what a representation draws is uploaded, for views that cache their image
    // The plane the gizmo in the view edits: of the representation of this id (its node id), at this index
    struct {
        uint32_t rep_id = 0;
        int      plane = -1;
    } gizmo;
};

// THE OPERATIONS ON THE COORDINATES OF A SYSTEM (SystemOperation), each turned on or off by its own flag, which keeps
// its settings while it is off. Of the system: every instance of it has its coordinates made the same way
// (SystemInstance), which is what superposes them.
struct SystemOperations {
    bool recenter = false;
    bool fixate_orientation = false;

    // Extra transform applied to the target in its own centred frame, between fixating the
    // orientation and placing it at the centre of the cell:
    //     T = translate(cell_centre) * alignment_mat * R * translate(-com)
    // Identity by default. This is the hook for aligning the structure to something once it is
    // centred - principal axes, a chosen bond, a user supplied orientation - without having to
    // touch the transform itself. Owned by whoever sets it: recentering never writes to it, and
    // it survives a change of recenter target.
    md_float4x4_t alignment_mat = MD_MAT4_IDENT_INIT;


    struct {
        AtomFilter filter;    // whether it depends on the frame (and is evaluated again in every one) is its result's
        bool enabled = false;
        uint64_t version = 1;
        uint64_t evaluated_version = 0;
        uint64_t ir_fingerprint = 0;
        md_bitset_t mask = { 0 };
    } recenter_query;

    // External state which user controls
    bool apply_pbc = false;
    bool unwrap_structures = false;

    // The bonds formed and broken in the frame shown, against the bonds of the system (md_bond_delta_compute), when
    // tracking them. The system is not changed by it: the bonds formed are drawn as sticks and picked, those broken
    // are not drawn, and Recompute Bonds is what changes the topology. Of the atoms the scope picks, all without one.
    struct {
        bool enabled = false;
        AtomFilter scope = {};                      // Empty for every atom
        float tolerance = MD_BOND_DELTA_TOLERANCE;
        uint32_t breakable = MD_BOND_DELTA_BREAKABLE;

        md_bitset_t scope_mask = {0};               // What the scope picks, in the persistent allocator
        bool scope_dirty = true;                    // To be evaluated again
        uint64_t scope_ir_fingerprint = 0;          // Of the script whose names it may use, when it was evaluated

        md_bond_delta_t delta = {};                 // Of the frame shown, in the persistent allocator
        // The most bonds formed in a frame since the delta was last empty. The picking indices of the bonds formed
        // are reserved for this many when a frame begins, before its delta is known (picked_bond): a frame which
        // forms more than any before it has the indices of the ones beyond run into those of the next domain.
        size_t formed_capacity = 0;
        // The chemistry of the frame where its bonds changed (md_chem_perceive_delta): the orders, aromaticity and
        // rings drawn there, and what the tooltips say of its atoms and bonds. Empty when it is not perceived (a
        // region too large). Its rings pick after those of the system, for which ring_capacity are reserved as
        // formed_capacity are for the bonds.
        md_chem_delta_t chem = {};
        size_t ring_capacity = 0;
        uint64_t key = 0;                           // What it was computed from, 0 when it is empty
        // The topology of the system it was found against: the bonds its broken ones are indices of, which are of
        // the system as inferred then. Another inference (an edit of the bonds, Recompute Bonds) puts other bonds
        // in their place, and the delta is not of them until it is computed again (bond_tracking_update).
        uint64_t topology_version = 0;
        double time_ms = 0.0;                       // How long that took
    } bond_tracking;

    // Manual selection mask for recentering / orientating, which atoms to consider for calculating the center of mass and principal axes
    // This is populated by assigning a user defined selection.
    md_bitset_t selection_mask = {0};

    struct {
        // The atoms the reference was built from, compared exactly against the active target. A query
        // which depends on the frame is evaluated again in every frame but mostly picks the same atoms,
        // and the reference only has to be rebuilt (from the first frame, read from disk) when they change.
        md_bitset_t target_mask = {0};    // in the persistent allocator
        bool valid = false;

		// Need to store the initial relative coordinate vectors and center of mass of reference frame for recentering and orientating the structure.
        md_float4_t* rel_xyzw = nullptr;
        md_float3_t  com = {};
    } initial_frame;
};

// AN INSTANCE OF A SYSTEM, as the scene shows it: where it is placed, the copies of it drawn, and whether it is drawn
// at all (SystemInstance::scene). The figure's: in the workspace ([Instance]).
struct SystemScene {
    bool visible = true;    // Its representations are drawn only while it is (and its system is: SystemContext::visible)

    // The periodic images shown: the system again in the neighboring cells, this many each way along a, b and c, of
    // the representations which have periodic_images (all, unless they opt out)
    int periodic_images[3] = {0, 0, 0};

    // The biological assembly shown (md_system_t assembly, of a deposited structure), -1 for the coordinates as
    // they are. Shown in place of them and of the periodic images: each group of it the instances of the group, in the
    // copies its operators make, drawn by instancing. The atoms of a copy are those of the system, which they pick
    // and select as. The operators place the coordinates as the file gives them, so the operations which move them
    // (centering, wrapping, making whole) are not applied while one is shown.
    int assembly = -1;

    // Where the instance is placed in the world: a rigid transform of its coordinates as they are drawn, identity for
    // as its source has them (overlaid on the others). Of the view only: the coordinates themselves (SystemInstance::
    // state), which the script measures and the operations write, are not moved by it - unlike centering. Side by side,
    // and superposed on another by a fit of a selection of each, are placements. Applied wherever the instance's
    // coordinates meet the world: its draws (as a model matrix, under its copies), the box and dipoles drawn of it,
    // region selection, framing; its electronic structure (of the quantum chemistry system alone: qm_placement).
    md_float4x4_t placement = MD_MAT4_IDENT_INIT;

    // DERIVED, NOT SAVED: the copies drawn of what is above, made from it every frame (update_view_param)

    // The placement as drawn in this frame and in the one before, so that a system being moved shows as moving
    md_float4x4_t placement_curr = MD_MAT4_IDENT_INIT;
    md_float4x4_t placement_prev = MD_MAT4_IDENT_INIT;

    // The copies of the system the representations with periodic_images are drawn in, by instancing (md_gl_copies_t):
    // its own first, then its periodic images (periodic_images), none without a cell or images, or with an assembly
    // shown. As of this frame and as of the last, so that an image which moves with the cell (a trajectory of changing
    // cells) shows as moving.
    md_array(md_float4x4_t) periodic_transforms      = 0;
    md_array(md_float4x4_t) periodic_transforms_prev = 0;

    // Of the assembly shown: the transforms of its copies, group by group, those of group g
    // [assembly_copy_offset[g], assembly_copy_offset[g + 1]). Empty when none is shown.
    md_array(md_float4x4_t)   assembly_transforms  = 0;
    md_array(uint32_t) assembly_copy_offset = 0;

    // Of the assembly shown, the instances of each group as a subset of the molecule (md_gl_subset_t), which its
    // copies draw, and the subset of the cartoon glyphs which that makes. A group of all instances has none (0): the
    // whole molecule. Made when another assembly is shown, or the system is loaded again (assembly_subsets_update).
    md_array(md_gl_subset_t) assembly_subsets       = 0;
    md_array(md_gl_subset_t) assembly_glyph_subsets = 0;
    int      assembly_subsets_of  = -1;     // The assembly they are of, of the molecule
    uint32_t assembly_subsets_mol = 0;      // (md_gl_mol_t id)


};


// THE DATA OVER THE FRAMES OF A RUN.
// Views onto the temporal attributes, plus the one array still owned here.
//
// 'stride' and 'count' are redundant with the attributes' own shape {frames, segments} and are
// kept only so the existing consumers keep their loops. Whoever next touches those loops should
// take the shape from the attribute and delete them - the point of the move was that a hand
// rolled shape beside table owned storage is a second answer waiting to disagree.
//
// The fingerprints are gone for the same reason: md_attributes_version is the one answer to
// "did this change", and a stamp next to the data can be forgotten by whoever writes it.
struct SystemTrajectoryData {
    struct {
        size_t stride = 0;                          // == shape[1] of 'backbone/secondary_structure'
        size_t count = 0;                           // == its element count
        md_secondary_structure_t* data = nullptr;   // VIEW into sys.attributes, not owned
    } secondary_structure;
    struct Render {
        size_t stride = 0; // = mol.backbone.count. Multiply frame idx with this to get the data
        size_t count = 0;  // = mol.backbone.count * num_frames. Defines the end of the data for assertions
        md_secondary_structure_t* data = nullptr;   // owned here: a presentation smoothing, not a quantity
        uint64_t fingerprint = 0;
    } secondary_structure_render;
    struct {
        size_t stride = 0;                          // == shape[1] of 'backbone/angle'
        size_t count = 0;                           // == its element count
        md_backbone_angles_t* data = nullptr;       // VIEW into sys.attributes, not owned
    } backbone_angles;
};

// A RUN OF A SYSTEM: a trajectory opened on it, and what is read and derived along its frames. Several runs of one
// system (replicas, a continuation, the candidates of a fit) share its topology and its attribute table, each under a
// prefix of its own (SystemContext::runs).
struct SystemRun {
    // "run/<name>": where the trajectory, and everything sampled along it, lives in the system's
    // sys.attributes (see RUNS in md_system.h). Empty while no trajectory is loaded. The run's
    // frame axis is "<path>/time"; freeing the trajectory removes everything under the prefix.
    // Unique within its system (run_path_from_file).
    char                path[128] = {0};

    // Playback's extraction contexts, one per pool thread and indexed by the task system's thread
    // number, so each keeps its own files open across the frames it loads. Made on first use,
    // ended with the trajectory (free_trajectory_data), since the files they hold are its.
    md_system_extract_t* frame_extract[64] = {};

    // Computing the backbone angles and secondary structure over the frames of the run (viamd.cpp)
    task_system::Job    backbone_computations;

    SystemTrajectoryData trajectory_data;
};

// WHAT IS SELECTED AND HIGHLIGHTED OF A SYSTEM, over its atoms. What the user selects spans the viewport, so the
// selection is the union of every system's; how it is selected (granularity, the windows) is the application's
// (ApplicationState::selection). Edited in place, every frame for the highlight: over any atom index, in the persistent
// allocator (system_context_init), so that they need no resizing when the system is loaded again.
struct SystemSelection {
    md_bitmask_t selection_mask {};
    md_bitmask_t highlight_mask {};
};

// THE PICKING RANGES OF A SYSTEM, reserved for it in every frame (EventType_ViamdPickingRangeReserve), one per domain
// its representations write: an index within one is the object of the system it names
struct SystemPicking {
    PickingRange atom {};       // atoms
    PickingRange bond {};       // bonds: of the graph, then those formed in the frame shown (picked_bond)
    PickingRange backbone {};   // protein backbone segments (cartoon, ribbons)
    PickingRange ring {};       // rings (circles of aromatic rings)
    PickingRange glyph {};      // the points of the cartoon of lipids and nucleic acids
    PickingRange gap {};        // the gaps of the chains, one index per residue
};

// The most systems held at once. A hard limit, so that what is kept per system and atom (picking, selection) has a
// fixed budget.
constexpr size_t SYSTEM_MAX = 8;
constexpr size_t RUN_MAX = 8;        // Of one system
constexpr size_t INSTANCE_MAX = 8;   // Of one system

struct SystemContext;

// AN INSTANCE OF A SYSTEM: the system made concrete in one of the several ways it may be shown at once. Its coordinates
// come from one place - a run, at the time shown or at a frame of its own, or the system's own coordinates - and are
// made the same way as every other instance's (the system's operations), so that instances of several runs superpose.
// It is placed in the world by a transform of its own (scene.placement), drawn as what it is shown as (periodic
// images, an assembly), with the representations which are of it (Representation::instance).
//
// A system has one at least. The first is the one the timeline, the script, the selection and the tools read
// (system_instance): what a single instance of a system was before there were several.
struct SystemInstance {
    SystemContext*      sys = nullptr;          // The system it is an instance of
    uint32_t            id = 0;                 // Within its system, never reused while it is loaded: what its
                                                // representations and the workspace name it by
    char                name[64] = "";          // Empty: called after where its coordinates come from

    // Where its coordinates come from: run (an index into its system's runs) at the time shown, or at frame when that
    // is not negative; the system's own (SystemData::state) for a run of -1
    int                 run = -1;
    int64_t             frame = -1;

    SystemScene         scene;

    // Its coordinates as drawn, made from the source by interpolate_system_state, then the operations. Its attributes
    // carry the backbone angles and secondary structure of the frame (md_system_state_backbone_angles,
    // md_system_state_secondary_structure).
    md_system_state_t   state = {};
    // Moves whenever state is made again (interpolate_system_state): what is computed of its coordinates is of a version
    uint64_t            state_version = 0;
    // The turn (about the cell centre) that the coordinates in state carry relative to the lattice of state.unitcell.
    // Identity unless keeping the orientation (or alignment_mat) turned the system. md_unitcell_t is lower triangular
    // and cannot describe a turned lattice, so this is the only record of it: wrapping and making whole only work in
    // the lattice frame (see apply_state_operations), and the box is drawn through it to stay in register with the
    // atoms. Back to identity whenever state is made again from its source.
    md_float4x4_t       state_rotation = MD_MAT4_IDENT_INIT;

    // On the GPU. The cartoon's secondary structure (blend weights between coil, helix and sheet) lives there only: it
    // is computed where it is produced and uploaded straight away (upload_secondary_structure_weights,
    // interpolate_system_state), never read back and kept nowhere else.
    md_gl_mol_t         gl_mol = {};
    // The spheres, tubes and capsules of the cartoon of the lipids and nucleic acids, built from the system beside
    // gl_mol and derived from its atoms on the GPU (see cartoon_glyphs.h)
    CartoonGlyphs       cartoon_glyphs = {};
    md_float4x4_t       unitcell_transform = MD_MAT4_IDENT_INIT;

    // The longest bond drawn without a cell (bond_length_limit in main.cpp), and the topology version it is of
    float               bond_length_limit = 0.0f;
    uint64_t            bond_length_limit_version = 0;

    bool                interpolate_system_state = false;   // state is to be made again
    uint32_t            dirty_gpu_buffers = 0;
    uint64_t            flags_hash = 0;             // what the atom flags uploaded were made of (main): its visibility, selection, highlight
    uint64_t            gpu_buffers_version = 0;    // bumped whenever gl_mol is updated, for views that cache what they render from it
    uint64_t            gpu_scene_version = 0;      // bumped as gpu_buffers_version is, but for an update of the atom flags
                                                    // alone: the image of the scene depends on them only through the
                                                    // visible and the selected atoms (representation.visibility_mask_hash,
                                                    // render.selection_hash), not through the highlighted

    // The frames of its run read around the time shown, and what was last made of them
    FrameCache          frame_cache;
    // The nearest frame last pushed to the renderer in Nearest interpolation mode - lets
    // interpolate_system_state skip redundant work when playback ticks without the displayed
    // frame changing. Reset to -1 whenever the underlying per-frame data can no longer be
    // trusted to match what was last pushed (new trajectory, secondary structure recomputed),
    // so the next call is never skipped by a stale match.
    int64_t             last_interpolated_nearest_frame = -1;
    // The (fractional) frame last interpolated, negative when there is none. A jump away from it
    // resets the temporal coherence of the backbone orientation (MolBit_ResetBackboneHistory).
    double              last_interpolated_frame = -1.0;
    // The stride of the frames the displayed state was last interpolated between (playback_stride).
    // Starting or stopping playback changes it without changing the time, so it is compared
    // against to know when the displayed state has to be made again.
    int                 last_interpolated_stride = 1;

    // The picking ranges reserved for it in this frame: the same atoms in another instance are other objects to pick
    SystemPicking       picking;
};

// A SYSTEM, as viamd holds it: what is had of each one loaded, apart from what the application has once.
//
// ApplicationState holds the systems (systems, at most SYSTEM_MAX). The first is held in place (system_context), with
// references to its parts under the names the code has always used for them (mold, system, representation, run,
// trajectory_data, file_list). Those names are the transition: a place that works on a system is given the system it
// works on (a SystemContext, or the part it needs) instead of reaching it through the application, one place at a
// time, and the references go when nothing uses them. What has not been given its system yet works on the first,
// through those names or through system_context.
//
// The others are allocated when added (system_add) and freed when removed (system_release), so that what captures one
// (a task, a cache) can hold its address. The first stays in place, empty when it is removed while others are loaded,
// and is the one the next system added is loaded into: no system is listed, drawn or saved while it is empty
// (system_loaded).
//
// What is not here yet, though it is of a system, and why:
//   timeline.x_values     the times of a run, in the unit the timeline shows: moved with the global time axis
//   script                one script for all systems, evaluated per system: decided with the script rewrite
struct ApplicationState;
struct SystemContext;

// Whether the system is drawn: it has atoms and coordinates. A topology alone has none until a file of them, or a run,
// joins it; it is a system all the same, listed, inspected and saved.
bool system_drawn(const SystemContext& sys);

struct SystemContext {
    // Its id, by which a workspace and what refers to it name it ([System] Id, System=), never reused while viamd
    // runs; and its name, which the script and the windows show it by: the stem of the file that built it, until the
    // user names it. Empty while nothing is loaded.
    uint32_t              id = 0;
    char                  name[64] = "";
    // Its directory: the outermost layer of the paths of what is of it, as the Attributes window lists them (system_dir)
    char                  dir[32] = "";

    // The application it is of: what the events about it (SystemEventPayload) reach the rest of viamd through
    ApplicationState*     app = nullptr;

    // The files the system is built from, in the order they were loaded (see file_list.h). Kept by
    // load_data_from_file, and by what frees the data they loaded (free_system_data, free_trajectory_data). A
    // workspace stores this list and loads it again in the same order.
    FileList              files;
    SystemData            data;
    // The runs opened on it, in the order they were opened, which is that of the file list's (FileList::runs): run i
    // of the one is run i of the other. Allocated when opened (open_run) and freed when closed (close_run), so that
    // what captures one (a task) can hold its address. The one shown is the run it is drawn at, and the timeline's
    // (system_run).
    SystemRun*            runs[RUN_MAX] = {};
    size_t                num_runs = 0;
    SystemOperations      operations;
    // While what loads its files places the runs it opens in its instances itself (a workspace, a system built again):
    // open_run adds none for a further run
    bool                  instances_held = false;
    // Its instances, the first always there (system_context_init): allocated, so that what captures one can hold its
    // address
    SystemInstance*       instances[INSTANCE_MAX] = {};
    size_t                num_instances = 0;
    uint32_t              next_instance_id = 1;
    bool                  visible = true;           // None of its instances is drawn while it is not
    SystemRepresentations representation;
    SystemSelection       selection;
};

// The colours of the palette (ApplicationState::palette, palette_builtin)
struct Palette {
    md_float4_t element[MD_Z_Count];
    md_float4_t monomer[MD_MONOMER_COUNT];
    // The secondary structures: a colour of each group (SecondaryStructureGroup), which its classes have unless they
    // have one of their own (ss_class_own: palette_ss_color), and of what is not known
    md_float4_t ss_unknown;
    md_float4_t ss_group[SS_GROUP_COUNT];
    md_float4_t ss_class[SS_CLASS_COUNT];
    bool        ss_class_own[SS_CLASS_COUNT];
};

struct ApplicationState {
    // --- APPLICATION ---
    application::Context app {};

#if MD_ENABLE_GPU
    md_gpu_device_t gpu_device = nullptr;

    // Evaluation scratch. Device scoped and deliberately NOT per dataset: the volume is a fixed
    // GPU_VOLUME_DIM^3 R32F texture - half a gigabyte - and every evaluation reads it straight back into its
    // own destination texture, so one is enough no matter how many datasets or representations are
    // asking. The coefficient buffer is grown to fit the widest basis loaded so far, for the same
    // reason. Whoever evaluates borrows these; nobody else frees them.
    md_gpu_stream_t gpu_stream  = nullptr;   // the device's default compute stream
    md_gpu_texture_t gpu_volume = nullptr;   // 3d R32F scratch for an evaluated orbital / density
    md_gpu_texture_t gpu_plane  = nullptr;   // GPU_PLANE_DIM^2 x 1 R32F scratch for the image of a plane
    // Where the evaluations write and their readbacks read: gpu_volume when null. An evaluation into another
    // scratch sets it for its call (gpu_eval_target_scope, viamd.cpp); a readback queued records its own.
    md_gpu_texture_t gpu_eval_target = nullptr;
    // Normals of what is evaluated (for the raycast isosurfaces, es_normals_of): when gpu_eval_normals names a GL
    // texture (RGB10_A2, the grid's dims) for the call, the orbitals' and densities' evaluations write their analytic
    // gradient into gpu_grad (RGBA32F) and set gpu_grad_written, and the readback packs it into gpu_normals (R32UI,
    // md_gto_gpu_gradient_pack) and fills that texture too. Both made for the largest grid asked for so far
    // (gpu_grad_dim), up to GPU_GRAD_MAX_VOXELS; past it no normals.
    uint32_t         gpu_eval_normals = 0;
    bool             gpu_grad_written = false;
    md_gpu_texture_t gpu_grad    = nullptr;
    md_gpu_texture_t gpu_normals = nullptr;
    int              gpu_grad_dim[3] = {};
    md_gpu_addr_t   gpu_coeff   = 0;         // AO coefficients, sized to the widest basis
    size_t          gpu_coeff_capacity = 0;

    // One in-flight readback of gpu_volume into a GL texture. A slot must outlive the call that
    // queued it, so these live here rather than in the frame arena - and here specifically because
    // every one of them references gpu_volume and gpu_stream above. Whoever destroys
    // those drains these first.
    struct GpuVolumeJob {
        bool              in_flight = false;
        ApplicationState* owner     = nullptr;
        uint32_t          tex_id    = 0;        // GL texture to receive the data
        md_gpu_texture_t  src       = nullptr;  // the scratch it is read from (gpu_volume or gpu_plane)
        md_gpu_mem_t      rb        = {};       // HOST_READ block; the data is read through rb.cpu
        size_t            size      = 0;
        uint32_t          normal_tex = 0;       // GL texture to receive the normals (gpu_eval_normals), 0: none
        md_gpu_mem_t      rb_normals = {};      // their HOST_READ block
        // A density evaluation not yet all issued (md_gto_gpu_job_t: issued in pieces over the frames): the
        // readback is queued when it is (gpu_volume_jobs_poll), of this grid
        void*             density_job = nullptr;
        md_grid_t         grid      = {};
        // The same once it is all issued, kept for its progress until the readback lands (its activity)
        void*             issued_job = nullptr;
        int64_t           t0 = 0;               // md_tick_t: when the evaluation began
        const char*       label = nullptr;      // what it evaluates, for the overlay (a static string)
    };
    // Readbacks are issued at most one per representation per change, and a change cannot be
    // requested again until the previous frame has been drawn, so a handful of slots is ample.
    // Running out falls back to blocking.
    static constexpr int GPU_VOLUME_JOB_SLOTS = 8;
    GpuVolumeJob gpu_volume_jobs[GPU_VOLUME_JOB_SLOTS] = {};
#endif

    struct {
        md_allocator_i* frame = nullptr;
        md_allocator_i* persistent = nullptr;
    } allocator;

    LoadDatasetWindowState load_dataset;

    // --- FILES ---
    struct {
        // The files of the system are the system's (SystemContext::files, file_list)
        char workspace[1024] = {0};
    } files;

    // The idea for the file load queue is to fill it with files that are dropped onto the application
    // Or passed via the commandline. The files need to be processed in a certain order (topology before trajectory)
    // It also provides a way to chain the file load dialogue.
    FileQueue file_queue = {};

    // The files of a drop which are of the atoms of more than one system: the user picks where each group of them goes
    // (draw_placement_prompt). What else the drop held is loaded; these wait, and what is queued after them, until the
    // question is answered (file_queue_process).
    struct {
        DropFiles*     drop = nullptr;                  // Asked while there is one: the drop, its files as matched
        PlacementPlan* plan = nullptr;                  // Its systems of PlacementAction_Ask are the questions
        int            choice[PLACEMENT_MAX_FILES] = {};// Per system of the plan: the candidate picked, num_candidates
                                                        // for a new system, -1 for none
        bool           answered = false;
    } placement_prompt;

    // A run of a system to be closed (the Scene window): closed by file_queue_process
    struct {
        uint32_t system = SYSTEM_ID_NONE;               // None asked while SYSTEM_ID_NONE
        size_t   run = 0;
    } run_close;

    // A file of a system to go into another system, or into a new one (the Scene window): moved by file_queue_process
    struct {
        uint32_t from = SYSTEM_ID_NONE;                 // None asked while SYSTEM_ID_NONE
        size_t   file = 0;
        uint32_t to = SYSTEM_ID_NONE;                   // SYSTEM_ID_NONE: a new system
    } file_move;

    // --- CAMERA ---
    struct {
        Camera camera{};
        TrackballControllerParam trackball_param;
        ViewParam param{};
        CameraMode mode = CameraMode::Perspective;

        struct {
            md_float2_t sequence[JITTER_SEQUENCE_MAX] {};
        } jitter;

		ViewTransform target = {};

        // Where it is drawn: the dockspace's central node, what the windows docked about it leave (CreateDockspace),
        // in points, as ImGui places things on the screen. Its size is that of the image made of it (at the framebuffer's
        // pixels per point); windows that float over it cover it.
        md_float2_t screen_pos  = {0, 0};
        md_float2_t screen_size = {0, 0};

        // The copies each system is drawn in are its own (SystemScene)
    } view;

    struct {
        bool  hide_gui = true;
        ScreenshotResolution resolution = ScreenshotResolution::View;
        int   res_x = 1920;
        int   res_y = 1080;
        md_str_t path_to_file = {};    // Being made: saved once the image has settled (render())
    } screenshot;

    struct {
        md_gl_shaders_t shaders = {};
        md_gl_shaders_t shaders_marks = {};    // The marks of the selection and the highlight (the mark pass)
    } gl;

    // --- THE SYSTEMS ---
    // The first system (SystemContext), held in place, and its parts under the names the code knows them by (see
    // SystemContext). It is there whether or not anything is loaded into it.
    SystemContext system_context = {};
    SystemData&            mold            = system_context.data;
    SystemRepresentations& representation  = system_context.representation;
    FileList&              file_list       = system_context.files;

    // Every system, in the order the Scene window lists them, the first being system_context. Not owned past the
    // first: what adds one allocates it (none does yet).
    SystemContext* systems[SYSTEM_MAX] = { &system_context };
    size_t         num_systems = 1;
    uint32_t       next_system_id = 1;      // The first is 0, as a workspace before several systems names it

    // The quantum chemistry system (qm_system): its id, and whether it was chosen (by the user, a workspace) rather
    // than taken for having electronic structure
    struct {
        uint32_t system_id = 0;
        bool     chosen = false;
    } qm;

    // The script property hovered in any view of it (a plot, the script editor), by identifier. A copy,
    // not a view: the identifiers live in the script IR, which is rebuilt when the script recompiles,
    // and the hover outlives that.
    char  hovered_property_label[64] = "";
    int   hovered_property_pop_idx = -1;

    // What was under the cursor when the context menu of the viewport opened, which its entries are about: the menu
    // covers the cursor once it is open
    PickingHit context_menu_hit = {};

    // --- ASYNC TASKS HANDLES ---
    struct {
        task_system::Job evaluate;     // of the script

    } tasks;

    // --- ATOM SELECTION ---
    struct {
        SelectionGranularity granularity = SelectionGranularity::Atom;
        SingleSelectionSequence single_selection_sequence;

        // The masks of what is selected and highlighted are each system's (SystemContext::selection), and the sets kept of
        // it are groups of the user's of each system (store_selection_as_group). Everything here holds sets, replaced
        // whole, but for the masks of the windows (query, grow), edited in place.

        // The same over the frames of the run: the frames selected (a temporal selection, made in the Shapespace and
        // shown in the timeline) and those highlighted, which whoever sets clears again
        md_bitmask_t frame_mask {};
        md_bitmask_t frame_highlight_mask {};

        // How the selected and the highlighted atoms are marked over the image, in display colours (see render()): a
        // tint where they are in view, a fainter fill where they are behind something else, and the outline of their
        // silhouette, in view or not. The alpha is how much of the colour. While something is selected, the atoms that
        // are not keep 'saturation' of their colour (AtomBit_Dimmed): the scene, not what is drawn over it.
        struct {
            struct {
                md_float4_t visible = {0.15f, 0.40f, 1.00f, 0.35f};
                md_float4_t hidden  = {0.15f, 0.40f, 1.00f, 0.15f};
                md_float4_t outline = {0.25f, 0.55f, 1.00f, 0.90f};
            } selection;

            struct {
                md_float4_t visible = {1.00f, 0.85f, 0.10f, 0.35f};
                md_float4_t hidden  = {1.00f, 0.85f, 0.10f, 0.15f};
                md_float4_t outline = {1.00f, 0.85f, 0.10f, 0.90f};
            } highlight;

            // The colours of the groups of the highlight after the first, which has the highlight's own; with the
            // alphas of its colours. In an order in which each differs much from the one before it, and the last
            // from the highlight's yellow.
            md_float3_t highlight_groups[HIGHLIGHT_GROUP_COLORS - 1] = {
                {0.25f, 0.85f, 0.35f},  // green
                {0.95f, 0.30f, 0.80f},  // magenta
                {0.15f, 0.85f, 0.95f},  // cyan
                {1.00f, 0.50f, 0.05f},  // orange
                {0.60f, 0.40f, 1.00f},  // violet
                {0.70f, 0.95f, 0.15f},  // lime
                {0.95f, 0.20f, 0.20f},  // red
            };

            float saturation = 0.3f;        // Taken into the G-buffer's shaders as they are made (output_snippet)
        } color;

        double highlight_time = 0.0;        // when a highlight last appeared where there was none (ImGui::GetTime()): it fades in from there

        // The groups of the highlight, when what set it has them: the elements of the expression hovered, residue(:) a
        // group per residue. Group i is marked in colour i % HIGHLIGHT_GROUP_COLORS (AtomBit_GroupMask): 0 is the
        // highlight's own, k > 0 color.highlight_groups[k - 1]. Not owned: set with the highlight every frame
        // (script_vis_update), cleared at its start. key: what they are, 0 for none.
        struct {
            const md_bitset_t* sets = nullptr;
            size_t count = 0;
            uint64_t key = 0;
        } highlight_groups;

        // The highlight of the groups of an energy term (highlight_energy_term_groups) is set every frame its term is
        // hovered, wherever that is; the first frame it is not, script_vis_update clears it, so that it never outlives
        // the hover in a window which does not clear the highlight itself
        bool term_groups_highlighted = false;       // this frame
        bool term_groups_highlighted_prev = false;  // the frame before

        struct {
            AtomFilter filter;
            md_bitmask_t mask = {};
            bool query_invalid = true;
            bool show_window = false;
        } query;

        struct {
            md_bitmask_t mask = {};
            SelectionGrowthMode mode = SelectionGrowthMode::CovalentBond;
            float extent = 1;
            bool mask_invalid = true;
            bool show_window = false;
        } grow;
    } selection;

    // The palette (Palette): the user's, of every system and workspace. Moves palette_version with every change.
    Palette  palette;
    uint64_t palette_version = 1;
    bool     show_palette_window = false;   // Visuals > Palette

    // --- GBUFFER ---
    // Primary gbuffer
    GBuffer gbuffer {};

    struct {
		immediate::Queue* world = nullptr;
		immediate::Queue* overlay = nullptr;
    } gfx;

    // --- PICKING ---
    PickingSurface picking_surface {};      // Surface for main viewport picking

    // The ranges reserved for each system are its own (SystemContext::picking)

    // The glyph under the cursor (PickingDomain_Glyph): its representation's id and its point. Marked in the marks of the
    // highlight (glyphs_draw_marks). rep_id 0: none.
    struct {
        uint32_t rep_id = 0;
        uint32_t system = SYSTEM_ID_NONE;   // whose representation it is: a representation's id is its system's own
        uint32_t point  = 0;
    } glyph_hover;
    // The same, hovered elsewhere than in the view (a row of the GOSTSHYP table): set by what draws it, every frame it
    // is hovered. Marked in place of the one in the view.
    struct {
        uint32_t rep_id = 0;
        uint32_t system = SYSTEM_ID_NONE;
        uint32_t point  = 0;
    } glyph_hover_ui;

    PickingHandler picking_handler {};      // Handler for managing picking interactions

    // --- ANIMATION ---
    struct {
        // The time shown, on the clock (clock.h): in the unit of the reference's run (clock_reference), time_unit, in
        // which it is kept when another run becomes the reference. Every instance which follows the time shows its run
        // at it. What "the current frame" is to whatever reads one run is that run's frame at it (clock_frame).
        double time = 0.0;
        md_unit_t time_unit = {};
        float fps = 10.f;   // Frames of the reference per second
        float tension = 0.0f;
        // While playing, the displayed state is interpolated between every stride:th frame only (0, stride,
        // 2*stride, ... and the last), which smooths the motion out at the cost of fidelity to the frames
        // skipped. A stopped animation always shows the trajectory itself. See playback_stride.
        int stride = 1;
        InterpolationMode interpolation = InterpolationMode::CubicSpline;
        PlaybackMode mode = PlaybackMode::Stopped;

        // What playing does at the ends (clock_advance), of the span it goes over: the clock's, or of a fixed timeline
        // period (by id, 0 for none: playback_range). direction is that of a bounce, +1 or -1, times the speed's sign.
        PlaybackLoop loop = PlaybackLoop_Once;
        uint32_t loop_period = 0;
        int direction = 1;

        // The transport, the slim bar of the time over the view (draw_transport, main.cpp): where it is
        TransportAnchor transport = TransportAnchor_Top;
    } animation;

    // --- TIMELINE---
    struct {
        // The named periods and markers along the timeline (see timeline_periods.h)
        TimelinePeriods periods;

        // The visible span of the run's time ("<run>/time"), in its own unit like everything on the
        // timeline: the display unit only changes how a time is printed (frame_axis.h)
        struct {
            double beg_x = 0;
            double end_x = 1;
        } view_range;

        // What each subplot draws, by attribute path (see plot_series.h)
        PlotSubplot subplots[PLOT_MAX_SUBPLOTS];
        int num_subplots = 1;

        bool show_window = false;
    } timeline;

    // --- DISTRIBUTIONS ---
    struct {
        // What each subplot draws, by attribute path (see plot_series.h)
        PlotSubplot subplots[PLOT_MAX_SUBPLOTS];
        int num_subplots = 1;

        bool show_window = false;
    } distributions;

    // --- VISUALS ---
    struct {
        struct {
            md_float3_t color = {1, 1, 1};   // shown as picked: the scene is put over it after the display transform
        } background;

        struct {
            postprocess_pipeline::ShadingMode mode = postprocess_pipeline::ShadingMode_PhysicallyBased;
            postprocess_pipeline::ColorSpace colors = postprocess_pipeline::ColorSpace_Linear;   // the look the colours had before
            float light = 1.0f;         // strength of the directional light
            float environment = 1.0f;   // strength of the light from the surroundings (physically based)
            float ambient = 0.3f;       // basic
            float specular = 0.25f;     // basic
        } shading;

        // The visibilities of the two lights, found in screen space: the shadows of the directional light and the
        // ambient occlusion of the surroundings. They share a temporal filter.
        struct {
            bool enabled = true;
            float softness = 10.0f;     // the light's angular radius, degrees
        } shadows;

        struct { 
            bool enabled = true;
            float intensity = 5.0f;
            postprocess_pipeline::SsaoMode mode = postprocess_pipeline::SsaoMode_Performance;
        } ssao;

        bool occlusion_temporal = true; // accumulate the shadows and the ambient occlusion over frames

#if EXPERIMENTAL_CONE_TRACED_AO == 1
        struct {
            bool enabled = true;
            float intensity = 1.0f;
            float step_scale = 1.0f;
        } cone_traced_ao;
#endif

        struct {
            bool enabled = true;
            postprocess_pipeline::Tonemapper tonemapper = postprocess_pipeline::Tonemapper_ACES;   // the look the shading had before
            float exposure = 1.0f;
            float gamma = 2.2f;
        } tonemapping;

        struct {
            bool enabled = false;
            float focus_depth = 10.0f;
            float aperture = 1.0f;      // percent of the view height (CoC of an object at infinity)
        } dof;

        // Whether jagged edges are smoothed, and how: the mode is kept while off. Every mode smooths, so there is no
        // being on with nothing applied. fxaa_active() and taa_active() say which pass runs.
        struct {
            bool enabled = true;
            AntiAliasingMode mode = AntiAliasingMode::TAA;
        } anti_aliasing;

        // The temporal resolve's settings, used in TAA mode
        struct {
            bool jitter = true;
            // The jitter's positions at the full resolution (a screenshot's): one of jitter_samples_options (main.cpp)
            int jitter_samples = 32;
            float feedback_min = 0.80f;
            float feedback_max = 0.95f;

            struct {
                bool enabled = true;
                float motion_scale = 1.0f;
            } motion_blur;
        } temporal_aa;

        struct {
            bool enabled = true;
            float weight = 1.0f;
        } sharpen;

        struct {
            bool draw_control_points = false;
            bool draw_spline = false;
        } spline;
    } visuals;

    bool show_representation_window = false;

    // Persisted in the ImGui .ini; see app_settings.h. Bound in main().
    struct {
        bool keep_representations = false;
        bool exact_isosurfaces = false;     // volume::IsoRenderDesc::iso.exact for every isosurface view
        // The look a session starts with (Visuals > Save as Default), as the factory sets it until saved. A workspace
        // has its own, and loading one leaves these be. What depends on the computer rather than on the figure
        // (ambient occlusion, anti-aliasing, VSync) is bound to the .ini as it is, in visuals.
        struct {
            decltype(visuals.background)  background;
            decltype(visuals.shading)     shading;
            decltype(visuals.tonemapping) tonemapping;
            decltype(visuals.dof)         dof;
        } look;
        char palette[4096] = "";            // The palette as the .ini keeps it (palette_write)
        int  isosurface_rendering = 0;      // IsosurfaceRendering of the Orbital Grid and Transition Analysis panels (a
                                            // representation has its own: ElectronicStructureRepresentation::rendering)
        float font_size = 18.0f;    // Matches the size the default font is baked at.
        // The scene rendered at this fraction of the display's resolution (each side) and upscaled to it: by the TAA,
        // which rebuilds the detail over frames, else a frame on its own. Screenshots are at the full resolution.
        float render_scale = 1.0f;
        // The scene rendered only while something in it changes, and settles; else the image made last is shown again
        // (render: lazy rendering)
        bool lazy_rendering = true;
        // The script language a session starts in, as last chosen under Settings. A workspace is written in its own,
        // and loading one switches to it without changing this
        script::Language script_language = script::Language::V1;
    } settings;


    struct {
        md_float4_t point_color      = {1,0,0,0.8f};
        md_float4_t line_color       = {0,0,0,0.6f};
        md_float4_t triangle_color   = {0.55f,0.55f,1.0f,0.5f};
        md_float4_t text_color       = {1,1,1,1};
        md_float4_t text_bg_color    = {0,0,0,0.5f};

        uint64_t text_hash = 0;     // Of the editor's text in the last frame: a change compiles it again (main.cpp)

        // The workspace's language: what the editor's text is compiled as, and what the filters are written in
        script::Language language = script::Language::V1;
        // The script in the other language, kept while this one is used: the editor shows one buffer per language,
        // and switching swaps them (script_switch_language)
        std::string other_text;

        // ir: the last compilation of the editor's text, what its markers refer to. eval_ir: the compilation the
        // evaluation (eval) was made from, what the plots show. Often the same IR. Assigned directly: an IR is
        // never freed while either points at it (see script_ir_collect).
        script::Program*  ir = nullptr;
        script::Program*  eval_ir = nullptr;

        // Every IR made (script_ir_create), freed by script_ir_collect once nothing uses it (no field, no task)
        md_array(script::Program*) all_irs = nullptr;

        // The evaluation of the run the first instance of the first system showed when it was made (eval_run): every
        // property over every frame, and the distributions and volumes shown over timeline periods in slots over their
        // frames (SCRIPT_SLOT_PERIOD). What the tools, the hover and a series without a run read.
        script::Eval*     eval = nullptr;
        char              eval_run[128] = "";
        // One more of eval_ir per other run of the first system an instance shows (script_run_evals_sync): what a
        // series of that run shows (SeriesKey::run). Evaluated with eval, after it.
        md_array(ScriptRunEval) run_evals = nullptr;
        md_array(uint32_t) eval_frames = nullptr;   // what the running evaluation task works through
        const script::Eval* evaluating = nullptr;   // the evaluation the task runs over, NULL when none runs
        bool eval_queued = false;                   // eval is to be evaluated, once what runs now has ended
        // The visualization: the target is set anew every frame, the job computes it off the main thread, and
        // what is shown is the last result while it is still for the target (see script_vis_update)
        ScriptVisTarget vis_target = {};
        ScriptVisJob* vis_running = nullptr;
        ScriptVisJob* vis_result = nullptr;
        task_system::Job vis_task;
        const script::Drawing* vis_shown = nullptr;
        // The plot of a value over the trajectory hovered, read from eval on the main thread (script_vis_value_plot)
        script::ValuePlot trajectory_plot = {};
        uint64_t trajectory_plot_key = 0;
        md_allocator_i* trajectory_plot_arena = nullptr;

        // Controls current subindex being visualized (if array)
        int sub_idx = -1;

        bool compile_ir = false;
        bool eval_init = false;
        bool evaluate = false;     // a new evaluation, to be started
        double time_since_last_change = 0.0;
        size_t   num_compile_markers = 0;   // the editor's markers from compiling; the evaluation's come after them
        size_t   first_mark_marker = 0;     // where the marks' markers start among them, in the order of script::marks()
        uint64_t marks_units_version = 0;   // display_units::version() the marks' texts were written in
        uint64_t eval_errors_version = 0;   // of the evaluation's errors marked (script::eval_errors)
        uint64_t ir_fingerprint = 0;
    } script;

    struct {
        bool show_window = false;
        // Of which system (by id, system_or_first) and which of its instances (by id, instance_resolve): its run's
        // frames, at the frame it shows
        uint32_t system_id = 0;
        uint32_t instance_id = 0;
        int selected_atom_filter = 0;
        int selected_traj_filter = 0;       // all frames, the current frame, or the frames of selected_traj_period
        uint32_t selected_traj_period = 0;   // a timeline period, by id
        int selected_file_format = 0;
        struct {
            AtomFilter filter;
            md_bitset_t mask = {0};
            bool requires_evaluation = false;
        } query;
    } structure_export;

    bool show_script_window = true;
    bool show_script_reference_window = false;
    bool show_debug_window = false;
    bool show_profiler_window = false;
    bool profiler_keep_isomeshes = false;   // Windows > Profiler: the isosurface meshes drawn every frame (isomesh::RenderDesc::keep_drawing)
    bool show_property_export_window = false;

    // --- LAZY RENDERING ---
    // The main view's image is made again only when what it is made from changes (scene_key), and for as many frames
    // after as it takes to settle (the TAA's history, the temporal ambient occlusion, the meshed isosurfaces); else the
    // image made last (gbuffer.tex.image) is shown again. Input (a click, a drag, a key, the wheel) counts as a change:
    // what it edits need not be in the key. The marks are drawn again when they change, the scene or not.
    struct {
        uint64_t scene_key = 0;         // of the image kept
        uint64_t marks_key = 0;         // of the marks kept, 0 for none
        uint32_t settle_frames = 0;     // still to make after the last change
        bool image_valid = false;       // gbuffer.tex.image holds the scene (not after the G-buffer is made anew)
        uint64_t version = 0;           // render_invalidate()

        // Of the main loop's: for the keys
        uint64_t selection_hash = 0;
        uint64_t highlight_hash = 0;

        // Of the last frames, whether the scene was made in each (bit i: i frames ago), for the profiler
        uint64_t history = 0;
    } render;

    TextEditor editor = {};
    script_editor::Markers editor_markers = {};
    bool editor_focused = false;    // the script editor had keyboard focus when it was last drawn
};

struct ViamdEventHandler : viamd::EventHandler {
    ApplicationState* state = nullptr;

    explicit ViamdEventHandler(ApplicationState* s) : state(s) {
        MD_ASSERT(state);
        viamd::event_system_register_handler(*this);
    }

    void process_events(const viamd::Event* events, size_t num_events) final;
};

struct LoadDataPayload {
    ApplicationState* app_state;
    SystemContext*    system;       // What was loaded into
    loader::LoaderState loader_state;
    md_str_t path_to_file;
};

// The system a system event (its payload EventPayloadType_System, or EventPayloadType_LoadData) is of, null for
// another event
static inline const SystemContext* event_system(const viamd::Event& e) {
    if (e.payload_type == viamd::EventPayloadType_System) return (const SystemContext*)e.payload;
    if (e.payload_type == viamd::EventPayloadType_LoadData && e.payload) return ((const LoadDataPayload*)e.payload)->system;
    return nullptr;
}

// Whether a system event is of the first system (ApplicationState::system_context). What is not per system yet (some
// components, the System window) acts on that one alone, and leaves the events of the others be.
static inline bool event_of_first_system(const viamd::Event& e) {
    const SystemContext* sys = event_system(e);
    return sys && sys->app && sys == &sys->app->system_context;
}

// A picking tooltip, gathered from every handler that recognises the hit and drawn as one window: a title names what
// is under the cursor, then rows give one property each, drawn as two aligned columns (key, value). Handlers add to it
// with tooltip_title and tooltip_row; a title after rows begins a new section. What goes in one: PICKING TOOLTIP in
// viamd.cpp.
enum class TooltipLineKind : uint8_t { Title, Row };

struct TooltipLine {
    TooltipLineKind kind = TooltipLineKind::Title;
    md_str_t key  = {};    // Row only, and must live for the frame (a literal, a label of the system, frame memory)
    md_str_t text = {};
};

struct PickingTooltipTextRequest {
    const ApplicationState& app;
    const PickingHit& hit;
    md_allocator_i* alloc = nullptr;    // Frame memory, for the lines and their text
    md_array(TooltipLine) lines = 0;
};

// In UTF-8 escapes: not every compiler reads the source as UTF-8
#define TOOLTIP_SEP      " \xC2\xB7 "   // ' · ' between the parts of a title
#define TOOLTIP_ANGSTROM "\xC3\x85"     // Å

void tooltip_title(PickingTooltipTextRequest* req, const char* format, ...);
void tooltip_row(PickingTooltipTextRequest* req, md_str_t key, const char* format, ...);

enum ViewFitRound {
    ViewFitRound_Highlight = 0,
    ViewFitRound_Selection,
    ViewFitRound_Visible,
};

// The payload of EventType_ViamdRenderKey: what a handler draws into the main view from, folded into key
struct RenderKey {
    const ApplicationState* state = nullptr;
    uint64_t key = 0;
};

// This is the supplied event payload for requesting a view fit.
// The idea is to construct a weighted point cloud of whatever is 'highlighted'
// And then the view fit system can decide how to best fit the view based on that point cloud.
// This is dispatched in rounds starting with the top most priority being Highlight, then Selection and lastly Visible.
// If the point cloud is empty for a certain mode, the next mode will be dispatched, until a non empty point cloud is found or all modes are exhausted.

struct ViewFitRequest {
    const ApplicationState& app;
	uint64_t surface_id = 0;        // Source of the view fit request, can be used to filter out irrelevant requests or for debugging.
	md_array(md_float4_t) xyzw = nullptr;
	md_allocator_i* alloc = nullptr;
    ViewFitRound round;
};

enum class InteractionSelectionMode {
    None,
    Append,
    Remove,
};

struct InteractionSurfaceState {
    InteractionSurfaceID surface_id = 0;
    ImGuiID item_id = 0; // Internal ImGui ID for the *button* used to control interactions.

    bool hovered = false;
    bool active = false;
    bool activated = false;
    bool deactivated = false;

    InteractionSelectionMode selection_mode = InteractionSelectionMode::None;

    md_float2_t mouse_local  = {};
    md_float2_t surface_size = {};

    // Local coordinates
    md_float2_t region_min = {};
    md_float2_t region_max = {};
};

static inline const ImVec4& vec_cast(const md_float4_t& v) { return *(const ImVec4*)(&v); }
static inline const md_float4_t& vec_cast(const ImVec4& v) { return *(const md_float4_t*)(&v); }
static inline const ImVec2& vec_cast(const md_float2_t& v) { return *(const ImVec2*)(&v); }
static inline const md_float2_t& vec_cast(const ImVec2& v) { return *(const md_float2_t*)(&v); }

static inline ImVec4& vec_cast(md_float4_t& v) { return *(ImVec4*)(&v); }
static inline md_float4_t& vec_cast(ImVec4& v) { return *(md_float4_t*)(&v); }
static inline ImVec2& vec_cast(md_float2_t& v) { return *(ImVec2*)(&v); }
static inline md_float2_t& vec_cast(ImVec2& v) { return *(md_float2_t*)(&v); }

// Applies op to a mask (the selection, the highlight) with the set s as the operand
static inline void modify_field(md_bitmask_t* m, md_bitset_t s, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        md_bitmask_or(m, s);
        break;
    case SelectionOperator::And:
        md_bitmask_and(m, s);
        break;
    case SelectionOperator::AndNot:
    case SelectionOperator::Clear:
        md_bitmask_andnot(m, s);
        break;
    case SelectionOperator::Set:
        md_bitmask_assign(m, s);
        break;
    default:
        MD_ASSERT(false);
    }
}

static inline void modify_field(md_bitmask_t* m, md_urange_t range, SelectionOperator op) {
    modify_field(m, md_bitset_range(range.beg, range.end), op);
}

// Applies op to the frames selected (ApplicationState::selection.frame_mask)
static inline void frame_selection_modify(ApplicationState* app, md_bitset_t frames, SelectionOperator op) {
    modify_field(&app->selection.frame_mask, frames, op);
}

static inline void frame_selection_clear(ApplicationState* app) {
    md_bitmask_clear(&app->selection.frame_mask);
}

static inline void mask_set_molecule(md_bitmask_t* mask, size_t mol_idx, const md_system_t& sys) {
    const md_molecule_t mol = md_system_molecule(&sys, mol_idx);
    for (size_t k = 0; k < mol.count; ++k) md_bitmask_set_bit(mask, (uint32_t)mol.atom_idx[k]);
}

// A fragment is not a range of atoms (its hydrogens are interleaved with other atoms), so its atoms are set one by one
static inline void mask_set_fragment(md_bitmask_t* mask, md_fragment_idx_t frag, const md_system_t& sys) {
    for (md_fragment_atom_iter_t it = md_fragment_atom_iter(&sys, frag); md_fragment_atom_iter_next(&it);) {
        md_bitmask_set_bit(mask, (uint32_t)it.atom);
    }
}

static inline void mask_set_chain(md_bitmask_t* mask, size_t chain_idx, const md_system_t& sys) {
    const md_urange_t* end = md_system_chain_run_end(&sys, chain_idx);
    for (const md_urange_t* it = md_system_chain_run_beg(&sys, chain_idx); it && it < end; ++it) {
        const md_urange_t range = md_system_chain_run_atom_range(&sys, *it);
        if (range.beg < range.end) md_bitmask_set_range(mask, range.beg, range.end);
    }
}

// What Chain granularity takes of a residue: its chain (md_chain_data_t), else its part (the data gives it no chain
// id: a strand of a polymer is a part, as polymer(i) of the script has it, and so is a ligand or a water), else the
// residue itself
static inline void mask_set_chain_of_res(md_bitmask_t* mask, size_t res_idx, const md_system_t& sys) {
    const md_chain_idx_t chain = md_system_chain_find_by_residue_idx(&sys, res_idx);
    if (chain >= 0) {
        mask_set_chain(mask, (size_t)chain, sys);
        return;
    }
    const md_part_idx_t part = md_system_part_find_by_residue_idx(&sys, res_idx);
    const md_urange_t r = part >= 0 ? md_system_part_atom_range(&sys, (size_t)part) : md_system_residue_atom_range(&sys, res_idx);
    if (r.beg < r.end) md_bitmask_set_range(mask, r.beg, r.end);
}

// Grows the mask to whole fragments, residues, parts, chains or molecules: those it has any atom of. A chain is not one range of atoms (the waters
// filed under a chain come after every polymer), so what is touched is found first and set after.
static inline void grow_mask_by_selection_granularity(md_bitmask_t* mask, SelectionGranularity granularity, const md_system_t& sys) {
    MD_ASSERT(mask);
    if (granularity == SelectionGranularity::Atom || md_bitmask_empty(mask)) return;

    switch(granularity) {
    case SelectionGranularity::Fragment: {
        // Without fragments (a system whose molecules were not inferred) a fragment is taken to be its residue
        const size_t num_frags = md_system_fragment_count(&sys);
        if (num_frags == 0) {
            grow_mask_by_selection_granularity(mask, SelectionGranularity::Residue, sys);
            break;
        }
        md_allocator_i* alloc = md_get_heap_allocator();
        md_bitmask_t touched = md_bitmask_create(0, (uint32_t)num_frags, alloc);
        md_bitset_iter_t it = md_bitmask_iter(mask);
        while (md_bitset_iter_next(&it)) {
            const md_fragment_idx_t frag = md_system_fragment_find_by_atom_idx(&sys, it.idx);
            if (frag != MD_FRAGMENT_NONE) md_bitmask_set_bit(&touched, frag);
        }
        md_bitset_iter_t fit = md_bitmask_iter(&touched);
        while (md_bitset_iter_next(&fit)) {
            mask_set_fragment(mask, (md_fragment_idx_t)fit.idx, sys);
        }
        md_bitmask_free(&touched);
        break;
    }
    case SelectionGranularity::Residue:
        for (size_t i = 0; i < md_system_residue_count(&sys); ++i) {
            md_urange_t range = md_system_residue_atom_range(&sys, i);
            if (md_bitmask_test_any_range(mask, range.beg, range.end)) {
                md_bitmask_set_range(mask, range.beg, range.end);
            }
        }
        break;
    case SelectionGranularity::Part:
        // Without parts (a system whose parts were not inferred) a part is taken to be its residue
        if (md_system_part_count(&sys) == 0) {
            grow_mask_by_selection_granularity(mask, SelectionGranularity::Residue, sys);
            break;
        }
        for (size_t i = 0; i < md_system_part_count(&sys); ++i) {
            const md_urange_t range = md_system_part_atom_range(&sys, i);
            if (range.beg < range.end && md_bitmask_test_any_range(mask, range.beg, range.end)) {
                md_bitmask_set_range(mask, range.beg, range.end);
            }
        }
        break;
    case SelectionGranularity::Molecule: {
        // Without molecules a molecule is taken to be its part, as at Chain granularity outside a chain
        if (sys.molecule.count == 0) {
            grow_mask_by_selection_granularity(mask, SelectionGranularity::Part, sys);
            break;
        }
        // A molecule's atoms need not be one range, so each molecule is taken once whatever order its atoms come in
        md_allocator_i* alloc = md_get_heap_allocator();
        md_bitmask_t touched = md_bitmask_create(0, (uint32_t)sys.molecule.count, alloc);
        md_bitset_iter_t it = md_bitmask_iter(mask);
        while (md_bitset_iter_next(&it)) {
            const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, it.idx);
            if (mol >= 0) md_bitmask_set_bit(&touched, (uint32_t)mol);
        }
        md_bitset_iter_t mit = md_bitmask_iter(&touched);
        while (md_bitset_iter_next(&mit)) {
            mask_set_molecule(mask, mit.idx, sys);
        }
        md_bitmask_free(&touched);
        break;
    }
    default:
        MD_ASSERT(false);
    }
}

// Sets the bits of what holds atom_idx at the granularity: the atom, its fragment (its residue where it is in none),
// its residue, its part, its chain (its part where it is in no chain), or its molecule. Cheaper than setting a single bit followed by
// grow_mask_by_selection_granularity, as it only looks up the owner.
static inline void mask_set_atom_by_selection_granularity(md_bitmask_t* mask, size_t atom_idx, SelectionGranularity granularity, const md_system_t& sys) {
    MD_ASSERT(mask);
    if (atom_idx >= sys.atom.count) return;

    md_urange_t range = {(uint32_t)atom_idx, (uint32_t)atom_idx + 1};
    switch (granularity) {
    case SelectionGranularity::Atom:
        break;
    case SelectionGranularity::Fragment: {
        const md_fragment_idx_t frag = md_system_fragment_find_by_atom_idx(&sys, atom_idx);
        if (frag != MD_FRAGMENT_NONE) {
            mask_set_fragment(mask, frag, sys);
            return;
        }
        const md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(&sys, atom_idx);
        if (res_idx >= 0) range = md_system_residue_atom_range(&sys, res_idx);
        break;
    }
    case SelectionGranularity::Residue: {
        md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(&sys, atom_idx);
        if (res_idx >= 0) {
            range = md_system_residue_atom_range(&sys, res_idx);
        }
        break;
    }
    case SelectionGranularity::Part: {
        const int32_t part_idx = md_system_part_find_by_atom_idx(&sys, atom_idx);
        if (part_idx >= 0) {
            range = md_system_part_atom_range(&sys, (size_t)part_idx);
        } else {
            const md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(&sys, atom_idx);
            if (res_idx >= 0) range = md_system_residue_atom_range(&sys, res_idx);
        }
        break;
    }
    case SelectionGranularity::Molecule: {
        const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, atom_idx);
        if (mol >= 0) {
            mask_set_molecule(mask, (size_t)mol, sys);
            return;
        }
        const md_part_idx_t part_idx = md_system_part_find_by_atom_idx(&sys, atom_idx);
        if (part_idx >= 0) range = md_system_part_atom_range(&sys, (size_t)part_idx);
        break;
    }
    default:
        MD_ASSERT(false);
    }
    md_bitmask_set_range(mask, range.beg, range.end);
}

// Sets the bits of a residue, or of what holds it at a coarser granularity, for a hit that names a residue rather
// than an atom (a backbone segment of the cartoon or ribbons). There is no atom to pick there, so Atom granularity
// selects the residue as well.
static inline void mask_set_residue_by_selection_granularity(md_bitmask_t* mask, size_t res_idx, SelectionGranularity granularity, const md_system_t& sys) {
    MD_ASSERT(mask);
    if (res_idx >= sys.residue.count) return;

    const md_urange_t range = md_system_residue_atom_range(&sys, res_idx);
    if (range.beg == range.end) return;
    if (granularity == SelectionGranularity::Atom) {
        granularity = SelectionGranularity::Residue;
    }
    // Any atom of the residue finds the residue, and its chain
    mask_set_atom_by_selection_granularity(mask, range.beg, granularity, sys);
}

// Sets the bits of the atoms of a ring (count atoms), each grown by the granularity, for a hit on the circle of an
// aromatic ring (picked_ring): at Atom granularity the ring itself, at coarser ones what its atoms are in
static inline void mask_set_ring_by_selection_granularity(md_bitmask_t* mask, const int32_t* atoms, size_t count, SelectionGranularity granularity, const md_system_t& sys) {
    MD_ASSERT(mask);
    for (size_t i = 0; i < count; ++i) {
        if (atoms[i] >= 0 && (size_t)atoms[i] < sys.atom.count) {
            mask_set_atom_by_selection_granularity(mask, (size_t)atoms[i], granularity, sys);
        }
    }
}

static inline void single_selection_sequence_clear(SingleSelectionSequence* seq) {
    MD_ASSERT(seq);
    md_memset(seq->idx, -1, sizeof(seq->idx));
}

static inline void single_selection_sequence_push_idx(SingleSelectionSequence* seq, int32_t idx) {
    MD_ASSERT(seq);
    for (size_t i = 0; i < MD_COUNTOF(seq->idx); ++i) {
        // Already in the sequence: keep the position it was first picked at. Appending it again
        // would desync the sequence from the selection mask, which is an OR and does not grow.
        if (seq->idx[i] == idx) {
            break;
        }
        if (seq->idx[i] == -1) {
            seq->idx[i] = idx;
            break;
        }
    }
}

static inline void single_selection_sequence_pop_idx(SingleSelectionSequence* seq, int32_t idx) {
    MD_ASSERT(seq);
    for (size_t i = 0; i < MD_COUNTOF(seq->idx); ++i) {
        if (seq->idx[i] == idx) {
            for (size_t j = i; j < MD_COUNTOF(seq->idx) - 1; ++j) {
                seq->idx[j] = seq->idx[j+1];
            }
            seq->idx[MD_COUNTOF(seq->idx)-1] = -1;
            break;
        }
    }
}

static inline void single_selection_sequence_pop_back(SingleSelectionSequence* seq) {
    MD_ASSERT(seq);
    size_t i = 0;
    for (; i < MD_COUNTOF(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    if (i > 0) {
        seq->idx[i-1] = -1;
    }
}

static inline int32_t single_selection_sequence_last(const SingleSelectionSequence* seq) {
    MD_ASSERT(seq);
    size_t i = 0;
    for (; i < MD_COUNTOF(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    if (i > 0) {
        return seq->idx[i-1];
    }
    return -1;
}

static inline size_t single_selection_sequence_count(const SingleSelectionSequence* seq) {
    size_t i = 0;
    for (; i < MD_COUNTOF(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    return i;
}

static inline uint64_t generate_fingerprint() {
    return (uint64_t)md_tick_now();
}

void draw_picking_tooltip_window(const PickingHit& hit, const ApplicationState& state);

//void extract_picking_data(PickingData& out_picking, GBuffer& gbuffer, const md_float2_t& coord, const md_float4x4_t& inv_MVP);

void interrupt_async_tasks(ApplicationState* state);

// Dataset loading


// Loads a file into the system sys (see file_list.h), what it is to the system told by what it and the system are:
//   a system of its own (a structure, a topology, a restart)   builds sys when it has no atoms, else joins its atoms
//                                                              (md_system_merge), and fails when they are not its
//   a trajectory                                               opens the run (sys alone: a system of points first)
//   data of a system's atoms (an index file, QM data of them)  adds to sys
//   a series (role FileRole_Series)                            joins run 'along' of sys, the one shown for -1
// opens_own_run: a system file of several frames opens them as the run too. A file that loads is recorded in the
// system's files. loaded: the system the file is, loaded already (to match it: DropFiles), which is taken rather than
// read again; NULL to read it.
LoadResult load_data(ApplicationState* app, SystemContext& sys, md_str_t filepath, const loader::LoaderState& load_state, FileRole role, bool opens_own_run = false,
                     const md_system_t* loaded = nullptr, const md_system_state_t* loaded_state = nullptr, int along = -1);

// What follows a system just loaded into sys (load_data), and what goes with it before another replaces it. The
// application's part of it that is not yet per system (the script's evaluation) is of the first system.
void init_system_data(ApplicationState* app, SystemContext& sys);
void free_system_data(ApplicationState* app, SystemContext& sys);

// ## The systems (ApplicationState::systems)

// What a system holds from the start, empty or not: its arena and its masks. Once, before it is first loaded into; and
// what that made freed again (system_context_destroy) once it is unloaded and nothing uses it.
void system_context_init(SystemContext* sys, ApplicationState* app);
void system_context_destroy(SystemContext* sys);

// Whether anything is loaded into it: atoms, or files (which may be missing). An empty system is not listed or drawn.
bool system_loaded(const SystemContext& sys);

// A new, empty system to load into, listed last: the first's place when that is empty, else one allocated with a new
// id. NULL when SYSTEM_MAX are loaded.
SystemContext* system_add(ApplicationState* app);

// A new, empty system allocated with a new id, listed last, whether the first is empty or not (a workspace's systems):
// NULL when SYSTEM_MAX are held
SystemContext* system_alloc(ApplicationState* app);

// Takes the system, which has been unloaded (system_remove), out of the list and frees it; the first only stays empty
void system_release(ApplicationState* app, SystemContext* sys);

// The system with the id, NULL when there is none
SystemContext* system_find(ApplicationState* app, uint32_t id);
const SystemContext* system_find(const ApplicationState& app, uint32_t id);
// The system with the id, the first when there is none: what a tool of a system reads, its system gone
static inline SystemContext& system_or_first(ApplicationState* app, uint32_t id) {
    SystemContext* sys = system_find(app, id);
    return sys ? *sys : app->system_context;
}
static inline const SystemContext& system_or_first(const ApplicationState& app, uint32_t id) {
    const SystemContext* sys = system_find(app, id);
    return sys ? *sys : app.system_context;
}

// How many are loaded (system_loaded)
size_t system_count(const ApplicationState& app);

// THE QUANTUM CHEMISTRY SYSTEM: the one system whose quantum chemistry is studied, one at a time. The Quantum Chemistry
// windows are of it, and its electronic structure representations are the ones evaluated (those of another system are
// not). It is the one chosen (qm_system_choose: the windows' system, a workspace's [QuantumChemistry] System=) while it
// has electronic structure, else the first loaded system that has (qm_system_update, as systems are loaded and freed).
SystemContext& qm_system(ApplicationState& app);
const SystemContext& qm_system(const ApplicationState& app);
// Whether it could be: loaded, with electronic structure (es_source_mask)
bool qm_system_candidate(const SystemContext& sc);
// Where its coordinates are in the world, as drawn this frame: the placement of its first instance, whose coordinates its
// quantum chemistry is of (the electronic structure evaluated, the windows). What is drawn of it in the view (the
// isosurfaces, volumes and planes of its representations, the critical points) is placed by it.
const md_float4x4_t& qm_placement(const ApplicationState& app);
// It is the system with the id from now on: what is under way of the one before is taken first, and the representations
// of both are made again
void qm_system_set(ApplicationState* app, uint32_t id);
void qm_system_choose(ApplicationState* app, uint32_t id);
// After a system was loaded, or before one goes (going): the one chosen while it has electronic structure or is still
// to be loaded, else the first loaded one with electronic structure
void qm_system_update(ApplicationState* app, const SystemContext* going = nullptr);

// THE PALETTE: the colours things are given unless what they are of says otherwise (a type with a colour of its own, a
// representation's): of the elements, of the standard monomers (the amino acids, the nucleotides) and of the secondary
// structures. The user's, as they like to see things: kept as they set it (the .ini, settings.palette), never with a
// workspace, which keeps what overrides it alone. Edited under Visuals > Palette. The types of every system that follow
// it (UseDefaults: System > Atom Types, Residue Types) take its colours.
void palette_builtin(Palette* palette);
// As the .ini keeps it: the entries other than the built in ones, "e:C=909090ff m:HIS=... ss:helix=...". Read onto the
// built in ones; what it does not say is built in.
void palette_write(char* buf, size_t cap, const Palette& palette);
void palette_read(Palette* palette, md_str_t text);
// After an edit of it: its version moved, and what the .ini keeps of it written
void palette_changed(ApplicationState* app);
// The colour of a class of secondary structure: its own, else its group's; that of what is not known for UNKNOWN
md_float4_t palette_ss_color(const Palette& palette, md_secondary_structure_t ss);
// The names they are shown by: "Helix", and "α helix", "3₁₀ helix", "strand"... (UTF-8)
const char* ss_group_label(SecondaryStructureGroup group);
const char* ss_class_label(md_secondary_structure_t ss);

// A system's directory (SystemContext::dir): the outermost layer of the paths of what is of it, as the Attributes window
// lists them ("0/run/md/time", "mutant/atom/charge"). Unique among the systems: a number, the lowest free, given as it is
// made, or a name the user gives it (the System window). Kept with the workspace ([SystemScene] Dir=).
// Why the user cannot give sys the directory, NULL when they can: empty, longer than it holds, a character other than a
// letter, a digit, '_', '-' or '.', all digits (the numbers are given), or another system's
const char* system_dir_invalid(const ApplicationState& app, const SystemContext& sys, md_str_t dir);
// Gives sys the directory the user gives it, false (and nothing changed) where it cannot (system_dir_invalid)
bool system_set_dir(SystemContext& sys, md_str_t dir);
// Gives sys the directory a workspace says it had: a number too. Another system's number gives way, as it was given and
// the workspace says what that one's is; where another has it as a name, sys keeps its own.
void system_dir_take(SystemContext& sys, md_str_t dir);
// The lowest free number, as a system is made or emptied
void system_dir_number(SystemContext& sys);

// The system a hit is of (PickingHit::system), NULL for a hit of none or of a system since removed
SystemContext* picking_hit_system(ApplicationState* app, const PickingHit& hit);
const SystemContext* picking_hit_system(const ApplicationState& app, const PickingHit& hit);
// Every run of the system closed: its files, what was read and derived along them, the series loaded along them
void free_trajectory_data(ApplicationState* app, SystemContext& sys);

// The first instance of the system: the one the timeline, the script, the selection and the tools read
static inline SystemInstance&       system_instance(SystemContext& sys)       { return *sys.instances[0]; }
static inline const SystemInstance& system_instance(const SystemContext& sys) { return *sys.instances[0]; }

// The run an instance shows, NULL for the system's own coordinates
SystemRun*       instance_run(SystemInstance& inst);
const SystemRun* instance_run(const SystemInstance& inst);

// The run the first instance of the system shows: the timeline's, the script's. NULL without one.
SystemRun*       system_run(SystemContext& sys);
const SystemRun* system_run(const SystemContext& sys);

// The instance shows run i (at the time shown, or at frame when it is not negative), or the system's own coordinates for
// -1, from now on
void instance_show(SystemInstance& inst, int run, int64_t frame = -1);

// An instance added to the system, showing its own coordinates; NULL when it has as many as it can (INSTANCE_MAX)
SystemInstance* instance_add(SystemContext& sys);
// Instance i removed, and the representations of it; the last is not
void instance_remove(SystemContext& sys, size_t i);
SystemInstance* instance_find(SystemContext& sys, uint32_t id);
const SystemInstance* instance_find(const SystemContext& sys, uint32_t id);

// The instance with the id, the first where there is none: what a representation or a group naming it is of (0, which
// no instance has, names the first)
static inline SystemInstance& instance_resolve(SystemContext& sys, uint32_t id) {
    SystemInstance* inst = instance_find(sys, id);
    return inst ? *inst : system_instance(sys);
}
static inline const SystemInstance& instance_resolve(const SystemContext& sys, uint32_t id) {
    const SystemInstance* inst = instance_find(sys, id);
    return inst ? *inst : system_instance(sys);
}
// The instance a representation is of
static inline SystemInstance&       representation_instance(SystemContext& sys, const Representation& rep)             { return instance_resolve(sys, rep.instance); }
static inline const SystemInstance& representation_instance(const SystemContext& sys, const Representation& rep) { return instance_resolve(sys, rep.instance); }
// Whether the representation is of the instance (as representation_instance has it)
static inline bool representation_of(const SystemContext& sys, const Representation& rep, const SystemInstance& inst) {
    return &representation_instance(sys, rep) == &inst;
}

// Every instance of the system to make its coordinates again (the system's own changed, the operations), and to upload
// what bits names of them
void system_instances_update(SystemContext& sys, uint32_t dirty_bits = 0);
// Every instance of the system to upload what bits names of it (MolBit_*)
static inline void system_instances_flag(SystemContext& sys, uint32_t dirty_bits) {
    for (size_t k = 0; k < sys.num_instances; ++k) sys.instances[k]->dirty_gpu_buffers |= dirty_bits;
}

// Run i closed: its file and the series along it go out of the file list (file_list_remove_run), and what was read
// and derived along it with them. The instances which showed it go, and their representations, but for the first, which
// shows the system's own coordinates.
void close_run(ApplicationState* app, SystemContext& sys, size_t i);

// A run's frame axis ("<run>/time"): its frame count, its times (resident, F64, never decreasing) and their unit,
// which is none when the file does not know time and the values are frame ordinals. NULL, zero and empty for none.
// Without the run, of the run shown.
const md_attribute_t* run_time_axis(const SystemContext& sys, const SystemRun* run);
size_t                run_num_frames(const SystemContext& sys, const SystemRun* run);
const double*         run_frame_times(const SystemContext& sys, const SystemRun* run);
md_unit_t             run_time_unit(const SystemContext& sys, const SystemRun* run);
static inline const md_attribute_t* run_time_axis(const SystemContext& sys)  { return run_time_axis(sys, system_run(sys)); }
static inline size_t                run_num_frames(const SystemContext& sys) { return run_num_frames(sys, system_run(sys)); }
static inline const double*         run_frame_times(const SystemContext& sys){ return run_frame_times(sys, system_run(sys)); }
static inline md_unit_t             run_time_unit(const SystemContext& sys)  { return run_time_unit(sys, system_run(sys)); }

// "run/<name>" of the run shown, empty without one: what the script and the series are read along
md_str_t run_path(const SystemContext& sys);

// The time at a fractional frame along the run's frame axis, in its unit (frame_axis.h). Zero without a trajectory.
double run_time_at_frame(const SystemContext& sys, double frame);

// THE CLOCK (clock.h): the time shown, ApplicationState::animation.time, which every instance that follows the time shows
// its run at.

// The instance the clock counts in: the first system's first instance when it shows a run, else the first one of the
// systems, in their order, which does. The frames of its run are what stepping, the frame slider and the speed count,
// what the timeline's periods are in, and its unit of time is the clock's. NULL when no instance shows a run.
const SystemInstance* clock_reference(const ApplicationState* app);
// The unit of the clock: of the reference's run, none (frames numbered) without one
md_unit_t clock_unit(const ApplicationState* app);
// A run's frame axis as it is placed on the clock: its times, its offset and the factor into the clock's unit, which is 0
// when they cannot be taken into it (it is placed frame by frame as the reference then)
ClockAxis run_clock_axis(const ApplicationState* app, const SystemContext& sys, const SystemRun* run);
// The offset a run is moved by along the clock, in the unit of its times: kept with its files (RunFiles), in the
// workspace with its [Run]. Setting it, every instance of the run shows its frame of the time again.
double run_time_offset(const SystemContext& sys, const SystemRun* run);
void   run_set_time_offset(SystemContext& sys, const SystemRun* run, double offset);
// The run's offset set so that its first frame is at the reference's first frame on the clock: the two start together.
// False, and nothing changed, for the reference's own run or a run not placed by its times.
bool   run_align_start(SystemContext& sys, const SystemRun* run);
// The fractional frame of a run at a time on the clock, clamped to its frames; and the time on the clock of a frame of it
double run_frame_at_clock(const ApplicationState* app, const SystemContext& sys, const SystemRun* run, double time);
double run_clock_time(const ApplicationState* app, const SystemContext& sys, const SystemRun* run, double frame);
// The fractional frame an instance shows now: the one it keeps to, or its run's at the clock's time. -1 without a run.
double instance_current_frame(const SystemInstance& inst);
// The reference's fractional frame at the clock's time: "the current frame" of whatever reads the first system's run
// (the script, the tools, the periods). 0 without a reference.
double clock_frame(const ApplicationState* app);
// The clock set to the time of a fractional frame of the reference
void   clock_set_frame(ApplicationState* app, double frame);
// The times the clock goes over: from the first frame to the last of every run an instance shows. Both 0 without one.
void   clock_span(const ApplicationState* app, double* beg, double* end);
// The clock moved by frames of the reference (clock_step_time): ctrl and the arrow keys
void   clock_step(ApplicationState* app, double frames);
// The mean time from frame to frame of the reference, on the clock: what one frame per second plays at
double clock_frame_duration(const ApplicationState* app);
// PLAYBACK, of the clock
// The span playing goes over: of the fixed period it loops over (animation.loop_period), within the clock's span, else
// the clock's span
void   playback_range(const ApplicationState* app, double* beg, double* end);
// Plays, or pauses where it is. Played from an end it stopped at (or from outside the span it goes over), from where
// it begins in the direction played.
void   playback_toggle(ApplicationState* app);
// Once a frame, while playing: the clock moved on by dt seconds of fps frames of the reference each, as the loop has
// it at the ends (clock_advance)
void   playback_update(ApplicationState* app, double dt);

// Once a frame: the time kept in the unit of the clock when the reference changes to a run in another (converted, or to
// the clock's first time where it cannot be), and within the clock's span. The view of the timeline and the periods,
// which are on the clock too, are converted with it.
void   clock_update(ApplicationState* app);
// The frames [beg, end] of a run a timeline period covers at the clock's time (timeline_period_frames): those whose
// times on the clock are within it; of a run not placed by its times, the reference's there, as far as it has them.
// False when it covers none of the run.
bool   run_period_frames(int* beg, int* end, const ApplicationState* app, const SystemContext& sys, const SystemRun* run, const TimelinePeriod& period);
// Of the first system's run, which the tools read
static inline bool period_frames(int* beg, int* end, const ApplicationState* app, const TimelinePeriod& period) {
    return run_period_frames(beg, end, app, app->system_context, system_run(app->system_context), period);
}

// One frame of a run (the one shown, without it) into out: coordinates when out has storage for them, and the cell.
// For a single frame now and then; anything that walks many frames keeps an extraction context of its own
// (md_system_extract_begin) so the files stay open across them.
bool extract_frame(const SystemContext& sys, const SystemRun* run, int64_t frame, md_system_state_t* out);
static inline bool extract_frame(const SystemContext& sys, int64_t frame, md_system_state_t* out) { return extract_frame(sys, system_run(sys), frame, out); }

// Infers the system again (in the background, system_inference_start) from the given whole frame of the run, or from
// the coordinates shown when there is no run: the bonds inferred are found anew, those of a file, a topology or a user
// are kept, and what follows from them (rings, structures, bond orders) is of the same coordinates, which become the
// reference of the topology (md_system_t::reference).
bool recompute_bonds(SystemContext& sys, int64_t frame);

// An identifier the script does not use yet: base, or base with a number after it (main.cpp)
md_str_t create_unique_identifier(const script::Program* ir, md_str_t base, md_allocator_i* alloc);

// INFERENCE IN THE BACKGROUND (SystemInference), of each system on its own
// Starts inferring the system as it is now, from state (NULL: the coordinates shown): all of what is derived, put in
// place of what it has once done
void system_inference_start(SystemContext& sys, const md_system_state_t* state);
// Once per frame: starts again what was cancelled (the inference, the runs' backbone angles and secondary structure),
// or was inferred from a topology which changed since
void system_inference_update(SystemContext& sys);
// Whether the system shows its core, what inference derives not there yet
bool system_inference_pending(const SystemContext& sys);
// Waits for it and puts it in place: before what changes the topology, or reads what is derived from it, at once
// (a supplement, a workspace, an edit of the bonds)
void system_inference_finish(SystemContext& sys);

// To be called before the core of the system is edited (the bonds: md_system_bond_insert, md_system_bond_remove):
// what reads it on other threads stops. The script evaluation is interrupted (the script is compiled again once the
// system is inferred again, and starts over), the other tasks of the pool finish, and an inference under way, whose work
// shares the core, is dropped: the system is to be inferred again (system_core_changed, or system_inference_finish).
// What the system shows stays what it was inferred as until then, and so do the indices of its bonds.
void system_core_edit_begin(SystemContext& sys);

// After the core was edited: the system is inferred again in the background, from the coordinates shown, and what is
// drawn of it, the script (whose properties are of the topology) and the representations follow once that is in place.
void system_core_changed(SystemContext& sys);

// The kind of a bond by its word (md_bond_kind_name: "coordinate", "hydrogen"), as a workspace stores it.
// bond_kind_parse reads it back (flags of the kind alone), false for a word which is not of one kind; "interaction", of
// workspaces from before it named every interaction, is a contact.
const char* bond_kind_str(md_bond_flags_t flags);
bool        bond_kind_parse(md_bond_flags_t* out, md_str_t str);
// Origins of bonds, a mask of 1 << md_bond_origin_t, by name: "file,inferred"
const char* bond_origins_str(char* buf, size_t cap, uint32_t origins);
bool        bond_origins_parse(uint32_t* out, md_str_t str);

// Brings the bond changes tracked (operations.bond_tracking) up to date with the frame shown, the topology, the scope
// and the parameters. Once a frame, before the bonds are uploaded.
void bond_tracking_update(SystemContext& sys);

// The bond changes drawn, those of the frame shown: NULL when they are not tracked
const md_bond_delta_t* bond_tracking_delta(const SystemContext& sys);
// The chemistry of the frame shown where its bonds changed: NULL when the changes are not tracked
const md_chem_delta_t* bond_tracking_chem(const SystemContext& sys);

// What an atom is in the frame shown, its chemistry tracked (bond_tracking_chem) or that of the system: the flags,
// formal charge and hydrogens of the atom
md_atom_flags_t frame_atom_flags(const SystemContext& sys, size_t atom_idx, int* formal_charge, int* hydrogen_count);

// The atoms and flags of ring idx of the ring picking domain: a ring of the system, or one of the frame's chemistry
// after them. False when there is none.
bool picked_ring(const SystemContext& sys, size_t idx, const int32_t** atoms, size_t* count, md_ring_flags_t* flags);

// What a picked bond is (picked_bond)
enum class PickedBondSource {
    System,     // A bond of the system
    Formed,     // Formed in the frame shown, the bonds tracked (bond_tracking): not a bond of the system
    Found,      // Found for the frame shown by the filter of a Bonds representation (a hydrogen, ionic or halogen bond, a bond formed): not a bond of the system
    Unbonded,   // No bond: two atoms a representation draws a bond between (a pair its filter gave) which none joins
};

// The bond a hit of the bond picking domain names: of the system, a bond of the graph or one formed in the frame shown,
// which come after them as they are drawn; or a weak bond of a representation, that of the system (or of the frame)
// between its atoms, as of the frame shown. Its atoms, its flags (none when unbonded) and what it is.
bool picked_bond(const ApplicationState& app, const PickingHit& hit, md_atom_pair_t* pair, md_bond_flags_t* flags, PickedBondSource* source);

// Sets the bits of what a hit stands for, grown by the granularity: what hovering it highlights, and what clicking it
// selects (the atoms of a bond, the residue of a cartoon segment, ...)
void picking_hit_mask_set(md_bitmask_t* mask, const SystemContext& sys, const PickingHit& hit, SelectionGranularity granularity);

// "<run>/<leaf>" in a run, e.g. "run/md/backbone/angle" for "backbone/angle". Empty without the run, or when it does
// not fit in buf. Without the run, of the run shown.
md_str_t run_attribute_path(char* buf, size_t cap, const SystemRun* run, md_str_t leaf);
static inline md_str_t run_attribute_path(char* buf, size_t cap, const SystemContext& sys, md_str_t leaf) { return run_attribute_path(buf, cap, system_run(sys), leaf); }
// The attribute at "<run>/<leaf>" of a run; NULL without the run or without such an attribute
const md_attribute_t* run_attribute(const SystemContext& sys, const SystemRun* run, md_str_t leaf);
static inline const md_attribute_t* run_attribute(const SystemContext& sys, md_str_t leaf) { return run_attribute(sys, system_run(sys), leaf); }

// The secondary structure to show for the instance's state, one per segment of sys.protein_backbone: the denoised copy
// of its run at the frame it shows (secondary_structure_render), else the state's own assignment
// (md_system_state_secondary_structure). For presentation; NULL when there is neither.
const md_secondary_structure_t* displayed_secondary_structure(const SystemInstance& inst);

// The instance a hit is of: by its key for an instance's own objects (picking_instance_key), by the representation for
// what is keyed by one (a weak bond, a glyph); the first for anything else
const SystemInstance& picking_hit_instance(const SystemContext& sc, const PickingHit& hit);

// Frame cache operations
void clear_system_frame_cache(SystemContext& sys);

// Interpolate state
void interpolate_system_state(SystemInstance& inst);

// The stride of the frames the displayed state is interpolated between: animation.stride while playing,
// every frame otherwise.
int playback_stride(const ApplicationState* app);

// Workspace
// See #workspace in viamd.cpp for what a workspace holds and the order it is loaded in

// The system a section of a workspace being read is of, by the id the workspace names it by (System=; -1 where the
// section names none, which is of the first system). NULL for one it does not load: the section is left out.
typedef SystemContext* (*WorkspaceSystemFn)(void* user, int id);
// The same for a section a component reads (EventType_ViamdDeserialize): the system of the workspace being read that
// it names by id, NULL for one it does not load or outside a workspace
SystemContext* workspace_system_of(const viamd::deserialization_state_t& state, int id);

void load_workspace(ApplicationState* app, md_str_t file);
// True when the whole file was written; it is then the current workspace (files.workspace)
bool save_workspace(ApplicationState* app, md_str_t file);

// Whether this window is open belongs to the workspace: stored under [Windows] as name=0/1. name
// is the window's identity in the file, so it stays the same when the window is renamed on screen.
// Call once, at initialization; show has to outlive the application.
void workspace_register_window(const char* name, bool* show);
// The flag that shows the window registered by that name, NULL for none: to open one from another
bool* workspace_window_flag(const char* name);

// Selections
// What filters are evaluated and completed against: the system as it is shown, its run (which attr() reads), and the
// names of the script in the editor. The script is compiled against the first system (system_context) alone until it
// is per system, so its names are of that one: another system's filters see none of them.
script::FilterCtx filter_context(const SystemContext& sys);
// The same of an instance: its coordinates and its run, as a tool of it evaluates its filters in. The first instance's
// of a system is filter_context of the system.
script::FilterCtx filter_context(const SystemInstance& inst);

// The user's groups (MD_GROUP_ORIGIN_USER), what were stored selections: named sets of atoms of the system, kept in
// the workspace, which VIAL selects as group('name') and v1 names as identifiers. The System window lists, renames and
// removes them. A name is the user's handle on a set, so no two groups of the system share one of a user's.
//
// The system's selection as a new group of the user's, named sel1, sel2, ...: its index, or -1 for an empty selection
int store_selection_as_group(SystemContext& sys);
// Whether a user's group may be named so: not empty, and no other group of the system's name (self: the group renamed)
bool user_group_name_free(const md_system_t& sys, md_str_t name, size_t self = SIZE_MAX);

// Representations and their groups: representation.h

// Enumerates the complete dipole GROUPS in the system's attribute table, in name order - a handful
// of entries, and no values are read. Returns the total number found and writes at most cap, so
// pass cap 0 to count. Cheap enough to call per frame: it is a prefix query over a small sorted
// array plus two lookups per group, and there is nothing to cache or invalidate because the
// attribute table IS the state.
// Groups missing either half, or whose shapes are not one 3-residue value, are skipped.
size_t dipole_groups_gather(DipoleGroup out[], size_t cap, const md_system_t& sys);

// The same description for ONE group, straight from its vector attribute - for a representation or
// a picking hit which already holds the key and has no reason to walk the list. False when the key
// names nothing, which is what a workspace pointed at a file that no longer has that group looks like.
bool dipole_group_from_key(DipoleGroup* out, const md_system_t& sys, md_attribute_id_t key);

// One moment out of a group: the vector as stored and the origin in Angstrom, which extraction
// ENFORCES is a length. False when the key names nothing, the group has lost its origin, or the
// index is past the end. This is the whole read path - drawing and the picking tooltip both take
// the (key, index) they already hold and come straight here.
bool dipole_moment_read(md_float3_t* out_vec, md_float3_t* out_origin, const md_system_t& sys, md_attribute_id_t key, uint32_t index);

// Every density property the system carries, in the attribute table's own path order. Returns how
// many there are, which may exceed 'cap' - only the first 'cap' are written. Pass (NULL, 0) to count.
size_t density_properties_gather(DensityProperty out[], size_t cap, const md_system_t& sys);

// Which electronic structure sources this system can actually show. Asked of the attribute table on
// the spot rather than cached: each source needs particular attributes, so their presence IS the
// answer, and a second reader publishing them satisfies it without anyone being told.
ElectronicStructureSourceFlags es_source_mask(const md_system_t& sys);

// Points an electronic structure representation at something this system carries. A source the
// system has is left as it is; otherwise the first one it has, in the enum's order, is taken with that
// source's defaults - so orbitals stay the first choice wherever there are any, and a file with only
// density properties gets a density property rather than an orbital source it cannot evaluate. A
// density property source is also pointed at an actual property, as it otherwise evaluates nothing
// until the representation's window has been drawn. Returns false when the system carries no
// electronic structure at all, and leaves the representation untouched.
bool electronic_structure_select_available_source(ElectronicStructureRepresentation* es, const md_system_t& sys);

// Takes what the field a representation is coloured by has finished evaluating (surface_field_poll):
// once a frame, whatever else updates
void electronic_structure_field_poll(Representation* rep);

// Brings what an electronic structure representation draws up to date: its volume, evaluated again only
// when what it is evaluated from changed (the frame among it), and the colours of its surfaces - from
// atom_colors, the colours of the atoms as the representation has them, or from the field it is coloured by.
// False when something it waits for without blocking (the volume's last readback, a density matrix made on a
// worker thread, the density a field update reads) is not there yet: it is to be called again next frame.
// The planes display evaluated on its planes (on the GPU): each slice's image, made again where
// the plane, the source or the frame changed (not while the plane is dragged), one at a time per plane. A plane whose
// image is not current is drawn from the volume, which this display keeps too. Called every frame; false while one is
// under way.
bool electronic_structure_plane_images_update(ApplicationState* state, Representation* rep);
// For measurements (the timings shown with it): bumped with each change to how the images are evaluated
#define ES_PLANE_IMAGES_REVISION "plane-images r5"
// Whether it is: the planes display, and the source on the GPU
bool es_plane_images_supported(ApplicationState* state, const ElectronicStructureRepresentation& es);
// Whether the representation's isosurfaces are drawn as meshes: its setting, for what meshes can be made of
bool es_mesh_mode(const ApplicationState* state, const ElectronicStructureRepresentation& es);

bool electronic_structure_representation_update(ApplicationState* state, Representation* rep, const uint32_t* atom_colors);



// Turns a dipole group name into something presentable: "ground_state" -> "Ground State".
// Writes at most cap-1 characters plus a terminator, returns the length written.
int dipole_label_pretty(char* buf, size_t cap, md_str_t group);

// The same, plus a 1 based element number when the group holds more than one moment:
// "Electric Transition 3". A single moment gets no number, because there is nothing to tell apart.
int dipole_entry_label(char* buf, size_t cap, const DipoleGroup& group, uint32_t index);

// Per atom scalar fields are not a list anybody builds or keeps: they are whatever the system's
// attribute table holds under atom/, and these read it directly. Nothing to rebuild when the data
// changes, and nothing that can go stale against it.

// Ids of the attributes under atom/ which are per atom scalar fields, in path order.
// Returns the total number found and writes at most cap, so pass cap 0 to count.
// An attribute qualifies when its values are single residues and its LAST index axis is this
// system's atom count; an optional leading axis is the variant axis. Everything else under atom/ -
// a position, a velocity, anything several residues wide - is data this colouring cannot express
// and is skipped rather than mangled. Cheap enough to call per frame.
size_t atom_attribute_query(md_attribute_id_t out_ids[], size_t cap, const md_system_t& sys);

// What to show for it: the attribute's label, or its leaf path segment when it has none. A view
// into the table's storage, null terminated, so it can go straight to ImGui.
md_str_t atom_attribute_label(const md_attribute_t* attr);

// Number of variants: the leading index axis, 1 when the field has none.
int atom_attribute_variant_count(const md_attribute_t* attr);

// Span of the values, over EVERY variant so a colour ramp does not jump as the variant changes, and
// over the atoms of 'mask' (all of them when NULL). attr is read from table, the one it was found in.
// Derived rather than stored: it belongs to whoever is drawing the ramp, not to the table. Scans the
// whole attribute, so not something to call per frame without a reason. Atoms without a value play
// no part; invalid when there are none with one.
ColorScaleSpan atom_attribute_span(const md_attributes_t* table, const md_attribute_t* attr, const md_bitset_t* mask);

// What a legend calls it: the label, and which variant when there are several ("Charge [2/5]")
void atom_attribute_legend_label(char* buf, size_t cap, const md_attribute_t* attr, int variant_idx);

// Whether an atom has NO value in a per atom field. A producer marks that with NAN
// (md_attributes_publish_atom_column): the QM atoms in an embedding's charges, the embedding's sites
// in a column of the QM calculation, a blank in an mmCIF column. Tested on the bits, because under
// fast math every float spelling of the test (v != v, isnan) may be folded away.
static inline bool atom_attribute_value_absent(float v) {
    uint32_t u;
    md_memcpy(&u, &v, sizeof(u));
    return (u & 0x7fffffffu) > 0x7f800000u;
}

// Builds the per dataset GPU data - the uploaded GTO basis and the atom buffer - from the system's
// own basis/ attributes, and grows the device coefficient scratch to fit. Returns false when the
// system publishes no basis, which is the normal case rather than a failure.
bool system_gpu_data_update(SystemContext& sys, double cutoff);
void system_gpu_data_free(SystemContext& sys);

// One orbital's AO coefficient row, out of the attribute at coefficient_path narrowed by 'slice',
// in the double precision md_gto asks for. Allocated from temp; NULL when the path is absent or the
// slice does not narrow it to a single row. out_num_ao receives the row length.
//
// md_attribute_slice_1(m) for a {M,A} matrix of molecular orbitals, md_attribute_slice_2(s,l) for
// the {S,L,A} natural transition orbitals - the caller says which row, and nothing here has to know
// which of the two it was addressing.
double* orbital_coefficients_extract(size_t* out_num_ao, md_temp_scope_t temp, const md_system_t& sys, md_str_t coefficient_path, const md_attribute_slice_t* slice);

// Positions for a system's BASIS atoms, in the bohr and interleaved layout md_gto works in, taken
// out of `state` through the QM -> system map. Allocated from temp; null with a zero count when the
// system publishes no basis or the state cannot serve it.
//
// This is what an evaluation takes INSTEAD of a state. Which geometry an orbital is drawn at - the
// current frame, a reference, an optimisation step, a geometry that never was - is the call site's
// decision, and an evaluation is in no position to make it.
md_float3_t* basis_atom_positions_extract(size_t* out_count, md_temp_scope_t temp, const md_system_t& sys, const md_system_state_t& state);

// Blocks until every queued volume readback has landed in its GL texture. Call before destroying
// the device scratch, the pools or a GL texture one of them writes to, and wherever a caller needs
// the texture's contents NOW rather than next frame - an export, say. A no-op without the GPU
// backend, where evaluation writes the texture directly.
void gpu_volume_jobs_drain(ApplicationState* state);

// Issues what the density jobs have room for (their pieces, a few in flight), and queues the readbacks of those
// all issued: once a frame. Each evaluation of a volume on the GPU (a density, ELF or LOL) is listed in the overlay
// while it runs, as a task_system activity, from its first piece until its readback lands.
void gpu_volume_jobs_poll(ApplicationState* state);

// Drops the density evaluation into vol_tex that is not yet all issued, if there is one: what it has in flight
// completes, nothing more is issued, and the texture keeps what it held. True when one was dropped.
bool gpu_volume_job_cancel(ApplicationState* state, uint32_t vol_tex);

// A line of the task overlay from a component, given in answer to EventType_ViamdTaskProgressRequest: what is
// under way, how far it is (0 to 1), and optionally a flag the overlay's close button raises to stop it. Read in
// the frame it is given in, so the flag only has to live that long. This is how a component shows its work
// there: main.cpp asks rather than calls, as a component may not be built (the QM module is a build option).
struct TaskProgress {
    const char*       label    = "";
    float             fraction = 0.0f;
    volatile int32_t* cancel   = nullptr;
};

struct TaskProgressRequest {
    TaskProgress items[8];
    size_t       count = 0;
};

static inline void task_progress_add(TaskProgressRequest* req, const char* label, float fraction, volatile int32_t* cancel) {
    if (req->count < MD_COUNTOF(req->items)) {
        req->items[req->count++] = TaskProgress{ label, MD_CLAMP(fraction, 0.0f, 1.0f), cancel };
    }
}


// The QM ATOM DOMAIN: the atoms a quantum chemistry calculation covered, in its own order and at
// its own geometry. NOT the system's atoms - a calculation can cover part of a loaded system, so
// the two spaces differ in length and in order, and system_index is the only bridge.
//
// Published under qm/atom/. basis/shell/atom_index is an index INTO this space - the basis is one
// thing defined over these atoms, not the other way round - and 'count' is what "how many atoms
// does the QM cover" means; there is nothing else to ask.
//
// GATHER IT WHERE YOU USE IT; DO NOT STORE IT. The pointers are the attribute table's own storage,
// and md_attributes_replace frees the old buffer before writing the new one - so anything that
// republishes qm/atom/*, a supplemental QM file loaded onto an already loaded system above all,
// leaves a stored copy pointing at freed memory. A cache would have to be invalidated on a signal
// that does not exist: the supplement need not touch the orbital coefficients, so their version
// does not move, and nothing else announces that the domain was rewritten.
//
// The cost of not caching is three binary searches over a sorted array (md_attributes_find), which
// is nothing against what a caller does with the result. Gather once per function and hoist it out
// of your own loop, the way every caller already did with the stored copy.
struct QmAtoms {
    size_t          count         = 0;
    const uint8_t*  atomic_number = nullptr;  // [N]
    const md_double3_t*  coordinate    = nullptr;  // [N], Angstrom, the geometry the calculation ran at
    const uint32_t* system_index  = nullptr;  // [N], or null for the identity
};

// False when the system carries no QM atom domain at all. An individual member may still be null
// when that column was not published; a null system_index is not a failure, it is the standalone
// case and reads as the identity.
bool es_qm_atoms(QmAtoms* out, const md_system_t& sys);

// Which system atom QM atom i is. The one place the two index spaces meet.
inline size_t qm_to_system_atom(const QmAtoms& qm, size_t i) {
    return qm.system_index ? (size_t)qm.system_index[i] : i;
}

// The extent of the MO coefficient matrix: M molecular orbitals over A atomic orbitals, read from
// the attribute's own SHAPE. Deliberately not two more published scalars - a matrix cannot disagree
// with its own dimensions, and a reader that publishes coefficients has already said this. Either
// out pointer may be NULL. False when no coefficients are present, which is "this system is not a
// QM calculation" rather than a failure.
bool es_orbital_extent(const md_system_t& sys, size_t* out_num_mo, size_t* out_num_ao);

// The frontier orbitals of one spin channel, from the occupations at occupation_path
// (es_path::alpha_occupation or its beta sibling). DERIVED and not published: it is one scan of a
// column that is already in the table, and a stored copy is a second thing to keep in agreement.
//
// HOMO is the highest index carrying occupation, LUMO the next one up. Both are -1 when the column
// is absent; lumo_idx is num_orbitals when every orbital is occupied, which is out of range on
// purpose - there is no LUMO to point at. NOTE this differs from md_vlx's own homo_idx/lumo_idx,
// which take the FIRST zero occupancy and leave both at 0 when there is none; the two agree on
// every aufbau-ordered canonical set, and this one also answers for the fully occupied case.
// Natural orbitals, whose occupations are fractional and never exactly zero, have no frontier in
// this sense and will report the last orbital.
struct OrbitalFrontier {
    int homo_idx      = -1;
    int lumo_idx      = -1;
    int num_orbitals  = 0;
};
bool es_orbital_frontier(OrbitalFrontier* out, const md_system_t& sys, md_str_t occupation_path);

// The orbital energies of one spin channel (es_path::alpha_energy or its beta sibling), in Hartree,
// or NULL when that column is absent. Hands back the stored buffer rather than a copy: the column is
// resident F64 and every caller only reads it. out_count receives its length, which is the orbital
// count the column itself states - compare it against es_orbital_extent's before indexing the two
// together rather than assuming they agree.
const double* es_orbital_energies(size_t* out_count, const md_system_t& sys, md_str_t energy_path);

// True when the system holds a SECOND, DISTINCT set of orbital coefficients - an unrestricted
// calculation. Restricted and restricted-openshell publish beta as an alias of alpha: one datum
// under two names, so a spin selector over them would offer a choice that changes nothing.
bool es_has_distinct_beta_orbitals(const md_system_t& sys);

// True when the two spin channels carry DIFFERENT occupations. Together with the above this is the
// whole of the SCF type distinction, taken from the data rather than from a reader's enum:
//
//   restricted            same coefficients, same occupations
//   restricted open shell same coefficients, different occupations   (the singly occupied case)
//   unrestricted          different coefficients
//
// Which is the point - a second QM reader publishing the same columns says the same thing about
// itself without anyone having to map its own enum onto md_vlx's.
bool es_has_distinct_beta_occupations(const md_system_t& sys);

// Electrons in one spin channel: the sum of its occupations, which is what an occupation IS. Not a
// published number, because a published number is a second thing that can disagree with the column
// it summarises. Returns 0.0 when the column is absent.
//
// A double rather than a count: natural orbitals occupy fractionally, and rounding here would hide
// that from a caller who cares. Canonical SCF sums to an integer within float noise.
double es_electron_count(const md_system_t& sys, md_str_t occupation_path);

// Excited states carrying natural transition orbitals, from the weight table's leading axis.
size_t es_excited_state_count(const md_system_t& sys);

// One excited state's NTO weights, largest first, stopping at the first below 'cutoff'. The ragged
// lambda axis is zero padded, so the cutoff ends the run in exactly the place the pairs do.
// Returns how many were written.
size_t es_nto_lambdas(double* out_values, size_t cap, const md_system_t& sys, size_t state_idx, double cutoff);

// Whether one excited state's NTOs can be read without waiting. A reader may compute them the first time they
// are read (VeloxChem's, from the response vectors: a tenth of a second for a few hundred occupied orbitals,
// seconds for a thousand) and keep them after, so what the frame reads - the weights, the volume, the flow
// diagram - asks this first: false while they are being read on a worker thread, which this begins when they
// have not been. What waits on it keeps showing what it had. True when there is nothing to wait for (resident
// or absent).
bool es_nto_ready(ApplicationState* state, size_t state_idx);
// What es_nto_ready keeps, freed: with the system, after its pool tasks have ended
void es_nto_warm_free(SystemContext& sys);

// The application's density cache (density_cache.h), made when first asked for. Main thread.
struct DensityCache* density_cache_of(ApplicationState* state);

// One AO x AO density matrix, out of the attribute at density_path narrowed by 'slice', whole and in
// double precision. Allocated from temp; NULL when the path is absent or the slice does not narrow
// it to one symmetric matrix. Either storage is read: square ({A,A}, or {S,A,A} with
// md_attribute_slice_1(state)), or packed (MD_ATTRIBUTE_FLAG_PACKED_SYMMETRIC, as the density
// properties are). out_dim receives A.
double* density_matrix_extract(size_t* out_dim, md_temp_scope_t temp, const md_system_t& sys, md_str_t density_path, const md_attribute_slice_t* slice);

// ---------------------------------------------------------------------------
// The attribute paths an electronic structure representation draws from
// ---------------------------------------------------------------------------
// Named once here rather than spelled at each use, so the whole set is visible in one place.
//
// The 'orbital/' paths are FORMAT NEUTRAL: md_gto's own normalised representation, which any QM
// reader fills with the same meaning. The 'vlx/rsp/' two are still VeloxChem's own - they are what
// is waiting for a neutral name, and the moment a second reader publishes something comparable is
// the moment to give them one and alias these onto it.
namespace es_path {
    inline const md_str_t alpha_coefficient   = MD_STR_INIT("orbital/alpha/coefficient");
    inline const md_str_t beta_coefficient    = MD_STR_INIT("orbital/beta/coefficient");

    inline const md_str_t alpha_energy        = MD_STR_INIT("orbital/alpha/energy");
    inline const md_str_t beta_energy         = MD_STR_INIT("orbital/beta/energy");

    inline const md_str_t alpha_occupation    = MD_STR_INIT("orbital/alpha/occupation");
    inline const md_str_t beta_occupation     = MD_STR_INIT("orbital/beta/occupation");

    // The AO metric the coefficients and densities above are expressed against. Singular for a
    // basis that was stored spherically - see the note on md_vlx_publish_attributes.
    inline const md_str_t overlap             = MD_STR_INIT("basis/overlap");

    inline const md_str_t alpha_density       = MD_STR_INIT("orbital/alpha/density");
    inline const md_str_t beta_density        = MD_STR_INIT("orbital/beta/density");
    inline const md_str_t total_density       = MD_STR_INIT("orbital/total/density");
    inline const md_str_t difference_density  = MD_STR_INIT("orbital/difference/density");

    inline const md_str_t nto_particle        = MD_STR_INIT("vlx/rsp/nto/particle/coefficient");
    inline const md_str_t nto_hole            = MD_STR_INIT("vlx/rsp/nto/hole/coefficient");

    // {S,Lmax}, the weight of each NTO pair, padded with zeros on the ragged lambda axis and
    // indexing the same space as the two coefficient attributes above.
    inline const md_str_t nto_lambda          = MD_STR_INIT("vlx/rsp/nto/lambda");

    inline const md_str_t attachment_density  = MD_STR_INIT("vlx/rsp/transition_density/attachment");
    inline const md_str_t detachment_density  = MD_STR_INIT("vlx/rsp/transition_density/detachment");
    inline const md_str_t transition_diff     = MD_STR_INIT("vlx/rsp/transition_density/difference");

    inline const md_str_t density_property    = MD_STR_INIT("vlx/density_property");
}

// The only interpretation of use_magnitude there is. ABS is a MODIFIER bit on an accumulate op, not
// an op of its own, which is why it is or'd onto SET rather than replacing it.
inline md_gto_op_t es_gto_op(bool use_magnitude) {
    return MD_GTO_OP_SET | (use_magnitude ? MD_GTO_OP_ABS : 0);
}

// Which density path a spin selection names. The four are siblings in the table precisely so that
// this is a lookup and not four cases of arithmetic at the point of use.
inline md_str_t es_electron_density_path(ElectronicStructureSpin spin) {
    switch (spin) {
    case ElectronicStructureSpin::Alpha:      return es_path::alpha_density;
    case ElectronicStructureSpin::Beta:       return es_path::beta_density;
    case ElectronicStructureSpin::Difference: return es_path::difference_density;
    case ElectronicStructureSpin::None:
    case ElectronicStructureSpin::Total:
    default:                                  return es_path::total_density;
    }
}

// The grid an electronic structure representation is evaluated on: an object aligned box around the
// basis atoms, padded, at samples_per_unit_length samples per BOHR. False when the system publishes
// no basis, which is the normal "nothing to draw" case rather than a failure.
bool electronic_structure_grid_init(md_grid_t* grid, const md_system_t& sys, const md_system_state_t& state, double samples_per_unit_length);

// The entry points a representation calls: name the attribute, get a filled 3D texture, evaluated on the GPU
// (md_gpu) at the system's current state. One function per kind of thing being evaluated. coefficient_path is a rank
// 2 {M,A} attribute of doubles, "orbital/alpha/coefficient" or its beta sibling, and slice selects the row;
// density_path and slice name a density matrix - "orbital/total/density" with no slice,
// "vlx/rsp/transition_density/attachment" with the excited state index. False, and said once, without a GPU device.
//
// The texture is filled from the frame loop rather than before the call returns; gpu_volume_jobs_drain() forces it
// where the contents are needed immediately.
bool orbital_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, md_str_t coefficient_path,
                      const md_attribute_slice_t* slice, md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff);
// The same for an orbital the caller made itself (a block NTO, a combination of orbitals): its AO coefficients,
// num_ao of them, which must match the system's basis. Copied before the call returns.
bool orbital_coefficients_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, const double* ao_coeffs, size_t num_ao,
                                   md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff);
// async (the representations): on the GPU, the density is issued as a job (md_gto_gpu_job_t) the frame loop
// completes (gpu_volume_jobs_poll), where the launch would wait for its counting pass
bool density_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, md_str_t density_path,
                      const md_attribute_slice_t* slice, md_gto_op_t op, bool async = false);

// density_evaluate over a matrix the caller already has in hand rather than a path to extract it from. The
// transition densities are VIRTUAL and are rebuilt from the response solution on every read, so a caller which needs
// one twice holds on to it and calls this instead.
bool density_matrix_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid,
                             const double* density_matrix, size_t dim, md_gto_op_t op);

// Whether the evaluations above can run: the GPU is there, and the main thread only queues them. False when they
// cannot, and what would be evaluated ahead of being asked for should not be asked for.
bool gpu_volume_evaluation_available(const ApplicationState* state);

// The analytic normals of a density volume's texels (RGB10_A2: ElectronicStructureRepresentation::normal_tex), where its
// last evaluation made them along with it; 0 where it did not (the texels came without: from the CPU, of ELF, or before)
uint32_t es_normals_of(uint32_t vol_tex);

#if MD_ENABLE_GPU
// The first half of density_evaluate: evaluates into the DEVICE scratch volume (ApplicationState::
// gpu_volume) and leaves it there, for a consumer which records another kernel over it on the same
// stream instead of reading it back. False when the GPU backend is unavailable, in which case there
// is no scratch volume to run anything over.
bool density_evaluate_to_gpu_volume(ApplicationState* state, const md_grid_t& grid, md_str_t density_path,
                                    const md_attribute_slice_t* slice, md_gto_op_t op);

// The same, over a matrix the caller already holds - see density_matrix_evaluate.
bool density_matrix_evaluate_to_gpu_volume(ApplicationState* state, const md_grid_t& grid,
                                           const double* density_matrix, size_t dim, md_gto_op_t op);
#endif

// Points a representation at an attribute, and starts its colour scale out as fits the values: a
// symmetric range on a diverging map for values on both sides of zero (a charge), the range of the
// values on a sequential map otherwise. A colour map the user chose is kept; only the two defaults
// trade places. The range follows the values from there unless the user takes it over.
void atom_attribute_select(AtomAttributeColoring* coloring, md_attribute_id_t key, const md_system_t& sys);

// Recentering operations (low level)

void recenter_mark_query_dirty(SystemContext& sys);
const md_bitset_t& recenter_get_active_target_mask(const SystemContext& sys);
bool recenter_update_query_mask(SystemContext& sys);
void recenter_update(SystemContext& sys);

// Builds the reference that keeping the orientation fits against, from the first frame of the run, when
// the target's atoms differ from the ones it was built for. Does nothing while the orientation is not kept.
void recenter_update_target_data(SystemContext& sys);

// The recentering transform T = rotation * translation, split where the periodic images have to be
// settled (see apply_state_operations):
//     translation = translate(cell_centre - com)                                  lattice preserving
//     rotation    = translate(cell_centre) * alignment_mat * R * translate(-cell_centre)
// rotation is exactly the identity unless the orientation is kept or alignment_mat is set, and the
// return value says whether it is not.
bool recenter_calculate_transform(md_float4x4_t* translation, md_float4x4_t* rotation, const SystemInstance& inst);

// Applies the coordinate operations to mold.state in the one order that keeps them valid together:
//     take off operations.state_rotation -> translate -> wrap -> make whole -> turn
// Wrapping and making whole are only defined in the lattice frame of the cell, and turning before the
// periodic images are settled turns the arbitrary image each atom arrived in into a real displacement
// of R times a lattice vector. So the turn comes last, and a fresh turn is always preceded by a wrap.
// Returns true if the coordinates were touched.
bool apply_state_operations(SystemInstance& inst, bool recenter, bool pbc, bool unwrap);

// Picking

void picking_handler_new_frame(PickingHandler* handler);
PickingSpace* picking_handler_current_space(PickingHandler* handler);
const PickingSpace* picking_handler_find_space(const PickingHandler& handler, uint32_t submitted_frame_idx);

// Reserves a range within the picking space for a specific domain (atoms, bonds, etc).
// Returns true if the range was successfully reserved, false if there was not enough space. If successful, out_range will be filled with the reserved range.
// out_range is optional and can be null if the caller does not need the details of the reserved range (e.g. just needs to know if the reservation was successful or not).
// Reserves 'count' consecutive indices for a domain. 'key' names what they index for a domain whose
// objects are not the system's own arrays; see PickingRange. 'system' is the id of the system whose objects they are,
// which a hit in the range is of. Fails once the space is full.
bool picking_range_reserve(PickingRange* out_range, PickingSpace* space, PickingDomainID domain, size_t count, uint64_t key = 0, uint32_t system = SYSTEM_ID_NONE);

// The range a keyed object of the system was given this frame, or NULL when it was not reserved one.
const PickingRange* picking_space_find_range(const PickingSpace& space, PickingDomainID domain, uint64_t key, uint32_t system = SYSTEM_ID_NONE);

// The key the ranges of an instance's own objects (its atoms, bonds, ...) are reserved under: 0 for the first, as
// before there were several, and past the ids of the representations for the others, whose weak bonds are keyed by
// those ids in the same domain
static inline uint64_t picking_instance_key(size_t i) { return i == 0 ? 0 : ((1ull << 32) | (uint64_t)i); }
// The place in its system of the instance a hit on its own objects is of, -1 when the key names none
static inline int picking_key_instance(uint64_t key) { return key == 0 ? 0 : (key >> 32) == 1 ? (int)(key & 0xFFFFFFFF) : -1; }

// The picking offsets for drawing the system with md_gl, one per reserved range. Every md_gl_draw of the system takes
// them from here: a base left at 0 makes that kind write into the atom range and pick as the wrong atoms.
static inline md_gl_picking_offset_t gl_picking_offset(const SystemInstance& inst) {
    return {
        .atom_base     = inst.picking.atom.beg,
        .bond_base     = inst.picking.bond.beg,
        .backbone_base = inst.picking.backbone.beg,
        .ring_base     = inst.picking.ring.beg,
        .gap_base      = inst.picking.gap.beg,
    };
}

// The picking offsets for drawing the glyphs of the cartoon of the lipids and nucleic acids (CartoonGlyphs), a molecule of
// their own: their spheres and capsules write the index of a point, their tubes that of a control point after the points,
// their slabs their index after those (cartoon_glyphs_picked_point)
static inline md_gl_picking_offset_t gl_picking_offset_glyphs(const SystemInstance& inst) {
    return {
        .atom_base     = inst.picking.glyph.beg,
        .bond_base     = inst.picking.glyph.beg,
        .backbone_base = inst.picking.glyph.beg + inst.cartoon_glyphs.num_points,
        .ring_base     = inst.picking.glyph.beg + inst.cartoon_glyphs.num_points + inst.cartoon_glyphs.num_tube_points,
    };
}

// The marks of the selection and the highlight as post-processing puts them over the image (Settings::marks, from
// selection.color). highlight: how much of the highlight shows; outline_step: the outline's width, about half of it
static inline void marks_settings(postprocess_pipeline::Settings* settings, const ApplicationState& state, float highlight, int outline_step) {
    const auto& c = state.selection.color;
    settings->marks.visible[0] = c.selection.visible;
    settings->marks.hidden[0]  = c.selection.hidden;
    settings->marks.outline[0] = c.selection.outline;
    settings->marks.visible[1] = c.highlight.visible;
    settings->marks.hidden[1]  = c.highlight.hidden;
    settings->marks.outline[1] = c.highlight.outline;
    static_assert(sizeof(settings->marks.groups) == sizeof(c.highlight_groups), "A colour for each group after the first");
    for (int i = 0; i < HIGHLIGHT_GROUP_COLORS - 1; ++i) {
        settings->marks.groups[i] = c.highlight_groups[i];
    }
    settings->marks.highlight = highlight;
    settings->marks.outline_step = outline_step;
}

void picking_surface_init(PickingSurface* surface, PickingSourceID source);
void picking_surface_free(PickingSurface* surface);

// Submits a picking readback request for the given surface. The readback will be performed asynchronously, and the result can be polled using picking_surface_poll_hit.
bool picking_surface_submit_readback(
    PickingSurface* surface,
    uint32_t fbo,
    uint32_t width,
    uint32_t height,
    uint32_t submitted_frame_idx,
    md_float2_t surface_coord,
    md_float2_t screen_coord,
    const md_float4x4_t& inv_mvp
);

// Polls for the result of a picking readback request. If a hit is detected, out_hit will be filled with the details of the hit and the function will return true.
// If no hit is detected the function will return false.
bool picking_surface_poll_hit(
    PickingHit* out_hit,
    PickingSurface* surface,
    const PickingHandler& handler
);

// Convenience wrapper for the common frame loop path.
// Preserves the existing pipeline ordering by submitting the current frame readback and then polling the previous completed one.
bool picking_surface_submit_readback_and_poll_hit(
    PickingHit* out_hit,
    PickingSurface* surface,
    const PickingHandler& handler,
    const PickingReadbackRequest& request
);

enum InteractionSurfaceFlags : uint32_t {
    InteractionSurfaceFlags_None = 0,
    InteractionSurfaceFlags_NoRegionSelect  = 1 << 0,
};

// Creates an invisible interactive surface which forms the basis for picking and interaction.
InteractionSurfaceState interaction_surface(InteractionSurfaceID id, const md_float2_t& size, InteractionSurfaceFlags flags = InteractionSurfaceFlags_None);

struct InteractionSurfaceHitArgs {
    PickingSurface* picking_surface;
    const PickingHandler& picking_handler;
    uint32_t fbo;
    uint32_t width;
    uint32_t height;
    const md_float4x4_t& clip_to_world;
    float pixel_scale = 1.0f;   // the fbo's pixels per framebuffer pixel: the render scale, where it is upscaled
};

bool interaction_surface_hit_extract(PickingHit* out_hit, const InteractionSurfaceState& state, const InteractionSurfaceHitArgs& args);

struct InteractionSurfaceViewTransformArgs {
    const Camera& camera;
    const TrackballControllerParam& trackball_param = {};
};

struct InteractionSurfaceViewTransformResult {
    bool reset_requested = false;
};

// Uses the interaction surface state (e.g. mouse position, region selection) to calculate a view transform based on the provided camera and trackball parameters.
// Modifies the target view transform in place.
// Returns interaction outcomes which the caller may use to apply domain-specific behavior such as view reset.
InteractionSurfaceViewTransformResult interaction_surface_view_transform_apply(ViewTransform* target, const InteractionSurfaceState& state, const InteractionSurfaceViewTransformArgs& args);

enum class InteractionSurfaceEventKind {
    None,
    Hover,
    Click,
    ContextMenu,
    RegionSelect,
};

enum class InteractionSurfaceEventPhase {
    None,
    Update,
    Commit,
};

struct InteractionSurfaceEvent {
    InteractionSurfaceID surface_id = 0;
    ImGuiID item_id = 0;

    InteractionSurfaceEventKind kind = InteractionSurfaceEventKind::None;
    InteractionSelectionMode selection_mode = InteractionSelectionMode::None;
    InteractionSurfaceEventPhase region_phase = InteractionSurfaceEventPhase::None;

    md_float2_t mouse_local = {};
    md_float2_t surface_size = {};
    md_float2_t region_min = {};
    md_float2_t region_max = {};

    md_float4x4_t clip_to_world = MD_MAT4_IDENT_INIT;
    md_float4x4_t world_to_clip = MD_MAT4_IDENT_INIT;

    PickingHit hit = {};
};

void interaction_surface_event_extract(InteractionSurfaceEvent* out_event, const InteractionSurfaceState& state, const PickingHit& hit = {});

// Helper function for projecting world coordinates to surface coordinates, used for interaction surfaces and picking.
static inline md_float2_t world_to_surface_project(
    const md_float3_t& world_coord,
    const md_float4x4_t& world_to_clip,
    const md_float2_t& surface_size
) {
    const md_float4_t c = world_to_clip * md_float4_set(world_coord.x, world_coord.y, world_coord.z, 1.0f);
    md_float2_t out_coord = {
        ( c.x / c.w * 0.5f + 0.5f) * surface_size.x,
        (-c.y / c.w * 0.5f + 0.5f) * surface_size.y,
    };
    return out_coord;
}

// Sets mask to the points (of candidate_mask, all of them when NULL) which project into the region
void point_set_region_mask_compute(
    md_bitmask_t* mask,
    const md_float3_t xyz[],
    size_t count,
    const md_bitset_t* candidate_mask,
    const md_float4x4_t& world_to_clip,
    const md_float2_t& region_min,
    const md_float2_t& region_max,
    const md_float2_t& surface_size
);

// File Queue
bool file_queue_empty(const FileQueue* queue);
bool file_queue_full(const FileQueue* queue);
// Queues a file to be loaded: into the system with the id 'system', or where file_queue_process puts it
void file_queue_push(FileQueue* queue, md_str_t path, FileFlags flags = FileFlags_None, uint32_t system = SYSTEM_ID_NONE);

// Loads a file the user opened (the load dialog) into the system target, or into a new one when target is NULL
// (load_data); and what follows (load_into_done). A new system that nothing was loaded into is taken out again.
LoadResult load_file_into(ApplicationState* app, SystemContext* target, md_str_t path, const loader::LoaderState& load_state);

FileQueue::Entry file_queue_front(const FileQueue* queue);
FileQueue::Entry file_queue_pop(FileQueue* queue);

// Loads what is queued: a workspace, a file to ask about (the load dialog) on their own, else what is queued together as
// one drop, placed as a whole (drop_files.h). Waits while the load dialog is open or a placement is asked about, and
// moves a file between systems when one is asked to be (ApplicationState::file_move).
void file_queue_process(ApplicationState* state);

// The question of where files go which are of the atoms of more than one system (ApplicationState::placement_prompt)
void draw_placement_prompt(ApplicationState* state);

// Moves file i of the system from into the system to, or into a new one for SYSTEM_ID_NONE: loaded there if it is of its
// atoms (file_placement's explicit target), a trajectory with the series along its run, and the system it was in built
// again without them. False, and nothing changed, when it is not of to's atoms.
bool system_move_file(ApplicationState* app, uint32_t from, size_t i, uint32_t to);

// view
// A good view of the atoms in mask (all of them when null), framed in a view of aspect_ratio, its width over its height
// (view_aspect_ratio: the main view's; camera_compute_default_view). fov_y is that of the camera that looks through it.
// With an assembly of sys (an instance's scene.assembly, -1 for none), of those atoms in the instances of each of its
// groups, in all of the copies of the group.
void reset_view(ViewTransform* transform, const md_system_state_t& state, const md_bitset_t* mask = nullptr, float fov_y = Camera().fov_y, float aspect_ratio = 1.0f, const md_system_t* sys = nullptr, int assembly = -1);
static inline void reset_view(ViewTransform* transform, const md_system_state_t& state, const md_bitmask_t* mask, float fov_y = Camera().fov_y, float aspect_ratio = 1.0f, const md_system_t* sys = nullptr, int assembly = -1) {
    const md_bitset_t members = md_bitmask_view(mask);
    reset_view(transform, state, &members, fov_y, aspect_ratio, sys, assembly);
}

// The main view's width over its height (view.screen_size), 1 before it has a size
static inline float view_aspect_ratio(const ApplicationState& state) {
    const md_float2_t s = state.view.screen_size;
    return (s.x > 0.0f && s.y > 0.0f) ? s.x / s.y : 1.0f;
}

// The system whose cell the world is drawn in: the first loaded (the world is translated so that the center of its cell
// is at the origin, SystemData::unitcell_transform), the first held when none is. The others are drawn in the same
// frame, in their own coordinates: overlaid as the files have them.
const SystemContext& scene_frame_system(const ApplicationState& app);

// A good view of every system shown (their visible atoms), as reset_view of the one when only one is
void reset_view_scene(ViewTransform* transform, const ApplicationState& app, float fov_y, float aspect_ratio = 1.0f);

// ASSEMBLIES (system.assembly). The groups of assembly asm_idx are [a.group_offset[asm_idx], a.group_offset[asm_idx + 1]),
// and the instances and operators of group g those assembly_group gives, pointers into the arrays of a. Both read from the
// offsets as they are, not through the md_urange_t of md_assembly_group_range and its kin: MSVC 19.52 (/O2) compiled
// loops from the beg of such a range to index with the 64 bits that held both of its ends, the end in the upper half,
// and read far out of bounds (assembly_transforms, reset_view).
struct AssemblyGroup {
    const md_instance_idx_t* inst;
    const md_instance_idx_t* inst_end;
    const uint32_t* oper;
    const uint32_t* oper_end;
};

static inline AssemblyGroup assembly_group(const md_assembly_data_t& a, size_t g) {
    MD_ASSERT(g < a.group_count);
    return {a.inst + a.group_inst_offset[g], a.inst + a.group_inst_offset[g + 1], a.oper + a.group_oper_offset[g], a.oper + a.group_oper_offset[g + 1]};
}

// The transforms of the copies an assembly of the system makes, group by group, and with copy_offset (optional) where the
// copies of each part begin, its end last. Empty for an assembly the system has not (-1: none).
void assembly_transforms(md_array(md_float4x4_t)* transforms, md_array(uint32_t)* copy_offset, const md_system_t& sys, int assembly, md_allocator_i* alloc);

// The subsets of the parts of the assembly shown (view.assembly_subsets), made again if it or the molecule changed, and
// freed: before the molecule goes
void assembly_subsets_update(SystemInstance& inst);
void assembly_subsets_free(SystemInstance& inst);

// Script IR lifetime
// Every IR is made by script_ir_create, compiling src in the workspace's language, and nothing frees one directly. A
// task that reads an IR off the main thread retains it when it is made (script::program_retain, on the main thread)
// and releases it when it is done (from any thread): the evaluation task the IR it evaluates, the visualization job
// the IR it draws until the job is freed. script_ir_collect frees an IR that neither field (script.ir,
// script.eval_ir) points at and no task holds. Main thread, once per frame.
script::Program* script_ir_create(ApplicationState* state, md_str_t src, const script::CompileDesc& desc);

// Sets the workspace's script language: the editor's highlighting and completion follow, the script is compiled
// again, and every filter is evaluated again in it. The text is kept as it is: for a workspace being loaded, whose
// texts are set beside it.
void script_set_language(ApplicationState* state, script::Language lang);

// Switches to the other language's script: the editor's text is kept as that language's, the other's shown in its
// place, and the language set (script_set_language). What was evaluated is of the language left, and its properties
// show as not evaluated until Evaluate is pressed again.
void script_switch_language(ApplicationState* state, script::Language lang);
void script_ir_collect(ApplicationState* state);

// Frees state->script.eval, once what reads it has been told (EventType_ViamdScriptEvalFree), and leaves it NULL; and
// the evaluations of the other runs. Not while the evaluation task runs.
void script_eval_free(ApplicationState* state);
// Stops every evaluation of the script between frames
void script_evals_interrupt(ApplicationState* state);
// The evaluation of a run of the first system: script.eval for an empty run or its own, else the one of run_evals.
// NULL when there is none.
script::Eval* script_eval_of_run(const ApplicationState* state, md_str_t run);
// Keeps run_evals to the runs of the first system its instances show (other than eval_run): an evaluation of eval_ir
// made for a run shown anew, to be evaluated at the next Evaluate, and the one of a run no longer shown freed. Nothing
// while script.eval is not of eval_ir, nor while the evaluation task runs on one to be freed. Main thread, every frame.
void script_run_evals_sync(ApplicationState* state);

// Script visualization
// The script IR which made ref: the one being compiled (whose tokens the editor marks) or the one being
// evaluated (whose properties the plots show). NULL if neither, e.g. once it has been recompiled.
const script::Program* script_ir_of(const ApplicationState* state, script::Ref ref);
// What is visualized is set anew every frame by whatever is hovered, and the last one set in a frame is what
// is shown. It is computed on a worker thread, never stalling the frame, and shown once ready: from then on for
// as long as it is still what is set. When the atoms move it is computed again, while the previous result
// stays on screen.
//
// Sets ref, visualized in the IR which made it. Nothing if it no longer resolves.
void script_visualize_ref(ApplicationState* state, script::Ref ref, int subidx, uint32_t flags = 0);
// What the value visualized by reference (script_visualize_ref) is at the frame shown, in a few lines: what a hover shows
// beside the type. Empty while it is computed, and for a visualization of an expression.
md_str_t script_vis_value_text(ApplicationState* state);
// How that value is plotted beside its text (script::value_plot), its numbers at the frame shown; one over the trajectory
// once evaluated, while nothing evaluates. NULL when it is not plotted, or while it is computed.
const script::ValuePlot* script_vis_value_plot(ApplicationState* state);
// Sets an expression, compiled against the script being edited (its identifiers)
void script_visualize_str(ApplicationState* state, md_str_t str, uint32_t flags = 0);
// Clears what is set. Start of every frame, before anything sets one.
void script_vis_begin_frame(ApplicationState* state);
// Takes a finished computation, starts one for what is set if it is not computed yet, and decides what is
// shown (script.vis_shown), highlighting its atoms. After everything which sets one, before drawing.
void script_vis_update(ApplicationState* state);
// What is shown of the script (script.vis_shown), for the key of what the main view draws: 0 when nothing is
uint64_t script_vis_shown_key(const ApplicationState* state);
// The main view is made again (lazy rendering): for a change to what it draws that is not in its key
static inline void render_invalidate(ApplicationState* state) { state->render.version += 1; }
// Drops every result, waiting for a running computation. Before the system goes.
void script_vis_reset(ApplicationState* state);

// GROUPS OF AN ENERGY TERM. The groups of the system a term of an energy file is of (md_edr_term_groups): the two energy
// groups of Coul-SR:Protein-SOL, the T-coupling group of T-Protein. attr is a series of the system, whose label keeps the
// term's name. False, and out untouched, for any other series. Its text: "Protein and SOL", "Protein".
bool   energy_term_groups(md_edr_term_groups_t* out, const ApplicationState* state, const md_attribute_t* attr);
size_t energy_term_groups_text(char* buf, size_t cap, const ApplicationState* state, const md_edr_term_groups_t& tg);
// Whether the series is an energy term of groups by its name, whether or not the system has them: for saying so when it
// does not (md_edr_term_is_of_groups)
bool   energy_term_names_groups(const md_attribute_t* attr);
// Highlights the groups of an energy term, each in a colour of its own (selection.highlight_groups): hovering
// Coul-SR:Protein-SOL marks Protein and SOL apart. Every frame it is hovered, as any highlight. False when the series is
// of no groups.
bool   highlight_energy_term_groups(ApplicationState* state, const md_attribute_t* attr);
void script_set_hovered_property(ApplicationState* state, md_str_t label, int population_idx = -1);
