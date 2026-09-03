#pragma once

#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_lru_cache.inl>
#include <core/md_hash.h>
#include <core/md_str_builder.h>

#include <md_system.h>
#include <md_gto.h>
#include <core/md_grid.h>
#include <md_trajectory.h>
#include <md_script.h>
#include <md_gl.h>

#if MD_ENABLE_GPU
#include <core/md_gpu.h>
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
#include <loader.h>
#include <event.h>

#define IMGUI_DEFINE_MATH_OPERATORS

#include <TextEditor.h>
#include <implot.h>
#include <imgui_notify.h>

#include <bitset>
#include <stdint.h>
#include <stddef.h>

#define JITTER_SEQUENCE_SIZE 8  // Number of samples for temporal AA
#define DEFAULT_COLORMAP 5      // This corresponds to plasma colormap (Do not want to include implot.h just for this)
#define FRAME_CACHE_SIZE 4      // Number of frames used in the frame cache

#if FRAME_CACHE_SIZE == 4
#define FRAME_CACHE_LRU_TYPE md_lru_cache4_t
#elif FRAME_CACHE_SIZE == 8
#define FRAME_CACHE_LRU_TYPE md_lru_cache8_t
#else
#error "Unsupported frame cache size"
#endif

// For cpu profiling
#define PUSH_CPU_SECTION(lbl) {};
#define POP_CPU_SECTION() {};

// For gpu profiling
#define PUSH_GPU_SECTION(lbl) { if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); }
#define POP_GPU_SECTION()     { if (glPopDebugGroup) glPopDebugGroup(); }

// For logging
#define VIAMD_LOG_INFO  MD_LOG_INFO
#define VIAMD_LOG_DEBUG MD_LOG_DEBUG
#define VIAMD_LOG_ERROR MD_LOG_ERROR
#define VIAMD_LOG_SUCCESS(...) ImGui::InsertNotification(ImGuiToast(ImGuiToastType_Success, 6000, __VA_ARGS__))

#define DISPLAY_PROPERTY_MAX_POPULATION_SIZE 256
#define DISPLAY_PROPERTY_MAX_TEMPORAL_SUBPLOTS 10
#define DISPLAY_PROPERTY_MAX_DISTRIBUTION_SUBPLOTS 10

#define HIGHLIGHT_PULSE_TIME_SCALE  5.0
#define HIGHLIGHT_PULSE_ALPHA_SCALE 0.1

#define INVALID_PICKING_IDX (~0U)

constexpr ImGuiKey KEY_PLAY_PAUSE               = ImGuiKey_Space;
constexpr ImGuiKey KEY_SKIP_TO_PREV_FRAME       = ImGuiKey_LeftArrow;
constexpr ImGuiKey KEY_SKIP_TO_NEXT_FRAME       = ImGuiKey_RightArrow;
constexpr ImGuiKey KEY_RECOMPILE_SHADERS        = ImGuiKey_F5;
constexpr ImGuiKey KEY_SHOW_DEBUG_WINDOW        = ImGuiKey_F11;
constexpr ImGuiKey KEY_SCRIPT_EVALUATE          = ImGuiKey_Enter;
constexpr ImGuiKey KEY_SCRIPT_EVALUATE_MOD      = ImGuiMod_Shift;
constexpr ImGuiKey KEY_RECENTER_ON_HIGHLIGHT    = ImGuiKey_F1;

constexpr str_t WORKSPACE_FILE_EXTENSION = STR_LIT("via");
constexpr str_t SCRIPT_IMPORT_FILE_EXTENSIONS[] = { STR_LIT("edr"), STR_LIT("xvg"), STR_LIT("csv") };

typedef uint64_t PickingDomainID;
typedef uint64_t PickingSourceID;
typedef uint64_t InteractionSurfaceID;

constexpr PickingDomainID PickingDomain_Atom   = HASH_STR_LIT64("picking domain atom");
constexpr PickingDomainID PickingDomain_Bond   = HASH_STR_LIT64("picking domain bond");
constexpr PickingDomainID PickingDomain_Dipole = HASH_STR_LIT64("picking domain dipole");

constexpr uint64_t interaction_surface_main = HASH_STR_LIT64("interaction surface main"); // This is the main interaction surface which corresponds to the main interaction window, but we want to keep it separate from the picking source and domain ids as we may want to have different picking sources/domains for different interaction surfaces in the future

enum class PlaybackMode { Stopped, Playing, Count };
enum class SelectionGranularity { Atom, Component, Instance, Count };
enum class SelectionOperator { Or, And, AndNot, Set, Clear, Count };
enum class SelectionGrowthMode { CovalentBond, Radial, Count };
enum class CameraMode { Perspective, Orthographic, Count };
enum class InterpolationMode { Nearest, Linear, CubicSpline, Count };

enum class ColorMapping {
    Uniform,
    Type,
    Serial,
    CompName,
    CompSeqId,
    CompIndex,
    InstId,
    InstIndex,
    SecondaryStructure,
    Property,
    Count
};

enum class ElectronicStructureSource {
    MolecularOrbital,
    NaturalTransitionOrbital,
    TransitionDensity,
    ElectronDensity,
    DensityProperty,
    Count
};

enum class ElectronicStructureField {
    Amplitude,
    Density,
    Count
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
};

enum MolBit_ {
    MolBit_None                     = 0,
    MolBit_DirtyPosition            = 1u << 0,
    MolBit_DirtyRadius              = 1u << 1,
    MolBit_DirtySecondaryStructure  = 1u << 2,
    MolBit_DirtyFlags               = 1u << 3,
    MolBit_DirtyBonds               = 1u << 4,
    MolBit_ClearVelocity            = 1u << 5,
};

enum class RepresentationType {
    SpaceFill = MD_GL_REP_SPACE_FILL,
    Licorice = MD_GL_REP_LICORICE,
    BallAndStick = MD_GL_REP_BALL_AND_STICK,
    Ribbons = MD_GL_REP_RIBBONS,
    Cartoon = MD_GL_REP_CARTOON,
    ElectronicStructure,
    DipoleMoment,
    Count
};

inline const char* selection_granularity_str[(int)SelectionGranularity::Count] = {
    "Atom",
    "Component",
    "Instance",
};

inline const char* interpolation_mode_str[(int)InterpolationMode::Count] = {
    "Nearest",
    "Linear",
    "Cubic Spline",
};

inline const char* representation_type_str[(int)RepresentationType::Count] = {
    "Spacefill",
    "Licorice",
    "Ball And Stick",
    "Ribbons",
    "Cartoon",
    "Electronic Structure",
    "Dipole Moment"
};

inline const char* color_mapping_str[(int)ColorMapping::Count] = {
    "Uniform Color",
    "Type",
    "Serial",
    "Res Name",
    "Seq Id",
    "Res Idx",
    "Chain Id",
    "Chain Idx",
    "Secondary Structure",
    "Property",
};


inline const char* electronic_structure_source_str[(int)ElectronicStructureSource::Count] = {
    "Molecular Orbital",
    "Natural Transition Orbital",
    "Transition Density",
    "Electron Density",
    "Density Property",
};

enum ElectronicStructureSourceFlag_ : uint32_t {
    ElectronicStructureSourceFlag_MolecularOrbital          = 1u << (int)ElectronicStructureSource::MolecularOrbital,
    ElectronicStructureSourceFlag_NaturalTransitionOrbital  = 1u << (int)ElectronicStructureSource::NaturalTransitionOrbital,
    ElectronicStructureSourceFlag_TransitionDensity         = 1u << (int)ElectronicStructureSource::TransitionDensity,
    ElectronicStructureSourceFlag_ElectronDensity           = 1u << (int)ElectronicStructureSource::ElectronDensity,
    ElectronicStructureSourceFlag_DensityProperty           = 1u << (int)ElectronicStructureSource::DensityProperty,
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

// This is viamd's representation of a property
struct DisplayProperty {
    enum Type {
        Type_Temporal,
        Type_Distribution,
        Type_Volume,
        Type_Count
    };

    enum PlotType {
        PlotType_Line,      // Single line
        PlotType_Area,      // Shaded area
        PlotType_Bars,      // Bar chart
        PlotType_Scatter,   // Scatter plot
        PlotType_Count
    };

    enum ColorType {
        ColorType_Solid,
        ColorType_Colormap,
        ColorType_Count
    };

    // This is the payload passed to getters for display properties
    struct Payload {
        DisplayProperty* display_prop;
        int dim_idx;
    };

    // Callback signature for printing out the value (when hovering with mouse for example)
    typedef int (*PrintValue)(char* buf, size_t buf_cap, int sample_idx, Payload* data);

    struct Histogram {
        int num_bins;
        int dim = 0;
        // Can be multidimensional
        // Total number of entries will be dim * num_bins
        md_array(float) bins = 0;
        double x_min;
        double x_max;
        double y_min;
        double y_max;
        md_allocator_i* alloc;
    };

    Type type = Type_Temporal;

    char label[32] = "";

    ColorType color_type = ColorType_Solid;
    ImVec4 color = {1,1,1,1};
    ImPlotColormap colormap = ImPlotColormap_Plasma;
    float colormap_alpha = 1.0f;

    PlotType plot_type = PlotType_Line;
    ImPlotMarker marker_type = ImPlotMarker_Square;
    float marker_size = 1.0f;
    double bar_width_scale = 1.0;

    // We need two getters to support areas (min / max)
    ImPlotGetter getter[2] = {0,0};
    PrintValue   print_value = 0;

    bool aggregate_histogram = false;

    int dim = 1;                // Number of values per sample
    int num_samples = 0;        // Number of samples (length of x)
    const float* y_values = 0;  // Values (y)
    const float* x_values = 0;  // Corresponding x values

    int num_bins = 128;         // Requested number of bins for histogram

    md_unit_t unit[2] = {md_unit_none(), md_unit_none()};
    char unit_str[2][32] = {"",""};

    const md_script_eval_t* eval = NULL;

    md_script_property_flags_t prop_flags = MD_SCRIPT_PROPERTY_FLAG_NONE;
    const md_script_property_data_t* prop_data = NULL;
    const md_script_vis_payload_o* vis_payload = NULL;

    uint64_t prop_fingerprint = 0;

    // Encodes which temporal subplots this property is visible in
    uint32_t temporal_subplot_mask = 0;

    // Encodes which distribution subplots this property is visible in
    uint32_t distribution_subplot_mask = 0;

    bool show_in_volume = false;
    bool partial_evaluation = false;

    // Encodes which indices of the population to show (if applicable, i.e. dim > 1)
    std::bitset<DISPLAY_PROPERTY_MAX_POPULATION_SIZE> population_mask = {};

    STATIC_ASSERT(DISPLAY_PROPERTY_MAX_TEMPORAL_SUBPLOTS     <= sizeof(temporal_subplot_mask) * 8,     "Cannot fit temporal subplot mask");
    STATIC_ASSERT(DISPLAY_PROPERTY_MAX_DISTRIBUTION_SUBPLOTS <= sizeof(distribution_subplot_mask) * 8, "Cannot fit distribution subplot mask");

    Histogram hist = {};
};

struct LoadDatasetWindowState {
    char path_buf[1024] = "";
    char atom_format_buf[128] = "";
    char err_buf[128] = "";
    bool path_is_valid = false;
    bool path_changed = false;
    bool load_topology = false;
    bool load_trajectory = false;
    bool coarse_grained = false;
    bool show_window = false;
    bool show_file_dialog = false;
    bool atom_format_valid = false;
    int  loader_idx = 0;
    int  atom_format_idx = 0;
};

// Hint flags for operations to be performed for the files
enum {
    FileFlags_None = 0,
    FileFlags_ShowDialogue = 1,
    FileFlags_CoarseGrained = 2,
    FileFlags_DisableCacheWrite = 4,
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
    Window,
    FHD,
    QHD,
    UHD_4K,
    UHD_8K,
    Custom,
    Count,
};

inline const char* screenshot_resolution_str[(int)ScreenshotResolution::Count] = {
    "Window",
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

inline const char* bond_color_mode_str[(int)BondColorMode::Count] = {
    "Nearest Atom",
    "Smooth Atom",
    "Uniform Color",
};

struct FileQueue {
    struct Entry {
        str_t path;
        FileFlags flags;
        int prio;
    };
    Entry arr[8] = {};
    uint32_t head = 0;
    uint32_t tail = 0;

    // We use a ring alloc here to do a deep copy of the paths within the added entries
    // Because its a ring alloc, there is no need to free entries
    md_allocator_i* ring = NULL;
};

struct SingleSelectionSequence {
    int32_t idx[4] = {-1, -1, -1, -1};
};

struct Selection {
    char name[64] = "sel";
    md_bitfield_t atom_mask {};
};

// One dipole gathered from the system's attribute table, not accumulated by anyone.
// A group is two attributes which the producer published, sharing a shape:
//     dipole/<group>/vector    the moment, in `unit`
//     dipole/<group>/origin    where to draw it from, Angstrom
// A group holds ONE moment (a ground state) or one per excited state, and gathering flattens that:
// each element of each group is its own entry. key alone therefore does not address a dipole -
// (key, index) does, which is why a representation stores both.
// Build these with dipole_moments_gather; nothing owns or caches them.
struct DipoleMoment {
    md_attribute_id_t key = MD_ATTRIBUTE_INVALID;   // id of the vector attribute, stable across a reload
    uint32_t index = 0;                             // which element within the group
    uint32_t count = 0;                             // elements in the group; 1 means index is decorative
    str_t  label  = { 0 };                          // group name, prettified for display
    vec3_t vec    = { 0, 0, 0 };
    vec3_t origin = { 0, 0, 0 };
    md_unit_t unit = md_unit_none();
};

struct DensityProperty {
    uint64_t key    = 0;
    str_t label     = { 0 };
};

// Struct to fill in for the different components
// Which provides information of what representations are available for the currently loaded datasets
// What is LEFT of a fan-in that used to be filled by whichever component happened to have parsed
// the file. Both remaining members are answers about what the SYSTEM holds, computed in
// update_representation_info by reading its attribute table - which is the exact test, and one a
// second reader satisfies without any component being involved.
//
// The orbital and excited-state members went the same way: main.cpp asks the table (es_orbital_*,
// es_excited_state_count, es_nto_lambdas) rather than being handed a precomputed copy.
struct RepresentationInfo {
    ElectronicStructureSourceFlags electronic_structure_source_mask = 0;

    md_array(DensityProperty) density_properties = nullptr;

    md_allocator_i* alloc = nullptr;
};

// ---------------------------------------------------------------------------
// Volume geometry
// ---------------------------------------------------------------------------
// Fitting an object aligned box around a set of points and turning it into a grid and a volume.
// None of it knows what is being evaluated in that volume or who asked, which is why it sits here
// rather than in whichever component first needed it.

// Object aligned bounding box
struct OABB {
    mat3_t orientation = mat3_ident();
    vec3_t min_ext = { 0 };
    vec3_t max_ext = { 0 };
};

struct AABB {
    vec3_t min_ext = { 0 };
    vec3_t max_ext = { 0 };
};

static mat3_t mat3_PCA(const vec4_t* xyzw, size_t count) {
    vec3_t acc = vec3_zero();
    for (size_t i = 0; i < count; ++i) {
        acc = vec3_add(acc, vec3_from_vec4(xyzw[i]));
    }
    vec3_t mean = acc / (float)count;

    mat3_t cov = mat3_covariance_matrix_vec4(xyzw, nullptr, count, mean);
    mat3_eigen_t eigen = mat3_eigen(cov);
    mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));
    return PCA;
}

static void calculate_bounds(float out_min[3], float out_max[3], const vec4_t* xyzw, size_t count, const mat3_t& orientation = mat3_ident()) {
    vec4_t min_v = vec4_set1( FLT_MAX);
    vec4_t max_v = vec4_set1(-FLT_MAX);

    mat4_t rot = mat4_from_mat3(mat3_transpose(orientation));

    for (size_t i = 0; i < count; ++i) {
        vec4_t v = mat4_mul_vec4(rot, xyzw[i]);
        min_v = vec4_min(min_v, v);
        max_v = vec4_max(max_v, v);
    }

    // Padding
    const float pad = 6.0f;
    min_v -= pad;
    max_v += pad; 

	MEMCPY(out_min, &min_v, sizeof(float) * 3);
	MEMCPY(out_max, &max_v, sizeof(float) * 3);
}

// Construct texture to world transformation matrix for Volume
// extent is the extent of the volume (dim * voxel_size)
static inline mat4_t compute_texture_to_world_mat(const mat3_t& orientation, const vec3_t& origin, const vec3_t& extent) {
    mat4_t T = mat4_translate_vec3(origin);
    mat4_t R = mat4_from_mat3(orientation);
    mat4_t S = mat4_scale_vec3(extent);
    return T * R * S;
}

static inline mat4_t compute_world_to_model_mat(const mat3_t& orientation, const vec3_t& origin) {
    mat4_t world_to_model = mat4_from_mat3(mat3_transpose(orientation)) * mat4_translate_vec3(-origin);
    return world_to_model;
}

static inline mat4_t compute_index_to_world_mat(const mat3_t& orientation, const vec3_t& in_origin, const vec3_t& stepsize) {
    vec3_t step_x = orientation.col[0] * stepsize.x;
    vec3_t step_y = orientation.col[1] * stepsize.y;
    vec3_t step_z = orientation.col[2] * stepsize.z;
    // Shift origin by half voxel
    vec3_t origin = in_origin + orientation * (stepsize * 0.5f);

    mat4_t index_to_world = {
        step_x.x, step_x.y, step_x.z, 0.0f,
        step_y.x, step_y.y, step_y.z, 0.0f,
        step_z.x, step_z.y, step_z.z, 0.0f,
        origin.x, origin.y, origin.z, 1.0f,
    };

    return index_to_world;
}

// Attempts to compute fitting volume dimensions given an input extent and a suggested number of samples per length unit
static inline void compute_dim(int out_dim[3], const vec3_t& in_ext, double samples_per_unit_length) {
    out_dim[0] = CLAMP(ALIGN_TO((int)(in_ext.x * samples_per_unit_length), 8), 8, 512);
    out_dim[1] = CLAMP(ALIGN_TO((int)(in_ext.y * samples_per_unit_length), 8), 8, 512);
    out_dim[2] = CLAMP(ALIGN_TO((int)(in_ext.z * samples_per_unit_length), 8), 8, 512);
}

// Grid units are BOHR, matching what md_gto evaluates in; a Volume's transforms are Angstrom,
// matching the world the camera lives in. init_volume is where the two meet, and it is the only
// place that conversion belongs.
static inline void init_grid(md_grid_t* grid, const mat3_t& orientation, const vec3_t& min_ext, const vec3_t& max_ext, double samples_per_unit_length) {
    ASSERT(grid);
    vec3_t extent = max_ext - min_ext;
    compute_dim(grid->dim, extent, samples_per_unit_length);
    vec3_t voxel_size = vec3_div(extent, vec3_set((float)grid->dim[0], (float)grid->dim[1], (float)grid->dim[2]));
    grid->orientation = orientation;
    grid->origin = orientation * min_ext;
    grid->spacing = voxel_size;
}

struct Volume {
    mat4_t world_to_model   = {};   // Roto-translation into volume local axes, no scaling applied (preserves world length units)
    mat4_t texture_to_world = {};   // Texture space [0,1] to world coordinates
    vec3_t voxel_size  = {1,1,1};   // Size of each voxel in world units
    int dim[3] = {128, 128, 128};
    uint32_t tex_id = 0;
};

// Sizes and allocates the GL 3D texture for a grid, and fills in the transforms that place it in
// the world. Defined in viamd.cpp because it is the half that touches GL.
void init_volume(Volume* vol, const md_grid_t& grid, GLenum format = GL_R16F);

// Descriptor for handling iso surfaces
struct IsoDesc {
    size_t count;
    float values[8];
    vec4_t colors[8];
    float optical_densities[8];
};

struct ElectronicStructureRepresentation {
    Volume density_vol = {};
    Volume color_vol   = {};
    VolumeResolution resolution = VolumeResolution::Mid;

    // Shared for all electronic structure representations
    double iso_value = 0.05;

    // These are default values for different volume types.
    vec4_t col_psi_pos = {0.f/255.f,75.f/255.f,135.f/255.f,0.75f};
    vec4_t col_psi_neg = {255.f/255.f,205.f/255.f,0.f/255.f,0.75f};
    vec4_t col_den     = {255.f/255.f,255.f/255.f,255.f/255.f,0.75f};
    vec4_t col_att     = {0, 162.0f/255.0f, 135.0f/255.0f, 0.75f};
    vec4_t col_det     = {162.0f/255.0f, 35.0f/255.0f, 135.0f/255.0f, 0.75f};

    vec4_t tint_psi_pos = { 1.0f, 1.0f, 1.0f, 0.75f };
    vec4_t tint_psi_neg = { 1.0f, 1.0f, 1.0f, 0.75f };
    vec4_t tint_den     = { 1.0f, 1.0f, 1.0f, 0.75f };
    vec4_t tint_att     = { 1.0f, 1.0f, 1.0f, 0.75f };
    vec4_t tint_det     = { 1.0f, 1.0f, 1.0f, 0.75f };

    struct {
        int num_isos = 0;
        double values[8];
        vec4_t colors[8];
    } density_property;

    // Scaling factor of *power* in the gaussians to splat the color volume (when using atom colors for volumes)
    double gaussian_splatting_power = 10.0;

    // Optical scaling factor which controls attenuation of light within iso surfaces.
    double iso_optical_density = 0.005;

    bool use_atom_colors = false;

    struct {
        bool enabled = false;
        uint32_t tf_tex = 0;
        int colormap = DEFAULT_COLORMAP;
    } dvr;

    ElectronicStructureSource source = ElectronicStructureSource::MolecularOrbital;
    bool use_magnitude = false;
    ElectronicStructureSpin spin = ElectronicStructureSpin::Alpha;
    ElectronicStructureNtoComponent nto_component = ElectronicStructureNtoComponent::Particle;
    ElectronicStructureTransitionDensityComponent transition_density_component = ElectronicStructureTransitionDensityComponent::Attachment;

    int orbital_idx = 0;
    int excited_state_idx = 0;
    int nto_lambda_idx = 0;
    // The density property this representation draws, named by the id of its attribute rather than
    // by a position in a list. An index silently re-points at a different property whenever the set
    // changes across a reload; an id is derived from the path, so it names the same thing or
    // nothing at all.
    md_attribute_id_t density_property_key = MD_ATTRIBUTE_INVALID;

	uint64_t col_hash = 0;
	uint64_t vol_hash = 0;
    uint64_t tf_hash = 0;
};

static inline ElectronicStructureSourceFlags electronic_structure_source_flag(ElectronicStructureSource source) {
    return 1u << (uint32_t)source;
}

static inline bool electronic_structure_source_supported(ElectronicStructureSourceFlags mask, ElectronicStructureSource source) {
    return (mask & electronic_structure_source_flag(source)) != 0;
}

static inline bool electronic_structure_is_density_property(const ElectronicStructureRepresentation& rep) {
    return rep.source == ElectronicStructureSource::DensityProperty;
}

static inline bool electronic_structure_uses_orbital_idx(const ElectronicStructureRepresentation& rep) {
    return rep.source == ElectronicStructureSource::MolecularOrbital;
}

static inline bool electronic_structure_uses_excited_state_idx(const ElectronicStructureRepresentation& rep) {
    return rep.source == ElectronicStructureSource::NaturalTransitionOrbital ||
           rep.source == ElectronicStructureSource::TransitionDensity;
}

static inline bool electronic_structure_uses_nto_lambda_idx(const ElectronicStructureRepresentation& rep) {
    return rep.source == ElectronicStructureSource::NaturalTransitionOrbital;
}

static inline bool electronic_structure_uses_magnitude_toggle(const ElectronicStructureRepresentation& rep) {
    return rep.source == ElectronicStructureSource::MolecularOrbital ||
           rep.source == ElectronicStructureSource::NaturalTransitionOrbital ||
           (rep.source == ElectronicStructureSource::ElectronDensity && rep.spin == ElectronicStructureSpin::Difference);
}

static inline bool electronic_structure_uses_spin(const ElectronicStructureRepresentation& rep) {
    return rep.source == ElectronicStructureSource::MolecularOrbital ||
           rep.source == ElectronicStructureSource::ElectronDensity;
}

static inline bool electronic_structure_is_signed(const ElectronicStructureRepresentation& rep) {
    if (electronic_structure_is_density_property(rep)) {
        return false;
    }
    if (rep.source == ElectronicStructureSource::TransitionDensity &&
        rep.transition_density_component == ElectronicStructureTransitionDensityComponent::Difference) {
        return true;
    }
    return electronic_structure_uses_magnitude_toggle(rep) && !rep.use_magnitude;
}

static inline bool electronic_structure_uses_density_iso_squared(const ElectronicStructureRepresentation& rep) {
    return rep.source == ElectronicStructureSource::TransitionDensity ||
           rep.source == ElectronicStructureSource::ElectronDensity;
}

static inline vec4_t electronic_structure_density_property_default_color(const ElectronicStructureRepresentation& rep, int idx) {
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
    ASSERT(rep);
    rep->density_property.num_isos = 2;
    rep->density_property.values[0] =  0.05;
    rep->density_property.values[1] = -0.05;
    rep->density_property.colors[0] = rep->col_psi_pos;
    rep->density_property.colors[1] = rep->col_psi_neg;
    for (int i = 2; i < (int)ARRAY_SIZE(rep->density_property.values); ++i) {
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
    ASSERT(iso);
    *iso = {};    

    if (electronic_structure_is_density_property(rep)) {
        const int num_isos = CLAMP(rep.density_property.num_isos, 0, (int)ARRAY_SIZE(rep.density_property.values));
        iso->count = (size_t)num_isos;
        for (int i = 0; i < num_isos; ++i) {
            iso->values[i] = (float)rep.density_property.values[i];
            iso->colors[i] = rep.density_property.colors[i];
            iso->optical_densities[i] = (float)rep.iso_optical_density;
        }
        return;
    }

    const float iso_value = (float)rep.iso_value;
    const float iso_threshold = electronic_structure_uses_density_iso_squared(rep) ? iso_value * iso_value : iso_value;

    if (electronic_structure_is_signed(rep)) {
        iso->count = 2;
        iso->values[0] =  iso_threshold;
        iso->values[1] = -iso_threshold;
        if (rep.use_atom_colors) {
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
        if (rep.use_atom_colors) {
            if (rep.source == ElectronicStructureSource::TransitionDensity && rep.transition_density_component == ElectronicStructureTransitionDensityComponent::Attachment) {
                iso->colors[0] = rep.tint_att;
            } else if (rep.source == ElectronicStructureSource::TransitionDensity && rep.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment) {
                iso->colors[0] = rep.tint_det;
            } else {
                iso->colors[0] = rep.tint_den;
            }
        } else {
            if (rep.source == ElectronicStructureSource::TransitionDensity && rep.transition_density_component == ElectronicStructureTransitionDensityComponent::Attachment) {
                iso->colors[0] = rep.col_att;
            } else if (rep.source == ElectronicStructureSource::TransitionDensity && rep.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment) {
                iso->colors[0] = rep.col_det;
            } else {
                iso->colors[0] = rep.col_den;
            }
        }
        iso->optical_densities[0] = (float)rep.iso_optical_density;
    }
}

static inline ElectronicStructureField electronic_structure_legacy_field(const ElectronicStructureRepresentation& rep) {
    return rep.use_magnitude ? ElectronicStructureField::Density : ElectronicStructureField::Amplitude;
}

static inline void electronic_structure_set_legacy_field(ElectronicStructureRepresentation* rep, int field) {
    ASSERT(rep);
    rep->use_magnitude = (rep->source == ElectronicStructureSource::MolecularOrbital ||
                          rep->source == ElectronicStructureSource::NaturalTransitionOrbital) &&
                         (ElectronicStructureField)field == ElectronicStructureField::Density;
}

static inline void electronic_structure_set_source_defaults(ElectronicStructureRepresentation* rep) {
    ASSERT(rep);
    switch (rep->source) {
    case ElectronicStructureSource::MolecularOrbital:
        if (rep->spin != ElectronicStructureSpin::Alpha && rep->spin != ElectronicStructureSpin::Beta) {
            rep->spin = ElectronicStructureSpin::Alpha;
        }
        rep->use_magnitude = false;
        break;
    case ElectronicStructureSource::NaturalTransitionOrbital:
        rep->spin = ElectronicStructureSpin::None;
        rep->use_magnitude = false;
        break;
    case ElectronicStructureSource::TransitionDensity:
        rep->spin = ElectronicStructureSpin::None;
        rep->use_magnitude = false;
        break;
    case ElectronicStructureSource::ElectronDensity:
        rep->spin = ElectronicStructureSpin::Total;
        rep->use_magnitude = false;
        break;
    case ElectronicStructureSource::DensityProperty:
        rep->spin = ElectronicStructureSpin::None;
        rep->use_magnitude = false;
        rep->use_atom_colors = false;
        if (rep->density_property.num_isos <= 0) {
            electronic_structure_density_property_init_defaults(rep);
        }
        break;
    default:
        ASSERT(false);
        break;
    }
}

static inline void electronic_structure_set_legacy_type(ElectronicStructureRepresentation* rep, int type) {
    ASSERT(rep);
    switch ((ElectronicStructureLegacyType)type) {
    case ElectronicStructureLegacyType::MolecularOrbital:
        rep->source = ElectronicStructureSource::MolecularOrbital;
        rep->use_magnitude = false;
        break;
    case ElectronicStructureLegacyType::MolecularOrbitalDensity:
        rep->source = ElectronicStructureSource::MolecularOrbital;
        rep->use_magnitude = true;
        break;
    case ElectronicStructureLegacyType::NaturalTransitionOrbitalParticle:
        rep->source = ElectronicStructureSource::NaturalTransitionOrbital;
        rep->use_magnitude = false;
        rep->nto_component = ElectronicStructureNtoComponent::Particle;
        break;
    case ElectronicStructureLegacyType::NaturalTransitionOrbitalHole:
        rep->source = ElectronicStructureSource::NaturalTransitionOrbital;
        rep->use_magnitude = false;
        rep->nto_component = ElectronicStructureNtoComponent::Hole;
        break;
    case ElectronicStructureLegacyType::NaturalTransitionOrbitalDensityParticle:
        rep->source = ElectronicStructureSource::NaturalTransitionOrbital;
        rep->use_magnitude = true;
        rep->nto_component = ElectronicStructureNtoComponent::Particle;
        break;
    case ElectronicStructureLegacyType::NaturalTransitionOrbitalDensityHole:
        rep->source = ElectronicStructureSource::NaturalTransitionOrbital;
        rep->use_magnitude = true;
        rep->nto_component = ElectronicStructureNtoComponent::Hole;
        break;
    case ElectronicStructureLegacyType::AttachmentDensity:
        rep->source = ElectronicStructureSource::TransitionDensity;
        rep->use_magnitude = false;
        rep->transition_density_component = ElectronicStructureTransitionDensityComponent::Attachment;
        break;
    case ElectronicStructureLegacyType::DetachmentDensity:
        rep->source = ElectronicStructureSource::TransitionDensity;
        rep->use_magnitude = false;
        rep->transition_density_component = ElectronicStructureTransitionDensityComponent::Detachment;
        break;
    case ElectronicStructureLegacyType::ElectronDensity:
        rep->source = ElectronicStructureSource::ElectronDensity;
        rep->use_magnitude = false;
        rep->spin = ElectronicStructureSpin::Total;
        break;
    default:
        ASSERT(false);
        break;
    }
}


static inline int electronic_structure_legacy_type(const ElectronicStructureRepresentation& rep) {
    switch (rep.source) {
    case ElectronicStructureSource::MolecularOrbital:
        return (int)(!rep.use_magnitude ?
            ElectronicStructureLegacyType::MolecularOrbital :
            ElectronicStructureLegacyType::MolecularOrbitalDensity);
    case ElectronicStructureSource::NaturalTransitionOrbital:
        if (!rep.use_magnitude) {
            return (int)(rep.nto_component == ElectronicStructureNtoComponent::Particle ?
                ElectronicStructureLegacyType::NaturalTransitionOrbitalParticle :
                ElectronicStructureLegacyType::NaturalTransitionOrbitalHole);
        }
        return (int)(rep.nto_component == ElectronicStructureNtoComponent::Particle ?
            ElectronicStructureLegacyType::NaturalTransitionOrbitalDensityParticle :
            ElectronicStructureLegacyType::NaturalTransitionOrbitalDensityHole);
    case ElectronicStructureSource::TransitionDensity:
        return (int)(rep.transition_density_component == ElectronicStructureTransitionDensityComponent::Detachment ?
            ElectronicStructureLegacyType::DetachmentDensity :
            ElectronicStructureLegacyType::AttachmentDensity);
    case ElectronicStructureSource::ElectronDensity:
        return (int)ElectronicStructureLegacyType::ElectronDensity;
    case ElectronicStructureSource::DensityProperty:
        return (int)ElectronicStructureLegacyType::ElectronDensity;
    default:
        ASSERT(false);
        return 0;
    }
}

// Colouring by a per atom scalar field. The field itself is an attribute on the system, so what
// is stored here is its id and nothing else about it - no copied label, no cached list position.
// The id is a hash of the path and stable across a reload, so a representation keeps pointing at
// the same quantity when the data is reloaded; an index into a gathered list would not.
struct AtomicPropertyRepresentation {
    int colormap = DEFAULT_COLORMAP;
    md_attribute_id_t key = MD_ATTRIBUTE_INVALID;
    int variant_idx = 0;                // position along the leading axis, when the field has one
    float value_min = 0.0f;             // span of the DATA, refreshed when key changes
    float value_max = 1.0f;
    float range_beg = 0.0f;             // span the colour ramp is mapped over, user adjustable
    float range_end = 1.0f;
    bool  range_symmetric_zero = true;  // Use a symmetric min and max value around zero
};

struct DipoleRepresentation {
    md_attribute_id_t dipole_key = MD_ATTRIBUTE_INVALID;   // the vector attribute, so it survives a reload
    uint32_t dipole_index = 0;                             // which element of that group; see DipoleMoment
    vec4_t color = { 0, 0, 0, 1 };
    vec3_t offset = { 0, 0, 0 };
    double scale = 1.0;
	float radius = 0.05f;
};

struct Representation {
    char name[64] = "rep";
    char filt[256] = "all";
    char filt_error[256] = "";

    RepresentationType type = RepresentationType::BallAndStick;
    ColorMapping color_mapping = ColorMapping::Type;
    md_bitfield_t atom_mask = {};
    md_gl_rep_t md_rep = {};
#if EXPERIMENTAL_GFX_API
    md_gfx_handle_t gfx_rep = {};
#endif

    bool enabled = true;
	bool needs_update = true;
    bool type_is_valid = false;
    bool filt_is_dirty = true;
    bool filt_is_valid = false;
    bool filt_is_dynamic = false;
    bool dynamic_evaluation = true;

    // User defined base color used in uniform mode
    vec4_t base_color = {1.0f, 1.0f, 1.0f, 1.0f};

    struct {
        vec4_t color_unknown = {0.50f, 0.50f, 0.50f, 1.0f};
        vec4_t color_coil    = {0.86f, 0.86f, 0.86f, 1.0f};
        vec4_t color_helix   = {0.12f, 0.86f, 0.12f, 1.0f};
        vec4_t color_sheet   = {0.12f, 0.12f, 0.86f, 1.0f};
    } secondary_structure;

    // Global post processing parameters applied to the final colors.
    vec4_t tint_color = {1.0f, 0.0f, 0.0f, 1.0f};
    float tint_scale = 0.0f;
    float saturation = 1.0f;

    // scaling parameter (radius, width, height, etc depending on type)
    vec4_t scale = {1.0f, 1.0f, 1.0f, 1.0f};

	BondColorMode bond_color = BondColorMode::NearestAtom;
    float  bond_sharpness = 0.5f; // 0 = sharper, 1 = smoother
	vec4_t bond_base_color = { 1.0f, 1.0f, 1.0f, 1.0f };

    ElectronicStructureRepresentation electronic_structure = {};
    AtomicPropertyRepresentation atomic_property = {};
	DipoleRepresentation dipole = {};
};

struct FrameCache {
    FRAME_CACHE_LRU_TYPE lru = {};
    md_system_state_t states[FRAME_CACHE_SIZE] = {};
    int32_t frame_idx[FRAME_CACHE_SIZE] = {};
};

struct PickingRange {
    PickingDomainID domain = 0;
    uint32_t beg = 0;
    uint32_t end = 0;
};

struct PickingSpace {
    size_t num_ranges = 0;
    PickingRange ranges[8] = {};
};

struct PickingReadbackSlot {
    uint32_t color_pbo = 0;
    uint32_t depth_pbo = 0;
    
    uint32_t submitted_frame_idx = 0;
    bool pending = false;
    
    uint32_t viewport_width = 0;
    uint32_t viewport_height = 0;
    
    vec2_t surface_coord = {0};
    vec2_t screen_coord  = {0};
    mat4_t clip_to_world = MD_MAT4_IDENT_INIT;
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

    uint32_t frame_idx = 0;    uint32_t raw_idx = INVALID_PICKING_IDX;
    uint32_t local_idx = 0;

    vec2_t surface_coord = {0};
    vec2_t screen_coord = {0};
    vec3_t world_pos = {0};
    float depth = 1.0f;
};

struct PickingReadbackRequest {
    uint32_t fbo = 0;
    uint32_t width = 0;
    uint32_t height = 0;
    vec2_t surface_coord = {0};
    vec2_t screen_coord = {0};
    mat4_t clip_to_world = MD_MAT4_IDENT_INIT;
};

struct ApplicationState {
    // --- APPLICATION ---
    application::Context app {};

#if MD_ENABLE_GPU
    md_gpu_device_t gpu_device = nullptr;

    // Evaluation scratch. Device scoped and deliberately NOT per dataset: the volume is a fixed
    // 512^3 R32F texture - half a gigabyte - and every evaluation reads it straight back into its
    // own destination texture, so one is enough no matter how many datasets or representations are
    // asking. The coefficient buffer is grown to fit the widest basis loaded so far, for the same
    // reason. Whoever evaluates borrows these; nobody else frees them.
    md_gpu_stream_t gpu_stream  = nullptr;   // the device's default compute stream
    md_gpu_pool_t   gpu_pool    = nullptr;   // device-local: basis, atoms, coefficients
    md_gpu_pool_t   gpu_rb_pool = nullptr;   // host-readable: volume readback staging
    md_gpu_tex_t    gpu_volume  = 0;         // 3d R32F scratch for an evaluated orbital / density
    md_gpu_ptr_t    gpu_coeff   = nullptr;   // AO coefficient staging, sized to the widest basis
    size_t          gpu_coeff_capacity = 0;

    // One in-flight readback of gpu_volume into a GL texture. A slot must outlive the call that
    // queued it, so these live here rather than in the frame arena - and here specifically because
    // every one of them references gpu_volume, gpu_rb_pool and gpu_stream above. Whoever destroys
    // those drains these first.
    struct GpuVolumeJob {
        bool              in_flight = false;
        ApplicationState* owner     = nullptr;
        uint32_t          tex_id    = 0;        // GL texture to receive the data
        md_gpu_ptr_t      rb        = nullptr;  // HOST_READ staging block
        size_t            size      = 0;
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
    // for keeping track of open files
    struct {
        char molecule[1024]   = {0};
        char trajectory[1024] = {0};
        char workspace[1024]  = {0};

        bool coarse_grained = false;
    } files;

    // The idea for the file load queue is to fill it with files that are dropped onto the application
    // Or passed via the commandline. The files need to be processed in a certain order (topology before trajectory)
    // It also provides a way to chain the file load dialogue.
    FileQueue file_queue = {};

    // --- CAMERA ---
    struct {
        Camera camera{};
        TrackballControllerParam trackball_param;
        ViewParam param{};
        CameraMode mode = CameraMode::Perspective;

        struct {
            vec2_t sequence[JITTER_SEQUENCE_SIZE] {};
        } jitter;

		ViewTransform target = {};
    } view;

    struct {
        bool  hide_gui = true;
        ScreenshotResolution resolution = ScreenshotResolution::Window;
        int   res_x = 1920;
        int   res_y = 1080;
        int   sample_count = 0;
        int   sample_target = 0;
        str_t path_to_file = {};
    } screenshot;

    struct {
        md_gl_shaders_t shaders = {};
        md_gl_shaders_t shaders_lean_and_mean = {};
    } gl;

    // --- MDLIB DATA ---
    struct {
        // The arena backing this dataset is sys.alloc - there is no second handle to it. It is
        // created once at startup and survives every load: free_system_data rewinds it and puts the
        // handle back, because zeroing a system also zeroes the allocator it was using.
        md_gl_mol_t         gl_mol = {};

        // Derived FOR DRAWING, in the renderer's own format, never read back - the same rule the
        // GPU block below states, which is why this belongs beside gl_mol rather than in one of the
        // attribute tables. md_gl_secondary_structure_t is a pair of blend weights, not a physical
        // quantity: nothing would ask "what is the helix weight of segment 12" as a question about
        // the molecule, and its only consumer is md_gl_mol_set_backbone_secondary_structure.
        //
        // Per DATASET, so it replicates with the rest of this block when several can be loaded.
        struct {
            md_array(md_gl_secondary_structure_t) secondary_structure = nullptr;
        } interpolated_properties;
#if EXPERIMENTAL_GFX_API
        md_gfx_handle_t     gfx_structure = {};
#endif
        md_system_t         sys = {};

        // The actual interpolated state of the system.
        // Is only used for rendering and visualization of properties.
        md_system_state_t   state = {};

		mat4_t 			    unitcell_transform = MD_MAT4_IDENT_INIT;

        FrameCache          frame_cache;

        vec3_t              sys_aabb_min = {};
        vec3_t              sys_aabb_max = {};

        bool                interpolate_system_state = false;
        uint32_t            dirty_gpu_buffers = 0;

#if MD_ENABLE_GPU
        // GPU side data derived from sys, sitting beside gl_mol and for the same reason: it is
        // built FROM the system so that something can draw or evaluate it, and mdlib never reads it
        // back. The test is simple - no mdlib function takes an md_system_t* and touches these -
        // which is what keeps them out of md_system_t itself, unlike the trajectory it does read.
        //
        // This whole block is per DATASET, so it is what gets replicated when several systems can
        // be loaded at once. That is also why the evaluation scratch is not here but on the device.
        md_gto_gpu_basis_t  gpu_basis = nullptr;   // built from the basis/ attributes, not from a loader
        md_gpu_ptr_t        gpu_atoms = nullptr;   // packed float4 positions, xyz in Bohr
        // Hash of the positions currently in gpu_atoms, 0 when nothing has been uploaded. The
        // positions come from the system STATE, so they move with the trajectory; comparing what is
        // uploaded against what is wanted is the only test that cannot go stale, and it costs a hash
        // over the QM atoms against an upload of the same array.
        uint64_t            gpu_atoms_hash = 0;
#endif
    } mold;

    DisplayProperty* display_properties = nullptr;
    str_t hovered_display_property_label = STR_LIT("");
    int   hovered_display_property_pop_idx = -1;

    // --- ASYNC TASKS HANDLES ---
    struct {
        task_system::ID backbone_computations = task_system::INVALID_ID;
        task_system::ID evaluate_full = task_system::INVALID_ID;
        task_system::ID evaluate_filt = task_system::INVALID_ID;
    } tasks;

    // --- ATOM SELECTION ---
    struct {
        SelectionGranularity granularity = SelectionGranularity::Atom;
        SingleSelectionSequence single_selection_sequence;

        md_bitfield_t selection_mask {};
        md_bitfield_t highlight_mask {};
        Selection* stored_selections = NULL;

        struct {
            struct {
                vec4_t visible = {0.0f, 0.0f, 1.0f,  0.3f};
                vec4_t hidden  = {0.0f, 0.0f, 0.25f, 0.4f};
            } selection;

            struct {
                vec4_t visible = {1.0f, 1.0f, 0.0f, 0.3f};
                vec4_t hidden  = {0.5f, 0.5f, 0.0f, 0.4f};
            } highlight;

            float saturation = 0.3f;
        } color;

        struct {
            char buf[256] = "";
            char error[256] = "";
            md_bitfield_t mask = {0};
            bool query_ok = false;
            bool query_invalid = true;
            bool show_window = false;
        } query;

        struct {
            md_bitfield_t mask = {0};
            SelectionGrowthMode mode = SelectionGrowthMode::CovalentBond;
            float extent = 1;
            bool mask_invalid = true;
            bool show_window = false;
        } grow;
    } selection;

    // --- GBUFFER ---
    // Primary gbuffer
    GBuffer gbuffer {};

    struct {
		immediate::Queue* world = nullptr;
		immediate::Queue* overlay = nullptr;
    } gfx;

    // --- PICKING ---
    PickingSurface picking_surface {};      // Surface for main viewport picking

    PickingRange   picking_range_atom {};   // Reserved picking range for atoms
    PickingRange   picking_range_bond {};   // Reserved picking range for bonds
    PickingRange   picking_range_dipole {}; // Reserved picking range for dipoles

    PickingHandler picking_handler {};      // Handler for managing picking interactions

    // --- ANIMATION ---
    struct {
        double frame = 0.0;  // double precision for long trajectories
        float fps = 10.f;
        float tension = 0.0f;
        InterpolationMode interpolation = InterpolationMode::CubicSpline;
        PlaybackMode mode = PlaybackMode::Stopped;

        bool show_window = true;
    } animation;

    // --- TIMELINE---
    struct {
        struct {
            bool enabled = false;
            double beg_frame = 0;
            double end_frame = 1;

            struct {
                bool enabled = false;
                double extent_in_frames = 10;
            } temporal_window;

            uint64_t fingerprint = 0;
        } filter;

        struct {
            double beg_x = 0;
            double end_x = 1;
        } view_range;

        // Holds the timestamps for each frame
        md_array(float) x_values = 0;

        bool show_window = false;
    } timeline;

    // --- DISTRIBUTIONS ---
    struct {
        struct {
            bool enabled = false;
        } filter;
        bool show_window = false;
    } distributions;

    // --- VISUALS ---
    struct {
        struct {
            vec3_t color = {1, 1, 1};
            float intensity = 24.0f;
        } background;

        struct { 
            bool enabled = true;
            float intensity = 6.0f;
            float radius = 6.0f;
            float bias = 0.1f;
        } ssao;

#if EXPERIMENTAL_CONE_TRACED_AO == 1
        struct {
            bool enabled = true;
            float intensity = 1.0f;
            float step_scale = 1.0f;
        } cone_traced_ao;
#endif

        struct {
            bool enabled = true;
            postprocess_pipeline::Tonemapper tonemapper = postprocess_pipeline::Tonemapper_ACES;
            float exposure = 1.0f;
            float gamma = 2.2f;
        } tonemapping;

        struct {
            bool enabled = false;
            float focus_depth = 10.0f;
            float focus_scale = 10.0f;
        } dof;

        struct {
            bool enabled = true;
        } fxaa;

        struct {
            bool enabled = true;
            bool jitter = true;
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

    struct {
        bool enabled = false;
        vec4_t color = {0, 0, 0, 0.5f};
    } simulation_box;

    // --- REPRESENTATIONS ---
    struct {
        RepresentationInfo info = {};
        md_array(Representation) reps = 0;
        md_bitfield_t visibility_mask = {0};
        uint64_t visibility_mask_hash = 0;
        bool atom_visibility_mask_dirty = false;
        bool show_window = false;
        bool needs_update = false;
        bool advanced_mode = false;
    } representation;

    struct {
        bool recenter = false;
        bool fixate_orientation = false;

        // Extra transform applied to the target in its own centred frame, between fixating the
        // orientation and placing it at the centre of the cell:
        //     T = translate(cell_centre) * alignment_mat * R * translate(-com)
        // Identity by default. This is the hook for aligning the structure to something once it is
        // centred - principal axes, a chosen bond, a user supplied orientation - without having to
        // touch the transform itself. Owned by whoever sets it: recentering never writes to it, and
        // it survives a change of recenter target.
        mat4_t alignment_mat = MD_MAT4_IDENT_INIT;

        struct {
            char query[256] = "";
            char error[256] = "";
            bool valid = false;
            bool enabled = false;
            bool dynamic = false; // represents whether the query needs to be re-evaluated per frame.
            uint64_t version = 1;
            uint64_t evaluated_version = 0;
            uint64_t ir_fingerprint = 0;
            md_bitfield_t mask = { 0 };
        } recenter_query;

        // External state which user controls
        bool apply_pbc = false;
        bool unwrap_structures = false;
        bool recalc_bonds = false;

        // Manual selection mask for recentering / orientating, which atoms to consider for calculating the center of mass and principal axes
        // This is populated by assigning a user defined selection.
        uint64_t selection_version = 1;
        md_bitfield_t selection_mask = {0};

        struct {
            uint64_t target_version = 0;

			// Need to store the initial relative coordinate vectors and center of mass of reference frame for recentering and orientating the structure.
            vec4_t* rel_xyzw = nullptr;
            vec3_t  com = {};
        } initial_frame;

    } operations;

    struct {
        bool keep_representations = false;
        bool prefetch_frames = true;
    } settings;

    // Views onto the temporal attributes, plus the one array still owned here.
    //
    // 'stride' and 'count' are redundant with the attributes' own shape {frames, segments} and are
    // kept only so the existing consumers keep their loops. Whoever next touches those loops should
    // take the shape from the attribute and delete them - the point of the move was that a hand
    // rolled shape beside table owned storage is a second answer waiting to disagree.
    //
    // The fingerprints are gone for the same reason: md_attributes_version is the one answer to
    // "did this change", and a stamp next to the data can be forgotten by whoever writes it.
    struct {
        struct {
            size_t stride = 0;                          // == shape[1] of 'backbone/secondary_structure'
            size_t count = 0;                           // == its element count
            md_secondary_structure_t* data = nullptr;   // VIEW into sys.attributes, not owned
        } secondary_structure;
        struct {
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
    } trajectory_data;

    struct {
        vec4_t point_color      = {1,0,0,0.8f};
        vec4_t line_color       = {0,0,0,0.6f};
        vec4_t triangle_color   = {0.55f,0.55f,1.0f,0.5f};
        vec4_t text_color       = {1,1,1,1};
        vec4_t text_bg_color    = {0,0,0,0.5f};

        uint64_t text_hash;

        // A bit confusing and a bit of a hack,
        // But we want to preserve the ir while evaluating it (visualizing it etc)
        // So we only commit the 'new' ir to eval_ir upon starting evaluation
        md_script_ir_t*   ir = nullptr;
        md_script_ir_t*   eval_ir = nullptr;

        md_script_eval_t* full_eval = nullptr;
        md_script_eval_t* filt_eval = nullptr;
        md_script_vis_t vis = {};

        // Controls current subindex being visualized (if array)
        int sub_idx = -1;

        bool compile_ir = false;
        bool eval_init = false;
        bool evaluate_full = false;
        bool evaluate_filt = false;
        double time_since_last_change = 0.0;
        uint64_t ir_fingerprint = 0;
    } script;

    struct {
        bool show_window = false;
        int selected_atom_filter = 0;
        int selected_traj_filter = 0;
        int selected_file_format = 0;
        struct {
            char buf[256] = "";
            char error[256] = "";
            md_bitfield_t mask = {0};
            bool is_valid = false;
            bool is_dynamic = false;
            bool requires_evaluation = false;
        } query;
    } structure_export;

    bool show_script_window = true;
    bool show_debug_window = false;
    bool show_property_export_window = false;

    TextEditor editor = {};
};

struct ViamdEventHandler : viamd::EventHandler {
    ApplicationState* state = nullptr;

    explicit ViamdEventHandler(ApplicationState* s) : state(s) {
        ASSERT(state);
        viamd::event_system_register_handler(*this);
    }

    void process_events(const viamd::Event* events, size_t num_events) final;
};

struct LoadDataPayload {
    ApplicationState* app_state;
    loader::LoaderState loader_state;
    str_t path_to_file;
};

struct PickingTooltipTextRequest {
    const ApplicationState& app;
    const PickingHit& hit;
    md_strb_t sb = {};
};

enum ViewFitRound {
    ViewFitRound_Highlight = 0,
    ViewFitRound_Selection,
    ViewFitRound_Visible,
};

// This is the supplied event payload for requesting a view fit.
// The idea is to construct a weighted point cloud of whatever is 'highlighted'
// And then the view fit system can decide how to best fit the view based on that point cloud.
// This is dispatched in rounds starting with the top most priority being Highlight, then Selection and lastly Visible.
// If the point cloud is empty for a certain mode, the next mode will be dispatched, until a non empty point cloud is found or all modes are exhausted.

struct ViewFitRequest {
    const ApplicationState& app;
	uint64_t surface_id = 0;        // Source of the view fit request, can be used to filter out irrelevant requests or for debugging.
	md_array(vec4_t) xyzw = nullptr;
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

    vec2_t mouse_local  = {};
    vec2_t surface_size = {};

    // Local coordinates
    vec2_t region_min = {};
    vec2_t region_max = {};
};

static inline const ImVec4& vec_cast(const vec4_t& v) { return *(const ImVec4*)(&v); }
static inline const vec4_t& vec_cast(const ImVec4& v) { return *(const vec4_t*)(&v); }
static inline const ImVec2& vec_cast(const vec2_t& v) { return *(const ImVec2*)(&v); }
static inline const vec2_t& vec_cast(const ImVec2& v) { return *(const vec2_t*)(&v); }

static inline ImVec4& vec_cast(vec4_t& v) { return *(ImVec4*)(&v); }
static inline vec4_t& vec_cast(ImVec4& v) { return *(vec4_t*)(&v); }
static inline ImVec2& vec_cast(vec2_t& v) { return *(ImVec2*)(&v); }
static inline vec2_t& vec_cast(ImVec2& v) { return *(vec2_t*)(&v); }

static inline void modify_field(md_bitfield_t* bf, const md_bitfield_t* mask, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        md_bitfield_or_inplace(bf, mask);
        break;
    case SelectionOperator::And:
        md_bitfield_and_inplace(bf, mask);
        break;
    case SelectionOperator::AndNot:
        md_bitfield_andnot_inplace(bf, mask);
        break;
    case SelectionOperator::Set:
        md_bitfield_copy(bf, mask);
        break;
    default:
        ASSERT(false);
    }
}

static inline void modify_field(md_bitfield_t* bf, md_urange_t range, SelectionOperator op) {
    switch(op) {
    case SelectionOperator::Or:
        md_bitfield_set_range(bf, range.beg, range.end);
        break;
    case SelectionOperator::And:
        md_bitfield_clear_range(bf, 0, range.beg);
        md_bitfield_clear_range(bf, range.end, bf->end_bit);
        break;
    case SelectionOperator::AndNot:
        md_bitfield_clear_range(bf, range.beg, range.end);
        break;
    case SelectionOperator::Set:
        md_bitfield_clear(bf);
        md_bitfield_set_range(bf, range.beg, range.end);
        break;
    case SelectionOperator::Clear:
        md_bitfield_clear_range(bf, range.beg, range.end);
        break;
    default:
        ASSERT(false);
    }
}

static inline void grow_mask_by_selection_granularity(md_bitfield_t* mask, SelectionGranularity granularity, const md_system_t& sys) {
    ASSERT(mask);
    switch(granularity) {
    case SelectionGranularity::Atom:
        break;
    case SelectionGranularity::Component:
        for (size_t i = 0; i < md_system_component_count(&sys); ++i) {
            md_urange_t range = md_system_component_atom_range(&sys, i);
            if (md_bitfield_popcount_range(mask, range.beg, range.end)) {
                md_bitfield_set_range(mask, range.beg, range.end);
            }
        }
        break;
    case SelectionGranularity::Instance:
        for (size_t i = 0; i < md_system_instance_count(&sys); ++i) {
            md_urange_t range = md_system_instance_atom_range(&sys, i);
            if (md_bitfield_popcount_range(mask, range.beg, range.end)) {
                md_bitfield_set_range(mask, range.beg, range.end);
            }
        }
        break;
    default:
        ASSERT(false);
    }
}

static inline void single_selection_sequence_clear(SingleSelectionSequence* seq) {
    ASSERT(seq);
    MEMSET(seq->idx, -1, sizeof(seq->idx));
}

static inline void single_selection_sequence_push_idx(SingleSelectionSequence* seq, int32_t idx) {
    ASSERT(seq);
    for (size_t i = 0; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) {
            seq->idx[i] = idx;
            break;
        }
    }
}

static inline void single_selection_sequence_pop_idx(SingleSelectionSequence* seq, int32_t idx) {
    ASSERT(seq);
    for (size_t i = 0; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == idx) {
            for (size_t j = i; j < ARRAY_SIZE(seq->idx) - 1; ++j) {
                seq->idx[j] = seq->idx[j+1];
            }
            seq->idx[ARRAY_SIZE(seq->idx)-1] = -1;
            break;
        }
    }
}

static inline void single_selection_sequence_pop_back(SingleSelectionSequence* seq) {
    ASSERT(seq);
    size_t i = 0;
    for (; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    if (i > 0) {
        seq->idx[i-1] = -1;
    }
}

static inline int32_t single_selection_sequence_last(const SingleSelectionSequence* seq) {
    ASSERT(seq);
    size_t i = 0;
    for (; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    if (i > 0) {
        return seq->idx[i-1];
    }
    return -1;
}

static inline size_t single_selection_sequence_count(const SingleSelectionSequence* seq) {
    size_t i = 0;
    for (; i < ARRAY_SIZE(seq->idx); ++i) {
        if (seq->idx[i] == -1) break;
    }
    return i;
}

static inline uint64_t generate_fingerprint() {
    return (uint64_t)md_time_now();
}

void draw_picking_tooltip_window(const PickingHit& hit, const ApplicationState& state);

//void extract_picking_data(PickingData& out_picking, GBuffer& gbuffer, const vec2_t& coord, const mat4_t& inv_MVP);

void interrupt_async_tasks(ApplicationState* state);

// Dataset loading
bool load_data_from_file(ApplicationState* app, str_t filepath, const loader::LoaderState& load_state);
void init_system_data(ApplicationState* app);
void init_trajectory_data(ApplicationState* app);

// Frame cache operations
void clear_system_frame_cache(ApplicationState* app);

// Interpolate state
void interpolate_system_state(ApplicationState* app);

// Workspace
void load_workspace(ApplicationState* app, str_t file);
void save_workspace(ApplicationState* app, str_t file);

// Selections
Selection* create_selection(ApplicationState* app, str_t name, md_bitfield_t* bf = 0);
void remove_selection(ApplicationState* app, size_t idx);
void remove_all_selections(ApplicationState* app);

// Representations
Representation* create_representation(ApplicationState* app, RepresentationType type = RepresentationType::SpaceFill, ColorMapping color_mapping = ColorMapping::Type, str_t filter = STR_LIT("all"));
Representation* clone_representation(ApplicationState* app, const Representation& rep);
void remove_representation(ApplicationState* app, size_t idx);
void update_representation(ApplicationState* app, Representation* rep);
void update_representation_info(ApplicationState* app);
void update_all_representations(ApplicationState* app);
bool representation_uses_atom_colors(const Representation& rep);

void flag_representation_as_dirty(Representation* rep);
void flag_all_representations_as_dirty(ApplicationState* app);

void remove_all_representations(ApplicationState* app);
void create_default_representations(ApplicationState* app);
void recompute_atom_visibility_mask(ApplicationState* app);

// Enumerates the complete dipole groups in the system's attribute table, in name order.
// Returns the total number found and writes at most cap, so pass cap 0 to count. Cheap enough to
// call per frame: it is a prefix query over a small sorted array, and there is nothing to cache
// or invalidate because the attribute table IS the state.
// Groups missing either half, or whose shapes are not one 3-component value, are skipped.
size_t dipole_moments_gather(DipoleMoment out[], size_t cap, const md_system_t& sys);

// Turns a dipole group name into something presentable: "ground_state" -> "Ground State".
// Writes at most cap-1 characters plus a terminator, returns the length written.
int dipole_label_pretty(char* buf, size_t cap, str_t group);

// The same, plus a 1 based element number when the group holds more than one moment:
// "Electric Transition 3". A single moment gets no number, because there is nothing to tell apart.
int dipole_entry_label(char* buf, size_t cap, const DipoleMoment& dipole);

// Per atom scalar fields are not a list anybody builds or keeps: they are whatever the system's
// attribute table holds under atom/, and these read it directly. Nothing to rebuild when the data
// changes, and nothing that can go stale against it.

// Ids of the attributes under atom/ which are per atom scalar fields, in path order.
// Returns the total number found and writes at most cap, so pass cap 0 to count.
// An attribute qualifies when its values are single components and its LAST index axis is this
// system's atom count; an optional leading axis is the variant axis. Everything else under atom/ -
// a position, a velocity, anything several components wide - is data this colouring cannot express
// and is skipped rather than mangled. Cheap enough to call per frame.
size_t atom_property_query(md_attribute_id_t out_ids[], size_t cap, const md_system_t& sys);

// What to show for it: the attribute's label, or its leaf path segment when it has none. A view
// into the table's storage, null terminated, so it can go straight to ImGui.
str_t atom_property_label(const md_attribute_t* attr);

// Number of variants: the leading index axis, 1 when the field has none.
int atom_property_variant_count(const md_attribute_t* attr);

// Span of the values, over EVERY variant so a colour ramp does not jump as the variant changes.
// Derived rather than stored: it belongs to whoever is drawing the ramp, not to the table.
// Scans the whole attribute, so call it when the selection changes, not per frame.
bool atom_property_value_range(float* out_min, float* out_max, const md_attribute_t* attr);

// Builds the per dataset GPU data - the uploaded GTO basis and the atom buffer - from the system's
// own basis/ attributes, and grows the device coefficient scratch to fit. Returns false when the
// system publishes no basis, which is the normal case rather than a failure.
bool system_gpu_data_update(ApplicationState* state, double cutoff);
void system_gpu_data_free(ApplicationState* state);

// One orbital's AO coefficient row, out of the attribute at coefficient_path narrowed by 'slice',
// in the double precision md_gto asks for. Allocated from temp; NULL when the path is absent or the
// slice does not narrow it to a single row. out_num_ao receives the row length.
//
// md_attribute_slice_1(m) for a {M,A} matrix of molecular orbitals, md_attribute_slice_2(s,l) for
// the {S,L,A} natural transition orbitals - the caller says which row, and nothing here has to know
// which of the two it was addressing.
double* orbital_coefficients_extract(size_t* out_num_ao, md_temp_scope_t temp, const md_system_t& sys, str_t coefficient_path, const md_attribute_slice_t* slice);

// Evaluates one molecular orbital into a 3D texture, sourcing everything it needs from the system:
// the GTO basis and the AO coefficients from the attribute table, the atom positions from the state.
// coefficient_path is a rank 2 {M,A} attribute of doubles, "orbital/alpha/coefficient" or its beta
// sibling; orbital_idx selects the row.
//
// The caller owns the grid, the destination texture and the evaluation parameters, which is why no
// registry of "fields" is needed to sit between: what varies per request is the request, and what
// is data is data.
bool orbital_evaluate_gl(uint32_t vol_tex, const md_grid_t& grid, const md_system_t& sys,
                         const vec3_t* atom_pos, size_t num_atom_pos,
                         str_t coefficient_path, const md_attribute_slice_t* slice, md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff);

// Positions for a system's BASIS atoms, in the bohr and interleaved layout md_gto works in, taken
// out of `state` through the QM -> system map. Allocated from temp; null with a zero count when the
// system publishes no basis or the state cannot serve it.
//
// This is what an evaluation takes INSTEAD of a state. Which geometry an orbital is drawn at - the
// current frame, a reference, an optimisation step, a geometry that never was - is the call site's
// decision, and an evaluation is in no position to make it.
vec3_t* basis_atom_positions_extract(size_t* out_count, md_temp_scope_t temp, const md_system_t& sys, const md_system_state_t& state);

// Blocks until every queued volume readback has landed in its GL texture. Call before destroying
// the device scratch, the pools or a GL texture one of them writes to, and wherever a caller needs
// the texture's contents NOW rather than next frame - an export, say. A no-op without the GPU
// backend, where evaluation writes the texture directly.
void gpu_volume_jobs_drain(ApplicationState* state);


// The QM ATOM DOMAIN: the atoms a quantum chemistry calculation covered, in its own order and at
// its own geometry. NOT the system's atoms - a calculation can cover part of a loaded system, so
// the two spaces differ in length and in order, and system_index is the only bridge.
//
// Published under qm/atom/. basis/shell/atom_index is an index INTO this space - the basis is one
// thing defined over these atoms, not the other way round - and 'count' is what "how many atoms
// does the QM cover" means; there is nothing else to ask.
//
// The pointers are into the table's own storage and are invalidated by the next attribute create
// or remove - which, for a loaded system, means the next load.
struct QmAtoms {
    size_t          count         = 0;
    const uint8_t*  atomic_number = nullptr;  // [N]
    const dvec3_t*  coordinate    = nullptr;  // [N], Angstrom, the geometry the calculation ran at
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
bool es_orbital_frontier(OrbitalFrontier* out, const md_system_t& sys, str_t occupation_path);

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
double es_electron_count(const md_system_t& sys, str_t occupation_path);

// Excited states carrying natural transition orbitals, from the weight table's leading axis.
size_t es_excited_state_count(const md_system_t& sys);

// One excited state's NTO weights, largest first, stopping at the first below 'cutoff'. The ragged
// lambda axis is zero padded, so the cutoff ends the run in exactly the place the pairs do.
// Returns how many were written.
size_t es_nto_lambdas(double* out_values, size_t cap, const md_system_t& sys, size_t state_idx, double cutoff);

// One AO x AO density matrix, out of the attribute at density_path narrowed by 'slice', in the
// double precision md_gto asks for. Allocated from temp; NULL when the path is absent or the slice
// does not narrow it to a square matrix. Pass a zero slice (or NULL) for an attribute which is
// already {A,A}, md_attribute_slice_1(state) for one indexed by state. out_dim receives A.
double* density_matrix_extract(size_t* out_dim, md_temp_scope_t temp, const md_system_t& sys, str_t density_path, const md_attribute_slice_t* slice);

// Evaluates one density into a 3D texture, from the same two sources as orbital_evaluate_gl: the
// basis/ attributes and the state's atom positions. density_path plus 'slice' name the matrix -
// "orbital/total/density" with no slice, "vlx/rsp/transition_density/attachment" with the excited
// state index - and whether that matrix is stored or computed on demand is the table's business.
bool density_evaluate_gl(uint32_t vol_tex, const md_grid_t& grid, const md_system_t& sys,
                         const vec3_t* atom_pos, size_t num_atom_pos,
                         str_t density_path, const md_attribute_slice_t* slice, md_gto_op_t op);

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
    inline const str_t alpha_coefficient   = STR_LIT("orbital/alpha/coefficient");
    inline const str_t beta_coefficient    = STR_LIT("orbital/beta/coefficient");

    inline const str_t alpha_energy        = STR_LIT("orbital/alpha/energy");
    inline const str_t beta_energy         = STR_LIT("orbital/beta/energy");

    inline const str_t alpha_occupation    = STR_LIT("orbital/alpha/occupation");
    inline const str_t beta_occupation     = STR_LIT("orbital/beta/occupation");

    // The AO metric the coefficients and densities above are expressed against. Singular for a
    // basis that was stored spherically - see the note on md_vlx_publish_attributes.
    inline const str_t overlap             = STR_LIT("basis/overlap");

    inline const str_t alpha_density       = STR_LIT("orbital/alpha/density");
    inline const str_t beta_density        = STR_LIT("orbital/beta/density");
    inline const str_t total_density       = STR_LIT("orbital/total/density");
    inline const str_t difference_density  = STR_LIT("orbital/difference/density");

    inline const str_t nto_particle        = STR_LIT("vlx/rsp/nto/particle/coefficient");
    inline const str_t nto_hole            = STR_LIT("vlx/rsp/nto/hole/coefficient");

    // {S,Lmax}, the weight of each NTO pair, padded with zeros on the ragged lambda axis and
    // indexing the same space as the two coefficient attributes above.
    inline const str_t nto_lambda          = STR_LIT("vlx/rsp/nto/lambda");

    inline const str_t attachment_density  = STR_LIT("vlx/rsp/transition_density/attachment");
    inline const str_t detachment_density  = STR_LIT("vlx/rsp/transition_density/detachment");
    inline const str_t transition_diff     = STR_LIT("vlx/rsp/transition_density/difference");

    inline const str_t density_property    = STR_LIT("vlx/density_property");
}

// The only interpretation of use_magnitude there is. ABS is a MODIFIER bit on an accumulate op, not
// an op of its own, which is why it is or'd onto SET rather than replacing it.
inline md_gto_op_t es_gto_op(bool use_magnitude) {
    return MD_GTO_OP_SET | (use_magnitude ? MD_GTO_OP_ABS : 0);
}

// Which density path a spin selection names. The four are siblings in the table precisely so that
// this is a lookup and not four cases of arithmetic at the point of use.
inline str_t es_electron_density_path(ElectronicStructureSpin spin) {
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

// The entry points a representation calls: name the attribute, get a filled 3D texture. One
// function per kind of thing being evaluated, with the choice of backend inside - so a caller never
// asks which device is available, and both backends share the basis, the atom gather and the map.
//
// With the GPU backend the texture is filled from the frame loop rather than before the call
// returns; gpu_volume_jobs_drain() forces it where the contents are needed immediately.
bool orbital_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, str_t coefficient_path,
                      const md_attribute_slice_t* slice, md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff);
bool density_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, str_t density_path,
                      const md_attribute_slice_t* slice, md_gto_op_t op);

#if MD_ENABLE_GPU
// The first half of density_evaluate: evaluates into the DEVICE scratch volume (ApplicationState::
// gpu_volume) and leaves it there, for a consumer which records another kernel over it on the same
// stream instead of reading it back. False when the GPU backend is unavailable, in which case there
// is no scratch volume to run anything over and the caller wants the GL path instead.
bool density_evaluate_to_gpu_volume(ApplicationState* state, const md_grid_t& grid, str_t density_path,
                                    const md_attribute_slice_t* slice, md_gto_op_t op);
#endif

// Points a representation at an attribute and seeds its drawing range from that attribute's own
// span. Selecting a field and choosing the range to draw it over are one action the first time and
// separate afterwards, which is why the range is seeded here and never recomputed behind the user.
void atom_property_select(AtomicPropertyRepresentation* prop, md_attribute_id_t key, const md_system_t& sys);

// Recentering operations (low level)

void recenter_mark_query_dirty(ApplicationState* app);
void recenter_mark_selection_dirty(ApplicationState* app);
const md_bitfield_t& recenter_get_active_target_mask(const ApplicationState* app);
uint64_t recenter_get_active_target_version(const ApplicationState* app);
bool recenter_update_query_mask(ApplicationState* app);
void recenter_update(ApplicationState* app);

// Update the required initial frame data for the recentering target (if needed)
void recenter_update_target_data(ApplicationState* app);
void recenter_calculate_transform(float M[4][4], const ApplicationState* app);

// Picking

void picking_handler_new_frame(PickingHandler* handler);
PickingSpace* picking_handler_current_space(PickingHandler* handler);
const PickingSpace* picking_handler_find_space(const PickingHandler& handler, uint32_t submitted_frame_idx);

// Reserves a range within the picking space for a specific domain (atoms, bonds, etc).
// Returns true if the range was successfully reserved, false if there was not enough space. If successful, out_range will be filled with the reserved range.
// out_range is optional and can be null if the caller does not need the details of the reserved range (e.g. just needs to know if the reservation was successful or not).
bool picking_range_reserve(PickingRange* out_range, PickingSpace* space, PickingDomainID domain, size_t count);

void picking_surface_init(PickingSurface* surface, PickingSourceID source);
void picking_surface_free(PickingSurface* surface);

// Submits a picking readback request for the given surface. The readback will be performed asynchronously, and the result can be polled using picking_surface_poll_hit.
bool picking_surface_submit_readback(
    PickingSurface* surface,
    uint32_t fbo,
    uint32_t width,
    uint32_t height,
    uint32_t submitted_frame_idx,
    vec2_t surface_coord,
    vec2_t screen_coord,
    const mat4_t& inv_mvp
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
InteractionSurfaceState interaction_surface(InteractionSurfaceID id, const vec2_t& size, InteractionSurfaceFlags flags = InteractionSurfaceFlags_None);

struct InteractionSurfaceHitArgs {
    PickingSurface* picking_surface;
    const PickingHandler& picking_handler;
    uint32_t fbo;
    uint32_t width;
    uint32_t height;
    const mat4_t& clip_to_world;
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

    vec2_t mouse_local = {};
    vec2_t surface_size = {};
    vec2_t region_min = {};
    vec2_t region_max = {};

    mat4_t clip_to_world = MD_MAT4_IDENT_INIT;
    mat4_t world_to_clip = MD_MAT4_IDENT_INIT;

    PickingHit hit = {};
};

void interaction_surface_event_extract(InteractionSurfaceEvent* out_event, const InteractionSurfaceState& state, const PickingHit& hit = {});

// Helper function for projecting world coordinates to surface coordinates, used for interaction surfaces and picking.
static inline vec2_t world_to_surface_project(
    const vec3_t& world_coord,
    const mat4_t& world_to_clip,
    const vec2_t& surface_size
) {
    const vec4_t c = world_to_clip * vec4_set(world_coord.x, world_coord.y, world_coord.z, 1.0f);
    vec2_t out_coord = {
        ( c.x / c.w * 0.5f + 0.5f) * surface_size.x,
        (-c.y / c.w * 0.5f + 0.5f) * surface_size.y,
    };
    return out_coord;
}

void point_set_region_mask_compute(
    md_bitfield_t* mask,
    const float x[],
    const float y[],
    const float z[],
    size_t count,
    const md_bitfield_t* candidate_mask,
    const mat4_t& world_to_clip,
    const vec2_t& region_min,
    const vec2_t& region_max,
    const vec2_t& surface_size
);

// File Queue
bool file_queue_empty(const FileQueue* queue);
bool file_queue_full(const FileQueue* queue);
void file_queue_push(FileQueue* queue, str_t path, FileFlags flags = FileFlags_None);

FileQueue::Entry file_queue_front(const FileQueue* queue);
FileQueue::Entry file_queue_pop(FileQueue* queue);

void file_queue_process(ApplicationState* state);

// view
void reset_view(ViewTransform* transform, const md_system_state_t& state, const md_bitfield_t* mask = nullptr);

// Script visualization
void script_visualize_payload(ApplicationState* state, const md_script_vis_payload_o* payload, int subidx, md_script_vis_flags_t flags = 0);
void script_visualize_str(ApplicationState* state, str_t str, md_script_vis_flags_t flags = 0);
void script_set_hovered_property(ApplicationState* state, str_t label, int population_idx = -1);
