#pragma once

#include <core/md_str.h>
#include <core/md_os.h>
#include <core/md_lru_cache.inl>
#include <core/md_hash.h>
#include <core/md_str_builder.h>

#include <md_system.h>
#include <md_trajectory.h>
#include <md_script.h>
#include <md_gl.h>

#include <app/IconsFontAwesome6.h>
#include <app/application.h>
#include <gfx/camera.h>
#include <gfx/camera_utils.h>
#include <gfx/view_param.h>
#include <gfx/postprocessing_utils.h>

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

constexpr PickingDomainID PickingDomain_Atom = HASH_STR_LIT64("picking domain atom");
constexpr PickingDomainID PickingDomain_Bond = HASH_STR_LIT64("picking domain bond");

constexpr uint64_t interaction_surface_main = HASH_STR_LIT64("interaction surface main"); // This is the main interaction surface which corresponds to the main interaction window, but we want to keep it separate from the picking source and domain ids as we may want to have different picking sources/domains for different interaction surfaces in the future

enum class PlaybackMode { Stopped, Playing };
enum class SelectionGranularity { Atom, Component, Instance };
enum class SelectionOperator { Or, And, AndNot, Set, Clear };
enum class SelectionGrowthMode { CovalentBond, Radial };
enum class CameraMode { Perspective, Orthographic };

enum class InterpolationMode {
    Nearest,
    Linear,
    CubicSpline,
    Count
};

static const char* interpolation_mode_str[(int)InterpolationMode::Count] = {
    "Nearest",
    "Linear",
    "Cubic Spline",
};

// These bits are a compressed form of flags which are passed onto rendering as the rendering only supports 8-bits
enum AtomBit_ {
    AtomBit_Highlighted = 1,
    AtomBit_Selected    = 2,
    AtomBit_Visible     = 4,
};

enum class RepresentationType {
    SpaceFill = MD_GL_REP_SPACE_FILL,
    Licorice = MD_GL_REP_LICORICE,
    BallAndStick = MD_GL_REP_BALL_AND_STICK,
    Ribbons = MD_GL_REP_RIBBONS,
    Cartoon = MD_GL_REP_CARTOON,
    ElectronicStructure,
//    DipoleMoment,
    Count
};

static const char* representation_type_str[(int)RepresentationType::Count] = {
    "Spacefill",
    "Licorice",
    "Ball And Stick",
    "Ribbons",
    "Cartoon",
    "Electronic Structure",
//    "Dipole Moment"
};

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

static const char* color_mapping_str[(int)ColorMapping::Count] = {
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

enum class ElectronicStructureType {
    MolecularOrbital,
    MolecularOrbitalDensity,
    NaturalTransitionOrbitalParticle,
    NaturalTransitionOrbitalHole,
    NaturalTransitionOrbitalDensityParticle,
    NaturalTransitionOrbitalDensityHole,
    AttachmentDensity,
    DetachmentDensity,
    ElectronDensity,
    Count
};

static const char* electronic_structure_type_str[(int)ElectronicStructureType::Count] = {
    (const char*)u8"Molecular Orbital (Ψ)",
    (const char*)u8"Molecular Orbital Density (Ψ²)",
    (const char*)u8"Natural Transition Orbital (NTO) Particle",
    (const char*)u8"Natural Transition Orbital (NTO) Hole",
    (const char*)u8"Natural Transition Orbital Density (NTO²) Particle",
    (const char*)u8"Natural Transition Orbital Density (NTO²) Hole",
    (const char*)u8"Attachment Density (A)",
    (const char*)u8"Detachment Density (D)",
    (const char*)u8"Electron Density (ρ)",
};

enum MolBit_ {
    MolBit_DirtyPosition            = 1u << 0,
    MolBit_DirtyRadius              = 1u << 1,
    MolBit_DirtySecondaryStructure  = 1u << 2,
    MolBit_DirtyFlags               = 1u << 3,
    MolBit_DirtyBonds               = 1u << 4,
    MolBit_ClearVelocity            = 1u << 5,
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

enum class VolumeResolution {
    Low,
    Mid,
    High,
    Count,
};

static const char* volume_resolution_str[(int)VolumeResolution::Count] = {
    "Low",
    "Mid",
    "High",
};

enum class ScreenshotResolution {
    Window,
    FHD,
    QHD,
    UHD_4K,
    UHD_8K,
    Custom,
    Count,
};

static const char* screenshot_resolution_str[(int)ScreenshotResolution::Count] = {
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

static const char* bond_color_mode_str[(int)BondColorMode::Count] = {
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
    md_bitfield_t atom_mask{};
};

struct DipoleMoment {
    size_t num_dipoles = 0;
    str_t* label = nullptr;
    vec3_t* vec  = nullptr;
};

struct NaturalTransitionOrbitalLambda {
    size_t num_lambdas = 0;
    str_t* label = nullptr;
    double* value = nullptr;
};

struct NaturalTransitionOrbital {
    size_t num_orbitals = 0;
    str_t* label = nullptr;
    NaturalTransitionOrbitalLambda* lambda = nullptr;
};

struct MolecularOrbital {
    size_t homo_idx = 0;
    size_t lumo_idx = 0;
    size_t num_orbitals = 0;
    str_t* label = nullptr;
    double* occupation = nullptr;
    double* energy = nullptr;
};

struct AtomProperty {
    uint64_t key    = 0;
    str_t label     = { 0 };
    int num_idx     = 0;
    float value_min = 0;
    float value_max = 0;
};

// Struct to fill in for the different components
// Which provides information of what representations are available for the currently loaded datasets
struct RepresentationInfo {
    MolecularOrbital alpha;
    MolecularOrbital beta;
    NaturalTransitionOrbital nto;
    DipoleMoment electric_dipoles;
    DipoleMoment magnetic_dipoles;
    DipoleMoment velocity_dipoles;

    uint32_t electronic_structure_type_mask;

    md_array(AtomProperty) atom_properties = nullptr;

    md_allocator_i* alloc = nullptr;
};

struct Volume {
    mat4_t world_to_model   = {};   // Roto-translation into volume local axes, no scaling applied (preserves world length units)
    mat4_t texture_to_world = {};   // Texture space [0,1] to world coordinates
    vec3_t voxel_size  = {1,1,1};   // Size of each voxel in world units
    int dim[3] = {128, 128, 128};
    uint32_t tex_id = 0;
};

// Descriptor for handling iso surfaces
struct IsoDesc {
    bool enabled;
    size_t count;
    float values[8];
    vec4_t colors[8];
    float optical_densities[8];
};

struct EvalAtomProperty {
    uint64_t key = 0;
	int idx = 0;    // Represents the index if the property is multidimensional.
    size_t num_values = 0;
    float* dst_values = nullptr;
    bool output_written = false;
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

    ElectronicStructureType type = ElectronicStructureType::MolecularOrbital;
    int mo_idx = 0;
    int nto_idx = 0;
    int nto_lambda_idx = 0;

	uint64_t col_hash = 0;
	uint64_t vol_hash = 0;
    uint64_t tf_hash = 0;
};

struct AtomicPropertyRepresentation {
    int colormap = DEFAULT_COLORMAP;
    float range_beg = 0.0f;
    float range_end = 1.0f;
    bool  range_symmetric_zero = true; // Use a symmetric min and max value around zero
    int idx = 0;
	int sub_idx = 0;
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
};

// Event Payload when an electronic structure is to be evaluated
struct EvalElectronicStructure {
    const md_system_t* sys = 0;
    double frame = 0.0;
    Representation* rep = 0;
    uint32_t* atom_colors = 0;
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
    mat4_t clip_to_world = mat4_ident();
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
    mat4_t clip_to_world = mat4_ident();
};

struct ApplicationState {
    // --- APPLICATION ---
    application::Context app {};

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
        md_allocator_i*     sys_alloc = nullptr;

        md_gl_mol_t         gl_mol = {};
#if EXPERIMENTAL_GFX_API
        md_gfx_handle_t     gfx_structure = {};
#endif
        md_system_t         sys = {};

        // The actual interpolated state of the system.
        // Is only used for rendering and visualization of properties.
        md_system_state_t   state = {};

		mat4_t 			    unitcell_transform = mat4_ident();

        FrameCache          frame_cache;

        vec3_t              sys_aabb_min = {};
        vec3_t              sys_aabb_max = {};
        ViewTransform       default_view = {};

        bool                interpolate_system_state = false;
        uint32_t            dirty_gpu_buffers = 0;
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

    // --- FRAMEBUFFER ---
    GBuffer gbuffer {};

    // --- PICKING ---
    PickingSurface picking_surface {};      // Surface for main viewport picking

    PickingRange   picking_range_atom {};   // Reserved picking range for atoms
    PickingRange   picking_range_bond {};   // Reserved picking range for bonds

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
            postprocessing::Tonemapping tonemapper = postprocessing::Tonemapping_ACES;
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
        bool rotate   = false;

        bool apply_pbc = false;
        bool unwrap_structures = false;
        bool recalc_bonds = false;

        // For recentering / orientating, which atoms to consider for calculating the center of mass and principal axes
        md_bitfield_t target_mask = {0};

        struct {
            uint64_t hash = 0;

            // Need to store the initial frame position for recentering and orienting to work properly when applying on trajectories
            vec4_t* xyzw = nullptr;
            vec3_t  com = {};

            // Alignment matrix for orienting the structure based on principal axes. This is calculated based on the initial frame and applied to all frames for consistent orientation.
            mat4_t alignment_mat = mat4_ident();
        } initial_frame;

    } operations;

    struct {
        bool keep_representations = false;
        bool prefetch_frames = true;
    } settings;

    struct {
        md_array(md_gl_secondary_structure_t) secondary_structure = nullptr;
    } interpolated_properties;

    struct { 
        struct {
            size_t stride = 0; // = mol.backbone.count. Multiply frame idx with this to get the data
            size_t count = 0;  // = mol.backbone.count * num_frames. Defines the end of the data for assertions
            md_secondary_structure_t* data = nullptr;
            uint64_t fingerprint = 0;
        } secondary_structure;
        struct {
            size_t stride = 0; // = mol.backbone.count. Multiply frame idx with this to get the data
            size_t count = 0;  // = mol.backbone.count * num_frames. Defines the end of the data for assertions
            md_backbone_angles_t* data = nullptr;
            uint64_t fingerprint = 0;
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
    loader::State loader_state;
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
    for (size_t i = 0; i < ARRAY_SIZE(seq->idx); ++i) {
        seq->idx[i] = -1;
    }
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
bool load_data_from_file(ApplicationState* state, str_t filepath, const loader::State& load_state);
void init_system_data(ApplicationState* state);
void init_trajectory_data(ApplicationState* state);

// Frame cache operations
void clear_system_frame_cache(ApplicationState* state);

// Interpolate state
void interpolate_system_state(ApplicationState* state);

// Workspace
void load_workspace(ApplicationState* state, str_t file);
void save_workspace(ApplicationState* state, str_t file);

// Selections
Selection* create_selection(ApplicationState* state, str_t name, md_bitfield_t* bf = 0);
void remove_selection(ApplicationState* state, int idx);
void remove_all_selections(ApplicationState* state);

// Representations
Representation* create_representation(ApplicationState* state, RepresentationType type = RepresentationType::SpaceFill, ColorMapping color_mapping = ColorMapping::Type, str_t filter = STR_LIT("all"));
Representation* clone_representation(ApplicationState* state, const Representation& rep);
void remove_representation(ApplicationState* state, int idx);
void update_representation(ApplicationState* state, Representation* rep);
void update_representation_info(ApplicationState* state);
void update_all_representations(ApplicationState* state);
bool representation_uses_atom_colors(const Representation& rep);

void flag_representation_as_dirty(Representation* rep);
void flag_all_representations_as_dirty(ApplicationState* state);

void remove_all_representations(ApplicationState* state);
void create_default_representations(ApplicationState* state);
void recompute_atom_visibility_mask(ApplicationState* state);

// Recentering operations (low level)

// Update the required initial frame data for the recentering target (if needed)
void recenter_update_target_data(ApplicationState* state);
void recenter_calculate_transform(float M[4][4], const ApplicationState* state);

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
    const ViewTransform& reset_transform = {};
};

// Uses the interaction surface state (e.g. mouse position, region selection) to calculate a view transform based on the provided camera and trackball parameters.
// Modifies the target view transform in place.
// Reset transform supplied in args represents the *reset target* transform which is optionally applied when the user double clicks the surface.
void interaction_surface_view_transform_apply(ViewTransform* target, const InteractionSurfaceState& state, const InteractionSurfaceViewTransformArgs& args);

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

    InteractionSurfaceEventKind kind = InteractionSurfaceEventKind::Hover;
    InteractionSelectionMode selection_mode = InteractionSelectionMode::None;
    InteractionSurfaceEventPhase region_phase = InteractionSurfaceEventPhase::None;

    vec2_t mouse_local = {};
    vec2_t surface_size = {};
    vec2_t region_min = {};
    vec2_t region_max = {};

    mat4_t clip_to_world = mat4_ident();
    mat4_t world_to_clip = mat4_ident();

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
void reset_view(ViewTransform* transform, const md_system_t& sys, const md_bitfield_t* mask = nullptr);

// Script visualization
void script_visualize_payload(ApplicationState* state, const md_script_vis_payload_o* payload, int subidx, md_script_vis_flags_t flags = 0);
void script_visualize_str(ApplicationState* state, str_t str, md_script_vis_flags_t flags = 0);
void script_set_hovered_property(ApplicationState* state, str_t label, int population_idx = -1);