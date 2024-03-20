#pragma once

#include <core/md_str.h>
#include <core/md_os.h>
#include <md_molecule.h>
#include <md_trajectory.h>
#include <md_script.h>
#include <md_gl.h>

#include <app/application.h>
#include <gfx/camera.h>
#include <gfx/camera_utils.h>
#include <gfx/view_param.h>
#include <gfx/postprocessing_utils.h>
#include <task_system.h>

#include <implot.h>

#include <stdint.h>
#include <stddef.h>

#define JITTER_SEQUENCE_SIZE 32

// For cpu profiling
#define PUSH_CPU_SECTION(lbl) {};
#define POP_CPU_SECTION() {};

// For gpu profiling
#define PUSH_GPU_SECTION(lbl) { if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl); }
#define POP_GPU_SECTION()     { if (glPopDebugGroup) glPopDebugGroup(); }

enum class PlaybackMode { Stopped, Playing };
enum class InterpolationMode { Nearest, Linear, CubicSpline };
enum class SelectionLevel { Atom, Residue, Chain };
enum class SelectionOperator { Or, And, AndNot, Set, Clear };
enum class SelectionGrowth { CovalentBond, Radial };
enum class TrackingMode { Absolute, Relative };
enum class CameraMode { Perspective, Orthographic };

enum class RepresentationType {
    SpaceFill = MD_GL_REP_SPACE_FILL,
    Licorice = MD_GL_REP_LICORICE,
    BallAndStick = MD_GL_REP_BALL_AND_STICK,
    Ribbons = MD_GL_REP_RIBBONS,
    Cartoon = MD_GL_REP_CARTOON,
    Orbital,
    DipoleMoment,
    Count
};

static const char* representation_type_str[(int)RepresentationType::Count] = {
    "Spacefill",
    "Licorice",
    "Ball And Stick",
    "Ribbons",
    "Cartoon",
    "Orbital",
    "Dipole Moment"
};

enum class ColorMapping {
    Uniform,
    Cpk,
    AtomLabel,
    AtomIndex,
    ResName,
    ResId,
    ChainId,
    ChainIndex,
    SecondaryStructure,
    Property,
    Count
};

static const char* color_mapping_str[(int)ColorMapping::Count] = {
    "Uniform Color",
    "CPK",
    "Atom Label",
    "Atom Idx"
    "Res Id",
    "Res Idx",
    "Chain Id",
    "Chain Idx",
    "Secondary Structure",
    "Property"
};

enum class OrbitalType {
    Psi,
    PsiSquared,
    Nto,
    Count
};

static const char* orbital_type_str[(int)OrbitalType::Count] = {
    (const char*)u8"Orbital (Ψ)",
    (const char*)u8"Density (Ψ²)",
    "Natural Transition Orbital",
};

enum MolBit_ {
    MolBit_DirtyPosition            = 0x01,
    MolBit_DirtyRadius              = 0x02,
    MolBit_DirtySecondaryStructure  = 0x04,
    MolBit_DirtyFlags               = 0x08,
    MolBit_DirtyBonds               = 0x10,
};

struct DisplayProperty;

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
    int  loader_idx = -1;
    int  atom_format_idx = -1;
};

// Hint flags for operations to be performed for the files
enum {
    FileFlags_None = 0,
    FileFlags_ShowDialogue = 1,
    FileFlags_CoarseGrained = 2,
    FileFlags_DisableCacheWrite = 4,
};

typedef uint32_t FileFlags;

struct FileQueue {
    struct Entry {
        str_t path;
        FileFlags flags;
        int prio;
    };
    Entry arr[8] = {};
    int head = 0;
    int tail = 0;

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

struct AtomElementMapping {
    char lbl[31] = "";
    md_element_t elem = 0;
};

// We use this to represent a single entity within the loaded system, e.g. a residue type
struct DatasetItem {
    char label[32] = "";
    char query[32] = "";
    uint32_t count = 0;
    float fraction = 0;
};

struct DipoleMoment {
    str_t  label;
    vec3_t vector;
};

struct MolecularOrbital {
    int idx;
    float occupation;
    float energy;
};

// Struct to fill in for the different components
// Which provides information of what representations are available for the currently loaded datasets
struct RepresentationInfo {
    int mo_homo_idx = 0;
    int mo_lumo_idx = 0;

    md_array(MolecularOrbital) molecular_orbitals = 0;
    md_array(DipoleMoment) dipole_moments = 0;

    md_allocator_i* alloc;
};

// Event Payload when an orbital is to be evaluated
struct ComputeOrbital {
    // Input information
    OrbitalType type = OrbitalType::Psi;
    int orbital_idx = 0;

    // Output information
    bool output_written = false;
    mat4_t mdl_mat = {};
    mat4_t tex_mat = {};
    vec3_t voxel_spacing = {};
    uint32_t *dst_texture = 0;
};

struct RepresentationVolume {
    uint32_t vol_tex = 0;
    mat4_t mdl_mat = mat4_ident();
    mat4_t tex_mat = mat4_ident();
    vec3_t voxel_spacing = {};

    float density_scale = 1.0f;

    struct {
        bool enabled = true;
        int count = 2;
        float  values[8] = {0.05f, -0.05};
        vec4_t colors[8] = {{215.f/255.f,25.f/255.f,28.f/255.f,0.75f}, {44.f/255.f,123.f/255.f,182.f/255.f,0.75f}};
    } iso;

    struct {
        bool enabled = false;
        uint32_t tf_tex = 0;
        ImPlotColormap colormap = ImPlotColormap_Plasma;
    } dvr;
};

struct Representation {
    char name[64] = "rep";
    char filt[256] = "all";
    char filt_error[256] = "";

    RepresentationType type = RepresentationType::SpaceFill;
    ColorMapping color_mapping = ColorMapping::Cpk;
    md_bitfield_t atom_mask = {};
    md_gl_representation_t md_rep = {};
#if EXPERIMENTAL_GFX_API
    md_gfx_handle_t gfx_rep = {};
#endif

    bool enabled = true;
    bool type_is_valid = false;
    bool filt_is_dirty = true;
    bool filt_is_valid = false;
    bool filt_is_dynamic = false;
    bool dynamic_evaluation = false;

    // User defined color used in uniform mode
    vec4_t uniform_color = {1.0f, 1.0f, 1.0f, 1.0f};

    // scaling parameter (radius, width, height, etc depending on type)
    vec4_t scale = {1.0f, 1.0f, 1.0f, 1.0f};

    struct {
        RepresentationVolume vol = {};
        OrbitalType type = OrbitalType::Psi;
        int orbital_idx = 0;
		uint64_t vol_hash = 0;
        uint64_t tf_hash = 0;
    } orbital;

    struct {
        ImPlotColormap color_map = ImPlotColormap_Plasma;
        float map_beg = 0.0f;
        float map_end = 1.0f;
        float map_min = 0.0f;
        float map_max = 1.0f;
        char ident[64] = "";
    } prop;
};

struct ApplicationState {
    // --- APPLICATION ---
    application::Context app {};

    struct {
        md_allocator_i* frame = 0;
        md_allocator_i* persistent = 0;
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

        struct {
            vec3_t target_position = {};
            quat_t target_orientation = {};
            float  target_distance = 0;
        } animation;
    } view;

    struct {
        bool  hide_gui = true;
        str_t path_to_file = {};
    } screenshot;

    // --- MDLIB DATA ---
    struct {
        md_allocator_i*     mol_alloc = nullptr;
        md_gl_shaders_t     gl_shaders = {};
        md_gl_shaders_t     gl_shaders_lean_and_mean = {};
        md_gl_molecule_t    gl_mol = {};
#if EXPERIMENTAL_GFX_API
        md_gfx_handle_t     gfx_structure = {};
#endif
        md_molecule_t       mol = {};
        md_trajectory_i*    traj = nullptr;

        vec3_t              mol_aabb_min = {};
        vec3_t              mol_aabb_max = {};

        uint32_t dirty_buffers = 0;
    } mold;

    DisplayProperty* display_properties = nullptr;
    str_t hovered_display_property_label = STR_LIT("");
    int   hovered_display_property_pop_idx = -1;

    // --- ASYNC TASKS HANDLES ---
    struct {
        task_system::ID backbone_computations = task_system::INVALID_ID;
        task_system::ID prefetch_frames = task_system::INVALID_ID;
        task_system::ID evaluate_full = task_system::INVALID_ID;
        task_system::ID evaluate_filt = task_system::INVALID_ID;
    } tasks;

    // --- ATOM SELECTION ---
    struct {
        SelectionLevel granularity = SelectionLevel::Atom;

        struct {
            int32_t hovered = -1;
            int32_t right_click = -1;
        } atom_idx;

        struct {
            int32_t hovered = -1;
            int32_t right_click = -1;
        } bond_idx;

        SingleSelectionSequence single_selection_sequence;

        md_bitfield_t selection_mask{};
        md_bitfield_t highlight_mask{};
        Selection* stored_selections = NULL;

        struct {
            struct {
                vec4_t visible = {0.0f, 0.0f, 1.0f, 0.25f};
                vec4_t hidden  = {0.0f, 0.0f, 0.25f, 0.4f};
            } selection;

            struct {
                vec4_t visible = {1.0f, 1.0f, 0.0f, 0.25f};
                vec4_t hidden  = {0.5f, 0.5f, 0.0f, 0.40f};
            } highlight;

            float saturation = 0.5f;
        } color;

        bool selecting = false;

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
            SelectionGrowth mode = SelectionGrowth::CovalentBond;
            float extent = 1;
            bool mask_invalid = true;
            bool show_window = false;
        } grow;
    } selection;

    // --- FRAMEBUFFER ---
    GBuffer gbuffer {};

    PickingData picking {};

    // --- ANIMATION ---
    struct {
        double frame = 0.f;  // double precision for long trajectories
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

    struct {
        bool show_window = false;
        bool enabled = false;

        struct {
            bool enabled = true;
            struct {
                uint32_t id = 0;
                float alpha_scale = 1.f;
                ImPlotColormap colormap = ImPlotColormap_Plasma;
                bool dirty = true;
            } tf;
        } dvr;

        struct {
            bool enabled = false;
            float values[8] = {};
            vec4_t colors[8] = {};
            size_t count = 0;
            //IsoSurfaces isosurfaces;
        } iso;

        struct {
            uint32_t id = 0;
            bool dirty = false;
            int  dim[3] = {0};
            float max_value = 1.f;
        } volume_texture;

        GBuffer fbo = {0};

        struct {
            vec3_t min = {0, 0, 0};
            vec3_t max = {1, 1, 1};
        } clip_volume;

        struct {
            bool enabled = true;
            bool checkerboard = true;
            int  colormap_mode = 2;
        } legend;

        float density_scale = 1.f;
        vec3_t voxel_spacing = {1.0f, 1.0f, 1.0f};
        float resolution_scale = 2.0f;

        vec4_t clip_volume_color = {1,0,0,1};
        vec4_t bounding_box_color = {0,0,0,1};

        bool show_bounding_box = true;
        bool show_reference_structures = true;
        bool show_reference_ensemble = false;
        bool show_density_volume = false;
        bool show_coordinate_system_widget = true;

        bool dirty_rep = false;
        bool dirty_vol = false;

        struct {
            RepresentationType type = RepresentationType::BallAndStick;
            ColorMapping colormap = ColorMapping::Cpk;
            float param[4] = {1,1,1,1};
            vec4_t color = {1,1,1,1};
        } rep;

        md_gl_representation_t* gl_reps = nullptr;
        mat4_t* rep_model_mats = nullptr;
        mat4_t model_mat = {0};        

        Camera camera = {};
        quat_t target_ori;
        vec3_t target_pos;
        float  target_dist;
    } density_volume;

    // --- VISUALS ---
    struct {
        struct {
            vec3_t color = {1, 1, 1};
            float intensity = 24.f;
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
            bool enabled = false;
            float focus_depth = 10.0f;
            float focus_scale = 10.0f;
        } dof;

        struct {
            bool enabled = true;
            bool jitter = true;
            float feedback_min = 0.88f;
            float feedback_max = 0.97f;

            struct {
                bool enabled = true;
                float motion_scale = 1.0f;
            } motion_blur;
        } temporal_reprojection;

        struct {
            bool enabled = true;
            postprocessing::Tonemapping tonemapper = postprocessing::Tonemapping_Filmic;
            float exposure = 1.f;
            float gamma = 2.2f;
        } tonemapping;

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
        bool atom_visibility_mask_dirty = false;
        bool show_window = false;
    } representation;

    struct {
        bool show_window = false;
        AtomElementMapping* atom_element_remappings = 0;

        md_array(DatasetItem) chains = 0;
        md_array(DatasetItem) residue_names = 0;
        md_array(DatasetItem) atom_types = 0;
    } dataset;

    struct {
        bool apply_pbc = false;
        bool unwrap_structures = false;
    } operations;

    struct {
        bool keep_representations = false;
    } settings;

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
        vec4_t triangle_color   = {1,1,0,0.5f};

        str_t text; // The current text in the texteditor
        uint64_t text_hash;

        // A bit confusing and a bit of a hack,
        // But we want to preserve the ir while evaluating it (visualizing it etc)
        // So we only commit the 'new' ir to eval_ir upon starting evaluation
        md_script_ir_t*   ir = nullptr;
        md_script_ir_t*   eval_ir = nullptr;

        md_script_eval_t* full_eval = nullptr;
        md_script_eval_t* filt_eval = nullptr;
        md_script_vis_t vis = {};

        // Semaphore to control access to IR
        md_semaphore_t ir_semaphore = {};

        bool compile_ir = false;
        bool eval_init = false;
        bool evaluate_full = false;
        bool evaluate_filt = false;
        double time_since_last_change = 0.0;
        uint64_t ir_fingerprint = 0;
    } script;

    bool show_script_window = true;
    bool show_debug_window = false;
    bool show_property_export_window = false;
};

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

static inline void modify_field(md_bitfield_t* bf, md_range_t range, SelectionOperator op) {
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

static inline void grow_mask_by_selection_granularity(md_bitfield_t* mask, SelectionLevel granularity, const md_molecule_t& mol) {
    ASSERT(mask);
    switch(granularity) {
    case SelectionLevel::Atom:
        break;
    case SelectionLevel::Residue:
        for (size_t i = 0; i < mol.residue.count; ++i) {
            md_range_t range = md_residue_atom_range(mol.residue, i);
            if (md_bitfield_popcount_range(mask, range.beg, range.end)) {
                md_bitfield_set_range(mask, range.beg, range.end);
            }
        }
        break;
    case SelectionLevel::Chain:
        for (size_t i = 0; i < mol.chain.count; ++i) {
            md_range_t range = md_chain_atom_range(mol.chain, i);
            if (md_bitfield_popcount_range(mask, range.beg, range.end)) {
                md_bitfield_set_range(mask, range.beg, range.end);
            }
        }
        break;
    default:
        ASSERT(false);
    }
}