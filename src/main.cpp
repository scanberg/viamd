#include <core/md_compiler.h>

#if MD_COMPILER_MSVC
#   ifndef _CRT_SECURE_NO_WARNINGS
#       define _CRT_SECURE_NO_WARNINGS
#   endif
#   pragma warning( disable : 26812 4244 )
#endif

#define IMGUI_DEFINE_MATH_OPERATORS

#include <md_util.h>
#include <md_gl.h>
#include <md_gfx.h>
#include <md_filter.h>
#include <md_script.h>
#include <md_system.h>
#include <md_xvg.h>
#include <md_csv.h>
#include <md_lammps.h>
#include <md_pdb.h>
#include <md_gto_int.h>

#include <core/md_log.h>
#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_intrinsics.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_ring_allocator.h>
#include <core/md_tracking_allocator.h>
#include <core/md_simd.h>
#include <core/md_os.h>
#include <core/md_unit.h>
#include <core/md_writer.h>
#include <core/md_parse.h>
#include <core/md_hash.h>

#include <gfx/gl.h>
#include <gfx/gl_utils.h>
#include <gfx/camera.h>
#include <gfx/camera_utils.h>
#include <gfx/immediate_draw_utils.h>
#include <gfx/postprocessing_utils.h>
#include <gfx/volumerender_utils.h>
#include <gfx/isomesh_utils.h>
#include <iso_meshes.h>
#include <surface_field.h>
#include <md_gto_iso.h>

#include <resource_path.h>
#include <app_settings.h>
#include <display_units.h>
#include <imgui_widgets.h>
#include <property_grid.h>
#include <implot_widgets.h>
#include <task_system.h>
#include <color_utils.h>
#include <image.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <implot.h>
#include <implot_internal.h>

#include <stdio.h>
#include <cmath>
#include <string>

#include <viamd.h>
#include <representation.h>
#include <timeline_periods_ui.h>
#include <plot_periodic.h>
#include <script_reference.h>
#include <script_suggestions.h>
#include <viamd_event.h>
#include <event.h>

#define EXPERIMENTAL_GFX_API 0
#define COMPILATION_TIME_DELAY_IN_SECONDS 1.0
#define NOTIFICATION_DISPLAY_TIME_IN_SECONDS 5.0
#define MEASURE_EVALUATION_TIME 1
#define VIAMD_RECOMPUTE_ORBITAL_PER_FRAME 0
#define VIS_FLAGS (MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY | MD_SCRIPT_VISUALIZE_TEXT)

// Global data for application
static md_allocator_i* frame_alloc = 0; // Linear allocator for scratch data which only is valid for the frame and then is reset
static md_allocator_i* persistent_alloc = 0;

static bool use_gfx = false;

// The main view's isosurface meshes (IsosurfaceRendering_Mesh): the running mean of their frames, and whether it has
// converged (the view is made until it has)
static isomesh::Target iso_mesh_target;
static bool iso_mesh_converged = true;

constexpr str_t shader_output_snippet = STR_LIT(R"(
layout(location = 0) out vec4 out_color;
layout(location = 1) out vec4 out_normal;
layout(location = 2) out vec4 out_velocity;
layout(location = 3) out vec4 out_atom_index;

vec2 encode_normal (vec3 n) {
   float p = sqrt(n.z * 8 + 8);
   return n.xy / p + 0.5;
}

vec4 encode_index(uint index) {
    return vec4(
        (index & 0x000000FFU) >> 0U,
        (index & 0x0000FF00U) >> 8U,
        (index & 0x00FF0000U) >> 16U,
        (index & 0xFF000000U) >> 24U) / 255.0;
}

vec2 compute_ss_vel(vec3 view_coord, vec3 view_vel) {
    vec3 prev_view_coord = view_coord - view_vel;
    vec4 prev_clip_coord = u_curr_view_to_prev_clip * vec4(prev_view_coord, 1);

    vec4 clip_coord = u_view_to_clip * vec4(view_coord, 1);
    vec2 curr_ndc = clip_coord.xy / clip_coord.w;
    vec2 prev_ndc = prev_clip_coord.xy / prev_clip_coord.w;
    vec2 ss_vel = (curr_ndc - prev_ndc) * 0.5 + (u_jitter_uv.xy - u_jitter_uv.zw);
    return ss_vel;
}

void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index, uint atom_flags) {
   // Not selected while something is (AtomBit_Dimmed): less of its colour
   if ((atom_flags & 16u) != 0u) {
       color.rgb = mix(vec3(dot(color.rgb, vec3(0.2126, 0.7152, 0.0722))), color.rgb, DIMMED_SATURATION);
   }
   out_color  = color;
   out_normal = vec4(encode_normal(view_normal), 0, 0);
   out_velocity = vec4(compute_ss_vel(view_coord, view_vel), 0, 0);
   out_atom_index = encode_index(atom_index);
}
)");

// The G-buffer's snippet as the shaders are made of it: with the saturation of the dimmed atoms in it
static_assert(AtomBit_Dimmed == 16u, "The G-buffer's snippet tests the flags by this bit");
static str_t output_snippet(const ApplicationState& state) {
    return str_printf(persistent_alloc, "#define DIMMED_SATURATION %.3f\n%.*s", CLAMP(state.selection.color.saturation, 0.0f, 1.0f), STR_ARG(shader_output_snippet));
}

// The marks of the selection and the highlight (the mark pass, draw_representations_marks): of the selected, 1 in r
// where the fragment is in view and 1 in g where it is there at all, in view or behind the scene; of the highlighted
// the same in b and a; and in the second target, of the highlighted, the group they are in (AtomBit_GroupMask) + 1,
// 8 more where in view, so that the group in view wins. Blended by max. In view is within a tolerance of the depth of
// the scene (the G-buffer's, md_user_texture_0): that is of the render scale and jittered, the marks are neither, so
// the depth sampled is up to about a pixel of the scene off, which on a surface at an angle to the view is a step in
// depth: about two pixels' worth of its slope.
static_assert(AtomBit_Highlighted == 1u && AtomBit_Selected == 2u && AtomBit_GroupShift == 5 && HIGHLIGHT_GROUP_COLORS == 8,
              "The marks snippet tests the flags by these bits");
constexpr str_t shader_output_snippet_marks = STR_LIT(R"(
layout(location = 0) out vec4 out_marks;
layout(location = 1) out vec4 out_group;
uniform sampler2D md_user_texture_0;

void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint atom_index, uint atom_flags) {
    vec4  clip  = u_view_to_clip * vec4(view_coord, 1.0);
    vec2  uv    = clip.xy / clip.w * 0.5 + 0.5;
    ivec2 size  = textureSize(md_user_texture_0, 0);
    ivec2 texel = clamp(ivec2(uv * vec2(size)), ivec2(0), size - 1);
    float depth = texelFetch(md_user_texture_0, texel, 0).r;
    vec4  scene = u_clip_to_view * vec4(uv * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
    float scene_z = -scene.z / scene.w;
    float z = -view_coord.z;

    // The size of a pixel of the scene where the fragment is, and the slope of its surface to the view
#if ORTHO
    float pixel  = 2.0 / (u_view_to_clip[1][1] * float(size.y));
    float facing = abs(view_normal.z);
#else
    float pixel  = 2.0 * z / (u_view_to_clip[1][1] * float(size.y));
    float facing = abs(dot(view_normal, normalize(view_coord)));
#endif
    float slope = sqrt(max(1.0 - facing * facing, 0.0)) / max(facing, 0.125);
    float tolerance = 2.0 * pixel * slope + 2e-3 * z;

    bool in_view = z <= scene_z + tolerance;
    vec2 mark = vec2(in_view ? 1.0 : 0.0, 1.0);
    bool highlighted = (atom_flags & 1u) != 0u;
    out_marks = vec4((atom_flags & 2u) != 0u ? mark : vec2(0.0), highlighted ? mark : vec2(0.0));
    uint group = 1u + ((atom_flags >> 5u) & 7u) + (in_view ? 8u : 0u);
    out_group = vec4(highlighted ? float(group) / 255.0 : 0.0, 0.0, 0.0, 0.0);
}
)");


// What Evaluate would do with the script as it is in the editor (md_script_eval_plan): how each of its properties
// stands against what is evaluated, and so how what is shown stands against what the script now says
struct ScriptPlan {
    bool valid = false;         // the editor's script compiles, and there are frames to evaluate
    bool running = false;       // an evaluation is running
    bool same = false;          // the editor's script is the one evaluated
    bool changes = false;       // Evaluate would change which properties there are, or what they are called
    const md_script_eval_plan_t* props = nullptr;   // one per property of the editor's script
    const str_t* names = nullptr;
    size_t   num_props = 0;
    uint32_t num_frames = 0;
    uint32_t num_pending = 0;   // the frames Evaluate would evaluate
};

static ScriptPlan script_plan(const ApplicationState& state, md_allocator_i* alloc) {
    ScriptPlan plan;
    const md_script_ir_t* ir = state.script.ir;
    plan.num_frames = (uint32_t)run_num_frames(&state);
    plan.valid = md_script_ir_valid(ir) && plan.num_frames > 0;
    plan.running = task_system::task_is_running(state.tasks.evaluate);
    if (!plan.valid) return plan;

    plan.same = state.script.eval && ir == state.script.eval_ir;
    plan.num_props = md_script_ir_property_count(ir);
    plan.names = md_script_ir_property_names(ir);
    md_script_eval_plan_t* props = (md_script_eval_plan_t*)md_alloc(alloc, MAX(plan.num_props, 1) * sizeof(md_script_eval_plan_t));
    md_bitset_t pending = md_bitset_empty_set();
    if (!md_script_eval_plan(props, plan.num_props, &pending, state.script.eval, ir, plan.num_frames, alloc)) {
        plan.valid = false;
        return plan;
    }
    plan.props = props;
    plan.num_pending = (uint32_t)md_bitset_count(pending);

    // A property gone, or renamed, changes what is shown without anything to evaluate
    const size_t num_evaluated = state.script.eval ? md_script_eval_property_count(state.script.eval) : 0;
    plan.changes = num_evaluated != plan.num_props;
    for (size_t i = 0; i < plan.num_props; ++i) {
        plan.changes |= props[i].status != MD_SCRIPT_EVAL_STATUS_HELD || !str_eq(props[i].from, plan.names[i]);
    }
    return plan;
}

// How property i of the plan stands, as a line of text
static void script_plan_describe(char* buf, size_t cap, const ScriptPlan& plan, size_t i) {
    const md_script_eval_plan_t& p = plan.props[i];
    const str_t name = plan.names[i];
    switch (p.status) {
    case MD_SCRIPT_EVAL_STATUS_NEW:
        snprintf(buf, cap, STR_FMT ": not evaluated yet, %u frames to evaluate", STR_ARG(name), p.num_pending);
        break;
    case MD_SCRIPT_EVAL_STATUS_CHANGED:
        snprintf(buf, cap, STR_FMT ": changed since it was evaluated, %u frames to evaluate", STR_ARG(name), p.num_pending);
        break;
    case MD_SCRIPT_EVAL_STATUS_HELD: {
        char from[96] = "";
        if (!str_eq(p.from, name)) snprintf(from, sizeof(from), ", computed as " STR_FMT " was", STR_ARG(p.from));
        if (p.num_pending == 0) {
            snprintf(buf, cap, STR_FMT ": up to date%s", STR_ARG(name), from);
        } else if (p.num_frames < plan.num_frames) {
            snprintf(buf, cap, STR_FMT ": %u of %u frames evaluated%s", STR_ARG(name), p.num_frames, plan.num_frames, from);
        } else {
            snprintf(buf, cap, STR_FMT ": %u frames to add or take out over periods%s", STR_ARG(name), p.num_pending, from);
        }
        break;
    }
    }
}

// A mark beside the line a property is named on: how it stands. It is also a handle on the property, which is dragged
// from it into the window which shows its kind.
struct ScriptGutterMark {
    int line;
    md_script_eval_status_t status;
    float done;                 // the fraction of what it needs which it holds
    char text[160];
    char ident[64];
    uint32_t kinds;             // MD_SCRIPT_PROPERTY_FLAG_*
};

static void draw_script_gutter_mark(ApplicationState* app, const TextEditor::Decorator& d, const ScriptGutterMark& m) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const ImVec2 c  = p0 + ImVec2(d.width * 0.5f, d.height * 0.5f);
    const float  r  = d.glyphSize.y * 0.26f;
    const ImU32 green = IM_COL32(90, 200, 110, 255);
    const ImU32 amber = IM_COL32(240, 170, 60, 255);
    const ImU32 gray  = IM_COL32(150, 150, 150, 255);
    switch (m.status) {
    case MD_SCRIPT_EVAL_STATUS_NEW:
        dl->AddCircle(c, r, gray, 0, 1.5f);
        break;
    case MD_SCRIPT_EVAL_STATUS_CHANGED:
        dl->AddCircleFilled(c, r, amber);
        break;
    case MD_SCRIPT_EVAL_STATUS_HELD:
        if (m.done >= 1.0f) {
            dl->AddCircleFilled(c, r, green);
        } else {
            // How much of it there is, as a slice
            dl->AddCircle(c, r, green, 0, 1.5f);
            if (m.done > 0.0f) {
                dl->PathLineTo(c);
                dl->PathArcTo(c, r, -IM_PI * 0.5f, -IM_PI * 0.5f + 2.0f * IM_PI * m.done);
                dl->PathFillConvex(green);
            }
        }
        break;
    }

    // A temporal property goes to a timeline (which a distribution takes as well), a distribution to a distribution,
    // a volume to the density volume
    const char* dnd_type = nullptr;
    const char* hint = "";
    if (m.kinds & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL) {
        dnd_type = TIMELINE_SERIES_DND;
        hint = "Drag into a Timelines or Distributions subplot";
    } else if (m.kinds & MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) {
        dnd_type = DISTRIBUTION_SERIES_DND;
        hint = "Drag into a Distributions subplot";
    } else if (m.kinds & MD_SCRIPT_PROPERTY_FLAG_VOLUME) {
        dnd_type = VOLUME_SERIES_DND;
        hint = "Drag into the Density Volume";
    }
    ImGui::InvisibleButton("##mark", ImVec2(d.width, d.height));
    if (dnd_type && ImGui::BeginDragDropSource()) {
        char path[SERIES_PATH_CAP];
        snprintf(path, sizeof(path), "script/%s", m.ident);
        const SeriesKey key = series_key(SeriesSource_Script, str_from_cstr(path));
        char label[96];
        series_label(label, sizeof(label), app, key);
        series_set_drag_payload(dnd_type, key, -1, label, series_default_color(app, key));
        ImGui::EndDragDropSource();
    } else if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s%s%s", m.text, hint[0] ? "\n" : "", hint);
    }
}

static double frame_to_time(double frame, const ApplicationState& data) {
    const int64_t num_frames = md_array_size(data.timeline.x_values);
    ASSERT(num_frames);
    const int64_t f0 = CLAMP((int64_t)frame + 0, 0, num_frames - 1);
    const int64_t f1 = CLAMP((int64_t)frame + 1, 0, num_frames - 1);
    return lerp(data.timeline.x_values[f0], data.timeline.x_values[f1], fract(frame));
}

// Try to map time t back into frame
static double time_to_frame(double time, const md_array(float) frame_times) {
    const int64_t num_frames = md_array_size(frame_times);
    if (!num_frames) return 0.0;

    const double beg = frame_times[0];
    const double end = frame_times[num_frames - 1];
    time = CLAMP(time, beg, end);

    // Estimate the frame
    const double frame_est = CLAMP(((time - beg) / (end-beg)) * (num_frames - 1), 0, num_frames - 1);

    int64_t prev_frame_idx = CLAMP((int64_t)frame_est,     0, num_frames - 1);
    int64_t next_frame_idx = CLAMP((int64_t)frame_est + 1, 0, num_frames - 1);

    if (time < (double)frame_times[prev_frame_idx]) {
        // Linear search down
        for (prev_frame_idx = prev_frame_idx - 1; prev_frame_idx >= 0; --prev_frame_idx) {
            next_frame_idx = prev_frame_idx + 1;
            if ((double)frame_times[prev_frame_idx] <= time && time <= (double)frame_times[next_frame_idx])
                break;
        }
    }
    else if (time > (double)frame_times[next_frame_idx]) {
        // Linear search up
        for (next_frame_idx = next_frame_idx + 1; next_frame_idx < num_frames; ++next_frame_idx) {
            prev_frame_idx = next_frame_idx - 1;
            if ((double)frame_times[prev_frame_idx] <= time && time <= (double)frame_times[next_frame_idx])
                break;
        }
    }

    // Compute true fraction between timestamps
    double t = (time - (double)frame_times[prev_frame_idx]) / ((double)frame_times[next_frame_idx] - (double)frame_times[prev_frame_idx]);
    t = CLAMP(t, 0.0, 1.0);

    // Compose frame value (base + fraction)
    return (double)prev_frame_idx + t;
}

static void update_view_param(ApplicationState* state);
// The transforms of the copies of the scene the periodic images show (view.periodic_transforms)
static void compute_periodic_transforms(md_array(mat4_t)* transforms, const ApplicationState& state, md_allocator_i* alloc);

// The lowest render scale offered: below half, too little of the detail is left to rebuild
static constexpr float RENDER_SCALE_MIN = 0.5f;

// The range of the shadows' softness, the light's angular radius (degrees): from about the sun's to a large softbox's
static constexpr float SHADOW_SOFTNESS_MIN = 0.5f;
static constexpr float SHADOW_SOFTNESS_MAX = 20.0f;

// The G-buffer's size for an image shown of out_w x out_h: the render scale of each side, at least a pixel
static void render_size(const ApplicationState* state, uint32_t out_w, uint32_t out_h, uint32_t* w, uint32_t* h) {
    const float s = CLAMP(state->settings.render_scale, RENDER_SCALE_MIN, 1.0f);
    *w = MAX((uint32_t)(out_w * s + 0.5f), 1U);
    *h = MAX((uint32_t)(out_h * s + 0.5f), 1U);
}

// The lengths of the jitter offered (Visuals > Anti-Aliasing): powers of two, so that each divides the isosurface
// meshes' mean of isomesh::MAX_FRAMES frames (every position weighs the same in it) and the sequence's x (Halton, base
// 2) is stratified. The TAA's history of a still view (feedback_max 0.95) weighs about the last 30 frames: of the
// positions of a cycle it averages about 8 at 8, 15 at 16, 26 at 32, 36 at 64, and 39 at most. The meshes' mean gains
// up to 128.
static const int jitter_samples_options[] = {8, 16, 32, 64, 128};
static const char* jitter_samples_names[] = {"8 samples", "16 samples", "32 samples", "64 samples", "128 samples"};
static_assert(isomesh::MAX_FRAMES % JITTER_SEQUENCE_MAX == 0, "Each length divides the meshes' mean");

// The option nearest to n (a file can hold anything), by its index
static int jitter_samples_index(int n) {
    int best = 0;
    for (int i = 1; i < (int)ARRAY_SIZE(jitter_samples_options); ++i) {
        if (abs(jitter_samples_options[i] - n) < abs(jitter_samples_options[best] - n)) best = i;
    }
    return best;
}

// The jitter's positions: temporal_aa.jitter_samples at the full resolution, and as many per pixel shown when a pixel
// of the G-buffer covers several (1 / scale^2 of them), rounded up to a power of two as the options are, up to
// JITTER_SEQUENCE_MAX
static uint32_t jitter_count(const ApplicationState* state) {
    const uint32_t base = (uint32_t)jitter_samples_options[jitter_samples_index(state->visuals.temporal_aa.jitter_samples)];
    const float sx = (float)state->gbuffer.width  / (float)MAX(state->gbuffer.out_width, 1U);
    const float sy = (float)state->gbuffer.height / (float)MAX(state->gbuffer.out_height, 1U);
    const float n  = (float)base / MAX(sx * sy, 1.0e-3f);
    uint32_t count = base;
    while ((float)count < n && count < JITTER_SEQUENCE_MAX) count *= 2;
    return MIN(count, (uint32_t)JITTER_SEQUENCE_MAX);
}

// Which anti-aliasing passes run: the switch and the mode. Temporal carries more than the smoothing (the jitter,
// motion blur, sharpening), all of it off with it.
static inline bool fxaa_active(const ApplicationState* state) {
    const auto& aa = state->visuals.anti_aliasing;
    return aa.enabled && aa.mode == AntiAliasingMode::FXAA;
}
static inline bool taa_active(const ApplicationState* state) {
    const auto& aa = state->visuals.anti_aliasing;
    return aa.enabled && aa.mode == AntiAliasingMode::TAA;
}

// The frames the shadows and the ambient occlusion take to settle when they are filtered over time: their histories,
// twice through their sample patterns (shaders/ssao/ssao_temporal.frag)
static constexpr uint32_t OCCLUSION_SETTLE_FRAMES = 32;

static inline bool occlusion_temporal_active(const ApplicationState* state) {
    const bool shadows = state->visuals.shadows.enabled && state->visuals.shading.light > 0.0f;
    return (shadows || state->visuals.ssao.enabled) && state->visuals.occlusion_temporal;
}

// A screenshot is the image as lazy rendering leaves it once it has settled (render()): with the GUI, the window as it
// is shown; without, the image alone at the screenshot's size (the window's or another), the view not moved over for
// the windows. Either way it is made as the view is, frame by frame, until it has settled.
static inline bool screenshot_pending(const ApplicationState* state) {
    return !str_empty(state->screenshot.path_to_file);
}
static inline bool screenshot_offscreen(const ApplicationState* state) {
    return screenshot_pending(state) && state->screenshot.hide_gui;
}

// Where the view's focus is seen (ViewRegion): where the windows leave it uncovered, moved over in clip space by this. Not
// in a screenshot without the windows, which shows it in the middle.
static inline vec2_t view_shift(const ApplicationState* state) {
    return screenshot_offscreen(state) ? vec2_t{0, 0} : state->view.shift;
}

static ViewRegion main_view_uncovered_region(const ApplicationState& state, const ImGuiWindow* canvas);

//static void update_density_volume_texture(ApplicationState* state);

static void render(ApplicationState* state);

// The geometry pass by the kind of what it draws: a profiler section each (gfx/profiler.h), around one md_gl_draw of
// all the ops of that kind. The cartoon, the ribbons and the tubes include the evaluation of their spline, which
// md_gl_draw does for each call that draws one.
enum GeometryStage {
    GeometryStage_SpaceFill,
    GeometryStage_Licorice,
    GeometryStage_BallAndStick,
    GeometryStage_Ribbons,
    GeometryStage_Cartoon,
    GeometryStage_GlyphSpheres,     // Of the cartoon of lipids and nucleic acids: heads
    GeometryStage_GlyphTubes,       // Tails, sterol skeletons, nucleic backbones
    GeometryStage_GlyphSlabs,       // Sugars, bases, the sticks between them
    GeometryStage_Other,
    GeometryStage_Count
};

static const char* geometry_stage_name[GeometryStage_Count] = {
    "Space fill", "Licorice", "Ball and stick", "Ribbons", "Cartoon", "Heads", "Tubes", "Slabs", "Other",
};

static GeometryStage geometry_stage(md_gl_rep_type_t type, bool glyphs) {
    if (glyphs) {
        switch (type) {
        case MD_GL_REP_SPACE_FILL: return GeometryStage_GlyphSpheres;
        case MD_GL_REP_TUBES:      return GeometryStage_GlyphTubes;
        case MD_GL_REP_SLABS:      return GeometryStage_GlyphSlabs;
        default:                   return GeometryStage_Other;
        }
    }
    switch (type) {
    case MD_GL_REP_SPACE_FILL:     return GeometryStage_SpaceFill;
    case MD_GL_REP_LICORICE:       return GeometryStage_Licorice;
    case MD_GL_REP_BALL_AND_STICK: return GeometryStage_BallAndStick;
    case MD_GL_REP_RIBBONS:        return GeometryStage_Ribbons;
    case MD_GL_REP_CARTOON:        return GeometryStage_Cartoon;
    default:                       return GeometryStage_Other;
    }
}

// Draws ops a kind at a time, in the order the kinds first appear, each a profiler section
static void draw_ops_timed(md_gl_draw_args_t args, const md_gl_draw_op_t* ops, size_t count, bool glyphs, md_allocator_i* alloc) {
    bool done[16] = {};
    md_array(md_gl_draw_op_t) group = 0;
    for (size_t i = 0; i < count; ++i) {
        const uint32_t type = (uint32_t)ops[i].type;
        if (type < ARRAY_SIZE(done) && done[type]) continue;
        md_array_shrink(group, 0);
        for (size_t k = i; k < count; ++k) {
            if (ops[k].type == ops[i].type) md_array_push(group, ops[k], alloc);
        }
        if (type < ARRAY_SIZE(done)) done[type] = true;
        args.draw_operations = { .count = (uint32_t)md_array_size(group), .ops = group };
        PUSH_GPU_SECTION(geometry_stage_name[geometry_stage(ops[i].type, glyphs)]);
        md_gl_draw(&args);
        POP_GPU_SECTION();
    }
}

static void draw_representations_opaque(ApplicationState* state);
// The marks of the selected and highlighted atoms into the marks target (gbuffer_marks_begin), unjittered
static void draw_representations_marks(ApplicationState* state);
static uint64_t opaque_scene_hash(const ApplicationState* state);
// Returns true if anything was drawn into gbuffer.tex.transparency_hdr this frame
static bool draw_representations_transparent(ApplicationState* state);

static void draw_load_dataset_window(ApplicationState* state);
static void draw_main_menu(ApplicationState* state);
static void draw_context_popup(ApplicationState* state, const PickingHit& hit);
static void draw_selection_query_window(ApplicationState* state);
static void draw_selection_grow_window(ApplicationState* state);
static void draw_animation_window(ApplicationState* state);
static void draw_timeline_window(ApplicationState* state);
static void draw_distribution_window(ApplicationState* state);
static void draw_async_task_window(ApplicationState* state);
static void draw_script_editor_window(ApplicationState* state);
static void draw_script_reference_window(ApplicationState* state);
static void open_script_reference(ApplicationState* state, str_t topic, bool take_focus = true);
static void draw_coordinate_system_widget_window(ViewTransform* target, const ViewTransform& current);

static void draw_debug_window(ApplicationState* state);
static void draw_profiler_window(ApplicationState* state);
static void draw_property_export_window(ApplicationState* state);
static void draw_structure_export_window(ApplicationState* state);
static void draw_notifications_window();

static void update_md_buffers(ApplicationState* state);
static str_t create_unique_identifier(const md_script_ir_t* ir, str_t base, md_allocator_i* alloc);

static bool export_xvg(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename);
static bool export_csv(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename);

static void create_screenshot(str_t path);

// The sizes offered in the Settings menu. The stored setting is the size itself, not an
// index into this table, so the table can change without invalidating anyone's .ini.
static const float font_sizes[] = { 10.0f, 12.0f, 14.0f, 16.0f, 18.0f, 20.0f, 24.0f, 30.0f, 36.0f, 48.0f, 64.0f, 72.0f };
static const char* font_size_names[] = { "10", "12", "14", "16", "18", "20", "24", "30", "36", "48", "64", "72" };

static int nearest_font_size_index(float size) {
    int   best_idx  = 0;
    float best_dist = FLT_MAX;
    for (int i = 0; i < (int)ARRAY_SIZE(font_sizes); ++i) {
        const float dist = fabsf(font_sizes[i] - size);
        if (dist < best_dist) {
            best_dist = dist;
            best_idx  = i;
        }
    }
    return best_idx;
}

// Also the app_settings apply hook, so a size read from the .ini takes effect on the first frame.
static void apply_font_size(void* user_data) {
    ApplicationState* state = (ApplicationState*)user_data;
    ImGuiStyle& style = ImGui::GetStyle();
    style.FontSizeBase = state->settings.font_size;
    style._NextFrameFontSizeBase = style.FontSizeBase;  // From the demo, seems like a temporary fix
}

// The look of the Visuals menu a session starts with (settings.look) and the session's own (visuals): a workspace
// sets the latter, Save as Default copies it to the former, which the .ini keeps
static void visual_look_save(ApplicationState* state) {
    state->settings.look.background  = state->visuals.background;
    state->settings.look.shading     = state->visuals.shading;
    state->settings.look.tonemapping = state->visuals.tonemapping;
    state->settings.look.dof         = state->visuals.dof;
    app_settings::mark_dirty();
}

static void visual_look_restore(ApplicationState* state) {
    state->visuals.background  = state->settings.look.background;
    state->visuals.shading     = state->settings.look.shading;
    state->visuals.tonemapping = state->settings.look.tonemapping;
    state->visuals.dof         = state->settings.look.dof;
}

// Once the .ini is read: its enums within their ranges (a file can hold anything), and the session in the look saved
static void apply_visual_settings(void* user_data) {
    ApplicationState* state = (ApplicationState*)user_data;
    const decltype(state->visuals) factory = {};
    auto& v    = state->visuals;
    auto& look = state->settings.look;
    if (!(0 <= (int)v.ssao.mode && (int)v.ssao.mode < postprocess_pipeline::SsaoMode_Count)) v.ssao.mode = factory.ssao.mode;
    if (!(0 <= (int)v.anti_aliasing.mode && (int)v.anti_aliasing.mode < (int)AntiAliasingMode::Count)) v.anti_aliasing.mode = factory.anti_aliasing.mode;
    v.temporal_aa.jitter_samples = jitter_samples_options[jitter_samples_index(v.temporal_aa.jitter_samples)];
    state->settings.render_scale = CLAMP(state->settings.render_scale, RENDER_SCALE_MIN, 1.0f);
    v.shadows.softness = CLAMP(v.shadows.softness, SHADOW_SOFTNESS_MIN, SHADOW_SOFTNESS_MAX);
    if (!(0 <= (int)look.shading.mode && (int)look.shading.mode < postprocess_pipeline::ShadingMode_Count)) look.shading.mode = factory.shading.mode;
    if (!(0 <= (int)look.shading.colors && (int)look.shading.colors < postprocess_pipeline::ColorSpace_Count)) look.shading.colors = factory.shading.colors;
    if (!(0 <= (int)look.tonemapping.tonemapper && (int)look.tonemapping.tonemapper < postprocess_pipeline::Tonemapper_Count)) look.tonemapping.tonemapper = factory.tonemapping.tonemapper;
    visual_look_restore(state);
}

// timeline.x_values and view_range are stored in display units, so a change to the time preference
// has to be carried into them. The raw frame times are still with the trajectory, so this re-derives
// rather than compounding factors, and the visible range is moved by the ratio so the user keeps
// looking at the same stretch of trajectory.
static void update_timeline_time_unit(ApplicationState* state) {
    ASSERT(state);

    if (state->timeline.units_version == display_units::version()) {
        return;
    }
    state->timeline.units_version = display_units::version();

    const double prev_scl = state->timeline.time_scale;
    const double time_scl = display_units::factor(&state->timeline.time_unit, run_time_unit(state));
    state->timeline.time_scale = time_scl;

    if (time_scl == prev_scl) {
        return;
    }

    const double* frame_times = run_frame_times(state);
    const size_t  num_frames  = md_array_size(state->timeline.x_values);
    if (frame_times) {
        for (size_t i = 0; i < num_frames; ++i) {
            state->timeline.x_values[i] = (float)(frame_times[i] * time_scl);
        }
    }

    const double ratio = prev_scl != 0.0 ? time_scl / prev_scl : 1.0;
    state->timeline.view_range.beg_x *= ratio;
    state->timeline.view_range.end_x *= ratio;
}

static void modify_selection(ApplicationState* state, md_bitset_t atom_mask, SelectionOperator op = SelectionOperator::Set) {
    ASSERT(state);
    modify_field(&state->selection.selection_mask, atom_mask, op);
}

int main(int argc, char** argv) {
#if DEBUG
    persistent_alloc = md_tracking_allocator_create(md_get_heap_allocator());
#elif RELEASE
    persistent_alloc = md_get_heap_allocator();
#else
#error "Must define DEBUG or RELEASE"
#endif
    frame_alloc = md_vm_arena_create(GIGABYTES(4));

    struct NotificationState {
        md_mutex_t lock;
        uint64_t hash;
        md_tick_t time;
    };

    NotificationState notification_state = {
        .lock = md_mutex_create(),
        .hash = 0,
        .time = 0
    };

    md_logger_i notification_logger = {
        (md_logger_o*)&notification_state,
        [](struct md_logger_o* inst, enum md_log_type_t log_type, const char* msg) {
            NotificationState& state = *(NotificationState*)inst;            

            // Prevent spamming the logger with the same message by comparing its hash
            const md_tick_t time = md_tick_now();
            const uint64_t hash = md_hash64(msg, strlen(msg), 0);

            if (md_tick_to_seconds(time - state.time) < 1.0 && hash == state.hash) {
                return;
            }
            state.hash = hash;
            state.time = time;
            
            ImGuiToastType toast_type = ImGuiToastType_None;
            switch (log_type) {
            case MD_LOG_TYPE_INFO:
                toast_type = ImGuiToastType_Info;
                break;
            case MD_LOG_TYPE_ERROR:
                toast_type = ImGuiToastType_Error;
                break;
            case MD_LOG_TYPE_DEBUG:
            default:
                break;
            }
            if (toast_type != ImGuiToastType_None) {
                // @NOTE: This needs to be protected with a mutex as it pushes to an internal vector
                md_mutex_lock(&state.lock);
                ImGui::InsertNotification(ImGuiToast(toast_type, (uint64_t)(NOTIFICATION_DISPLAY_TIME_IN_SECONDS * 1000), msg));
                md_mutex_unlock(&state.lock);
            }
        }
    };

    md_log_register(&notification_logger);

    ApplicationState state;
    ViamdEventHandler event_handler(&state);

    state.allocator.persistent = persistent_alloc;
    state.allocator.frame = frame_alloc;
    state.file_queue.ring = md_ring_allocator_create(md_alloc(persistent_alloc, MEGABYTES(1)), MEGABYTES(1));
    // One arena per dataset, held by the system and the state that share it. Both are rewound
    // rather than destroyed between loads, so this handle is set once and never reassigned.
    state.mold.sys.alloc   = md_arena_allocator_create(persistent_alloc, MEGABYTES(1));
    state.mold.state.alloc = state.mold.sys.alloc;

    md_temp_arena_system_init();

    // The masks are edited in place. Over any atom index: they hold only the words they write, and a mask over a
    // system's atoms needs no resizing when another is loaded.
    md_bitmask_init(&state.selection.selection_mask,       0, MD_BITSET_MAX_INDEX, persistent_alloc);
    md_bitmask_init(&state.selection.highlight_mask,       0, MD_BITSET_MAX_INDEX, persistent_alloc);
    md_bitmask_init(&state.selection.query.mask,           0, MD_BITSET_MAX_INDEX, persistent_alloc);
    md_bitmask_init(&state.selection.grow.mask,            0, MD_BITSET_MAX_INDEX, persistent_alloc);
    md_bitmask_init(&state.representation.visibility_mask, 0, MD_BITSET_MAX_INDEX, persistent_alloc);
    // Sets, replaced whole
    state.operations.selection_mask = md_bitset_empty_set();
    state.operations.recenter_query.mask = md_bitset_empty_set();

    // Init platform
    VIAMD_LOG_DEBUG("Initializing GL...");
    if (!application::initialize(&state.app, 0, 0, STR_LIT("VIAMD"))) {
        VIAMD_LOG_ERROR("Could not initialize application...\n");
        return -1;
    }

    // Application settings live in the ImGui .ini. Bind first, then initialize: that call
    // reads the file, so the bound values are current by the time the loop starts.
    app_settings::bind(STR_LIT("keep_representations"), &state.settings.keep_representations);
    app_settings::bind(STR_LIT("exact_isosurfaces"), &state.settings.exact_isosurfaces);
    app_settings::bind(STR_LIT("isosurface_rendering"), &state.settings.isosurface_rendering);
    app_settings::bind(STR_LIT("focus_in_free_space"), &state.settings.focus_in_free_space);
    app_settings::bind(STR_LIT("font_size"), &state.settings.font_size);
    app_settings::on_apply(apply_font_size, &state);
    // What the computer can handle, not what a figure looks like: kept per user as it is, never in a workspace
    state.app.window.vsync = true;
    app_settings::bind(STR_LIT("vsync"), &state.app.window.vsync);
    app_settings::bind(STR_LIT("render_scale"), &state.settings.render_scale);
    app_settings::bind(STR_LIT("lazy_rendering"), &state.settings.lazy_rendering);
    app_settings::bind(STR_LIT("shadows"), &state.visuals.shadows.enabled);
    app_settings::bind(STR_LIT("shadows_softness"), &state.visuals.shadows.softness);
    app_settings::bind(STR_LIT("ssao"), &state.visuals.ssao.enabled);
    app_settings::bind_enum(STR_LIT("ssao_mode"), &state.visuals.ssao.mode);
    app_settings::bind(STR_LIT("ssao_intensity"), &state.visuals.ssao.intensity);
    app_settings::bind(STR_LIT("occlusion_temporal"), &state.visuals.occlusion_temporal);
    app_settings::bind(STR_LIT("anti_aliasing"), &state.visuals.anti_aliasing.enabled);
    app_settings::bind_enum(STR_LIT("anti_aliasing_mode"), &state.visuals.anti_aliasing.mode);
    app_settings::bind(STR_LIT("taa_jitter_samples"), &state.visuals.temporal_aa.jitter_samples);
    app_settings::bind(STR_LIT("taa_history_min"), &state.visuals.temporal_aa.feedback_min);
    app_settings::bind(STR_LIT("taa_history_max"), &state.visuals.temporal_aa.feedback_max);
    app_settings::bind(STR_LIT("motion_blur"), &state.visuals.temporal_aa.motion_blur.enabled);
    app_settings::bind(STR_LIT("motion_blur_strength"), &state.visuals.temporal_aa.motion_blur.motion_scale);
    app_settings::bind(STR_LIT("sharpen"), &state.visuals.sharpen.enabled);
    app_settings::bind(STR_LIT("sharpen_strength"), &state.visuals.sharpen.weight);
    // The look a session starts with (Visuals > Save as Default)
    {
        auto& look = state.settings.look;
        app_settings::bind(STR_LIT("look_background"), look.background.color.elem, 3);
        app_settings::bind_enum(STR_LIT("look_shading"), &look.shading.mode);
        app_settings::bind_enum(STR_LIT("look_colors"), &look.shading.colors);
        app_settings::bind(STR_LIT("look_light"), &look.shading.light);
        app_settings::bind(STR_LIT("look_environment"), &look.shading.environment);
        app_settings::bind(STR_LIT("look_ambient"), &look.shading.ambient);
        app_settings::bind(STR_LIT("look_specular"), &look.shading.specular);
        app_settings::bind(STR_LIT("look_tonemapping"), &look.tonemapping.enabled);
        app_settings::bind_enum(STR_LIT("look_tonemapper"), &look.tonemapping.tonemapper);
        app_settings::bind(STR_LIT("look_exposure"), &look.tonemapping.exposure);
        app_settings::bind(STR_LIT("look_gamma"), &look.tonemapping.gamma);
        app_settings::bind(STR_LIT("look_dof"), &look.dof.enabled);
        app_settings::bind(STR_LIT("look_dof_strength"), &look.dof.aperture);
    }
    app_settings::on_apply(apply_visual_settings, &state);
    display_units::register_settings();
    app_settings::initialize();

#if MD_ENABLE_GPU
    VIAMD_LOG_DEBUG("Initializing GPU device...");
    state.gpu_device = md_gpu_device_create(nullptr);
    if (!state.gpu_device) {
        const char* reason = md_gpu_last_error();
        VIAMD_LOG_ERROR("Failed to create GPU device: %s", reason ? reason : "unknown");
    } else {
        // The scratch belongs to the device, not to whichever component happens to evaluate
        // first. Components borrow it; see the note on ApplicationState.
        state.gpu_stream = md_gpu_stream_default(state.gpu_device, MD_GPU_STREAM_COMPUTE);

        md_gpu_texture_desc_t vol_desc = {
            .type            = MD_GPU_TEX_3D,
            .format          = MD_GPU_FORMAT_R32_FLOAT,
            .usage           = MD_GPU_TEX_STORAGE,
            .width           = GPU_VOLUME_DIM,
            .height          = GPU_VOLUME_DIM,
            .depth_or_layers = GPU_VOLUME_DIM,
            .label           = "Evaluation volume",
        };
        state.gpu_volume = md_gpu_texture_create(state.gpu_stream, &vol_desc);
        if (!state.gpu_volume) {
            VIAMD_LOG_ERROR("Failed to create the GPU evaluation volume: %s", md_gpu_last_error());
        }

        // Potentials (surface fields). Kernels are made on first use, so this costs nothing until a
        // surface is coloured by one.
        md_gto_int_gpu_initialize(state.gpu_device);
        // Isosurface meshes (the mesh rendering mode): kernels made on the first extraction
        md_gto_iso_gpu_initialize(state.gpu_device);
    }
#endif

    state.app.file_drop.user_data = &state;
    state.app.file_drop.callback = [](size_t num_files, const str_t file_paths[], void* user_data) {
        ApplicationState* state = (ApplicationState*)user_data;
        ASSERT(state);

        for (size_t i = 0; i < num_files; ++i) {
            file_queue_push(&state->file_queue, file_paths[i], FileFlags_None);
        }
    };

    VIAMD_LOG_DEBUG("Initializing framebuffer...");
    gbuffer_init(&state.gbuffer, state.app.framebuffer.width, state.app.framebuffer.height);
    picking_surface_init(&state.picking_surface, interaction_surface_main);

    for (int i = 0; i < (int)ARRAY_SIZE(state.view.jitter.sequence); ++i) {
        state.view.jitter.sequence[i].x = md_halton(i + 1, 2);
        state.view.jitter.sequence[i].y = md_halton(i + 1, 3);
    }



    // Init subsystems
    VIAMD_LOG_DEBUG("Initializing immediate draw...");
    immediate::initialize();
    state.gfx.world   = immediate::queue_create("world");
	state.gfx.overlay = immediate::queue_create("overlay");

    VIAMD_LOG_DEBUG("Initializing post processing...");

    postprocess_pipeline::initialize(state.gbuffer.width, state.gbuffer.height, state.gbuffer.out_width, state.gbuffer.out_height);
    VIAMD_LOG_DEBUG("Initializing volume...");
    volume::initialize();
    isomesh::initialize();
    VIAMD_LOG_DEBUG("Initializing task system...");
    md_os_sys_info_t sys_info = {0};
	md_os_sys_info_query(&sys_info);

    int num_threads = VIAMD_NUM_WORKER_THREADS == 0 ? sys_info.num_physical_cores : VIAMD_NUM_WORKER_THREADS;
	num_threads = CLAMP(num_threads, 2, sys_info.num_physical_cores);
    task_system::initialize((size_t)num_threads);

    md_gl_initialize();
    state.gl.shaders                = md_gl_shaders_create(output_snippet(state));
    state.gl.shaders_marks          = md_gl_shaders_create(shader_output_snippet_marks);

    // Which of these are open is part of a workspace; the components register their own on initialize
    workspace_register_window("Timelines",       &state.timeline.show_window);
    workspace_register_window("Distributions",   &state.distributions.show_window);
    workspace_register_window("Representations", &state.representation.show_window);
    workspace_register_window("ScriptEditor",    &state.show_script_window);
    workspace_register_window("Animation",       &state.animation.show_window);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdInitialize, viamd::EventPayloadType_ApplicationState, &state);

#if EXPERIMENTAL_GFX_API
    md_gfx_initialize(state.gbuffer.width, state.gbuffer.height, 0);
#endif

    ImGui::init_theme();

    state.editor.SetLanguage(script_editor::language());
    state.editor.SetPalette(TextEditor::GetDarkPalette());
    state.editor.SetInsertSpacesOnTabs(true);
    // The marks of how each property stands, left of the line numbers (see draw_script_gutter_mark)
    state.editor.SetDecorationsBeforeLineNumbers(true);
    script_editor::enable_autocomplete(state.editor, &state.mold.sys);

    {
#ifdef VIAMD_DEFAULT_DATASET
        {
            md_str_writer_t sb = md_str_writer_create(frame_alloc);
            str_t path = viamd::resource_path(&sb, STR_LIT(VIAMD_DEFAULT_DATASET));
            if (md_path_is_valid(path)) {
                // @NOTE: We want explicitly to disable writing of cache files for the default dataset
                // The motivation is that the dataset may reside in a shared folder on the system that has no write access.
                file_queue_push(&state.file_queue, path, FileFlags_DisableCacheWrite);
                state.editor.SetText("s1 = resname(\"ALA\")[2:8];\nd1 = distance(10,30);\na1 = angle(2,1,3) in resname(\"ALA\");\nr = rdf(element('C'), element('H'), 10.0);\nv = sdf(s1, element('H'), 10.0);\n{lin,plan,iso} = shape_weights(all);");
            }
        }
#endif
        if (argc > 1) {
            // Assume argv[1..] are files to load
            // Currently we do not support any command line flags
            // So anything here which is a file path is assumed to be a file to load
            for (int i = 1; i < argc; ++i) {
                str_t path = str_from_cstr(argv[i]);
                if (md_path_is_valid(path)) {
                    file_queue_push(&state.file_queue, path);
                }
            }
        }
    }

#if EXPERIMENTAL_SDF == 1
    draw::scan::test_scan();
#endif
    bool time_changed = true;
    bool time_stopped = true;

    //bool demo_window = true;

    // Main loop
    while (!state.app.window.should_close) {
        // The frame for the profiler, in sections: the start (events, ImGui's new frame, files, components' frame
        // ticks), the windows, the main viewport's input, the update of the state, its upload, the rendering
        profiler::frame_begin();
        PUSH_GPU_SECTION("Start");
        application::update(&state.app);
        
        // This needs to happen first (in imgui events) to enable docking of imgui windows
#if VIAMD_IMGUI_ENABLE_DOCKSPACE
        ImGui::CreateDockspace();
#endif

        const size_t num_frames  = run_num_frames(&state);
        const size_t last_frame  = num_frames > 0 ? num_frames - 1 : 0;
        const double   max_frame = (double)last_frame;

        file_queue_process(&state);

        script_vis_begin_frame(&state);

        picking_handler_new_frame(&state.picking_handler);
        viamd::event_system_broadcast_event(viamd::EventType_ViamdPickingRangeReserve, viamd::EventPayloadType_PickingSpace, picking_handler_current_space(&state.picking_handler));

#if MD_ENABLE_GPU
        // Retires stream-ordered frees and fires completion callbacks for work
        // that finished since the last frame. Callbacks run on the calling
        // thread, so this must stay on the GL thread: they issue GL calls.
        // Not optional -- without it, freed device memory is never reclaimed.
        if (state.gpu_device) {
            md_gpu_device_poll(state.gpu_device);
        }
#endif
        // The density evaluations whose counts have come back: the rest of them issued, their readbacks queued
        gpu_volume_jobs_poll(&state);
        // What the isosurface meshes' pool tasks left behind that nobody wants any more
        iso_meshes_collect();

        viamd::event_system_broadcast_event(viamd::EventType_ViamdFrameTick, viamd::EventPayloadType_ApplicationState, &state);

        POP_GPU_SECTION();

        // GUI
        PUSH_GPU_SECTION("Windows");
        if (state.show_script_window) draw_script_editor_window(&state);
        if (state.show_script_reference_window) draw_script_reference_window(&state);
        if (state.load_dataset.show_window) draw_load_dataset_window(&state);
        if (state.representation.show_window) draw_representations_window(&state);
        if (state.distributions.show_window) draw_distribution_window(&state);
        if (state.timeline.show_window) draw_timeline_window(&state);
        series_cache_gc(&state);
        if (state.selection.query.show_window) draw_selection_query_window(&state);
        if (state.selection.grow.show_window) draw_selection_grow_window(&state);
        if (state.show_property_export_window) draw_property_export_window(&state);
        if (state.structure_export.show_window) draw_structure_export_window(&state);
        if (state.show_debug_window) draw_debug_window(&state);
        if (state.show_profiler_window) draw_profiler_window(&state);
        if (state.animation.show_window) draw_animation_window(&state);

        draw_async_task_window(&state);
        draw_main_menu(&state);
        draw_notifications_window();

        //ImGui::ShowDemoWindow();

        draw_coordinate_system_widget_window(&state.view.target, state.view.camera);
        draw_representation_legends(state);
        POP_GPU_SECTION();

        PUSH_GPU_SECTION("Viewport");
        ImGui::BeginCanvas("Main interaction window", true);
        const ImGuiWindow* canvas_window = ImGui::GetCurrentWindow();
        ImVec2 view_size = ImGui::GetContentRegionAvail();
        InteractionSurfaceState surface_state = interaction_surface(interaction_surface_main, vec_cast(view_size));
        ImGui::EndCanvas();

        const mat4_t clip_to_world = camera_view_to_world_matrix(state.view.camera) * state.view.param.matrix.inv.proj;
        const mat4_t world_to_clip = state.view.param.matrix.curr.proj * state.view.param.matrix.curr.view;

        PickingHit hit = {};

        if (surface_state.hovered) {
            InteractionSurfaceHitArgs args = {
                .picking_surface = &state.picking_surface,
                .picking_handler = state.picking_handler,
                .fbo = state.gbuffer.fbo,
                .width = state.gbuffer.width,
                .height = state.gbuffer.height,
                .clip_to_world = clip_to_world,
                .pixel_scale = (float)state.gbuffer.width / (float)MAX(state.app.framebuffer.width, 1),
            };

            interaction_surface_hit_extract(&hit, surface_state, args);

            InteractionSurfaceEvent event = {};
            interaction_surface_event_extract(&event, surface_state, hit);

            event.clip_to_world = clip_to_world;
            event.world_to_clip = world_to_clip;

            if (event.kind == InteractionSurfaceEventKind::RegionSelect) {
                md_bitset_t candidate_mask = md_bitmask_view(&state.representation.visibility_mask);
                if (event.selection_mode == InteractionSelectionMode::Remove) {
                    // When removing, only consider currently selected atoms as candidates for region selection
                    candidate_mask = md_bitmask_view(&state.selection.selection_mask);
                }
                point_set_region_mask_compute(&state.selection.highlight_mask,
                    state.mold.state.xyz, state.mold.state.num_atoms,
                    &candidate_mask, world_to_clip, event.region_min, event.region_max, event.surface_size);

                grow_mask_by_selection_granularity(&state.selection.highlight_mask, state.selection.granularity, state.mold.sys);
                if (event.region_phase == InteractionSurfaceEventPhase::Commit) {
                    // Merge highlight into selection
                    if (event.selection_mode == InteractionSelectionMode::Append) {
                        md_bitmask_or(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
                    }
                    else if (event.selection_mode == InteractionSelectionMode::Remove) {
                        md_bitmask_andnot(&state.selection.selection_mask, md_bitmask_view(&state.selection.highlight_mask));
                    }
                    md_bitmask_clear(&state.selection.highlight_mask);
                    // A region select is not a sequence of individual picks, so whatever order was
                    // recorded before no longer describes the selection.
                    single_selection_sequence_clear(&state.selection.single_selection_sequence);
                }
            } else if (event.kind == InteractionSurfaceEventKind::ContextMenu) {
                ImGui::OpenPopup("Context Popup");
            }
            
            // Since this is the main interaction view, we still broadcast all kinds of events, even if we have handled some explicitly here.
            if (event.kind != InteractionSurfaceEventKind::None) {
                viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
            }
        }

        InteractionSurfaceViewTransformArgs view_args = {
            .camera = state.view.camera,
            .trackball_param = state.view.trackball_param,
            .view_shift = state.view.shift,
        };

        InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&state.view.target, surface_state, view_args);
        if (view_result.reset_requested) {
            // What is focused on is seen in the middle of what the windows leave in view, rather than behind them
            state.view.region = state.settings.focus_in_free_space ? main_view_uncovered_region(state, canvas_window) : ViewRegion{};
            ViewTransform reset_transform = {};
            if (hit.depth < 1.0f) {
                reset_transform.distance = state.view.target.distance;
                reset_transform.orientation = state.view.camera.orientation;
                reset_transform.position = hit.world_pos + state.view.camera.orientation * vec3_set(0, 0, state.view.target.distance);
            } else {
                reset_view(&reset_transform, state.mold.state, &state.representation.visibility_mask, state.view.camera.fov_y, state.view.region,
                           &state.mold.sys, state.view.assembly);
            }
            state.view.target = reset_transform;
        }

        draw_context_popup(&state, hit);

        camera_animate(&state.view.camera, state.view.target, state.app.timing.delta_s);
        camera_animate_shift(&state.view.shift, view_region_center(state.view.region), state.app.timing.delta_s);
        state.visuals.dof.focus_depth = state.view.camera.distance;

        ImGuiWindow* win = ImGui::GetCurrentContext()->HoveredWindow;
        if (win && strcmp(win->Name, "Main interaction window") == 0) {
            script_set_hovered_property(&state,  STR_LIT(""));
        }
        POP_GPU_SECTION();

        PUSH_GPU_SECTION("Update");

        // Capture non-window specific keyboard events
        if (!ImGui::GetIO().WantCaptureKeyboard) {
            if (ImGui::IsKeyPressed(KEY_SHOW_DEBUG_WINDOW)) {
                state.show_debug_window = true;
            }

            if (ImGui::IsKeyPressed(KEY_RECOMPILE_SHADERS)) {
                VIAMD_LOG_INFO("Recompiling shaders and re-initializing volume");
                postprocess_pipeline::initialize(state.gbuffer.width, state.gbuffer.height, state.gbuffer.out_width, state.gbuffer.out_height);
                volume::initialize();
                isomesh::initialize();
                md_gl_shaders_destroy(state.gl.shaders);
                md_gl_shaders_destroy(state.gl.shaders_marks);
                state.gl.shaders = md_gl_shaders_create(output_snippet(state));
                state.gl.shaders_marks = md_gl_shaders_create(shader_output_snippet_marks);
                render_invalidate(&state);
            }

            if (ImGui::IsKeyPressed(KEY_PLAY_PAUSE)) {
                switch (state.animation.mode) {
                    case PlaybackMode::Playing:
                        state.animation.mode = PlaybackMode::Stopped;
                        state.mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                        break;
                    case PlaybackMode::Stopped:
                        state.animation.mode = PlaybackMode::Playing;
                        if (state.animation.frame == max_frame && state.animation.fps > 0) {
                            state.animation.frame = 0;
                        } else if (state.animation.frame == 0 && state.animation.fps < 0) {
                            state.animation.frame = max_frame;
                        }
                        break;
                    default:
                        ASSERT(false);
                }

            }

            if (ImGui::IsKeyPressed(KEY_SKIP_TO_PREV_FRAME) || ImGui::IsKeyPressed(KEY_SKIP_TO_NEXT_FRAME)) {
                double step = ImGui::IsKeyDown(ImGuiMod_Ctrl) ? 10.0 : 1.0;
                if (ImGui::IsKeyPressed(KEY_SKIP_TO_PREV_FRAME)) step = -step;
                state.animation.frame = CLAMP(state.animation.frame + step, 0.0, max_frame);
            }
        }

        // Not while a screenshot is pending: it is of the frame shown, made until it has settled
        if (state.animation.mode == PlaybackMode::Playing && !screenshot_pending(&state)) {
            state.animation.frame += state.app.timing.delta_s * state.animation.fps;
            state.animation.frame = CLAMP(state.animation.frame, 0.0, max_frame);
            if (state.animation.frame >= max_frame) {
                state.animation.mode = PlaybackMode::Stopped;
                state.animation.frame = max_frame;
            } else if (state.animation.frame <= 0) {
                state.animation.mode = PlaybackMode::Stopped;
                state.animation.frame = 0;
            }
        }

        {
            static auto prev_frame = state.animation.frame;
            if (state.animation.frame != prev_frame) {
                time_changed = true;
                prev_frame = state.animation.frame;
            }
            else {
                time_changed = false;
            }
        }

        if (time_changed) {
            state.mold.interpolate_system_state = true;
            time_stopped = false;

            PUSH_CPU_SECTION("Flag dynamic representations for update")
            flag_dynamic_representations_as_dirty(&state);
            POP_CPU_SECTION()
        } else {
            if (!time_stopped) {
                time_stopped = true;
                state.mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
            }
        }

        // The stride only applies while playing: starting or stopping playback, like changing the stride, changes
        // the frames the displayed state is made from without changing the time. Stopped, it is the trajectory's own.
        if (playback_stride(&state) != state.mold.last_interpolated_stride) {
            state.mold.interpolate_system_state = true;
            state.mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
        }

        if (state.mold.interpolate_system_state) {
			state.mold.interpolate_system_state = false;
            if (run_num_frames(&state) > 0) {
                PUSH_CPU_SECTION("Interpolate System State")
                interpolate_system_state(&state);
                POP_CPU_SECTION()
                viamd::event_system_enqueue_event(viamd::EventType_ViamdSystemStateChanged, viamd::EventPayloadType_ApplicationState, &state);
            }
        }

        {
            // The script is compiled again whatever changed its text. Typing is reported by the editor (its Render,
            // above in the frame, which asks for a compilation once typing pauses), but an edit through the editor's
            // API is not: Render takes its change as its own state before the keyboard's, and reports only that.
            // Such an edit is a line added from the context menu, an example inserted from the script reference, a
            // file loaded, an undo from the Edit menu; or any edit while the Script Editor window is closed, and its
            // Render not called at all. Those are whole edits, not typing under way: compiled in this frame.
            std::string text = state.editor.GetText();
            const uint64_t text_hash = md_hash64(text.c_str(), text.length(), 0);
            if (text_hash != state.script.text_hash && !state.script.compile_ir) {
                state.script.compile_ir = true;
                state.script.time_since_last_change = 1.0e9;
            }
            state.script.text_hash = text_hash;
        }

        if (state.script.compile_ir) {
            state.script.time_since_last_change += state.app.timing.delta_s;

            script_editor::markers_clear(&state.editor_markers);

            if (state.script.time_since_last_change > COMPILATION_TIME_DELAY_IN_SECONDS) {
                // The evaluations are of what is about to be replaced, so they are stopped. Their IR is not
                // affected: the compilation goes into a new one, and the old is freed once nothing holds it.
                if (state.script.eval) md_script_eval_interrupt(state.script.eval);

                // Try aquire all semaphores
                {
                    defer {
                        flag_all_representations_as_dirty(&state);
                    };

                    // Now we hold all semaphores for the script
                    state.script.compile_ir = false;
                    state.script.time_since_last_change = 0;

                    // The previous IR is freed once nothing uses it (script_ir_collect)
                    state.script.ir = script_ir_create(&state);

                    std::string src = state.editor.GetText();
                    str_t src_str {src.data(), src.length()};

                    char buf[1024];
                    size_t len = md_path_write_cwd(buf, sizeof(buf));
                    str_t old_cwd = {buf, len};
                    defer {
                        md_path_set_cwd(old_cwd);
                    };
                    
                    str_t cwd = {};
                    if (state.files.workspace[0] != '\0') {
                        extract_folder_path(&cwd, str_from_cstr(state.files.workspace));
                    } else if (const FileEntry* run_file = file_list_run_file(&state.files.list)) {
                        // The trajectory's folder, else the structure's
                        extract_folder_path(&cwd, str_from_cstr(run_file->path));
                    }
                    if (!str_empty(cwd)) {
                        md_path_set_cwd(cwd);
                    }

                    const size_t num_stored_selections = md_array_size(state.selection.stored_selections);                       
                    for (size_t i = 0; i < num_stored_selections; ++i) {
                        str_t name = str_from_cstr(state.selection.stored_selections[i].name);
                        const md_bitset_t* bf = &state.selection.stored_selections[i].atom_mask;
                        md_script_ir_add_identifier_bitset(state.script.ir, name, *bf);
                    }
                    
                    if (src_str) {
                        md_script_ir_compile_from_source(state.script.ir, src_str, &state.mold.sys, NULL);

                        // The markers copy what they need from the IR (it is freed below if it did not compile)
                        script_editor::Markers* markers = &state.editor_markers;
                        script_editor::markers_set_source(markers, src_str);

                        const size_t num_errors = md_script_ir_num_errors(state.script.ir);
                        const md_log_token_t* errors = md_script_ir_errors(state.script.ir);
                        for (size_t i = 0; i < num_errors; ++i) {
                            // Errors win over everything they overlap
                            script_editor::markers_add(markers, script_editor::MarkerType_Error, INT32_MAX, errors[i].range, errors[i].text, errors[i].context);
                        }

                        const size_t num_warnings = md_script_ir_num_warnings(state.script.ir);
                        const md_log_token_t* warnings = md_script_ir_warnings(state.script.ir);
                        for (size_t i = 0; i < num_warnings; ++i) {
                            script_editor::markers_add(markers, script_editor::MarkerType_Warning, INT32_MAX - 1, warnings[i].range, warnings[i].text, warnings[i].context);
                        }

                        const size_t num_tokens = md_script_ir_num_vis_tokens(state.script.ir);
                        const md_script_vis_token_t* vis_tokens = md_script_ir_vis_tokens(state.script.ir);
                        for (size_t i = 0; i < num_tokens; ++i) {
                            const md_script_vis_token_t& tok = vis_tokens[i];
                            script_editor::markers_add(markers, script_editor::MarkerType_Visualization, tok.depth, tok.range, tok.text, nullptr, tok.ref);
                        }

                        if (md_script_ir_valid(state.script.ir)) {
                            uint64_t ir_figerprint = md_script_ir_fingerprint(state.script.ir);
                            if (state.script.ir_fingerprint != ir_figerprint) {
                                state.script.ir_fingerprint = ir_figerprint;
                            }
                        } else {
                            state.script.ir = nullptr;
                        }
                    }
                }
            }
        }

        // IRs replaced, once no task may be on them. Before the evaluation is started again below, which could
        // otherwise keep a task running every frame
        script_ir_collect(&state);

        if (num_frames > 0) {
            if (state.script.eval_init) {
                if (task_system::task_is_running(state.tasks.evaluate)) {
                    md_script_eval_interrupt(state.script.eval);
                } else {
                    state.script.eval_init = false;

                    // The new evaluation takes over from the previous one: every property computed from the
                    // same thing keeps its values and is not evaluated again (md_script_eval_create_from).
                    // Not remade if the last compilation failed (in the same frame as the request).
                    md_script_eval_t* prev = state.script.eval;
                    state.script.eval = nullptr;

                    if (md_script_ir_valid(state.script.ir)) {
                        state.script.eval_ir = state.script.ir;
                        // The slots over periods come along, and what they hold where their property is taken over
                        const md_script_eval_desc_t desc = { .num_frames = num_frames };
                        state.script.eval = md_script_eval_create_from(prev, state.script.eval_ir, &desc, state.allocator.persistent);
                    }

                    if (prev) {
                        // The versions of the new evaluation carry on from the previous one's, so nothing derived
                        // from those can be taken for these: the caches stay
                        md_script_eval_free(prev);
                    } else {
                        // Its versions start over, in tables possibly allocated where earlier ones were
                        series_cache_free(&state);
                    }

                    state.script.evaluate = true;
                }
            }

            if (state.script.eval && md_script_ir_valid(state.script.eval_ir) &&
                md_script_eval_ir_fingerprint(state.script.eval) == md_script_ir_fingerprint(state.script.eval_ir))
            {
                // The slots over periods follow their periods at once, also while evaluating (series_resolve_histogram).
                // Nothing is evaluated for them until Evaluate is pressed: what needs evaluation only starts when
                // asked to. A distribution needs none, it is kept per frame.
                if (state.script.evaluate) {
                    if (task_system::task_is_running(state.tasks.evaluate)) {
                        // Started again once it has stopped, over the frames which remain then. Nothing it has
                        // done is lost.
                        md_script_eval_interrupt(state.script.eval);
                    } else {
                        state.script.evaluate = false;
                        md_script_eval_reset_interrupt(state.script.eval);

                        // What is left to do, for the properties and the slots
                        const md_bitset_t pending = md_script_eval_pending_frames(state.script.eval, frame_alloc);
                        const size_t count = md_bitset_count(pending);
                        if (count > 0) {
                            md_array_resize(state.script.eval_frames, count, state.allocator.persistent);
                            size_t n = 0;
                            md_bitset_iter_t it = md_bitset_iter(pending);
                            while (md_bitset_iter_next(&it) && n < count) {
                                state.script.eval_frames[n++] = (uint32_t)it.idx;
                            }

                            // The IR is taken now: eval_ir may be replaced while it runs, and the IR it runs on is
                            // not freed until it has finished. Nor are the frames touched.
                            const md_script_ir_t* ir = state.script.eval_ir;
                            md_script_eval_t* eval = state.script.eval;
                            const uint32_t* frames = state.script.eval_frames;
                            state.tasks.evaluate = task_system::create_pool_task(STR_LIT("Evaluate script"), (uint32_t)n, [&state, ir, eval, frames](uint32_t beg, uint32_t end, uint32_t thread_num) {
                                (void)thread_num;
                                md_script_eval_frames(eval, ir, &state.mold.sys, str_from_cstr(state.mold.run), frames + beg, end - beg);
                            });

#if MEASURE_EVALUATION_TIME
                            uint64_t time = (uint64_t)md_tick_now();
                            task_system::ID time_task = task_system::create_pool_task(STR_LIT("##Time Eval"), [t0 = time, n]() {
                                uint64_t t1 = md_tick_now();
                                double s = md_tick_to_seconds(t1 - t0);
                                VIAMD_LOG_INFO("Evaluation of %zu frames completed in: %.3fs", n, s);
                            });
                            task_system::set_task_dependency(time_task, state.tasks.evaluate);
#endif
                            task_system::enqueue_task(state.tasks.evaluate);
                        }
                    }
                }
            }
        }

        recenter_update(&state);

        // The bonds formed and broken in the frame shown, before they are uploaded
        bond_tracking_update(&state);

		// Perform once per-frame updates of representations (if required)
		update_all_representations(&state);

        if (state.representation.atom_visibility_mask_dirty) {
            recompute_atom_visibility_mask(&state);
            state.representation.atom_visibility_mask_dirty = false;
        }

        // Everything which sets what the script visualizes has run
        script_vis_update(&state);

        if (state.script.vis_shown && state.script.vis_shown->text) {
            PUSH_CPU_SECTION("Draw vis text");
            ImGuiWindow* window = ImGui::FindWindowByName("Main interaction window");
            if (window) {
                ImDrawList* dl = window->DrawList;
                ASSERT(dl);

                const vec2_t res = { (float)state.app.window.width, (float)state.app.window.height };
                const mat4_t mvp = state.view.param.matrix.curr.proj_no_jitter * state.view.param.matrix.curr.view;

                // Script text
                const ImU32 text_color = convert_color(state.script.text_color);
                const ImU32 rect_color = convert_color(state.script.text_bg_color);
                const float rect_rounding = 5.f;
                const ImVec2 rect_padding = ImVec2(4.f, 2.f);

                const md_script_vis_text_t* texts = state.script.vis_shown->text;
                size_t num_text = md_array_size(texts);
                for (size_t i = 0; i < num_text; ++i) {
                    const md_script_vis_text_t& vis_text = texts[i];

                    const vec4_t p = mat4_mul_vec4(mvp, vec4_from_vec3(vis_text.pos, 1.0f));
                    const vec4_t c = p / p.w;

                    // A label carrying a unit is a quantity, so it is shown in the same unit as the
                    // plots and re-formatted here rather than taken as mdlib spelled it. Anything
                    // else is drawn as it came.
                    char buf[64];
                    str_t str = vis_text.str;
                    if (!md_unit_is_none(vis_text.unit)) {
                        char unit_buf[32];
                        const double scl = display_units::factor_print(unit_buf, sizeof(unit_buf), vis_text.unit);
                        // The degree sign hugs the number the way the convention has it, everything else takes a space.
                        const char* sep = strcmp(unit_buf, "\xC2\xB0") == 0 ? "" : " ";
                        const int len = snprintf(buf, sizeof(buf), "%.2f%s%s", vis_text.value * scl, sep, unit_buf);
                        str = {buf, (size_t)CLAMP(len, 0, (int)sizeof(buf) - 1)};
                    }

                    const ImVec2 text_size = ImGui::CalcTextSize(str.beg(), str.end());

                    if (-1 < c.x && c.x < 1 && -1 < c.y && c.y < 1 && -1 < c.z && c.z < 1) {
                        ImVec2 tc = {(c.x * 0.5f + 0.5f) * res.x, (-c.y * 0.5f + 0.5f) * res.y};
                        ImVec2 p0 = tc - text_size * 0.5f;
                        ImVec2 p1 = tc + text_size * 0.5f;
                        dl->AddRectFilled(p0 - rect_padding, p1 + rect_padding, rect_color, rect_rounding);
                        dl->AddText(p0, text_color, str.beg(), str.end());
                    }
                }
            }
            POP_CPU_SECTION();
        }

        if (ImGui::IsKeyPressed(KEY_RECENTER_ON_HIGHLIGHT) && !state.editor_focused) {
			ViewFitRequest fit_request = {
				.app = state,
				.surface_id = interaction_surface_main,
                .xyzw = nullptr,
                .alloc = state.allocator.frame,
				.round = ViewFitRound_Highlight,
            };
			viamd::event_system_broadcast_event(viamd::EventType_ViamdViewFit, viamd::EventPayloadType_ViewFitRequest, &fit_request);

			if (fit_request.xyzw == nullptr) {
                // If no highlight, fallback to selection
                fit_request.round = ViewFitRound_Selection;
                viamd::event_system_broadcast_event(viamd::EventType_ViamdViewFit, viamd::EventPayloadType_ViewFitRequest, &fit_request);
			}
			if (fit_request.xyzw == nullptr) {
                // If no selection, fallback to all
                fit_request.round = ViewFitRound_Visible;
                viamd::event_system_broadcast_event(viamd::EventType_ViamdViewFit, viamd::EventPayloadType_ViewFitRequest, &fit_request);
            }

            size_t count = md_array_size(fit_request.xyzw);
            if (count > 0) {
                // Extract OBB from returned pointcloud
                vec3_t com = md_util_com_compute_vec4(fit_request.xyzw, 0, count, 0);
                mat3_t cov = mat3_covariance_matrix_vec4(fit_request.xyzw, 0, count, com);
                mat3_eigen_t eigen = mat3_eigen(cov);
                mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));

                // Compute min and maximum extent along the PCA axes
                mat3_t basis = mat3_transpose(PCA);
                vec3_t min_obb = vec3_set1(FLT_MAX);
                vec3_t max_obb = vec3_set1(-FLT_MAX);

                // Transform the gto (x,y,z,cutoff) into the PCA frame to find the min and max extend within it
                for (size_t i = 0; i < count; ++i) {
                    vec3_t xyz = vec3_from_vec4(fit_request.xyzw[i]);
                    xyz = mat3_mul_vec3(PCA, xyz);
                    min_obb = vec3_min(min_obb, xyz);
                    max_obb = vec3_max(max_obb, xyz);
                }

				vec3_t ext = max_obb - min_obb;
				float max_ext = MAX(ext.x, MAX(ext.y, ext.z));

                if (max_ext > 3.0f) {
                    vec3_t half_ext = (max_obb - min_obb) * 0.5f;
                    state.view.target = compute_optimal_view(com, half_ext, basis);
                } else {
                    state.view.target.position = com + state.view.target.orientation * vec3_set(0, 0, state.view.target.distance);
                }
                state.view.target.position = mat4_mul_vec3(state.mold.unitcell_transform, state.view.target.position, 1.0f);
            }
            else {
                reset_view(&state.view.target, state.mold.state, (const md_bitset_t*)nullptr, state.view.camera.fov_y, state.view.region, &state.mold.sys, state.view.assembly);
            }
        }

        // The motivation for doing this is to reduce the frequency at which we invalidate and upload the atom flag field to the GPU
        // For large systems, this can be a costly operation: Consider a system of 100'000'000 atoms
        // The size of each bitfield to represent the mask would be 12.5 MB
        // Given a throughput of 10GB/s results in a time of 1.25 ms per bitfield.

        uint64_t v_hash = state.representation.visibility_mask_hash;
        uint64_t h_hash = md_bitset_hash64(md_bitmask_view(&state.selection.highlight_mask), 0);
        uint64_t s_hash = md_bitset_hash64(md_bitmask_view(&state.selection.selection_mask), 0);
        uint64_t f_hash = md_hash64_combine(v_hash, md_hash64_combine(h_hash, s_hash));
        f_hash = md_hash64_combine(f_hash, state.selection.highlight_groups.key);
        state.render.selection_hash = s_hash;
        state.render.highlight_hash = h_hash;
        
        // What each representation is: not the bytes of the struct, which hold handles, pointers and flags that
        // change without the representation changing (its atom mask, for one, is replaced whenever it is evaluated)
        uint64_t r_hash = md_array_size(state.representation.reps);
        for (size_t i = 0; i < md_array_size(state.representation.reps); ++i) {
            const Representation& rep = state.representation.reps[i];
            const uint32_t kind[3] = { (uint32_t)rep.type, (uint32_t)rep.color_mapping, (uint32_t)representation_shown(state, rep) };
            r_hash = md_hash64(kind, sizeof(kind), r_hash);
            const str_t filt = str_from_cstrn(rep.filter.text, sizeof(rep.filter.text));
            r_hash = md_hash64(filt.ptr, filt.len, r_hash);
            r_hash = md_bitset_hash64(rep.atom_mask, r_hash);
        }

        // These represent the 'current' state so we can compare against it to see if they were modified
        static uint64_t highlight_hash = 0;
        static bool     highlight_empty = true;
        static uint64_t selection_hash = 0;
        static uint64_t representation_hash = 0;
        static uint64_t flag_hash = 0;

        if (h_hash != highlight_hash) {
            highlight_hash = h_hash;
            // The highlight fades in where there was none: from one atom to the next it stays
            const bool empty = md_bitmask_empty(&state.selection.highlight_mask);
            if (highlight_empty && !empty) {
                state.selection.highlight_time = ImGui::GetTime();
            }
            highlight_empty = empty;
            // enqueue the event here and do not broadcast (stall)
            viamd::event_system_enqueue_event(viamd::EventType_ViamdHighlightMaskChanged, viamd::EventPayloadType_ApplicationState, &state);
        }

        if (s_hash != selection_hash) {
            selection_hash = s_hash;
            // enqueue the event here and do not broadcast (stall)
            viamd::event_system_enqueue_event(viamd::EventType_ViamdSelectionMaskChanged, viamd::EventPayloadType_ApplicationState, &state);
        }

        if (f_hash != flag_hash) {
            flag_hash = f_hash;
            state.mold.dirty_gpu_buffers |= MolBit_DirtyFlags;
        }

        if (r_hash != representation_hash) {
            representation_hash = r_hash;
            // enqueue the event here and do not broadcast (stall)
            viamd::event_system_enqueue_event(viamd::EventType_ViamdRepresentationChanged, viamd::EventPayloadType_ApplicationState, &state);
        }

        // Process the queued events BEFORE uploading to the GPU and rendering.
        // The ViamdSystemStateChanged handler applies the toggled operations (recenter, apply
        // pbc, unwrap) directly to mold.state coordinates. update_md_buffers uploads those
        // coordinates and then clears dirty_gpu_buffers, so anything processed after it is not
        // merely rendered a frame late - it is discarded outright, because the next frame
        // re-interpolates mold.state from the trajectory before it is ever uploaded.
        viamd::event_system_process_event_queue();

        POP_GPU_SECTION();

        PUSH_GPU_SECTION("Upload");
        update_md_buffers(&state);
        update_timeline_time_unit(&state);
        POP_GPU_SECTION();

        PUSH_GPU_SECTION("Render");
        render(&state);
        POP_GPU_SECTION();

        PUSH_GPU_SECTION("Main thread tasks");
        task_system::execute_main_task_queue();
        POP_GPU_SECTION();
        profiler::frame_end();

        // Reset frame allocator
        md_vm_arena_reset(frame_alloc);

        // Swap buffers
        application::swap_buffers(&state.app);
    }

    interrupt_async_tasks(&state);
    script_vis_reset(&state);
    series_cache_free(&state);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdShutdown);

    // shutdown subsystems
    VIAMD_LOG_DEBUG("Shutting down immediate draw...");
    immediate::queue_destroy(state.gfx.world);
    immediate::queue_destroy(state.gfx.overlay);
    immediate::shutdown();

    state.gfx.world = nullptr;
    state.gfx.overlay = nullptr;

    profiler::shutdown();
    VIAMD_LOG_DEBUG("Shutting down post processing...");
    postprocess_pipeline::shutdown();
    VIAMD_LOG_DEBUG("Shutting down volume...");
    volume::shutdown();
    isomesh::target_free(&iso_mesh_target);
    for (size_t i = 0; i < md_array_size(state.representation.reps); ++i) {
        electronic_structure_mesh_free(&state, &state.representation.reps[i].electronic_structure);
        // what its field has under way holds the GPU stream and pool tasks: both still here
        surface_field_free(&state.representation.reps[i].electronic_structure.field_vol);
    }
    iso_meshes_shutdown(&state);    // what the isosurface meshes' pool tasks made, on the GPU stream
    isomesh::shutdown();
    VIAMD_LOG_DEBUG("Shutting down task system...");
    task_system::shutdown();

    gbuffer_free(&state.gbuffer);
#if MD_ENABLE_GPU
    if (state.gpu_device) {
        // A queued readback writes into a GL texture and frees its staging memory, so nothing
        // below may go until the queue has run out.
        gpu_volume_jobs_drain(&state);
        system_gpu_data_free(&state);
        md_gto_int_gpu_shutdown();
        md_gto_iso_gpu_shutdown();

        md_gpu_free(state.gpu_stream, state.gpu_coeff);
        md_gpu_texture_destroy(state.gpu_volume);
        state.gpu_coeff = 0;
        state.gpu_coeff_capacity = 0;
        state.gpu_volume = nullptr;
        state.gpu_stream = nullptr;

        md_gpu_device_destroy(state.gpu_device);
        state.gpu_device = nullptr;
    }
#endif
    application::shutdown(&state.app);

    return 0;
}

// #misc
static void update_view_param(ApplicationState* data) {
    ViewParam& param = data->view.param;
    param.matrix.prev = param.matrix.curr;
    param.jitter.prev = param.jitter.curr;

    // The copies of the scene, as of this frame and the last: as many as before, or they did not move
    {
        md_array(mat4_t)& curr = data->view.periodic_transforms;
        md_array(mat4_t)& prev = data->view.periodic_transforms_prev;
        md_array_shrink(prev, 0);
        md_array_push_array(prev, curr, md_array_size(curr), persistent_alloc);
        compute_periodic_transforms(&curr, *data, persistent_alloc);
        if (md_array_size(prev) != md_array_size(curr)) {
            md_array_shrink(prev, 0);
            md_array_push_array(prev, curr, md_array_size(curr), persistent_alloc);
        }
    }
    // The copies of the assembly shown, which stand still
    if (data->view.assembly >= 0 && (size_t)data->view.assembly >= data->mold.sys.assembly.count) {
        data->view.assembly = -1;
    }
    assembly_transforms(&data->view.assembly_transforms, &data->view.assembly_copy_offset, data->mold.sys, data->view.assembly, persistent_alloc);
    // And the subsets of its parts, which its copies draw: made when another is shown
    assembly_subsets_update(data);

    param.clip_planes.near = data->view.camera.near_plane;
    param.clip_planes.far = data->view.camera.far_plane;
    param.fov_y = data->view.camera.fov_y;
    param.resolution = {(float)data->gbuffer.width, (float)data->gbuffer.height};

    param.matrix.curr.view = camera_world_to_view_matrix(data->view.camera) * data->mold.unitcell_transform;
    param.matrix.inv.view  = mat4_inverse(data->mold.unitcell_transform) * camera_view_to_world_matrix(data->view.camera);

    const float n = data->view.camera.near_plane;
    const float f = data->view.camera.far_plane;
    const float aspect_ratio = (float)data->gbuffer.width / (float)data->gbuffer.height;

    if (taa_active(data) && data->visuals.temporal_aa.jitter) {
        static uint32_t i = 0;
        i = (i+1) % jitter_count(data);
        param.jitter.curr = data->view.jitter.sequence[i] - 0.5f;
        if (data->view.mode == CameraMode::Perspective) {
            const vec2_t j = param.jitter.curr;
            param.matrix.curr.proj = camera_view_to_clip_matrix_persp(data->view.camera, data->gbuffer.width, data->gbuffer.height, j.x, j.y);
            param.matrix.inv.proj  = camera_clip_to_view_matrix_persp(data->view.camera, data->gbuffer.width, data->gbuffer.height, j.x, j.y);
            param.matrix.curr.proj_no_jitter = camera_view_to_clip_matrix_persp(data->view.camera, aspect_ratio);
        } else {
            const float h = data->view.camera.distance * tanf(data->view.camera.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            const vec2_t scl = {w / data->gbuffer.width * 2.0f, h / data->gbuffer.height * 2.0f};
            const vec2_t j = param.jitter.curr * scl;
            param.matrix.curr.proj = camera_view_to_clip_matrix_ortho(-w + j.x, w + j.x, -h + j.y, h + j.y, n, f);
            param.matrix.inv.proj  = camera_clip_to_view_matrix_ortho(-w + j.x, w + j.x, -h + j.y, h + j.y, n, f);
            param.matrix.curr.proj_no_jitter = camera_view_to_clip_matrix_ortho(-w, w, -h, h, n, f);

        }
    } else {
        param.jitter.curr = {0,0};
        if (data->view.mode == CameraMode::Perspective) {
            param.matrix.curr.proj = camera_view_to_clip_matrix_persp(data->view.camera, aspect_ratio);
            param.matrix.inv.proj = camera_clip_to_view_matrix_persp(data->view.camera, (float)data->gbuffer.width / (float)data->gbuffer.height);
        } else {
            const float h = data->view.camera.distance * tanf(data->view.camera.fov_y * 0.5f);
            const float w = aspect_ratio * h;
            param.matrix.curr.proj = camera_view_to_clip_matrix_ortho(-w, w, -h, h, n, f);
            param.matrix.inv.proj = camera_clip_to_view_matrix_ortho(-w, w, -h, h, n, f);
        }
        param.matrix.curr.proj_no_jitter = param.matrix.curr.proj;
    }

    // What the view is focused on is seen where the windows leave it uncovered (view_shift)
    const vec2_t shift = view_shift(data);
    if (shift != vec2_t{0, 0}) {
        const mat4_t to_shifted   = mat4_translate( shift.x,  shift.y, 0.0f);
        const mat4_t from_shifted = mat4_translate(-shift.x, -shift.y, 0.0f);
        param.matrix.curr.proj           = to_shifted * param.matrix.curr.proj;
        param.matrix.curr.proj_no_jitter = to_shifted * param.matrix.curr.proj_no_jitter;
        param.matrix.inv.proj            = param.matrix.inv.proj * from_shifted;
    }

    param.matrix.curr.norm = mat4_transpose(param.matrix.inv.view);
}

// The part of the main view that the windows leave in view, as a region of it: the main view is drawn over all of
// the main viewport, behind the windows. All of it when they leave too little to show anything in.
static ViewRegion main_view_uncovered_region(const ApplicationState& state, const ImGuiWindow* canvas) {
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const ImRect bounds = ImRect(viewport->Pos, viewport->Pos + viewport->Size);
    const ImRect rect = ImGui::LargestUncoveredRect(bounds, canvas);

    ViewRegion region = {.aspect_ratio = (float)state.gbuffer.width / (float)state.gbuffer.height};

    // Smallest side, relative to that of the viewport, that is still worth fitting into
    const float MIN_SIDE = 0.25f;
    if (ImMin(rect.GetWidth(), rect.GetHeight()) >= MIN_SIDE * ImMin(bounds.GetWidth(), bounds.GetHeight())) {
        // Normalized device coordinates have y up
        const ImVec2 p0 = (rect.Min - bounds.Min) / bounds.GetSize();
        const ImVec2 p1 = (rect.Max - bounds.Min) / bounds.GetSize();
        region.min = {p0.x * 2.0f - 1.0f, 1.0f - p1.y * 2.0f};
        region.max = {p1.x * 2.0f - 1.0f, 1.0f - p0.y * 2.0f};
    }
    return region;
}

// ### DRAW WINDOWS ###
static void draw_main_menu(ApplicationState* data) {
    ASSERT(data);
    bool new_clicked = false;
    char path_buf[2048] = "";

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Open File...", "CTRL+L")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open)) {
                    file_queue_push(&data->file_queue, str_from_cstr(path_buf), FileFlags_ShowDialogue);
                }
            }
            if (ImGui::MenuItem("Open Workspace...", "CTRL+O")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open, WORKSPACE_FILE_EXTENSION)) {
                    load_workspace(data, str_from_cstr(path_buf));
                }
            }
            if (ImGui::MenuItem("Save Workspace", "CTRL+S")) {
                if (strnlen(data->files.workspace, sizeof(data->files.workspace)) == 0) {
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, WORKSPACE_FILE_EXTENSION)) {
                        save_workspace(data, {path_buf, strnlen(path_buf, sizeof(path_buf))});
                    }
                } else {
                    save_workspace(data, str_from_cstr(data->files.workspace));
                }
            }
            if (ImGui::MenuItem("Save Workspace As...")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, WORKSPACE_FILE_EXTENSION)) {
                    save_workspace(data, {path_buf, strnlen(path_buf, sizeof(path_buf))});
                }
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Quit", "ALT+F4")) {
                data->app.window.should_close = true;
            }
            ImGui::EndMenu();
        }
        /*
        if (ImGui::BeginMenu("Edit")) {
        if (ImGui::MenuItem("Undo", "CTRL+Z")) {
        }
        if (ImGui::MenuItem("Redo", "CTRL+Y", false, false)) {
        }  // Disabled item
        ImGui::Separator();
        if (ImGui::MenuItem("Cut", "CTRL+X")) {
        }
        if (ImGui::MenuItem("Copy", "CTRL+C")) {
        }
        if (ImGui::MenuItem("Paste", "CTRL+V")) {
        }
        ImGui::EndMenu();
        }
        */
        if (ImGui::BeginMenu("Visuals")) {
            const bool basic_shading = data->visuals.shading.mode == postprocess_pipeline::ShadingMode_Basic;
            const float gap = ImGui::GetStyle().ItemInnerSpacing.x;

            // A menu sizes itself to what is in it, and the rows fill what there is: the width is set here
            ImGui::Dummy(ImVec2(ImGui::GetFontSize() * 26.0f, 0.0f));

            // The display first, then sections in the order the image is made: what is seen, how it is lit, what the camera
            // does with the light, how the result is smoothed; then how isosurfaces are drawn. The effects that
            // can be turned off have the switch on their header, and their settings only while they are on.
            prop::begin("visuals");

            // The computer's settings (VSync, the ambient occlusion and the anti-aliasing here; Isosurfaces marks its own)
            // are kept as they are set: the .ini
            bool machine_dirty = false;

            if (prop::section("Display")) {
                prop::row("vsync", "Limit the frame rate to the display's refresh rate");
                machine_dirty |= ImGui::Checkbox("##vsync", &data->app.window.vsync);

                // A fraction of each side, in steps of 5 %, and the size it renders at beside it
                prop::row("resolution", "Render the scene at this fraction of the display's resolution and upscale it.\n"
                                        "TAA rebuilds the detail over frames: a still view comes close to the full resolution, a moving one is softer.\n"
                                        "FXAA, or no anti-aliasing, upscales each frame on its own: softer.\n"
                                        "Screenshots are rendered at the full resolution.");
                uint32_t rw = 0, rh = 0;
                render_size(data, data->app.framebuffer.width, data->app.framebuffer.height, &rw, &rh);
                char size_text[32];
                snprintf(size_text, sizeof(size_text), "%u x %u", rw, rh);
                ImGui::SetNextItemWidth(prop::control_width_before(size_text));
                float percent = data->settings.render_scale * 100.0f;
                if (ImGui::SliderFloat("##resolution", &percent, RENDER_SCALE_MIN * 100.0f, 100.0f, "%.0f %%")) {
                    data->settings.render_scale = CLAMP(roundf(percent / 5.0f) * 0.05f, RENDER_SCALE_MIN, 1.0f);
                    machine_dirty = true;
                }
                ImGui::SameLine(0.0f, 2.0f * gap);
                ImGui::TextDisabled("%s", size_text);

                prop::row("lazy rendering", "Render the scene only when something in it changes, and show the last image otherwise: "
                                            "it saves power, on a laptop the battery.\n"
                                            "Hovering and selecting only draw their marks again.");
                machine_dirty |= ImGui::Checkbox("##lazy_rendering", &data->settings.lazy_rendering);
                prop::section_end();
            }

            if (prop::section("View")) {
                prop::row("camera");
                const char* reset_label = "reset";
                const float reset_w = ImGui::CalcTextSize(reset_label).x + 2.0f * ImGui::GetStyle().FramePadding.x;
                ImGui::SetNextItemWidth(MAX(prop::control_width() - reset_w - gap, 1.0f));
                ImGui::Combo("##camera", (int*)(&data->view.mode), "Perspective\0Orthographic\0");
                ImGui::SameLine(0.0f, gap);
                if (ImGui::Button(reset_label)) {
                    reset_view(&data->view.target, data->mold.state, &data->representation.visibility_mask, data->view.camera.fov_y, data->view.region,
                               &data->mold.sys, data->view.assembly);
                }
                ImGui::SetItemTooltip("Frame what is shown again, from the default direction");
                if (data->view.mode == CameraMode::Perspective) {
                    float fov = RAD_TO_DEG(data->view.camera.fov_y);
                    prop::row("field of view");
                    if (ImGui::SliderFloat("##fov", &fov, 12.5f, 80.0f, (const char*)u8"%.1f°")) {
                        data->view.camera.fov_y = DEG_TO_RAD(fov);
                    }
                }
                const md_assembly_data_t& assemblies = data->mold.sys.assembly;
                if (assemblies.count > 0) {
                    prop::row("assembly", "The biological assembly the file describes (a dimer, the shell of a virus) in place of the\n"
                                          "coordinates as they are: the copies of its chains its operators make, drawn by instancing.\n"
                                          "The atoms of a copy are those of the file, and pick and select as them. The operators place\n"
                                          "the coordinates as the file gives them: centering, wrapping and making whole (Operations)\n"
                                          "are not applied while one is shown, and coordinates they moved before stay moved.");
                    // A chain of an assembly: a copy of a polymer strand, a part of the instances it copies
                    auto num_chains = [&](size_t asm_idx) {
                        const md_system_t& sys = data->mold.sys;
                        size_t n = 0;
                        for (size_t g = assemblies.group_offset[asm_idx]; g < assemblies.group_offset[asm_idx + 1]; ++g) {
                            const AssemblyGroup group = assembly_group(assemblies, g);
                            size_t polymers = 0;
                            for (const md_instance_idx_t* k = group.inst; k < group.inst_end; ++k) {
                                const md_urange_t res = md_system_instance_residue_range(&sys, (size_t)*k);
                                for (uint32_t r = res.beg; r < res.end;) {
                                    const md_part_idx_t part = md_system_part_find_by_residue_idx(&sys, r);
                                    if (part < 0) break;
                                    polymers += md_part_kind_is_polymer(md_system_part_kind(&sys, (size_t)part)) ? 1 : 0;
                                    r = md_system_part_residue_range(&sys, (size_t)part).end;
                                }
                            }
                            n += polymers * (size_t)(group.oper_end - group.oper);
                        }
                        return n;
                    };
                    char preview[64] = "As deposited";
                    if (data->view.assembly >= 0 && (size_t)data->view.assembly < assemblies.count) {
                        const str_t id = md_assembly_id(&assemblies, (size_t)data->view.assembly);
                        snprintf(preview, sizeof(preview), "Assembly " STR_FMT " (%zu chains)", STR_ARG(id), num_chains((size_t)data->view.assembly));
                    }
                    int chosen = data->view.assembly;
                    ImGui::SetNextItemWidth(prop::control_width());
                    if (ImGui::BeginCombo("##assembly", preview)) {
                        if (ImGui::Selectable("As deposited", data->view.assembly < 0)) chosen = -1;
                        for (size_t k = 0; k < assemblies.count; ++k) {
                            char label[64];
                            const str_t id = md_assembly_id(&assemblies, k);
                            snprintf(label, sizeof(label), "Assembly " STR_FMT " (%zu chains)##%zu", STR_ARG(id), num_chains(k), k);
                            if (ImGui::Selectable(label, data->view.assembly == (int)k)) chosen = (int)k;
                        }
                        ImGui::EndCombo();
                    }
                    if (chosen != data->view.assembly) {
                        // Its copies at once (they are set every frame, update_view_param), and what is shown changes
                        // its extent: the view frames it
                        data->view.assembly = chosen;
                        assembly_transforms(&data->view.assembly_transforms, &data->view.assembly_copy_offset, data->mold.sys, chosen, persistent_alloc);
                        reset_view(&data->view.target, data->mold.state, &data->representation.visibility_mask, data->view.camera.fov_y, data->view.region, &data->mold.sys, chosen);
                        // Room to zoom out to all of it
                        data->view.trackball_param.max_distance = MAX(data->view.trackball_param.max_distance, data->view.target.distance * 10.0f);
                    }
                }
                if (md_unitcell_flags(&data->mold.state.unitcell) != 0) {
                    prop::row("periodic images", "The scene again in the neighboring cells, this many each way along a, b and c:\n"
                                                 "1 1 1 shows the 26 cells around it. A representation can leave itself out\n"
                                                 "(its periodic images, under Representations). Not with an assembly shown.");
                    ImGui::BeginDisabled(data->view.assembly >= 0);
                    ImGui::SliderInt3("##periodic_images", data->view.periodic_images, 0, PERIODIC_IMAGES_MAX, "%d", ImGuiSliderFlags_AlwaysClamp);
                    ImGui::EndDisabled();
                }
                prop::section_end();
            }

            // How the scene is lit, how its colors are read, and what it is shown over
            if (prop::section("Shading")) {
                const char* modes[postprocess_pipeline::ShadingMode_Count] = {"Physically Based", "Basic"};
                int mode = (int)data->visuals.shading.mode;
                prop::row("mode", "Physically Based: lit as a plastic-like material by a light and its surroundings, then tone mapped.\n"
                                  "Basic: the colors as they are picked where a surface faces the light, darker away from it. No tone mapping.");
                if (ImGui::Combo("##mode", &mode, modes, IM_ARRAYSIZE(modes))) {
                    data->visuals.shading.mode = (postprocess_pipeline::ShadingMode)mode;
                }
                prop::row("light", "Strength of the directional light, relative to the default (1).\n"
                                   "More makes the shapes stand out more, its shadows deeper, and the scene brighter. 0 leaves the surroundings only.");
                ImGui::SliderFloat("##light", &data->visuals.shading.light, 0.0f, 4.0f, (const char*)u8"%.2f×");
                if (!basic_shading) {
                    prop::row("environment", "Strength of the light from all around (a sky a little brighter than the ground), relative to the default (1).\n"
                                             "More softens the shading and brightens the side away from the light; the ambient occlusion darkens it in creases.\n"
                                             "At 1 for both, a color facing the viewer shows close to as picked (with sRGB colors and the Neutral curve);\n"
                                             "Exposure (Tone Mapping) scales all of the light.");
                    ImGui::SliderFloat("##environment", &data->visuals.shading.environment, 0.0f, 4.0f, (const char*)u8"%.2f×");
                } else {
                    prop::row("ambient", "How much of the color is left on the side facing away from the light");
                    ImGui::SliderFloat("##ambient", &data->visuals.shading.ambient, 0.0f, 1.0f, "%.2f");
                    prop::row("specular", "Strength of the white highlight, 0 for none");
                    ImGui::SliderFloat("##specular", &data->visuals.shading.specular, 0.0f, 1.0f, "%.2f");
                }
                const char* color_spaces[postprocess_pipeline::ColorSpace_Count] = {"sRGB", "Linear"};
                int colors = (int)data->visuals.shading.colors;
                prop::row("colors", "How the colors of atoms, cartoons and surfaces are read.\n"
                                    "sRGB: as picked and seen; a color in full light shows as its swatch.\n"
                                    "Linear: as the fraction of the light reflected per channel; lighter and paler than the swatch.\n"
                                    "The background always shows as picked.");
                if (ImGui::Combo("##colors", &colors, color_spaces, IM_ARRAYSIZE(color_spaces))) {
                    data->visuals.shading.colors = (postprocess_pipeline::ColorSpace)colors;
                }
                prop::row("background", "Shown as picked, whatever the shading and tone mapping. It does not light the scene.");
                prop::color("##background", data->visuals.background.color.elem);
                prop::section_end();
            }

            // The visibilities of the two lights, found in what is on screen: the shadows (of the directional light) and
            // the ambient occlusion (of the surroundings). Each has its switch and its own settings beside it, the header's
            // switch is for both; they share the temporal filter and the time.
            bool* occlusion_switches[2] = {&data->visuals.shadows.enabled, &data->visuals.ssao.enabled};
            bool occlusion_toggled = false;
            if (prop::section_toggle_all("Shadows & Occlusion", occlusion_switches, 2, &occlusion_toggled, false,
                                         "Darken what less light reaches: the shadows of the directional light and the occlusion of the surroundings.\n"
                                         "Turns both off, or on; each has its own switch in the section.")) {
                const bool shadows_on = data->visuals.shadows.enabled;
                const bool ao_on = data->visuals.ssao.enabled;

                prop::row("shadows", "Shadows of the directional light: darken cavities and what lies behind other parts, where the light is blocked.\n"
                                     "How dark they are is how strong the light is against the surroundings (Shading: Light, Environment).\n"
                                     "Found in what is on screen: what is off screen, or hidden behind a surface, casts none.");
                machine_dirty |= ImGui::Checkbox("##shadows", &data->visuals.shadows.enabled);
                if (shadows_on) {
                    ImGui::SameLine(0.0f, gap);
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    machine_dirty |= ImGui::SliderFloat("##softness", &data->visuals.shadows.softness, SHADOW_SOFTNESS_MIN, SHADOW_SOFTNESS_MAX,
                                                        (const char*)u8"softness %.1f°", ImGuiSliderFlags_AlwaysClamp);
                    ImGui::SetItemTooltip("The size of the light, as its angular radius. A larger light casts softer shadows,\n"
                                          "whose edges blur more the further they fall from what casts them.");
                }

                prop::row("ambient occlusion", "Darken creases and cavities, where less of the surroundings' light reaches (SSAO).\n"
                                               "Performance: point samples.\n"
                                               "Quality: screen-space slices with a visibility bitmask, finds creases and cavities more reliably and avoids halos behind small atoms.");
                machine_dirty |= ImGui::Checkbox("##ssao", &data->visuals.ssao.enabled);
                if (ao_on) {
                    const char* modes[postprocess_pipeline::SsaoMode_Count] = {"Performance", "Quality"};
                    int mode = (int)data->visuals.ssao.mode;
                    ImGui::SameLine(0.0f, gap);
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::Combo("##ssao_mode", &mode, modes, IM_ARRAYSIZE(modes))) {
                        data->visuals.ssao.mode = (postprocess_pipeline::SsaoMode)mode;
                        machine_dirty = true;
                    }
                    prop::row("occlusion strength", "How much the ambient occlusion darkens. The shadows have no such setting: they are as dark as the light is strong.");
                    machine_dirty |= ImGui::SliderFloat("##ssao_intensity", &data->visuals.ssao.intensity, 0.0f, 8.f, "%.2f");
                }

                if (shadows_on || ao_on) {
                    prop::row("temporal filter", "Accumulate the shadows and the occlusion over frames: removes the flicker in crevasses and smooths the soft shadows.\n"
                                                 "Off gives a fixed sample pattern instead.");
                    machine_dirty |= ImGui::Checkbox("##occlusion_temporal", &data->visuals.occlusion_temporal);

                    // Their GPU time in the main view (the component windows' have their own), with filtering
                    const float shadow_ms = shadows_on ? profiler::gpu_ms("Shadows", "Render") : 0.0f;
                    const float ao_ms = ao_on ? profiler::gpu_ms("SSAO", "Render") : 0.0f;
                    if (shadow_ms + ao_ms > 0.0f) {
                        prop::row("GPU time");
                        ImGui::TextDisabled("%.2f ms", shadow_ms + ao_ms);
                        ImGui::SetItemTooltip("GPU time per frame, with filtering\n"
                                              "Shadows: %.3f ms (rays %.3f ms)\n"
                                              "Ambient occlusion: %.3f ms (estimator %.3f ms)",
                                              shadow_ms, shadows_on ? profiler::gpu_ms("Shadow Rays", "Render") : 0.0f,
                                              ao_ms, ao_on ? profiler::gpu_ms("AO", "Render") : 0.0f);
                    }
                }
                prop::section_end();
            }
            machine_dirty |= occlusion_toggled;
#if EXPERIMENTAL_CONE_TRACED_AO == 1
            if (prop::section_toggle("Cone Traced AO", &data->visuals.cone_traced_ao.enabled)) {
                prop::row("intensity");
                ImGui::SliderFloat("##intensity", &data->visuals.cone_traced_ao.intensity, 0.01f, 5.f);
                prop::row("step scale");
                ImGui::SliderFloat("##step_scale", &data->visuals.cone_traced_ao.step_scale, 0.25f, 8.f);
                prop::section_end();
            }
#endif

            // What the camera and the display do with the light. Basic shading has none of it: shown, but not to be set.
            ImGui::BeginDisabled(basic_shading);
            if (prop::section_toggle("Tone Mapping", &data->visuals.tonemapping.enabled, nullptr, false,
                                     basic_shading ? "Not used by Basic shading, which shows the colors as they are picked"
                                                   : "Map the lighting into the colors a display can show.\nOff: no curve, what is brighter than white is clipped.")) {
                // ACES first: it is the default, the look the shading had before. Exposure Gamma is the same as None
                // (exposure and gamma apply to every curve), and shows as it.
                static const postprocess_pipeline::Tonemapper tonemappers[] = {
                    postprocess_pipeline::Tonemapper_ACES,
                    postprocess_pipeline::Tonemapper_Neutral,
                    postprocess_pipeline::Tonemapper_Filmic,
                    postprocess_pipeline::Tonemapper_Passthrough,
                };
                static const char* tonemapper_names[] = {"ACES", "Neutral", "Filmic", "None"};
                int idx = (int)ARRAY_SIZE(tonemappers) - 1;
                for (int i = 0; i < (int)ARRAY_SIZE(tonemappers); ++i) {
                    if (tonemappers[i] == data->visuals.tonemapping.tonemapper) idx = i;
                }
                prop::row("curve", "ACES: a filmic curve with more contrast, which shifts and desaturates the colors.\n"
                                   "Neutral: the colors as they are, only the highlights are compressed (Khronos PBR Neutral).\n"
                                   "Filmic: John Hable's Uncharted 2 curve.\n"
                                   "None: no curve, what is brighter than white is clipped.");
                if (ImGui::Combo("##curve", &idx, tonemapper_names, IM_ARRAYSIZE(tonemapper_names))) {
                    data->visuals.tonemapping.tonemapper = tonemappers[idx];
                }
                prop::row("exposure", "Scales the light reaching the camera. At 1 a color on a surface facing the viewer shows close to as picked (Neutral)");
                ImGui::SliderFloat("##exposure", &data->visuals.tonemapping.exposure, 0.01f, 10.f, "%.2f");
                prop::row("gamma", "An adjustment of the sRGB curve: 2.2 is the sRGB curve itself");
                ImGui::SliderFloat("##gamma", &data->visuals.tonemapping.gamma, 1.0f, 3.0f, "%.2f");
                prop::section_end();
            }
            ImGui::EndDisabled();

            if (prop::section_toggle("Depth of Field", &data->visuals.dof.enabled, nullptr, false, "Blur what lies away from the focus")) {
                prop::row("blur", "Blur of distant objects, in percent of the view height. Independent of zoom level");
                ImGui::SliderFloat("##blur", &data->visuals.dof.aperture, 0.0f, 4.0f, "%.2f %%");
                prop::section_end();
            }

            // The switch on the header, as the other effects', and the mode within; folded, as the switch is what is
            // mostly wanted, the mode and the temporal resolve's settings seldom
            bool aa_toggled = false;
            if (prop::section_toggle("Anti-Aliasing", &data->visuals.anti_aliasing.enabled, &aa_toggled, false,
                                     "Smooth jagged edges.\nOff: screenshots are a single sample.")) {
                const char* modes[(int)AntiAliasingMode::Count] = {"TAA", "FXAA"};
                int mode = (int)data->visuals.anti_aliasing.mode;
                prop::row("mode", "TAA: each frame blended with the ones before, rendered from slightly shifted views. The smoothest still image.\n"
                                  "FXAA: one pass over each frame on its own. Nothing is carried over from frame to frame, but edges crawl in motion and detail finer than a pixel flickers.");
                if (ImGui::Combo("##mode", &mode, modes, IM_ARRAYSIZE(modes))) {
                    data->visuals.anti_aliasing.mode = (AntiAliasingMode)mode;
                    machine_dirty = true;
                }
                if (taa_active(data)) {
                    // The temporal resolve's settings
                    prop::row("jitter", "The positions the view is shifted through, a frame each, that the edges of a still view are smoothed over.\n"
                                        "The image keeps a history of about its last 30 frames, so it gains little past 32; the isosurface meshes and\n"
                                        "transparent representations average 128 frames and gain up to 128. A still view takes longer to settle past 32.");
                    int jitter = jitter_samples_index(data->visuals.temporal_aa.jitter_samples);
                    if (ImGui::Combo("##jitter_samples", &jitter, jitter_samples_names, IM_ARRAYSIZE(jitter_samples_names))) {
                        data->visuals.temporal_aa.jitter_samples = jitter_samples_options[jitter];
                        machine_dirty = true;
                    }
                    prop::row("history min", "How much of the previous frames is kept where the image changes: lower is sharper in motion, higher is smoother");
                    machine_dirty |= ImGui::SliderFloat("##history_min", &data->visuals.temporal_aa.feedback_min, 0.5f, 1.0f, "%.2f");
                    prop::row("history max", "How much of the previous frames is kept where the image is still");
                    machine_dirty |= ImGui::SliderFloat("##history_max", &data->visuals.temporal_aa.feedback_max, 0.5f, 1.0f, "%.2f");
                    // On, and how much, on one row each: the amount means nothing while off
                    prop::row("motion blur");
                    machine_dirty |= ImGui::Checkbox("##motion_blur", &data->visuals.temporal_aa.motion_blur.enabled);
                    if (data->visuals.temporal_aa.motion_blur.enabled) {
                        ImGui::SameLine(0.0f, gap);
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        machine_dirty |= ImGui::SliderFloat("##motion_blur_strength", &data->visuals.temporal_aa.motion_blur.motion_scale, 0.f, 2.0f, "strength %.2f");
                    }
                    prop::row("sharpen");
                    machine_dirty |= ImGui::Checkbox("##sharpen", &data->visuals.sharpen.enabled);
                    if (data->visuals.sharpen.enabled) {
                        ImGui::SameLine(0.0f, gap);
                        ImGui::SetNextItemWidth(-FLT_MIN);
                        machine_dirty |= ImGui::SliderFloat("##sharpen_strength", &data->visuals.sharpen.weight, 0.0f, 4.0f, "strength %.2f");
                    }
                }
                prop::section_end();
            }
            machine_dirty |= aa_toggled;

            if (prop::section("Isosurfaces", false)) {
                int rendering = CLAMP(data->settings.isosurface_rendering, 0, 1);
                prop::row("rendering", "Ray Casting: the surfaces are found in the evaluated volume, ray by ray, every frame.\n"
                                       "Mesh: the surfaces of orbitals and densities are extracted as triangle meshes from the\n"
                                       "field itself rather than from the volume, at a spacing set by the representation's resolution,\n"
                                       "and drawn with stochastic transparency that is refined over the frames while the view\n"
                                       "is still. Extraction runs when the surface changes (an isovalue, an orbital, a frame);\n"
                                       "drawing is cheap. Volumes from scripts are always ray cast.");
                if (ImGui::Combo("##rendering", &rendering, "Ray Casting\0Mesh\0")) {
                    data->settings.isosurface_rendering = rendering;
                    app_settings::mark_dirty();
                    // Their volumes are evaluated only when something draws from them
                    for (size_t i = 0; i < md_array_size(data->representation.reps); ++i) {
                        if (data->representation.reps[i].type == RepresentationType::ElectronicStructure) {
                            flag_representation_as_dirty(&data->representation.reps[i]);
                        }
                    }
                }
                // Also for the volumes that are ray cast in either mode (from scripts)
                int iso_mode = data->settings.exact_isosurfaces ? 1 : 0;
                prop::row("intersection", "How ray cast isosurfaces are found (with Mesh rendering: only those of volumes from scripts).\n"
                                          "Fast: samples the volume once per voxel along each ray. Parts of a surface thinner than\n"
                                          "about a voxel can be missed, and up close or at grazing angles the sampling can show;\n"
                                          "a finer volume resolution makes both smaller.\n"
                                          "Exact: intersects every cell a ray passes through, so nothing is missed and the surfaces\n"
                                          "do not depend on any sampling. Roughly twice the GPU time of Fast.");
                if (ImGui::Combo("##intersection", &iso_mode, "Fast\0Exact\0")) {
                    data->settings.exact_isosurfaces = iso_mode == 1;
                    app_settings::mark_dirty();
                }
                prop::section_end();
            }

            prop::end();
            if (machine_dirty) {
                app_settings::mark_dirty();
            }

            // The look (View to Depth of Field but the camera, and the ambient occlusion) is the workspace's: a new
            // session starts with the one saved here
            ImGui::Separator();
            ImGui::PushID("Defaults");
            const float half_w = 0.5f * (ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ItemSpacing.x);
            if (ImGui::Button("save as default", ImVec2(half_w, 0))) {
                visual_look_save(data);
            }
            ImGui::SetItemTooltip("Start new sessions with this look: the shading, the background, the tone mapping and\n"
                                  "the depth of field. A workspace keeps its own and loads it. The rest of this menu is\n"
                                  "this computer's, kept as it is set.");
            ImGui::SameLine();
            if (ImGui::Button("reset to default", ImVec2(-FLT_MIN, 0))) {
                visual_look_restore(data);
            }
            ImGui::SetItemTooltip("Back to the look saved as default");
            ImGui::PopID();
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Windows")) {
            ImGui::Checkbox("Animation", &data->animation.show_window);
            ImGui::Checkbox("Representations", &data->representation.show_window);
            ImGui::Checkbox("Script Editor", &data->show_script_window);
            ImGui::Checkbox("Script Reference", &data->show_script_reference_window);
            ImGui::Checkbox("Timelines", &data->timeline.show_window);
            ImGui::Checkbox("Distributions", &data->distributions.show_window);
            ImGui::Checkbox("Structure Export", &data->structure_export.show_window);
            ImGui::Checkbox("Profiler", &data->show_profiler_window);

            viamd::event_system_broadcast_event(viamd::EventType_ViamdWindowDrawMenu);

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Selection")) {
            ImGui::Combo("Select by", (int*)(&data->selection.granularity), selection_granularity_str, (int)SelectionGranularity::Count);
            ImGui::SetItemTooltip("What a click or a dragged region in the viewport selects: single atoms, whole residues,\nwhole parts (a covalent piece of a molecule: a strand of a polymer with its caps, a glycan, a ligand),\nwhole chains (all that the file gives the chain id, or the molecule where it gives none),\nor whole molecules (atoms joined by bonds, with their residues and polymer chains)");
            size_t num_selected_atoms = md_bitmask_count(&data->selection.selection_mask);
            if (ImGui::MenuItem("Invert")) {
                md_bitmask_invert_range(&data->selection.selection_mask, 0, data->mold.sys.atom.count);
            }
            if (ImGui::IsItemHovered()) {
                md_bitmask_assign(&data->selection.highlight_mask, md_bitmask_view(&data->selection.selection_mask));
                md_bitmask_invert_range(&data->selection.highlight_mask, 0, data->mold.sys.atom.count);
            }
            if (ImGui::MenuItem("Query")) data->selection.query.show_window = true;
            if (num_selected_atoms == 0) ImGui::PushDisabled();
            if (ImGui::MenuItem("Grow"))  data->selection.grow.show_window = true;
            if (num_selected_atoms == 0) ImGui::PopDisabled();
            if (ImGui::MenuItem("Clear")) {
                md_bitmask_clear(&data->selection.selection_mask);
                single_selection_sequence_clear(&data->selection.single_selection_sequence);
            }
            ImGui::Spacing();
            ImGui::Separator();

            // STORED SELECTIONS
            {
                md_bitmask_clear(&data->selection.highlight_mask);
                // @NOTE(Robin): This ImGui ItemFlag can be used to force the menu to remain open after buttons are pressed.
                // Leave it here as a comment if we feel that it is needed in the future
                //ImGui::PushItemFlag(ImGuiItemFlags_SelectableDontClosePopup, true);

                ImGui::Text("Stored Selections");
                for (int i = 0; i < (int)md_array_size(data->selection.stored_selections); i++) {
                    auto& sel = data->selection.stored_selections[i];
                    const str_t name_str = str_from_cstr(sel.name);
                    bool is_valid = md_script_identifier_name_valid(name_str);
                    char error[64] = "";
                    if (!is_valid) {
                        snprintf(error, sizeof(error), "'%s' is not a valid identifier.", sel.name);
                    }

                    for (int j = 0; j < i; ++j) {
                        if (str_eq_cstr(name_str, data->selection.stored_selections[j].name)) {
                            is_valid = false;
                            snprintf(error, sizeof(error), "identifier '%s' is already taken.", sel.name);
                            break;
                        }
                    }

                    ImGui::PushID(i);
                    ImGui::InputQuery("##label", sel.name, sizeof(sel.name), is_valid, error);
                    ImGui::SameLine();
                    if (ImGui::Button("Load")) {
                        md_bitmask_assign(&data->selection.selection_mask, sel.atom_mask);
                        flag_all_representations_as_dirty(data);
                    }
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("Load the stored selection as the active selection");
                        md_bitmask_assign(&data->selection.highlight_mask, sel.atom_mask);
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Store")) {
                        ImGui::SetTooltip("Store the active selection into this selection");
                        md_bitset_assign(&sel.atom_mask, md_bitmask_view(&data->selection.selection_mask), persistent_alloc);
                        data->script.compile_ir = true;
                        flag_all_representations_as_dirty(data);
                    }
                    ImGui::SameLine();
                    if (ImGui::DeleteButton("Remove")) {
                        ImGui::SetTooltip("Remove this selection");
                        remove_selection(data, i);
                    }
                    ImGui::PopID();
                }

                if (ImGui::Button("Create New")) {
                    char name_buf[64];
                    snprintf(name_buf, sizeof(name_buf), "sel%i", (int)md_array_size(data->selection.stored_selections) + 1);
                    create_selection(data, str_from_cstr(name_buf), md_bitmask_view(&data->selection.selection_mask));
                }

                //ImGui::PopItemFlag();
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Screenshot")) {
            ImGui::Checkbox("Hide GUI", &data->screenshot.hide_gui);
            if (data->screenshot.hide_gui) {
                data->screenshot.resolution = (ScreenshotResolution)MIN((int)data->screenshot.resolution, (int)ScreenshotResolution::Count - 1);
                if (ImGui::BeginCombo("Resolution", screenshot_resolution_str[(int)data->screenshot.resolution])) {
                    for (int i = 0; i < (int)ScreenshotResolution::Count; ++i) {
                        if (ImGui::Selectable(screenshot_resolution_str[i], (i == (int)data->screenshot.resolution))) {
                            data->screenshot.resolution = (ScreenshotResolution)i;
                        }
                    }
                    ImGui::EndCombo();
                }

                switch (data->screenshot.resolution) {
                case ScreenshotResolution::Window:
                    data->screenshot.res_x = data->app.framebuffer.width;
                    data->screenshot.res_y = data->app.framebuffer.height;
                    break;
                case ScreenshotResolution::FHD:
                    data->screenshot.res_x = 1920;
                    data->screenshot.res_y = 1080;
                    break;
                case ScreenshotResolution::QHD:
                    data->screenshot.res_x = 2560;
                    data->screenshot.res_y = 1440;
                    break;
                case ScreenshotResolution::UHD_4K:
                    data->screenshot.res_x = 3840;
                    data->screenshot.res_y = 2160;
                    break;
                case ScreenshotResolution::UHD_8K:
                    data->screenshot.res_x = 7680;
                    data->screenshot.res_y = 4320;
                    break;
                case ScreenshotResolution::Custom:
                    ImGui::InputInt("Res X", &data->screenshot.res_x);
                    ImGui::InputInt("Res Y", &data->screenshot.res_y);

                    data->screenshot.res_x = CLAMP(data->screenshot.res_x, 640, 16384);
                    data->screenshot.res_y = CLAMP(data->screenshot.res_y, 480, 16384);
                    break;
                default:
                    ASSERT(false);
                }
            } else {
                data->screenshot.res_x = data->app.framebuffer.width;
                data->screenshot.res_y = data->app.framebuffer.height;
            }
            if (ImGui::MenuItem("Take Screenshot")) {
                if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("jpg,jpeg,png,bmp"))) {
                    size_t path_len = strnlen(path_buf, sizeof(path_buf));
                    str_t ext;
                    if (!extract_ext(&ext, {path_buf, path_len})) {
                        path_len += snprintf(path_buf + path_len, sizeof(path_buf) - path_len, ".jpg");
                        ext = STR_LIT("jpg");
                    }
                    if (image_format_from_ext(ext) != ImageFormat_Unknown) {
                        data->screenshot.path_to_file = str_copy({path_buf, path_len}, persistent_alloc);
                    }
                    else {
                        VIAMD_LOG_ERROR("Supplied image extension is not supported");
                    }
                }
                ImGui::GetCurrentWindow()->Hidden = true;
            }
            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Operations")) {
            bool do_recenter = false;
            bool do_pbc = false;
            bool do_unwrap = false;
            bool do_bonds = false;
            bool redo_every_frame = false;
            auto& tracking = data->operations.bond_tracking;

            // Each operation can be applied to the frame shown now, or to every frame as it is shown.
            // The labels say what happens to the atoms, not what the operation is called internally.
            // Those which move the atoms are not, while an assembly is shown: its operators place the coordinates as
            // the file gives them (view.assembly)
            const bool assembly_shown = data->view.assembly >= 0;
            if (assembly_shown) {
                ImGui::TextDisabled("Centering, wrapping and making whole are not applied\nwhile an assembly is shown (View)");
            }
            ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg;
            if (ImGui::BeginTable("##table", 3, flags)) {
                if (ImGui::IsWindowHovered()) {
                    md_bitmask_clear(&data->selection.highlight_mask);
                }
                ImGui::TableSetupColumn("Now");
                ImGui::TableSetupColumn("Every frame");
                ImGui::TableSetupColumn("Options");
                ImGui::TableHeadersRow();
                ImGui::BeginDisabled(assembly_shown);

                const float button_width = (ImGui::GetFontSize() / 20.f) * 175.f;
                const md_bitset_t& target_mask = recenter_get_active_target_mask(data);
                const bool recenter_available = !md_bitset_empty(target_mask);
                const char* no_target_tip =
                    "There is nothing to center on yet. Select atoms and choose 'Set as Centering Target' in the\n"
                    "right-click menu of the viewport, or tick " ICON_FA_COMMENT_DOTS " and pick the target with a query.";

                // ## Centering
                ImGui::TableNextRow();

                ImGui::TableSetColumnIndex(0);
                if (!recenter_available) {
                    ImGui::PushDisabled();
                }
                if (ImGui::Button("Center Target", ImVec2(button_width,0))) {
                    do_recenter = true;
                }
                if (!recenter_available) {
                    ImGui::PopDisabled();
                }
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                    if (!recenter_available) {
                        ImGui::SetTooltip("%s", no_target_tip);
                    } else {
                        ImGui::SetTooltip("Move everything so the target (highlighted) is in the middle of the box,\nor at the origin when there is no box");
                        // Highlight the target atoms that the system will be recentered around
                        md_bitmask_assign(&data->selection.highlight_mask, target_mask);
                    }
                }

                ImGui::TableSetColumnIndex(1);
                if (!recenter_available) {
                    ImGui::PushDisabled();
                }
                if (ImGui::Checkbox("##recenter", &data->operations.recenter) && data->operations.recenter) {
                    do_recenter = true;
                }
                if (!recenter_available) {
                    ImGui::PopDisabled();
                }
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) {
                    ImGui::SetTooltip("%s", recenter_available ? "Keep the target in the middle of the box in every frame" : no_target_tip);
                }

                ImGui::TableSetColumnIndex(2);
                if (ImGui::Checkbox(ICON_FA_ANCHOR_LOCK "##keep-orientation", &data->operations.fixate_orientation) && data->operations.recenter) {
                    // Turning it on or off while centering every frame takes effect on the frame shown, not the next one
                    redo_every_frame = true;
                }
                ImGui::SetItemTooltip("Keep orientation: when centering, also turn everything so the target\nkeeps the orientation it has in the first frame.\nEverything is wrapped into the box around the target before it is turned,\nand the box turns with it.");

                ImGui::SameLine();
                ImGui::Checkbox(ICON_FA_COMMENT_DOTS "##target-query", &data->operations.recenter_query.enabled);
                ImGui::SetItemTooltip("Target by query: center on the atoms a query picks (evaluated every frame\nwhen it depends on it) instead of the target set from a selection");

                if (data->operations.recenter_query.enabled) {
                    auto& recenter_query = data->operations.recenter_query;
                    ImGui::SameLine();
                    if (filter_input("##recenter-query", &recenter_query.filter, filter_context(data))) {
                        recenter_mark_query_dirty(data);
                    }
                }

                // ## Periodic box
                ImGui::TableNextRow();

                ImGui::TableSetColumnIndex(0);
                if (ImGui::Button("Wrap into Box", ImVec2(button_width,0))) {
                    do_pbc = true;
                }
                ImGui::SetItemTooltip("Move every atom that is outside the periodic box back in through the opposite side\n(apply periodic boundary conditions). Molecules across the boundary are split.");

                ImGui::TableSetColumnIndex(1);
                if (ImGui::Checkbox("##pbc", &data->operations.apply_pbc) && data->operations.apply_pbc) {
                    do_pbc = true;
                }
                ImGui::SetItemTooltip("Wrap every atom into the periodic box in every frame");

                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                if (ImGui::Button("Make Whole", ImVec2(button_width, 0))) {
                    do_unwrap = true;
                }
                ImGui::SetItemTooltip("Join molecules that are split across the periodic box boundary,\nso each one is drawn in one piece (unwrap)");

                ImGui::TableSetColumnIndex(1);
                if (ImGui::Checkbox("##unwrap", &data->operations.unwrap_structures) && data->operations.unwrap_structures) {
                    do_unwrap = true;
                }
                ImGui::SetItemTooltip("Keep every molecule whole in every frame");

                ImGui::EndDisabled();

                // ## Bonds
                ImGui::TableNextRow();
                ImGui::TableSetColumnIndex(0);
                if (ImGui::Button("Recompute Bonds", ImVec2(button_width, 0))) {
                    do_bonds = true;
                }
                ImGui::SetItemTooltip("Guess the covalent bonds again from the distances between atoms in this frame, for when\n"
                                      "the file's bonds are wrong or a reaction should be the topology from here on.\n"
                                      "Changes the bonds of the system: the bonds a file, a topology or you gave stay.");

                ImGui::TableSetColumnIndex(1);
                if (ImGui::Checkbox("##track-bonds", &tracking.enabled)) {
                    tracking.scope_dirty = true;
                }
                ImGui::SetItemTooltip("Track bond changes: in every frame, show the bonds that form and break against the bonds\n"
                                      "of the system. The bonds of the system stay as they are.");

                ImGui::TableSetColumnIndex(2);
                if (tracking.enabled) {
                    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 12.0f);
                    if (filter_input("##track-bonds-scope", &tracking.scope, filter_context(data), 0, nullptr, "all atoms")) {
                        tracking.scope_dirty = true;
                    }
                    if (ImGui::IsItemHovered() && !filter_completions_open() && tracking.scope.result.success) {
                        ImGui::SetTooltip("Where bonds are tracked: the bonds with an atom among those it picks.\nEmpty for every atom; a smaller scope is faster.");
                    }
                }
                ImGui::EndTable();
            }

            if (tracking.enabled) {
                const md_bond_delta_t* delta = bond_tracking_delta(data);
                if (delta) {
                    const size_t num_broken = md_bitset_count(delta->broken);
                    ImGui::Text("%zu formed, %zu broken in this frame (%.1f ms)", delta->num_formed, num_broken, tracking.time_ms);
                    if (ImGui::IsItemHovered()) {
                        // The atoms of the bonds that changed
                        md_bitmask_clear(&data->selection.highlight_mask);
                        for (size_t i = 0; i < delta->num_formed; ++i) {
                            md_bitmask_set_bit(&data->selection.highlight_mask, (uint32_t)delta->formed_pair[i].idx[0]);
                            md_bitmask_set_bit(&data->selection.highlight_mask, (uint32_t)delta->formed_pair[i].idx[1]);
                        }
                        md_bitset_iter_t it = md_bitset_iter(delta->broken);
                        while (md_bitset_iter_next(&it)) {
                            if (it.idx >= data->mold.sys.bond.count) break;
                            const md_atom_pair_t p = data->mold.sys.bond.pairs[it.idx];
                            md_bitmask_set_bit(&data->selection.highlight_mask, (uint32_t)p.idx[0]);
                            md_bitmask_set_bit(&data->selection.highlight_mask, (uint32_t)p.idx[1]);
                        }
                        ImGui::SetTooltip("The bonds formed are drawn as sticks and the broken ones are not drawn.\n"
                                          "Their atoms are highlighted.");
                    }
                } else if (!tracking.scope.result.success) {
                    ImGui::TextDisabled("The scope does not pick atoms: nothing is tracked");
                } else {
                    ImGui::TextDisabled("Nothing to track in this system");
                }

                ImGui::SameLine();
                if (ImGui::SmallButton("Copy to Script")) {
                    // The same as a script value, to count or plot: formed(d), broken(d), bond(d)
                    char args[512] = "";
                    int len = 0;
                    const str_t scope = str_trim(str_from_cstr(tracking.scope.text));
                    if (!str_empty(scope)) len += snprintf(args + len, sizeof(args) - len, STR_FMT, STR_ARG(scope));
                    if (tracking.tolerance != MD_BOND_DELTA_TOLERANCE && len < (int)sizeof(args)) {
                        len += snprintf(args + len, sizeof(args) - len, "%stolerance=%.3g", len ? ", " : "", tracking.tolerance);
                    }
                    if (tracking.breakable != MD_BOND_DELTA_BREAKABLE && len < (int)sizeof(args)) {
                        len += snprintf(args + len, sizeof(args) - len, "%sbreakable={", len ? ", " : "");
                        bool first = true;
                        for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED && len < (int)sizeof(args); ++o) {
                            if (!(tracking.breakable & (1u << o))) continue;
                            len += snprintf(args + len, sizeof(args) - len, "%s'%s'", first ? "" : ", ", md_bond_origin_name((md_bond_origin_t)o));
                            first = false;
                        }
                        if (len < (int)sizeof(args)) len += snprintf(args + len, sizeof(args) - len, "}");
                    }
                    const str_t ident = create_unique_identifier(data->script.ir, STR_LIT("bonds_delta"), data->allocator.frame);
                    char line[640];
                    snprintf(line, sizeof(line), STR_FMT " = bond_delta(%s);", STR_ARG(ident), args);
                    script_editor::append_line(data->editor, str_from_cstr(line));
                }
                ImGui::SetItemTooltip("Add the same bond changes to the script, to count them over the frames\n"
                                      "(count(formed(d))) or select them (bond(d) is the bonds of each frame)");

                if (ImGui::TreeNode("Advanced##track-bonds")) {
                    float tolerance = tracking.tolerance;
                    if (ImGui::SliderFloat("tolerance", &tolerance, 0.0f, 0.45f, "%.2f", ImGuiSliderFlags_AlwaysClamp)) {
                        tracking.tolerance = tolerance;
                    }
                    ImGui::SetItemTooltip("A bond breaks beyond (1 + tolerance) times its cutoff and forms within (1 - tolerance) times it,\n"
                                          "so that it does not flicker as it vibrates around the cutoff");
                    ImGui::TextUnformatted("breakable:");
                    for (int o = MD_BOND_ORIGIN_FILE; o <= MD_BOND_ORIGIN_INFERRED; ++o) {
                        ImGui::SameLine();
                        bool on = (tracking.breakable & (1u << o)) != 0;
                        if (ImGui::Checkbox(md_bond_origin_name((md_bond_origin_t)o), &on)) {
                            tracking.breakable = on ? (tracking.breakable | (1u << o)) : (tracking.breakable & ~(1u << o));
                        }
                    }
                    ImGui::SetItemTooltip("The bonds of which origins may break. A force field topology holds throughout a\n"
                                          "classical simulation, and bonds you added stay unless you say otherwise.");
                    if (ImGui::SmallButton("Defaults")) {
                        tracking.tolerance = MD_BOND_DELTA_TOLERANCE;
                        tracking.breakable = MD_BOND_DELTA_BREAKABLE;
                    }
                    ImGui::TreePop();
                }
            }

            if (redo_every_frame) {
                do_recenter = true;
                do_pbc     |= data->operations.apply_pbc;
                do_unwrap  |= data->operations.unwrap_structures;
            }

            // One pass, in the order that keeps them valid together, and in the lattice frame of the cell
            // even when the coordinates are turned (see apply_state_operations)
            if (!assembly_shown && apply_state_operations(data, do_recenter, do_pbc, do_unwrap)) {
                data->mold.dirty_gpu_buffers |= MolBit_DirtyPosition | MolBit_ClearVelocity;
            }

            if (do_bonds) {
                // The whole frame nearest the animation time, or the coordinates shown when there is no run
                recompute_bonds(data, (int64_t)(data->animation.frame + 0.5));
            }

            ImGui::EndMenu();
        }
        if (ImGui::BeginMenu("Settings")) {
            if (ImGui::Checkbox("Keep Representations", &data->settings.keep_representations)) {
                app_settings::mark_dirty();
            }
            ImGui::SetItemTooltip("Keep representations when loading new topology (does not apply for workspaces)\n");

            if (ImGui::Checkbox("Focus in Free Space", &data->settings.focus_in_free_space)) {
                if (!data->settings.focus_in_free_space) {
                    // Back to the middle of the whole view
                    data->view.region = {};
                }
                app_settings::mark_dirty();
            }
            ImGui::SetItemTooltip("Double clicking the view frames what it focuses on in the largest part of the view\nthat no window covers, rather than in the middle of the whole view");

            // Font
            int font_size_idx = nearest_font_size_index(data->settings.font_size);
            if (ImGui::Combo("Font Size", &font_size_idx, font_size_names, (int)ARRAY_SIZE(font_size_names))) {
                data->settings.font_size = font_sizes[font_size_idx];
                apply_font_size(data);
                app_settings::mark_dirty();
            }

			if (ImGui::TreeNode("Units")) {
                display_units::draw_settings_menu_items();
				ImGui::TreePop();
			}

            ImGui::EndMenu();
        }
        {
            // Fps counter
            static int num_frames = 0;
            static double acc_ms  = 0;
            static double avg_ms  = 0;

            double ms = (data->app.timing.delta_s * 1000);
            acc_ms += ms;
            num_frames += 1;

            if (acc_ms > 500) {
                avg_ms = acc_ms / num_frames;
                acc_ms = 0;
                num_frames = 0;
            }

            char fps_buf[64];
            snprintf(fps_buf, ARRAY_SIZE(fps_buf), "%.2f ms (%.1f fps)", avg_ms, 1000.f / avg_ms);
            const float w = ImGui::CalcTextSize(fps_buf).x;
            ImGui::SetCursorPosX(ImGui::GetWindowContentRegionMax().x - w);
            ImGui::Text("%s", fps_buf);
            if (ImGui::IsItemHovered()) {
                // Where the frame goes, on the GPU and the CPU (all of it: Windows > Profiler)
                ImGui::BeginTooltip();
                profiler::draw_summary();
                ImGui::EndTooltip();
            }
        }
        ImGui::EndMainMenuBar();
    }

    if (new_clicked) ImGui::OpenPopup("Warning New");
}

void draw_notifications_window() {
    // Render toasts on top of everything, at the end of your code!
    // You should push style vars here
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 5.f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(43.f / 255.f, 43.f / 255.f, 43.f / 255.f, 100.f / 255.f));
    ImGui::RenderNotifications();
    ImGui::PopStyleVar(1); // Don't forget to Pop()
    ImGui::PopStyleColor(1);
}

void draw_load_dataset_window(ApplicationState* data) {
    LoadDatasetWindowState& state = data->load_dataset;
    
    ImGui::SetNextWindowSize(ImVec2(300, 215), ImGuiCond_FirstUseEver);
    if (state.show_window) {
        ImGui::OpenPopup("Load Dataset");
    }

    if (ImGui::BeginPopupModal("Load Dataset", &state.show_window)) {
        const bool path_invalid = !state.path_is_valid && state.path_buf[0] != '\0';
        const int  loader_count = LoaderType_COUNT;

        if (path_invalid) ImGui::PushInvalid();
        if (ImGui::InputText("##path", state.path_buf, sizeof(state.path_buf))) {
            state.path_changed = true;
        }
        if (path_invalid) ImGui::PopInvalid();


        // @WORKAROUND(Robin): This show_file_dialog is only here to circumvent the issue that if you open a file dialog
        // Within the same frame as the button is clicked, the dialogue will open again after the closing the dialog.
        ImGui::SameLine();
        if (ImGui::Button("Browse") && !state.show_file_dialog) {
            state.show_file_dialog = true;
        }

        if (state.show_file_dialog) {
            state.show_file_dialog = false;
            if (application::file_dialog(state.path_buf, sizeof(state.path_buf), application::FileDialogFlag_Open)) {
                state.path_changed = true;
            }
        }

        str_t path = str_from_cstr(state.path_buf);

        if (state.path_changed) {
            state.path_changed = false;
            state.path_is_valid = md_path_is_valid(path) && !md_path_is_directory(path);

            // Try to assign loader_idx from extension
            state.loader_idx = 0;
            str_t ext;
            if (extract_ext(&ext, path)) {
                state.loader_idx = (int)loader::type_from_ext(ext);
            }
        }

        if (ImGui::BeginCombo("Loader", loader::type_name(state.loader_idx).ptr)) {
            for (int i = 1; i < loader_count; ++i) {
                const char* str = loader::type_name(i).ptr;
                if (!str) continue;
                if (ImGui::Selectable(str, state.loader_idx == i)) {
                    state.loader_idx = i;
                }
            }
            ImGui::EndCombo();
        }

        // True if the button should be enabled
        bool load_button_enabled = (state.path_is_valid && state.loader_idx > -1);

        LoaderType   type = (LoaderType)(state.loader_idx);
        LoaderFlags flags = loader::type_flags(type);

        // Draw Options
        bool show_lammps_atom_format = state.path_is_valid && (type == LoaderType_LAMMPSDATA);
        if (show_lammps_atom_format) {
            const char** atom_format_names = md_lammps_atom_format_names();
            const char** atom_format_strings = md_lammps_atom_format_strings();

            if (state.atom_format_idx == -1) {
                // Try to determine format from file
                md_lammps_atom_format_t format = md_lammps_atom_format_from_file(path);
                state.atom_format_idx = format;
				const char* format_str = atom_format_strings[format];
                snprintf(state.atom_format_buf, sizeof(state.atom_format_buf), "%s", format_str);
            }
            state.atom_format_idx = CLAMP(state.atom_format_idx, 0, MD_LAMMPS_ATOM_FORMAT_COUNT - 1);

            if (ImGui::BeginCombo("Atom Format", state.atom_format_idx > 0 ? atom_format_names[state.atom_format_idx] : "user defined")) {
                for (int i = 0; i < MD_LAMMPS_ATOM_FORMAT_COUNT; ++i) {
                    if (ImGui::Selectable(i > 0 ? atom_format_names[i] : "user defined", state.atom_format_idx == i)) {
                        state.atom_format_idx = i;
                        int source_idx = i > 0 ? i : MD_LAMMPS_ATOM_FORMAT_FULL;
                        const char* format_str = atom_format_strings[source_idx];
                        snprintf(state.atom_format_buf, sizeof(state.atom_format_buf), "%s", format_str);
                    }
                }
                ImGui::EndCombo();
            }

            if (state.atom_format_idx == MD_LAMMPS_ATOM_FORMAT_UNKNOWN) {
                bool valid = state.atom_format_valid;
                if (!valid) ImGui::PushInvalid();
                if (ImGui::InputText("##atom_format", state.atom_format_buf, sizeof(state.atom_format_buf))) {
                    state.atom_format_valid = md_lammps_validate_atom_format(state.err_buf, sizeof(state.err_buf), state.atom_format_buf);
                }
                if (!valid) ImGui::PopInvalid();
                if (ImGui::IsItemHovered() && !valid) {
                    ImGui::SetTooltip("%s", state.err_buf);
                }
            }
            else {
                ImGui::PushDisabled();
                ImGui::InputText("##atom_format", state.atom_format_buf, sizeof(state.atom_format_buf), ImGuiInputTextFlags_ReadOnly);
                ImGui::PopDisabled();
            }

            if (state.atom_format_idx < 0 || (state.atom_format_idx == MD_LAMMPS_ATOM_FORMAT_UNKNOWN && !state.atom_format_valid)) {
                load_button_enabled = false;
            }
        }

        enum Action {
            Action_None,
            Action_Cancel,
            Action_Load,
        };
        Action action = Action_None;

        if (!load_button_enabled) ImGui::PushDisabled();
        if (ImGui::Button("Load")) {
            action = Action_Load;
        }
        if (!load_button_enabled) ImGui::PopDisabled();

        ImGui::SameLine();
        if (ImGui::Button("Cancel")) {
            action = Action_Cancel;
        }

        switch (action) {
        case Action_Load: {
            loader::LoaderState load_state = {};
            load_state.type = type;
            load_state.flags = flags;
            if (type == LoaderType_LAMMPSDATA) {
                load_state.arg = state.atom_format_buf;
            }
            // A loader other than the extension would pick has to be stored with the file in a workspace
            str_t ext = {};
            extract_ext(&ext, path);
            if (type != loader::type_from_ext(ext)) {
                load_state.flags |= LoaderFlag_ExplicitType;
            }

            const LoadResult result = load_data_from_file(data, path, load_state);
            if (result == LoadResult_System || result == LoadResult_Trajectory) {
                if (result == LoadResult_System && !data->settings.keep_representations) {
                    remove_all_representations(data);
                    create_default_representations(data);
                }
                data->animation = {};
                reset_view(&data->view.target, data->mold.state, &data->representation.visibility_mask, data->view.camera.fov_y, data->view.region);
            }
            [[fallthrough]];
        }
        case Action_Cancel:
            // Reset state
            // @NOTE(Robin): Don't change this to {}, it won't work on GCC 9
            state = LoadDatasetWindowState();
            [[fallthrough]];
        case Action_None:
        default:
            break;
        }

        ImGui::EndPopup();
    }
}

/*
// TODO: Move these functions to dataset component
void clear_atom_elem_mappings(ApplicationState* data) {
    md_array_shrink(data->dataset.atom_element_remappings, 0);
}

AtomElementMapping* add_atom_elem_mapping(ApplicationState* data, str_t lbl, md_element_t elem) {
    // Check if we already have a mapping for the label -> overwrite
    size_t i = 0;
    for (; i < md_array_size(data->dataset.atom_element_remappings); ++i) {
        if (str_eq_cstr(lbl, data->dataset.atom_element_remappings[i].lbl)) break;
    }
    if (i == md_array_size(data->dataset.atom_element_remappings)) {
        AtomElementMapping mapping = {
            .elem = elem,
        };
        str_copy_to_char_buf(mapping.lbl, sizeof(mapping.lbl), lbl);
        md_array_push(data->dataset.atom_element_remappings, mapping, persistent_alloc);
        return md_array_last(data->dataset.atom_element_remappings);
    } else {
        data->dataset.atom_element_remappings[i].elem = elem;
        return &data->dataset.atom_element_remappings[i];
    }
}

void apply_atom_elem_mappings(ApplicationState* data) {
    if (data->mold.sys.atom.count == 0 || !data->mold.sys.atom.element) {
        return;
    }

    for (size_t j = 0; j < md_array_size(data->dataset.atom_element_remappings); ++j) {
        str_t lbl = str_from_cstr(data->dataset.atom_element_remappings[j].lbl);
        md_element_t elem = data->dataset.atom_element_remappings[j].elem;
        float radius = md_util_element_vdw_radius(elem);
        float mass = md_util_element_atomic_mass(elem);

        for (size_t i = 0; i < data->mold.sys.atom.count; ++i) {
            if (str_eq(lbl, data->mold.sys.atom.type[i])) {
                data->mold.sys.atom.element[i] = elem;
                data->mold.sys.atom.radius[i] = radius;
                data->mold.sys.atom.mass[i] = mass;
                data->mold.dirty_buffers |= MolBit_DirtyRadius;
            }
        }
    }
    md_system_t* mol = &data->mold.sys;
    
    
    md_array_free(mol->bond.pairs, data->mold.sys.alloc);
    md_array_free(mol->bond.order, data->mold.sys.alloc);

    md_array_free(mol->bond.conn.atom_idx, data->mold.sys.alloc);
    md_array_free(mol->bond.conn.bond_idx, data->mold.sys.alloc);

    md_index_data_free(&mol->molecule);
    md_ring_data_free(&mol->ring, data->mold.sys.alloc);
    
    md_system_state_t mol_state = {0};
    md_util_system_infer(mol, &mol_state, data->mold.sys.alloc, MD_UTIL_INFER_BOND_BIT | MD_UTIL_INFER_MOLECULE_BIT);
    data->mold.dirty_buffers |= MolBit_DirtyBonds;

    flag_all_representations_as_dirty(data);
}
*/

static int64_t find_identifier(const md_script_ir_t* ir, str_t ident) {
    const int64_t num_ident = md_script_ir_num_identifiers(ir);
    const str_t* idents = md_script_ir_identifiers(ir);
    for (int64_t i = 0; i < num_ident; ++i) {
        if (str_eq(ident, idents[i])) return i;
    }
    return -1;
}

static str_t create_unique_identifier(const md_script_ir_t* ir, str_t base, md_allocator_i* alloc) {
    char buf[128];
    for (int64_t i = 1; i < 1000; ++i) {
        int res = snprintf(buf, sizeof(buf), "%.*s%i", (int)base.len, base.ptr, (int)i);
        if (res > 0) {
            str_t ident = {buf, (size_t)res};
            if (find_identifier(ir, ident) == -1) {
                return str_copy(ident, alloc);
            }
        }
    }
    return str_t();
}

// The lines the Script menu offers (script_suggestions.h): the measurement of the atoms picked one by one, and the
// selection. Every one of them is evaluated against the system before it is offered, which is not something to do
// every frame the menu is open, so they are kept for as long as what they are of stays the same.
struct ScriptSuggestions {
    uint64_t key = 0;
    md_allocator_i* arena = nullptr;
    md_array(ScriptSuggestion) measure = 0;
    md_array(ScriptSuggestion) select = 0;
};

static const ScriptSuggestions& script_suggestions_get(ApplicationState* state, bool reset) {
    static ScriptSuggestions cache;

    const md_system_t& sys = state->mold.sys;
    const md_bitset_t sel = md_bitmask_view(&state->selection.selection_mask);
    const SingleSelectionSequence& seq = state->selection.single_selection_sequence;

    // The atoms of a measurement, in the order they were picked: the picks have to be what is selected, as the order
    // is unknown for atoms selected otherwise. Two atoms need no order.
    int32_t atoms[4];
    size_t num_atoms = 0;
    const size_t num_selected = md_bitset_count(sel);
    const size_t num_picked = single_selection_sequence_count(&seq);
    if (2 <= num_selected && num_selected <= 4) {
        bool picks = num_picked == num_selected;
        for (size_t i = 0; i < num_picked && picks; ++i) picks = md_bitset_test(sel, (uint32_t)seq.idx[i]);
        if (picks) {
            MEMCPY(atoms, seq.idx, num_picked * sizeof(int32_t));
            num_atoms = num_picked;
        } else if (num_selected == 2) {
            md_bitset_iter_extract_indices(atoms, 2, md_bitset_iter(sel));
            num_atoms = 2;
        }
    }

    uint64_t key = md_bitset_hash64(sel, 0);
    key = md_hash64(atoms, num_atoms * sizeof(int32_t), key);
    const uint64_t counts[3] = {sys.atom.count, sys.residue.count, sys.part.count};
    key = md_hash64(counts, sizeof(counts), key);
    // The names the lines assign to are the first ones the script does not have
    key = md_hash64_combine(key, md_script_ir_fingerprint(state->script.ir));

    if (reset || key != cache.key || !cache.arena) {
        if (!cache.arena) cache.arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(64));
        md_arena_allocator_reset(cache.arena);
        cache.measure = 0;
        cache.select = 0;
        cache.key = key;

        if (num_atoms) {
            static const str_t base[5] = {{}, {}, STR_INIT("dist"), STR_INIT("ang"), STR_INIT("dih")};
            const str_t ident = create_unique_identifier(state->script.ir, base[num_atoms], cache.arena);
            cache.measure = script_suggest_measurement(ident, atoms, num_atoms, sys, cache.arena);
        }
        if (num_selected) {
            const str_t ident = create_unique_identifier(state->script.ir, STR_LIT("sel"), cache.arena);
            cache.select = script_suggest_selection(ident, sel, sys, cache.arena);
        }
    }
    return cache;
}

static void draw_script_suggestion(ApplicationState* state, const ScriptSuggestion& s) {
    // A long line is cut short in the menu, and shown whole when hovered
    constexpr int max_label_len = 100;
    char label[max_label_len + 8];
    const bool cut = s.text.len > (size_t)max_label_len;
    if (cut) snprintf(label, sizeof(label), "%.*s...", max_label_len - 3, s.text.ptr);
    else snprintf(label, sizeof(label), STR_FMT, STR_ARG(s.text));

    if (ImGui::MenuItem(label)) {
        script_editor::append_line(state->editor, s.text);
        ImGui::CloseCurrentPopup();
    }
    if (ImGui::IsItemHovered()) {
        // What the line picks or measures is shown in the viewport
        script_visualize_str(state, s.text, VIS_FLAGS);
        if (cut || !str_empty(s.hint)) {
            ImGui::BeginTooltip();
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 40.0f);
            if (cut) {
                constexpr size_t max_tooltip_len = 1024;
                const size_t len = MIN(s.text.len, max_tooltip_len);
                ImGui::TextUnformatted(s.text.ptr, s.text.ptr + len);
                if (len < s.text.len) ImGui::TextDisabled("... (%zu characters)", s.text.len);
            }
            if (!str_empty(s.hint)) ImGui::TextDisabled(STR_FMT, STR_ARG(s.hint));
            ImGui::PopTextWrapPos();
            ImGui::EndTooltip();
        }
    }
}

// # context_menu
void draw_context_popup(ApplicationState* state, const PickingHit& hit) {
    // @Robin: Will hit ever be used here? (Possibly remove argument)
    (void)hit;
    ASSERT(state);

    if (!state->mold.sys.atom.count) return;

    const size_t num_atoms_selected = md_bitmask_count(&state->selection.selection_mask);

    if (ImGui::BeginPopup("Context Popup")) {
        if (num_atoms_selected == 2) {
            // A bond between the two, of a kind, or the removal of the one a user made
            int idx[2];
            md_bitset_iter_extract_indices(idx, 2, md_bitmask_iter(&state->selection.selection_mask));

            md_system_t& sys = state->mold.sys;
            md_bond_idx_t bond_idx = md_system_bond_find(&sys, idx[0], idx[1]);
            if (bond_idx == -1) bond_idx = md_bond_find_noncovalent(&sys.bond, idx[0], idx[1]);
            if (bond_idx == -1) {
                char buf[256];
                snprintf(buf, sizeof(buf), "Create Bond (%i, %i)", idx[0] + 1, idx[1] + 1);
                if (ImGui::BeginMenu(buf)) {
                    // Of the bond graph, or a non-covalent bond: those are drawn dashed, and are not followed by the
                    // structures, the rings or the bonds tracked. A kind chosen is kept as given.
                    struct BondKind {
                        const char* label;
                        md_bond_flags_t flags;
                        const char* tip;
                    };
                    static const BondKind kinds[] = {
                        {"Bond",          MD_BOND_FLAG_NONE,
                         "A bond of the bond graph: covalent, or coordinate between a metal and a non-metal.\nIt joins the molecules of its atoms."},
                        {"Hydrogen Bond", md_bond_flags_set_kind(MD_BOND_FLAG_KIND_GIVEN, MD_BOND_KIND_HYDROGEN),
                         "Non-covalent, drawn dashed. bond(kind='hydrogen') in the script."},
                        {"Interaction",   md_bond_flags_set_kind(MD_BOND_FLAG_KIND_GIVEN, MD_BOND_KIND_INTERACTION),
                         "Non-covalent, of no particular kind (a salt bridge, a halogen bond), drawn dashed.\nbond(kind='interaction') in the script."},
                    };
                    for (size_t k = 0; k < ARRAY_SIZE(kinds); ++k) {
                        if (ImGui::MenuItem(kinds[k].label)) {
                            const bool graph = !md_bond_is_noncovalent(kinds[k].flags);
                            system_bonds_edit_begin(state);
                            md_system_bond_insert(&sys, idx[0], idx[1], md_bond_flags_set_origin(kinds[k].flags, MD_BOND_ORIGIN_USER));
                            if (graph) md_util_system_infer_coordination(&sys);
                            system_bonds_changed(state, graph);
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::SetItemTooltip("%s", kinds[k].tip);
                    }
                    ImGui::EndMenu();
                }
            } else {
				md_bond_flags_t flags = md_system_bond_flags(&sys, bond_idx);
                const md_atom_pair_t pair = md_bond_pair(&sys.bond, bond_idx);
                char buf[256];
                snprintf(buf, sizeof(buf), "Show Bond (%i, %i) In", idx[0] + 1, idx[1] + 1);
                if (ImGui::BeginMenu(buf)) {
                    // A representation of its own, to style the bond by, or one of the Bonds representations of the
                    // whole system
                    if (ImGui::MenuItem("New Bonds Representation")) {
                        representation_show_bond(state, pair, flags, nullptr);
                        ImGui::CloseCurrentPopup();
                    }
                    ImGui::SetItemTooltip("A representation of this bond alone, drawn over the sticks, to style it by.\n"
                                          "Turn off 'non-covalent bonds' in a licorice or ball and stick to have it drawn only here.");
                    bool separated = false;
                    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
                        Representation& rep = state->representation.reps[i];
                        if (!representation_can_show_bond(*state, rep)) continue;
                        if (!separated) {
                            ImGui::Separator();
                            separated = true;
                        }
                        ImGui::PushID((int)rep.node.id);
                        const bool shows = representation_shows_bond(rep, pair);
                        if (ImGui::MenuItem(rep.name, nullptr, shows, !shows)) {
                            if (!representation_show_bond(state, pair, flags, &rep)) {
                                VIAMD_LOG_INFO("The filter of '%s' has no room for the bond", rep.name);
                            }
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::SetItemTooltip(shows ? "Shows it already" : "Add the bond to what its filter picks");
                        ImGui::PopID();
                    }
                    ImGui::EndMenu();
                }
                if (md_bond_origin(flags) == MD_BOND_ORIGIN_USER) {
                    snprintf(buf, sizeof(buf), "Remove Bond (%i, %i)", idx[0] + 1, idx[1] + 1);
                    if (ImGui::MenuItem(buf)) {
                        const bool graph = (size_t)bond_idx < sys.bond.count;
                        system_bonds_edit_begin(state);
                        md_system_bond_remove(&sys, bond_idx);
                        system_bonds_changed(state, graph);
                        ImGui::CloseCurrentPopup();
                    }
                }
            }
        }
        if (num_atoms_selected >= 1 && run_num_frames(state) > 0) {
            // The bond changes of the selection, written as the scope of the tracking: "x = <expression>;" of the
            // suggestions of the script. Written when it is asked for, it can take long for a large selection.
            auto selection_scope = [state]() -> str_t {
                const md_bitset_t sel = md_bitmask_view(&state->selection.selection_mask);
                const md_array(ScriptSuggestion) suggestions = script_suggest_selection(STR_LIT("x"), sel, state->mold.sys, frame_alloc);
                for (size_t i = 0; i < md_array_size(suggestions); ++i) {
                    str_t expr = suggestions[i].text;
                    if (!str_begins_with(expr, STR_LIT("x = "))) continue;
                    expr = str_trim(str_substr(expr, 4));
                    if (str_ends_with(expr, STR_LIT(";"))) expr = str_trim(str_substr(expr, 0, expr.len - 1));
                    if (!str_empty(expr) && expr.len < sizeof(state->operations.bond_tracking.scope.text)) return expr;
                }
                return {};
            };
            if (ImGui::MenuItem("Track Bond Changes Here")) {
                const str_t scope = selection_scope();
                if (!str_empty(scope)) {
                    auto& tracking = state->operations.bond_tracking;
                    tracking.enabled = true;
                    str_copy_to_char_buf(tracking.scope.text, sizeof(tracking.scope.text), scope);
                    tracking.scope_dirty = true;
                } else {
                    VIAMD_LOG_INFO("The selection is too scattered to write as the scope of the bond changes tracked; write one in Operations");
                }
                ImGui::CloseCurrentPopup();
            }
            ImGui::SetItemTooltip("Show the bonds that form and break in every frame, for the bonds with an atom selected\n"
                                  "(Operations > Track bond changes, with the selection as its scope)");
        }
        if (ImGui::BeginMenu("Script")) {
            const ScriptSuggestions& sug = script_suggestions_get(state, ImGui::IsWindowAppearing());
            for (size_t i = 0; i < md_array_size(sug.measure); ++i) {
                ImGui::PushID((int)i);
                draw_script_suggestion(state, sug.measure[i]);
                ImGui::PopID();
            }
            if (md_array_size(sug.measure) && md_array_size(sug.select)) {
                ImGui::Separator();
            }
            for (size_t i = 0; i < md_array_size(sug.select); ++i) {
                ImGui::PushID((int)(1000 + i));
                draw_script_suggestion(state, sug.select[i]);
                ImGui::PopID();
            }
            if (!md_array_size(sug.measure) && !md_array_size(sug.select)) {
                ImGui::TextDisabled(num_atoms_selected ? "No suggestions for the current selection" : "Nothing is selected");
            }
            ImGui::EndMenu();
        }

        /*
        if (data->selection.atom_idx.right_click != -1 && data->mold.sys.atom.element) {
            int idx = data->selection.atom_idx.right_click;
            if (0 <= idx && idx < (int)data->mold.sys.atom.count) {
                char label[64] = "";
                str_t type = data->mold.sys.atom.type[idx];
                snprintf(label, sizeof(label), "Remap Element for '%.*s'", (int)type.len, type.ptr);
                if (ImGui::BeginMenu(label)) {
                    static char input_buf[32] = "";
                    md_element_t elem = data->mold.sys.atom.element[idx];
                    str_t name = md_util_element_name(elem);
                    str_t sym  = md_util_element_symbol(elem);

                    ImGui::Text("Current Element: %.*s (%.*s)", (int)name.len, name.ptr, (int)sym.len, sym.ptr);

                    str_t elem_str = {input_buf, strnlen(input_buf, sizeof(input_buf))};
                    md_element_t new_elem = md_util_element_lookup(elem_str);
                    const bool is_valid = new_elem != 0;

                    ImGui::InputQuery("##Symbol", input_buf, sizeof(input_buf), is_valid, "Cannot recognize Element symbol");
                    str_t new_name = md_util_element_name(new_elem);
                    str_t new_sym  = md_util_element_symbol(new_elem);
                    ImGui::Text("New Element: %.*s (%.*s)", (int)new_name.len, new_name.ptr, (int)new_sym.len, new_sym.ptr);
                    if (!is_valid) ImGui::PushDisabled();
                    if (ImGui::Button("Apply") && is_valid) {
                        add_atom_elem_mapping(data, type, new_elem);
                        apply_atom_elem_mappings(data);
                        ImGui::CloseCurrentPopup();
                    }
                    if (!is_valid) ImGui::PopDisabled();
                    ImGui::EndMenu();
                }
            }
        }
        */

        if (ImGui::BeginMenu("Selection")) {
            if (ImGui::MenuItem("Invert")) {
                md_bitmask_invert_range(&state->selection.selection_mask, 0, state->mold.sys.atom.count);
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsItemHovered()) {
                md_bitmask_assign(&state->selection.highlight_mask, md_bitmask_view(&state->selection.selection_mask));
                md_bitmask_invert_range(&state->selection.highlight_mask, 0, state->mold.sys.atom.count);
            }
            if (ImGui::MenuItem("Query")) {
                state->selection.query.show_window = true;
                ImGui::CloseCurrentPopup();
            }
            if (num_atoms_selected > 0) {
                if (ImGui::MenuItem("Grow")) {
                    state->selection.grow.show_window = true;
                    ImGui::CloseCurrentPopup();
                }
                if (ImGui::MenuItem("Clear")) {
                    md_bitmask_clear(&state->selection.selection_mask);
                    single_selection_sequence_clear(&state->selection.single_selection_sequence);
                }
            }
            ImGui::EndMenu();
        }
        
        if (num_atoms_selected > 0) {
            if (ImGui::MenuItem("Set as Centering Target")) {
                md_bitset_assign(&state->operations.selection_mask, md_bitmask_view(&state->selection.selection_mask), persistent_alloc);
                state->operations.recenter_query.enabled = false;
                recenter_update_target_data(state);
                ImGui::CloseCurrentPopup();
            }
            if (ImGui::IsItemHovered()) {
                md_bitmask_assign(&state->selection.highlight_mask, md_bitmask_view(&state->selection.selection_mask));
            }
        }

        ImGui::EndPopup();
    }
}

static void draw_selection_grow_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(300,150), ImGuiCond_Always);
    if (ImGui::Begin("Selection Grow", &data->selection.grow.show_window, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse)) {
        ImGui::PushItemWidth(-1);
        static uint64_t sel_popcount = 0;
        const uint64_t popcount = md_bitmask_count(&data->selection.selection_mask);
        const bool mode_changed = ImGui::Combo("##Mode", (int*)(&data->selection.grow.mode), "Covalent Bond\0Radial\0\0");
        const char* fmt = (data->selection.grow.mode == SelectionGrowthMode::CovalentBond) ? "%.0f" : "%.2f";
        const bool extent_changed = ImGui::SliderFloat("##Extent", &data->selection.grow.extent, 1.0f, 20.f, fmt);
        const bool appearing = ImGui::IsWindowAppearing();
        const bool sel_changed = popcount != sel_popcount;
        ImGui::PopItemWidth();

        const bool apply = ImGui::Button("Apply");

        // Need to invalidate when selection changes
        data->selection.grow.mask_invalid |= (mode_changed || extent_changed || appearing || sel_changed);

        if (data->selection.grow.mask_invalid) {
            sel_popcount = popcount;
            data->selection.grow.mask_invalid = false;
            md_bitmask_t* grow = &data->selection.grow.mask;
            md_bitmask_assign(grow, md_bitmask_view(&data->selection.selection_mask));
            const md_bitset_t viable = md_bitmask_view(&data->representation.visibility_mask);

            switch (data->selection.grow.mode) {
            case SelectionGrowthMode::CovalentBond:
                md_util_mask_grow_by_bonds(grow, &data->mold.sys, (int)data->selection.grow.extent, &viable);
                break;
            case SelectionGrowthMode::Radial: {
                md_util_mask_grow_by_radius(grow, &data->mold.state, data->selection.grow.extent, &viable);
                break;
            }
            default:
                ASSERT(false);
            }

            grow_mask_by_selection_granularity(grow, data->selection.granularity, data->mold.sys);
        }

        const bool show_preview =   (ImGui::GetHoveredID() == ImGui::GetID("##Extent")) ||
                                    (ImGui::GetActiveID()  == ImGui::GetID("##Extent")) ||
                                    (ImGui::GetHoveredID() == ImGui::GetID("Apply"));

        if (show_preview) {
            md_bitmask_assign(&data->selection.highlight_mask, md_bitmask_view(&data->selection.grow.mask));
        }
        if (apply) {
            md_bitmask_assign(&data->selection.selection_mask, md_bitmask_view(&data->selection.grow.mask));
            data->selection.grow.mask_invalid = true;
        }
    }
    ImGui::End();
}

static void draw_selection_query_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(300,100), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Selection Query", &data->selection.query.show_window, ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse)) {

        auto& query = data->selection.query;
        // Escape closes the window, unless it closes the suggestions of the query
        if (ImGui::IsKeyPressed(ImGuiKey_Escape, false) && !filter_completions_open()) {
            query.show_window = false;
        }

        static double query_frame = 0.0;

        ImGui::PushItemWidth(-1);
        bool apply = false;
        if (filter_input("##query", &query.filter, filter_context(data), ImGuiInputTextFlags_AutoSelectAll, &apply)) {
            query.query_invalid = true;
        }
        ImGui::PopItemWidth();

        // What depends on the coordinates is picked again in another frame
        if (query.filter.result.is_dynamic && data->animation.frame != query_frame) {
            query.query_invalid = true;
        }
        bool preview = ImGui::IsItemFocused() || ImGui::IsItemHovered();

        if (ImGui::IsWindowAppearing()) {
            ImGui::SetKeyboardFocusHere(-1);
        }

        if (!query.filter.result.success) ImGui::PushDisabled();
        apply |= ImGui::Button("Apply");
        if (!query.filter.result.success) ImGui::PopDisabled();

        preview |= ImGui::IsItemHovered();

        if (query.query_invalid) {
            query.query_invalid = false;
            md_bitset_t result = md_bitset_empty_set();
            filter_evaluate(&query.filter, &result, filter_context(data), frame_alloc);
            query_frame = data->animation.frame;

            md_bitmask_assign(&query.mask, query.filter.result.success ? result : md_bitset_empty_set());
            grow_mask_by_selection_granularity(&query.mask, data->selection.granularity, data->mold.sys);
        }

        if (preview) {
            md_bitmask_assign(&data->selection.highlight_mask, md_bitmask_view(&query.mask));
        }

        if (apply && query.filter.result.success) {
            md_bitmask_assign(&data->selection.selection_mask, md_bitmask_view(&data->selection.query.mask));
            data->selection.query.show_window = false;
        }
    }
    ImGui::End();
}

static void draw_animation_window(ApplicationState* data) {
    ASSERT(data);
    size_t num_frames = run_num_frames(data);
    if (num_frames == 0) return;

    ASSERT(data->timeline.x_values);
    ASSERT(md_array_size(data->timeline.x_values) == num_frames);

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Animation", &data->animation.show_window, ImGuiWindowFlags_NoFocusOnAppearing)) {
        ImGui::Text("Num Frames: %zu", num_frames);
        md_unit_t time_unit = data->timeline.time_unit;
        double t   = frame_to_time(data->animation.frame, *data);
        double min = data->timeline.x_values[0];
        double max = data->timeline.x_values[num_frames - 1];
        char time_label[64];
        if (md_unit_is_none(time_unit)) {
            snprintf(time_label, sizeof(time_label), "Time");
        } else {
            char unit_buf[32];
            md_unit_print(unit_buf, sizeof(unit_buf), time_unit);
            snprintf(time_label, sizeof(time_label), "Time (%s)", unit_buf);
        }
        const float w = ImGui::CalcTextSize("Time (ps)").x;
        const float item_width = MAX(ImGui::GetContentRegionAvail().x - w, 100.f);
        ImGui::PushItemWidth(item_width);
        if (ImGui::BeginCombo("Interp.", interpolation_mode_str[(int)data->animation.interpolation])) {
            static const char* mode_tooltip[(int)InterpolationMode::Count] = {
                "The nearest frame, as it is",
                "Straight from frame to frame, through each of them",
                "Smooth (C1) from frame to frame, through each of them",
                "Smoother (C2), and damps the jitter from frame to frame,\n"
                "but passes through none of the frames but the first and the last:\n"
                "the coordinates shown are never a frame's own, not even when stopped on one",
            };
            for (int i = 0; i < (int)InterpolationMode::Count; ++i) {
                if (ImGui::Selectable(interpolation_mode_str[i], (int)data->animation.interpolation == i)) {
                    data->animation.interpolation = (InterpolationMode)i;
					data->mold.interpolate_system_state = true;
                    data->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                }
                ImGui::SetItemTooltip("%s", mode_tooltip[i]);
            }
            ImGui::EndCombo();
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Interpolation Mode for Atom Positions");
        }
        if (ImGui::SliderScalar(time_label, ImGuiDataType_Double, &t, &min, &max, "%.2f")) {
            data->animation.frame = time_to_frame(t, data->timeline.x_values);
        }
        ImGui::SliderFloat("Speed", &data->animation.fps, -200.0f, 200.f, "%.2f", ImGuiSliderFlags_Logarithmic);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Animation Speed in Frames Per Second");
        }
        {
            // animation.stride: while playing, only every Nth frame is kept and the motion between them is interpolated.
            // Not called keyframes: that name is kept for the keyframes of animations to come.
            // The slider covers up to 10 (or the last frame, beyond which a larger stride changes nothing), a larger one is
            // typed in (Ctrl+click): only the lower bound is enforced.
            const int max_stride = (int)CLAMP((int64_t)num_frames - 1, (int64_t)2, (int64_t)10);
            const char* fmt = data->animation.stride == 1 ? "every frame" : "every %d frames";
            ImGui::SliderInt("Stride", &data->animation.stride, 1, max_stride, fmt);
            data->animation.stride = MAX(1, data->animation.stride);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("While playing, interpolate between every Nth frame only: smoother motion\n"
                                  "at the expense of the frames in between (and the usual interpolation artifacts).\n"
                                  "The speed is unchanged. A stopped animation shows the trajectory's own frames.\n"
                                  "Ctrl+click to type a value beyond the slider's range.");
            }
        }
        if (data->animation.interpolation == InterpolationMode::CubicSpline) {
            ImGui::SliderFloat("Tension", &data->animation.tension, 0.0f, 1.0f, "%.2f");
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Tension of the Cubic Spline (0 = Catmull-Rom, 1 = ease in and out of every frame)");
            }
        }
        switch (data->animation.mode) {
            case PlaybackMode::Playing:
                if (ImGui::Button((const char*)ICON_FA_PAUSE)) {
                    data->animation.mode = PlaybackMode::Stopped;
                }
                break;
            case PlaybackMode::Stopped:
                if (ImGui::Button((const char*)ICON_FA_PLAY)) {
                    data->animation.mode = PlaybackMode::Playing;
                }
                break;
            default:
                ASSERT(false);
        }
        ImGui::SameLine();
        if (ImGui::Button((const char*)ICON_FA_STOP)) {
            data->animation.mode = PlaybackMode::Stopped;
            data->animation.frame = 0.0;
        }
        ImGui::PopItemWidth();
    }
    ImGui::End();
}

static void draw_async_task_window(ApplicationState* data) {
    constexpr float WIDTH = 300.f;
    constexpr float MARGIN = 10.f;

    task_system::ID tasks[256];
    size_t num_tasks = task_system::pool_running_tasks(tasks, ARRAY_SIZE(tasks));
    bool any_task_label_visible = false;
    for (size_t i = 0; i < num_tasks; i++) {
        str_t label = task_system::task_label(tasks[i]);
        if (!label || label[0] == '\0' || (label[0] == '#' && label[1] == '#')) continue;
        any_task_label_visible = true;
    }
    // the isosurface meshes extracted on the GPU over frames (iso_meshes.h)
    float iso_fraction[8];
    const size_t num_iso = iso_meshes_progress(iso_fraction, ARRAY_SIZE(iso_fraction));
    any_task_label_visible |= num_iso > 0;
    // the SCF density matrix a volume waits for, made on a worker thread
    double dmat_seconds = 0.0;
    const bool dmat_preparing = density_dmat_preparing(data, &dmat_seconds) && dmat_seconds > 0.25;
    any_task_label_visible |= dmat_preparing;
    // and the fields surfaces are coloured by, evaluated over frames (surface_field.h)
    SurfaceFieldProgress field_progress[8];
    const size_t num_field = surface_field_progress(field_progress, ARRAY_SIZE(field_progress));
    any_task_label_visible |= num_field > 0;
    // and what the components have under way that the task system does not show, such as the analytic critical
    // point search of the QM module: asked for through an event, as a component may not be built
    TaskProgressRequest component_progress = {};
    viamd::event_system_broadcast_event(viamd::EventType_ViamdTaskProgressRequest, viamd::EventPayloadType_TaskProgressRequest, &component_progress);
    any_task_label_visible |= component_progress.count > 0;
    
    if (any_task_label_visible) {
        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->Pos + ImVec2(data->app.window.width - WIDTH - MARGIN,
                                                       ImGui::GetCurrentContext()->FontSize + ImGui::GetStyle().FramePadding.y * 2.f + MARGIN));
        ImGui::SetNextWindowSize(ImVec2(WIDTH, 0));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
        ImGui::Begin("##Async Info", 0,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                     ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing);

        const float pad = 3.0f;
        const float size = ImGui::GetFontSize() + pad * 2;

        char buf[64];
        for (size_t i = 0; i < MIN(num_tasks, 8); i++) {
            const auto id = tasks[i];
            str_t label = task_system::task_label(id);
            float fract = task_system::task_fraction_complete(id);

            if (!label || label[0] == '\0' || (label[0] == '#' && label[1] == '#')) continue;

            snprintf(buf, sizeof(buf), "%.*s %.1f%%", (int)label.len, label.ptr, fract * 100.f);
            ImGui::ProgressBar(fract, ImVec2(ImGui::GetContentRegionAvail().x - (size + pad),0), buf);
            ImGui::SameLine();
            if (ImGui::DeleteButton((const char*)ICON_FA_XMARK, ImVec2(size, size))) {
                task_system::task_interrupt(id);
                if (id == data->tasks.evaluate) {
                    md_script_eval_interrupt(data->script.eval);
                }
            }
        }
        for (size_t i = 0; i < num_iso; i++) {
            if (iso_fraction[i] < 0.0f) {
                // being prepared on a worker thread: an indeterminate bar
                ImGui::ProgressBar(-0.5f * (float)ImGui::GetTime() - 1.0e-3f, ImVec2(ImGui::GetContentRegionAvail().x - (size + pad), 0), "Isosurface mesh: preparing");
            } else {
                snprintf(buf, sizeof(buf), "Isosurface mesh %.0f%%", iso_fraction[i] * 100.f);
                ImGui::ProgressBar(iso_fraction[i], ImVec2(ImGui::GetContentRegionAvail().x - (size + pad), 0), buf);
            }
        }
        if (dmat_preparing) {
            ImGui::ProgressBar(-0.5f * (float)ImGui::GetTime() - 1.0e-3f, ImVec2(ImGui::GetContentRegionAvail().x - (size + pad), 0), "Density matrix: preparing");
        }
        for (size_t i = 0; i < num_field; i++) {
            const SurfaceFieldProgress& p = field_progress[i];
            if (p.preparing) {
                snprintf(buf, sizeof(buf), "%s: preparing", p.label);
                // no fraction to show: an indeterminate bar (a negative fraction moving with time), to say it is not stuck
                ImGui::ProgressBar(-0.5f * (float)ImGui::GetTime() - 1.0e-3f, ImVec2(ImGui::GetContentRegionAvail().x - (size + pad), 0), buf);
            } else {
                snprintf(buf, sizeof(buf), "%s %.0f%%", p.label, p.fraction * 100.f);
                ImGui::ProgressBar(p.fraction, ImVec2(ImGui::GetContentRegionAvail().x - (size + pad), 0), buf);
            }
        }
        for (size_t i = 0; i < component_progress.count; i++) {
            const TaskProgress& p = component_progress.items[i];
            snprintf(buf, sizeof(buf), "%s %.0f%%", p.label, p.fraction * 100.f);
            ImGui::ProgressBar(p.fraction, ImVec2(ImGui::GetContentRegionAvail().x - (size + pad), 0), buf);
            if (p.cancel) {
                ImGui::SameLine();
                ImGui::PushID((int)i);
                if (ImGui::DeleteButton((const char*)ICON_FA_XMARK "##component task", ImVec2(size, size))) {
                    *p.cancel = 1;
                }
                ImGui::PopID();
            }
        }

        ImGui::End();
        ImGui::PopStyleColor();
    }
}

struct TimelineArgs {
    const char* lbl;
    uint32_t col;
    int plot_height;

    struct {
        int count;
        int dim_y;

        const float* x;
        const float* y;
        const float* y_mean;
        const float* y_var;
        const float* y_min;
        const float* y_max;

        float min_y;
        float max_y;
        str_t unit;
    } values;

    struct {
        double* beg;
        double* end;
        double  min;
        double  max;
    } view_range;

    struct {
        bool* is_dragging;
        bool* is_selecting;
    } input;

    struct {
        bool show;
        bool enabled;

        double* beg;
        double* end;

        double min;
        double max;
    } filter;

    struct {
        bool enabled;
        double min;
        double max;
    } value_filter;

    double* time;
};

struct TimePayload {
    const TimelineArgs* args;
    int y_idx;
};

static ImPlotPoint get_time_point(int index, void* user_data) {
    const TimePayload* payload = (const TimePayload*)user_data;
    const TimelineArgs* args = payload->args;
    return ImPlotPoint(args->values.x[index], args->values.y[index * args->values.dim_y + payload->y_idx]);
}

bool draw_property_timeline(const ApplicationState& data, const TimelineArgs& args) {
    const ImPlotAxisFlags axis_flags = ImPlotAxisFlags_NoSideSwitch | ImPlotAxisFlags_NoHighlight;
    const ImPlotAxisFlags axis_flags_x = axis_flags;
    const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit | ImPlotAxisFlags_NoLabel |ImPlotAxisFlags_NoTickLabels;
    
    const ImPlotFlags flags = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;
    
    const float pad_x = ImPlot::GetStyle().PlotPadding.x;
    ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(pad_x, 0));
    defer { ImPlot::PopStyleVar(1); };

    if (ImPlot::BeginPlot("##Timeline", ImVec2(-1,args.plot_height), flags)) {
        ImPlot::SetupAxisLinks(ImAxis_X1, args.view_range.beg, args.view_range.end);
        ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, args.view_range.min, args.view_range.max);
        ImPlot::SetupAxes(0, 0, axis_flags_x, axis_flags_y);
        ImPlot::SetupFinish();

        bool active = ImGui::IsItemActive();
        bool print_timeline_tooltip = false;

        if (args.value_filter.enabled) {
            float* y_vals = (float*)md_alloc(frame_alloc, args.values.count * sizeof(float));
            for (int i = 0; i < args.values.count; ++i) {
                float val = args.values.y[i];
                y_vals[i] = (args.value_filter.min < val && val < args.value_filter.max) ? args.values.max_y : -FLT_MAX;
            }

            const ImVec4 filter_frame_color = ImVec4(1, 1, 1, 1);
            ImPlot::SetNextFillStyle(filter_frame_color, 0.15f);
            ImPlot::PlotShaded("##value_filter", args.values.x, y_vals, args.values.count, -FLT_MAX);
        }

        if (args.filter.show) {
            if (!args.filter.enabled) ImGui::PushDisabled();
            ImPlot::DragRangeX("Time Filter", args.filter.beg, args.filter.end, args.filter.min, args.filter.max);
            if (!args.filter.enabled) ImGui::PopDisabled();
            *args.filter.beg = CLAMP(*args.filter.beg, args.filter.min, args.filter.max);
            *args.filter.end = CLAMP(*args.filter.end, args.filter.min, args.filter.max);
        }

        ImVec4 line_col = ImGui::ColorConvertU32ToFloat4(args.col);
        ImPlot::SetNextLineStyle(line_col);

        if (args.values.count > 0) {
            ASSERT(args.values.x);
            ASSERT(args.values.y);
            if (args.values.y_var) {
                ASSERT(args.values.y_min);
                ASSERT(args.values.y_max);
                char lbl[32];

                snprintf(lbl, sizeof(lbl), "%s", args.lbl);
                TimePayload payload = {
                    .args = &args,
                    .y_idx = 0,
                };
                for (int i = 0; i < args.values.dim_y; ++i) {
                    payload.y_idx = i;
                    ImPlot::PlotLineG(lbl, get_time_point, &payload, args.values.count);
                }

                ImPlot::SetNextLineStyle(line_col);
                snprintf(lbl, sizeof(lbl), "%s (mean)", args.lbl);
                ImPlot::PlotLine(lbl, args.values.x, args.values.y_mean, args.values.count);

                ImPlot::SetNextFillStyle(line_col, 0.4f);
                snprintf(lbl, sizeof(lbl), "%s (var)", args.lbl);
                ImPlot::PlotShadedG(lbl,
                    [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y_mean[idx] - args->values.y_var[idx]);
                    },
                    (void*)&args,
                    [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y_mean[idx] + args->values.y_var[idx]);
                    },
                    (void*)&args,
                    args.values.count
                );

                ImPlot::SetNextFillStyle(line_col, 0.2f);
                snprintf(lbl, sizeof(lbl), "%s (min,max)", args.lbl);
                ImPlot::PlotShadedG(lbl,
                    [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y_min[idx]);
                    },
                    (void*)&args,
                        [](int idx, void* payload) -> ImPlotPoint {
                        TimelineArgs* args = (TimelineArgs*)payload;
                        return ImPlotPoint(args->values.x[idx], args->values.y_max[idx]);
                    },
                    (void*)&args,
                    args.values.count
                );
            } 
            else {
                ImPlot::PlotLine(args.lbl, args.values.x, args.values.y, args.values.count);
            }
        }
        
        if (*args.input.is_dragging) {
            *args.time = ImPlot::GetPlotMousePos().x;
        }
        else if (*args.input.is_selecting) {
            *args.filter.end = MAX(ImPlot::GetPlotMousePos().x, *args.filter.beg);
        }
        else if (ImPlot::IsPlotHovered()) {
            if (active && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {           
                if (ImGui::GetIO().KeyMods == ImGuiMod_Shift) {
                    if (args.filter.show && args.filter.enabled) {
                        *args.input.is_selecting = true;
                        *args.filter.beg = ImPlot::GetPlotMousePos().x;
                    }
                } else {
                    *args.input.is_dragging = true;
                }
            }
        }

        if (ImPlot::IsPlotHovered()) {
            print_timeline_tooltip = true;
        }
        
        if (ImPlot::DragLineX(0, args.time, ImVec4(1,1,0,1))) {
            *args.time = CLAMP(*args.time, args.filter.min, args.filter.max);
        }
        if (ImGui::IsItemHovered()) {
            print_timeline_tooltip = true;
        }

        if (print_timeline_tooltip) {
            ImPlotPoint plot_pos = ImPlot::GetPlotMousePos();
            ImVec2 screen_pos = ImPlot::PlotToPixels(plot_pos);
            ImVec2 p0 = {screen_pos.x, ImPlot::GetPlotPos().y};
            ImVec2 p1 = {screen_pos.x, ImPlot::GetPlotPos().y + ImPlot::GetPlotSize().y};
            ImPlot::PushPlotClipRect();
            ImPlot::GetPlotDrawList()->AddLine(p0, p1, IM_COL32(255, 255, 255, 120));
            ImPlot::PopPlotClipRect();

            char buf[128] = "";
            int len = 0;

            double time = plot_pos.x;
            int frame_idx = CLAMP((int)(time_to_frame(time, data.timeline.x_values) + 0.5), 0, (int)md_array_size(data.timeline.x_values)-1);
            len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), "time: %.2f", time);

            md_unit_t time_unit = data.timeline.time_unit;
            if (!md_unit_is_none(time_unit)) {
                char unit_buf[32];
                md_unit_print(unit_buf, sizeof(unit_buf), time_unit);
                len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), " (%s)", unit_buf);
            }

            if (0 <= frame_idx && frame_idx < args.values.count) {
                const char* value_lbl = args.values.y_var ? "mean" : "value";
                if (args.values.y) {
                    len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), ", %s: %.2f", value_lbl, args.values.y[frame_idx]);
                }
                if (args.values.y_var) {
                    ASSERT(args.values.y_min);
                    ASSERT(args.values.y_max);
                    len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), ", var: %.2f, min: %.2f, max: %.2f",
                        args.values.y_var[frame_idx],
                        args.values.y_min[frame_idx],
                        args.values.y_max[frame_idx]);
                }

                if (!str_empty(args.values.unit)) {
                    len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), " (%.*s)", (int)args.values.unit.len, args.values.unit.ptr);
                }
            }
            ImGui::SetTooltip("%.*s", len, buf);
        }

        ImPlot::EndPlot();
    }

    return true;
}

static double distance_to_linesegment(ImPlotPoint p0, ImPlotPoint p1, ImPlotPoint p) {
    double vx = p1.x - p0.x;
    double vy = p1.y - p0.y;

    double ux = p.x - p0.x;
    double uy = p.y - p0.y;

    double d_uv = ux*vx + uy*vy;
    double d_vv = vx*vx + vy*vy;

    if (d_vv < 1.0e-7) {
        double d_uu = ux*ux + uy*uy;
        return sqrt(d_uu);
    }

    double t = d_uv / d_vv;

    if (t < 0.0) {
        double d_uu = ux*ux + uy*uy;
        return sqrt(d_uu);
    } else if (t > 1.0) {
        double wx = p.x - p1.x;
        double wy = p.y - p1.y;
        double d_ww = wx*wx + wy*wy;
        return sqrt(d_ww);
    } else {
        double wx = p.x - (p0.x + vx * t);
        double wy = p.y - (p0.y + vy * t);
        double d_ww = wx*wx + wy*wy;
        return sqrt(d_ww);
    }
}

// #plots
//
// The Timelines and Distributions windows. Both draw series named by attribute path (see
// plot_series.h) into subplots that own them, and share the lists, the drag and drop and the
// styling of a legend entry.

// The colour of population member k of an entry drawn with a population of population_size
static ImVec4 plot_series_member_color(const PlotSeries& s, int k, int population_size) {
    ImVec4 color = s.color;
    if (s.use_colormap && population_size > 1) {
        if (ImPlot::ColormapQualitative(s.colormap)) {
            color = ImPlot::GetColormapColor(k, s.colormap);
        } else {
            color = ImPlot::SampleColormap((float)k / (float)(population_size - 1), s.colormap);
        }
        color.w *= s.colormap_alpha;
    }
    return color;
}

static ImVec4 plot_highlight(ImVec4 color) {
    const float scl = 1.5f;
    return ImVec4(ImSaturate(color.x * scl), ImSaturate(color.y * scl), ImSaturate(color.z * scl), color.w);
}

// The part of a legend entry's popup both windows share: colour or colormap, and which members of
// a population are drawn. A member hovered in the popup is written to hovered_pop_idx.
static void plot_series_style_popup(ApplicationState* data, PlotSeries& s, int dim, const char* script_ident, md_script_vis_ref_t vis_ref, int* hovered_pop_idx) {
    if (dim > 1) {
        const char* color_type_labels[] = {"Solid", "Colormap"};
        int color_type = s.use_colormap ? 1 : 0;
        if (ImGui::Combo("Color Type", &color_type, color_type_labels, IM_ARRAYSIZE(color_type_labels))) {
            s.use_colormap = color_type == 1;
        }
    }
    if (s.use_colormap && dim > 1) {
        ImPlot::ColormapSelection("##Colormap", &s.colormap);
        ImGui::SliderFloat("Alpha", &s.colormap_alpha, 0.0f, 1.0f);
    } else {
        ImGui::ColorEdit4("Color", &s.color.x);
    }

    if (dim > 1) {
        ImGui::Separator();
        if (ImGui::Button("Set All")) {
            s.population_mask.set();
        }
        ImGui::SameLine();
        if (ImGui::Button("Clear All")) {
            s.population_mask.reset();
        }

        const float sz = ImGui::GetFontSize() * 1.5f;
        ImGui::PushStyleVar(ImGuiStyleVar_SelectableTextAlign, ImVec2(0.5f, 0.5f));
        for (int k = 0; k < MIN(dim, PLOT_MAX_POPULATION); ++k) {
            char lbl[32];
            snprintf(lbl, sizeof(lbl), "%d", k+1);
            if (ImGui::Selectable(lbl, s.population_mask.test(k), ImGuiSelectableFlags_DontClosePopups, ImVec2(sz, sz))) {
                s.population_mask.flip(k);
            }
            const str_t member = series_member_label(data, s.key, k);
            if (!str_empty(member)) {
                ImGui::SetItemTooltip(STR_FMT, STR_ARG(member));
            }
            if (ImGui::IsItemHovered()) {
                if (script_ident[0] != '\0') {
                    script_visualize_ref(data, vis_ref, k, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                    script_set_hovered_property(data, str_from_cstr(script_ident), k);
                }
                *hovered_pop_idx = k;
            }
            if (!k || ((k+1) % 10)) {
                ImGui::SameLine();
            }
        }
        ImGui::PopStyleVar();
    }
}

// One entry in a list of series: a click toggles it in the first subplot, a drag puts it in any
static void plot_series_list_item(ApplicationState* data, const char* dnd_type, PlotSubplot& first, const SeriesKey& key) {
    char label[96];
    series_label(label, sizeof(label), data, key);
    const ImVec4 color = series_default_color(data, key);
    const int idx = plot_find_series(first, key);

    ImGui::PushID(key.path);
    ImGui::PushID((int)key.variant * SeriesSource_Count + (int)key.source);
    ImPlot::ItemIcon(color);
    ImGui::SameLine();
    if (ImGui::Selectable(label, idx != -1, ImGuiSelectableFlags_DontClosePopups)) {
        if (idx != -1) {
            plot_remove_series(first, idx);
        } else {
            plot_add_series(data, first, key);
        }
    }
    if (ImGui::IsItemHovered()) {
        const str_t ident = series_script_ident(key);
        if (!str_empty(ident)) {
            script_visualize_ref(data, md_script_ir_property_vis_ref(data->script.eval_ir, ident), -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
            script_set_hovered_property(data, ident);
        } else {
            ImGui::SetTooltip("%s\nClick to show in the first subplot, drag into any", key.path);
        }
    }
    if (ImGui::BeginDragDropSource()) {
        series_set_drag_payload(dnd_type, key, -1, label, color);
        ImGui::EndDragDropSource();
    }
    ImGui::PopID();
    ImGui::PopID();
}

// The series loaded along the run, one submenu per file
static void plot_system_series_menu(ApplicationState* data, const char* dnd_type, PlotSubplot& first) {
    str_t groups[64];
    const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), data), ARRAY_SIZE(groups));
    size_t num_listed = 0;
    for (size_t g = 0; g < num_groups; ++g) {
        md_temp_scope_t temp = md_temp_begin();
        const size_t num = system_series_members(nullptr, 0, data, groups[g]);
        str_t* paths = md_temp_alloc_array(temp, str_t, num + 1);
        system_series_members(paths, num, data, groups[g]);
        if (num > 0) {
            char group_label[256];
            system_series_group_label(group_label, sizeof(group_label), data, groups[g]);
            if (ImGui::BeginMenu(group_label)) {
                for (size_t i = 0; i < num; ++i) {
                    plot_series_list_item(data, dnd_type, first, series_key(SeriesSource_System, paths[i]));
                }
                ImGui::EndMenu();
            }
            num_listed += 1;
        }
        md_temp_end(temp);
    }
    if (num_listed == 0) {
        ImGui::TextDisabled("Nothing loaded");
    }
}

// The members of a script property's population, if it has one: the extent of its value axis
static size_t script_property_population(const ApplicationState* data, const SeriesKey& key) {
    const md_attributes_t* table = series_table(data, key.source, str_from_cstr(key.path));
    const md_attribute_t* attr = table ? md_attributes_find(table, str_from_cstr(key.path)) : nullptr;
    if (!attr || attr->format.rank < 1 || attr->format.shape[0] == 0) return 0;
    return md_attribute_element_count(&attr->format) / attr->format.shape[0];
}

// #timeline

// What a getter is handed: the resolved series, and which member of its population. Periodic values (len > 0) are placed a
// whole number of periods away from where they are: into the window [ref, ref + len[, or, for near, as near ref as they come
// (what is measured against the mouse).
struct TemporalGetterPayload {
    const SeriesTemporalView* view;
    int k;
    double ref = 0.0;
    double len = 0.0;
    bool near = false;
};

static inline double temporal_place(const TemporalGetterPayload* p, double y) {
    if (p->len <= 0.0) return y;
    return p->near ? periodic_near(y, p->ref, p->len) : periodic_wrap(y, p->ref, p->len);
}

static ImPlotPoint temporal_getter_line(int i, void* payload) {
    const TemporalGetterPayload* p = (const TemporalGetterPayload*)payload;
    ImPlotPoint pt = series_temporal_point(*p->view, i, p->k);
    pt.y = temporal_place(p, pt.y);
    return pt;
}

// A band is placed by its centre, and keeps its width
static ImPlotPoint temporal_getter_band_lo(int i, void* payload) {
    const TemporalGetterPayload* p = (const TemporalGetterPayload*)payload;
    double lo, hi;
    series_temporal_band(*p->view, i, &lo, &hi);
    const double c = 0.5 * (lo + hi);
    return ImPlotPoint(p->view->x[i], lo + (temporal_place(p, c) - c));
}

static ImPlotPoint temporal_getter_band_hi(int i, void* payload) {
    const TemporalGetterPayload* p = (const TemporalGetterPayload*)payload;
    double lo, hi;
    series_temporal_band(*p->view, i, &lo, &hi);
    const double c = 0.5 * (lo + hi);
    return ImPlotPoint(p->view->x[i], hi + (temporal_place(p, c) - c));
}

// Arrays of points moved by a whole number of periods (off), for drawing a run of a periodic band below and above itself
struct PeriodicBandPayload {
    const double* x;
    const double* y;
    double off;
};

static ImPlotPoint periodic_band_getter(int i, void* payload) {
    const PeriodicBandPayload* p = (const PeriodicBandPayload*)payload;
    return ImPlotPoint(p->x[i], p->y[i] + p->off);
}

// A periodic series over time, drawn into the window [base, base + len[ of its view (plot_periodic.h): as a polyline cut where it
// crosses the edges of the window, as points in it, or as a band about its centre. copies is how many periods below and above
// the band is drawn as well: on an axis of the period, where the window covers what is shown, 1; otherwise none, and the band
// may stick out of the window.
static void plot_periodic_temporal(const SeriesTemporalView& v, int k, PlotType plot_type, double base, double len, int copies) {
    md_temp_scope_t temp = md_temp_begin();
    const int n = v.num_samples;
    if (v.band) {
        double* x  = md_temp_alloc_array(temp, double, n);
        double* lo = md_temp_alloc_array(temp, double, n);
        double* hi = md_temp_alloc_array(temp, double, n);
        for (int i = 0; i < n; ++i) {
            x[i] = v.x[i];
            series_temporal_band(v, i, &lo[i], &hi[i]);
        }
        double* rx  = md_temp_alloc_array(temp, double, 2 * n);
        double* rlo = md_temp_alloc_array(temp, double, 2 * n);
        double* rhi = md_temp_alloc_array(temp, double, 2 * n);
        int* run_beg   = md_temp_alloc_array(temp, int, n);
        int* run_count = md_temp_alloc_array(temp, int, n);
        const int num_runs = periodic_band_runs(rx, rlo, rhi, run_beg, run_count, x, lo, hi, n, base, len, len * PERIODIC_MAX_STEP);
        for (int r = 0; r < num_runs; ++r) {
            if (run_count[r] < 1) continue;
            for (int c = -copies; c <= copies; ++c) {
                PeriodicBandPayload pl_lo = { rx + run_beg[r], rlo + run_beg[r], c * len };
                PeriodicBandPayload pl_hi = { rx + run_beg[r], rhi + run_beg[r], c * len };
                ImPlot::PlotShadedG(v.plot_id, periodic_band_getter, &pl_lo, periodic_band_getter, &pl_hi, run_count[r]);
            }
        }
    } else if (plot_type == PlotType_Scatter) {
        TemporalGetterPayload payload = { &v, k, base, len, false };
        ImPlot::PlotScatterG(v.plot_id, temporal_getter_line, &payload, n);
    } else {
        const int cap = 4 * n + 1;
        double* px = md_temp_alloc_array(temp, double, cap);
        double* py = md_temp_alloc_array(temp, double, cap);
        const int m = periodic_polyline(px, py, cap, v.x, v.y, n, v.stride, k, v.y_scale, base, len, len * PERIODIC_MAX_STEP);
        ImPlot::PlotLine(v.plot_id, px, py, m);
    }
    md_temp_end(temp);
}

// Cuts the circle of view where the distributions in hist are least populated: the density of every member of every one of them
// which is of the period of view, all counted alike, resampled to a common resolution
static void periodic_cut_least_populated(PeriodicView* view, const SeriesHistogramView* hist, int count) {
    const int nb = 360;
    float density[nb] = {};
    bool any = false;
    for (int j = 0; j < count; ++j) {
        const SeriesHistogramView& h = hist[j];
        if (!h.bins || h.num_bins <= 0 || !series_periodic_valid(h.periodic)) continue;
        if (!series_periodic_equal(h.periodic, view->period)) continue;
        const int members = h.num_periods > 0 ? h.num_periods : h.dim;
        for (int k = 0; k < members; ++k) {
            if (h.num_periods > 0 && !h.periods[k].valid) continue;
            for (int b = 0; b < nb; ++b) {
                const int src = MIN((int)((b + 0.5) * h.num_bins / nb), h.num_bins - 1);
                density[b] += h.bins[(size_t)k * h.num_bins + src];
            }
            any = true;
        }
    }
    if (any) periodic_view_cut_least_populated(view, density, nb);
}

// The same for the series of a timeline subplot: their values over every frame
static void timeline_cut_least_populated(ApplicationState* data, const PlotSubplot& sp, PeriodicView* view) {
    SeriesHistogramView hist[PLOT_MAX_SERIES_PER_SUBPLOT];
    int count = 0;
    for (int j = 0; j < sp.count; ++j) {
        SeriesKey key = sp.series[j].key;
        key.variant = SeriesVariant_Values;
        if (series_resolve_histogram(&hist[count], data, key, 128)) count += 1;
    }
    periodic_cut_least_populated(view, hist, count);
}

static void draw_timeline_properties_menu(ApplicationState* data) {
    PlotSubplot& first = data->timeline.subplots[0];

    ImGui::SeparatorText("Script");
    int num_listed = 0;
    series_for_each_script_property(data, SeriesSource_Script, MD_SCRIPT_PROPERTY_FLAG_TEMPORAL, [&](const SeriesKey& key) {
        plot_series_list_item(data, TIMELINE_SERIES_DND, first, key);
        num_listed += 1;

        // A property with a population also has its summary over it
        char mean_path[SERIES_PATH_CAP + 8];
        const int len = snprintf(mean_path, sizeof(mean_path), "%s/mean", key.path);
        const md_attributes_t* table = series_table(data, key.source);
        if (table && len > 0 && md_attributes_find(table, str_t{mean_path, (size_t)len})) {
            ImGui::Indent();
            const SeriesVariant variants[] = { SeriesVariant_Mean, SeriesVariant_Sigma, SeriesVariant_Extent };
            for (SeriesVariant v : variants) {
                SeriesKey summary = key;
                summary.variant = v;
                plot_series_list_item(data, TIMELINE_SERIES_DND, first, summary);
            }
            ImGui::Unindent();
        }
    });
    if (num_listed == 0) {
        ImGui::TextDisabled("No temporal properties, define and evaluate them in the script editor");
    }

    ImGui::SeparatorText("Loaded along the run");
    plot_system_series_menu(data, TIMELINE_SERIES_DND, first);
}

static void draw_timeline_window(ApplicationState* data) {
    ASSERT(data);
    ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("Timelines", &data->timeline.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
        int& num_subplots = data->timeline.num_subplots;
        num_subplots = CLAMP(num_subplots, 1, PLOT_MAX_SUBPLOTS);

        const float* x_values   = data->timeline.x_values;
        const int num_x_values  = (int)md_array_size(data->timeline.x_values);
        const float min_x_value = num_x_values > 0 ? x_values[0] : 0.0f;
        const float max_x_value = num_x_values > 0 ? x_values[num_x_values - 1] : 1.0f;

        ImPlot::PushStyleVar(ImPlotStyleVar_PlotPadding, ImVec2(ImPlot::GetStyle().PlotPadding.x, 2));
        defer { ImPlot::PopStyleVar(); };

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Properties")) {
                draw_timeline_properties_menu(data);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Periods")) {
                timeline_periods_menu(&data->timeline.periods);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Subplots")) {
                ImGui::SliderInt("Num Subplots", &num_subplots, 1, PLOT_MAX_SUBPLOTS);
                if (ImGui::Button("Add Subplot")) {
                    num_subplots = CLAMP(num_subplots + 1, 1, PLOT_MAX_SUBPLOTS);
                }
                if (ImGui::Button("Clear All")) {
                    plot_clear(data->timeline.subplots, PLOT_MAX_SUBPLOTS);
                }
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        if (ImGui::IsWindowFocused() && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false)) {
            data->animation.mode = data->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
        }

        if (num_x_values > 0) {
            ImPlotInputMap old_map = ImPlot::GetInputMap();

            static bool is_dragging  = false;

            if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                is_dragging = false;
            }
            if (!ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                is_dragging = false;
            }

            // Shift-drag marks a period and ctrl-drag moves the current frame: neither may pan the plots
            // meanwhile, and nor may a period being dragged in the lane
            if (ImGui::GetIO().KeyShift || ImGui::GetIO().KeyCtrl || timeline_periods_dragging()) {
                ImPlot::GetInputMap().Pan = ImGuiMouseButton_Middle;
            }

            char x_label[64] = "Frame";
            char x_unit_str[32] = "";
            md_unit_t x_unit = data->timeline.time_unit;
            if (!md_unit_is_none(x_unit)) {
                md_unit_print(x_unit_str, sizeof(x_unit_str), x_unit);
                snprintf(x_label, sizeof(x_label), "Time (%s)", x_unit_str);
            }

            // The periods, above the plots and on their time axis
            TimelinePeriodAxis range_axis;
            range_axis.frame_times = x_values;
            range_axis.num_frames  = num_x_values;
            range_axis.view_beg    = &data->timeline.view_range.beg_x;
            range_axis.view_end    = &data->timeline.view_range.end_x;
            range_axis.time_unit   = x_unit_str;
            timeline_periods_lane(&data->timeline.periods, range_axis, &data->animation.frame);

            // After the lane, which ctrl-drag in it may have moved
            double time = frame_to_time(data->animation.frame, *data);

            ImPlot::BeginSubplots("##Temporal", num_subplots, 1, ImVec2(-1,-1));

            const ImPlotFlags plot_flags = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;
            // The labels of the axes are drawn within the plots (DrawInsetAxes): the plots take all of the width and height
            const ImPlotAxisFlags axis_flags   = ImPlotAxisFlags_NoSideSwitch | ImPlot::InsetAxisFlags;
            const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_Opposite | ImPlotAxisFlags_AutoFit | ImPlotAxisFlags_RangeFit;

            for (int i = 0; i < num_subplots; ++i) {
                PlotSubplot& sp = data->timeline.subplots[i];

                // Every entry is resolved once per plot, and everything below reads the views
                SeriesTemporalView views[PLOT_MAX_SERIES_PER_SUBPLOT];
                bool resolved[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                for (int j = 0; j < sp.count; ++j) {
                    resolved[j] = series_resolve_temporal(&views[j], data, sp.series[j].key);
                }

                // Periodic values are drawn into the window of their view, where it is cut (plot_periodic.h). When every series
                // here is of one period, the y axis is an axis of it: at most a period, turned by panning it, and linked to
                // every other axis of the period, in the Distributions as well.
                PeriodicView* pviews[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                PeriodicView* axis_view = nullptr;
                {
                    bool all_one = true;
                    int num = 0;
                    for (int j = 0; j < sp.count; ++j) {
                        if (!resolved[j]) continue;
                        num += 1;
                        pviews[j] = periodic_view(views[j].periodic, views[j].unit_str);
                        if (!pviews[j] || (axis_view && pviews[j] != axis_view)) all_one = false;
                        if (!axis_view) axis_view = pviews[j];
                    }
                    if (!all_one || num == 0) axis_view = nullptr;
                }

                int remove_idx = -1;

                if (ImPlot::BeginPlot("", ImVec2(), plot_flags)) {
                    ImPlot::SetupAxisLinks(ImAxis_X1, &data->timeline.view_range.beg_x, &data->timeline.view_range.end_x);
                    ImPlot::SetupAxisLimitsConstraints(ImAxis_X1, min_x_value, max_x_value);

                    const ImPlotAxisFlags axis_flags_x = axis_flags;

                    // When every series in the subplot is shown in the same unit, that unit labels the y axis
                    char y_unit_str[32] = "";
                    for (int j = 0; j < sp.count; ++j) {
                        if (!resolved[j]) continue;
                        if (y_unit_str[0] == '\0') {
                            str_copy_to_char_buf(y_unit_str, sizeof(y_unit_str), str_from_cstr(views[j].unit_str));
                        } else if (strcmp(y_unit_str, views[j].unit_str) != 0) {
                            y_unit_str[0] = '\0';
                            break;
                        }
                    }

                    char y_label[64] = "";
                    if (y_unit_str[0] != '\0') {
                        snprintf(y_label, sizeof(y_label), "(%s)", y_unit_str);
                    }

                    ImPlot::SetupAxis(ImAxis_X1, x_label, axis_flags_x);
                    if (axis_view) {
                        // Not fitted to the values: they fill the period, wherever it is cut
                        periodic_axis_setup(ImAxis_Y1, y_label, axis_view, axis_flags | ImPlotAxisFlags_Opposite);
                    } else {
                        ImPlot::SetupAxis(ImAxis_Y1, y_label, axis_flags_y);
                    }
                    ImPlot::SetupFinish();

                    if (axis_view && periodic_axis_popup(ImAxis_Y1, axis_view)) {
                        timeline_cut_least_populated(data, sp, axis_view);
                    }

                    // The periods shaded below the series, and shift-press marking a new one
                    timeline_periods_plot(&data->timeline.periods, range_axis, data->animation.frame);

                    // The entry, and the member of its population, under the mouse
                    int  hovered_idx      = -1;
                    int  hovered_pop_idx  = -1;
                    char hovered_label[192] = "";

                    bool print_timeline_tooltip = false;

                    if (ImPlot::IsPlotHovered()) {
                        md_bitmask_clear(&data->selection.highlight_mask);
                        script_set_hovered_property(data,  STR_LIT(""));

                        print_timeline_tooltip = true;
                        const ImPlotPoint mouse_pos = ImPlot::GetPlotMousePos();
                        const ImVec2 mouse_coord = ImPlot::PlotToPixels(mouse_pos);
                        const float max_rad = 20; // 20 pixels
                        const float area_dist = max_rad * 0.2f;

                        float min_dist = max_rad;

                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j]) continue;
                            const PlotSeries& s = sp.series[j];
                            const SeriesTemporalView& v = views[j];

                            ImPlotItem* item = ImPlot::GetItem(v.plot_id);
                            if (!item || !item->Show) {
                                continue;
                            }

                            // Each series has an axis of its own, so the samples around the mouse are its own
                            const int n = v.num_samples;
                            const double fi = series_temporal_index_at(v, mouse_pos.x);
                            const int fn = CLAMP((int)(fi + 0.5), 0, n - 1);  // Nearest index
                            const int f[4] = {
                                CLAMP((int)fi - 1, 0, n - 1),
                                CLAMP((int)fi,     0, n - 1),
                                CLAMP((int)fi + 1, 0, n - 1),
                                CLAMP((int)fi + 2, 0, n - 1),
                            };

                            const int dim = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            for (int k = 0; k < dim; ++k) {
                                if (dim > 1 && !s.population_mask.test(k)) {
                                    continue;
                                }
                                // Periodic values: the place of each sample nearest the mouse, of all a period apart
                                TemporalGetterPayload payload = { &v, k, mouse_pos.y, pviews[j] ? periodic_view_length(*pviews[j]) : 0.0, true };
                                double d = DBL_MAX;

                                if (v.band) {
                                    const ImVec2 p_min[4] = {
                                        ImPlot::PlotToPixels(temporal_getter_band_lo(f[0], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_lo(f[1], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_lo(f[2], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_lo(f[3], &payload)),
                                    };
                                    const ImVec2 p_max[4] = {
                                        ImPlot::PlotToPixels(temporal_getter_band_hi(f[0], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_hi(f[1], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_hi(f[2], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_band_hi(f[3], &payload)),
                                    };
                                    // Each segment forms a trapezoid with its left and right sides parallel to the y axis:
                                    // clamp the mouse to it and measure the distance to the clamped point
                                    for (int l = 0; l < 2; ++l) {
                                        const float x_min = MIN(p_min[l].x, p_min[l+1].x);
                                        const float x_max = MAX(p_min[l].x, p_min[l+1].x);
                                        const float t = (x_max > x_min) ? CLAMP((mouse_coord.x - x_min) / (x_max - x_min), 0.0f, 1.0f) : 0.0f;
                                        const float y[2] = {
                                            lerp(p_min[l].y, p_min[l+1].y, t),
                                            lerp(p_max[l].y, p_max[l+1].y, t)
                                        };
                                        const ImVec2 p = {
                                            CLAMP(mouse_coord.x, x_min, x_max),
                                            CLAMP(mouse_coord.y, MIN(y[0], y[1]), MAX(y[0], y[1]))
                                        };
                                        d = MIN(d, sqrt(ImLengthSqr(mouse_coord - p)));
                                    }
                                    d += area_dist;
                                } else if (s.plot_type == PlotType_Scatter) {
                                    ImVec2 p = ImPlot::PlotToPixels(temporal_getter_line(fn, &payload));
                                    d = sqrt(ImLengthSqr(mouse_coord - p));
                                } else {
                                    // Distance to the segments before, at and after the mouse:
                                    // the current one alone is not enough near a vertex
                                    ImVec2 p[4] = {
                                        ImPlot::PlotToPixels(temporal_getter_line(f[0], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_line(f[1], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_line(f[2], &payload)),
                                        ImPlot::PlotToPixels(temporal_getter_line(f[3], &payload)),
                                    };
                                    d = distance_to_linesegment(p[0], p[1], mouse_coord);
                                    d = MIN(distance_to_linesegment(p[1], p[2], mouse_coord), d);
                                    d = MIN(distance_to_linesegment(p[2], p[3], mouse_coord), d);
                                }

                                if (d < min_dist) {
                                    min_dist = (float)d;
                                    char value_buf[64] = "";
                                    series_temporal_print_value(value_buf, sizeof(value_buf), v, fn, k);

                                    hovered_idx = j;
                                    hovered_pop_idx = k;
                                    const str_t member = v.dim > 1 ? series_member_label(data, s.key, k) : str_t{};
                                    if (!str_empty(member)) {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s[%i] (" STR_FMT "): %s %s", v.label, k + 1, STR_ARG(member), value_buf, v.unit_str);
                                    } else if (v.dim > 1) {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s[%i]: %s %s", v.label, k + 1, value_buf, v.unit_str);
                                    } else {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s: %s %s", v.label, value_buf, v.unit_str);
                                    }
                                }
                            }
                        }

                        if (hovered_idx != -1 && views[hovered_idx].script_ident[0] != '\0') {
                            script_set_hovered_property(data, str_from_cstr(views[hovered_idx].script_ident), hovered_pop_idx);
                        }
                    } else if (data->hovered_property_label[0] != '\0') {
                        // Hovered in another view of it (the other plot window, a list): light up its values here
                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j] || sp.series[j].key.variant != SeriesVariant_Values) continue;
                            if (strcmp(data->hovered_property_label, views[j].script_ident) == 0) {
                                hovered_idx = j;
                                hovered_pop_idx = data->hovered_property_pop_idx;
                                break;
                            }
                        }
                    }

                    if (is_dragging) {
                        time = ImPlot::GetPlotMousePos().x;
                    } else if (ImPlot::IsPlotHovered() && ImGui::IsMouseDown(ImGuiMouseButton_Left) && ImGui::IsKeyDown(ImGuiKey_LeftCtrl)) {
                        is_dragging = true;
                    }

                    for (int j = 0; j < sp.count; ++j) {
                        PlotSeries& s = sp.series[j];
                        const SeriesTemporalView& v = views[j];

                        if (ImPlot::IsLegendEntryHovered(v.plot_id)) {
                            if (v.script_ident[0] != '\0') {
                                script_visualize_ref(data, v.vis_ref, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                script_set_hovered_property(data, str_from_cstr(v.script_ident));
                            }
                            if (!resolved[j]) {
                                ImGui::SetTooltip("Nothing to plot at '%s' right now", s.key.path);
                            } else if (v.state.status != SeriesStatus_Current) {
                                char text[160];
                                series_state_describe(text, sizeof(text), v.state);
                                ImGui::SetTooltip("%s", text);
                            }
                            hovered_idx = j;
                            hovered_pop_idx = -1;
                        }

                        // legend context menu
                        if (ImPlot::BeginLegendPopup(v.plot_id)) {
                            if (ImGui::DeleteButton("Remove")) {
                                remove_idx = j;
                                ImGui::CloseCurrentPopup();
                            }

                            if (resolved[j] && !v.band) {
                                const char* plot_type_names[] = {"Line", "Scatter"};
                                int plot_type = s.plot_type == PlotType_Scatter ? 1 : 0;
                                if (ImGui::Combo("Plot Type", &plot_type, plot_type_names, IM_ARRAYSIZE(plot_type_names))) {
                                    s.plot_type = plot_type == 1 ? PlotType_Scatter : PlotType_Line;
                                }
                                if (s.plot_type == PlotType_Scatter) {
                                    if (ImGui::BeginCombo("Marker", ImPlot::GetMarkerName(s.marker))) {
                                        for (int k = 0; k < ImPlotMarker_COUNT; ++k) {
                                            if (ImGui::Selectable(ImPlot::GetMarkerName(k), s.marker == k)) {
                                                s.marker = (ImPlotMarker)k;
                                            }
                                        }
                                        ImGui::EndCombo();
                                    }
                                    ImGui::SliderFloat("Marker Size", &s.marker_size, 0.1f, 10.0f, "%.2f");
                                }
                            }

                            int popup_hovered_pop = -1;
                            plot_series_style_popup(data, s, resolved[j] ? v.dim : 1, v.script_ident, v.vis_ref, &popup_hovered_pop);
                            if (popup_hovered_pop != -1) {
                                hovered_idx = j;
                                hovered_pop_idx = popup_hovered_pop;
                            }
                            ImPlot::EndLegendPopup();
                        }

                        if (!resolved[j]) {
                            // Still in the legend, so it can be seen and removed
                            ImPlot::PlotDummy(v.plot_id);
                        } else {
                            const int population_size = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            auto plot = [&](int k) {
                                ImVec4 color = plot_series_member_color(s, k, population_size);
                                if (v.state.status >= SeriesStatus_Outdated) {
                                    // Not what the script says any more: shown, faintly, until it is evaluated
                                    color.w *= 0.35f;
                                }
                                ImVec4 marker_line_color = {};
                                float  marker_line_weight = 0;
                                float  fill_alpha = 1.0f;
                                float  weight = 1.0f;

                                if (hovered_idx == j) {
                                    if (hovered_pop_idx == -1 || hovered_pop_idx == k) {
                                        color = plot_highlight(color);
                                        fill_alpha = 1.25f;
                                        marker_line_color = {1,1,1,1};
                                        marker_line_weight = 1.0f;
                                    }
                                    if (hovered_pop_idx == k) {
                                        weight = 2.0f;
                                    }
                                }

                                if (v.band) {
                                    ImPlot::SetNextFillStyle(color, fill_alpha);
                                } else if (s.plot_type == PlotType_Scatter) {
                                    ImPlot::SetNextMarkerStyle(s.marker, s.marker_size, color, marker_line_weight, marker_line_color);
                                } else {
                                    ImPlot::SetNextLineStyle(color, weight);
                                }

                                if (pviews[j]) {
                                    // On an axis of the period, into the window of its view, the band below and above it as
                                    // well. Beside values of other kinds, on an axis which reads them as they are, into the
                                    // period they are given in.
                                    const double base = axis_view ? pviews[j]->min : pviews[j]->period.beg;
                                    plot_periodic_temporal(v, k, s.plot_type, base, periodic_view_length(*pviews[j]), axis_view ? 1 : 0);
                                    return;
                                }

                                TemporalGetterPayload payload = { &v, k };
                                if (v.band) {
                                    ImPlot::PlotShadedG(v.plot_id, temporal_getter_band_lo, &payload, temporal_getter_band_hi, &payload, v.num_samples);
                                } else if (s.plot_type == PlotType_Scatter) {
                                    ImPlot::PlotScatterG(v.plot_id, temporal_getter_line, &payload, v.num_samples);
                                } else {
                                    ImPlot::PlotLineG(v.plot_id, temporal_getter_line, &payload, v.num_samples);
                                }
                            };

                            for (int k = 0; k < population_size; ++k) {
                                if (population_size > 1 && !s.population_mask.test(k)) continue;
                                if (hovered_idx == j && hovered_pop_idx == k) continue;
                                plot(k);
                            }
                            // The hovered member last, on top
                            if (hovered_idx == j && hovered_pop_idx != -1 && hovered_pop_idx < population_size) {
                                plot(hovered_pop_idx);
                            }
                        }

                        if (ImPlot::BeginDragDropSourceItem(v.plot_id)) {
                            series_set_drag_payload(TIMELINE_SERIES_DND, s.key, i, v.label, s.color);
                            ImPlot::EndDragDropSource();
                        }
                    }

                    if (ImPlot::DragLineX(0, &time, ImVec4(1,1,0,1), 1.0f, ImPlotDragToolFlags_NoFit)) {
                        time = CLAMP(time, min_x_value, max_x_value);
                    }

                    if (ImPlot::IsPlotHovered() && hovered_idx != -1 && views[hovered_idx].script_ident[0] != '\0') {
                        const SeriesTemporalView& v = views[hovered_idx];
                        const int pop_idx = v.dim > 1 ? hovered_pop_idx : -1;
                        script_visualize_ref(data, v.vis_ref, pop_idx, ~MD_SCRIPT_VISUALIZE_SDF);
                        script_set_hovered_property(data, str_from_cstr(v.script_ident), hovered_pop_idx);
                    }

                    if (ImPlot::BeginDragDropTargetPlot()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_SERIES_DND)) {
                            ASSERT(payload->DataSize == sizeof(SeriesDragPayload));
                            const SeriesDragPayload* dnd = (const SeriesDragPayload*)(payload->Data);
                            if (dnd->src_subplot >= 0 && dnd->src_subplot < PLOT_MAX_SUBPLOTS) {
                                // Out of one plot and into another moves it, keeping its style
                                plot_move_series(data, data->timeline.subplots[dnd->src_subplot], sp, dnd->key);
                            } else {
                                plot_add_series(data, sp, dnd->key);
                            }
                        }
                        ImPlot::EndDragDropTarget();
                    }

                    // A period being marked or dragged has a tooltip of its own
                    if (print_timeline_tooltip && !timeline_periods_dragging()) {
                        ImPlotPoint plot_pos = ImPlot::GetPlotMousePos();
                        ImVec2 screen_pos = ImPlot::PlotToPixels(plot_pos);
                        ImVec2 p0 = {screen_pos.x, ImPlot::GetPlotPos().y};
                        ImVec2 p1 = {screen_pos.x, ImPlot::GetPlotPos().y + ImPlot::GetPlotSize().y};
                        ImPlot::PushPlotClipRect();
                        ImPlot::GetPlotDrawList()->AddLine(p0, p1, IM_COL32(255, 255, 255, 120));
                        ImPlot::PopPlotClipRect();

                        double t = plot_pos.x;
                        if (md_unit_is_none(x_unit)) {
                            int32_t frame_idx = CLAMP((int)(time_to_frame(t, x_values) + 0.5), 0, num_x_values-1);
                            ImGui::SetTooltip("Frame: %i\n%s", frame_idx, hovered_label);
                        } else {
                            ImGui::SetTooltip("Time: %.2f (%s)\n%s", t, x_unit_str, hovered_label);
                        }
                    }

                    // Over everything else. The x axis labelled only under the last plot: they all share it.
                    ImPlot::DrawInsetAxes(i < num_subplots - 1 ? ImPlotInsetFlags_NoX : 0);

                    ImPlot::EndPlot();
                }

                if (remove_idx != -1) {
                    plot_remove_series(sp, remove_idx);
                }
            }

            ImPlot::EndSubplots();

            time = CLAMP(time, (double)min_x_value, (double)max_x_value);
            data->animation.frame = time_to_frame(time, data->timeline.x_values);

            ImPlot::GetInputMap() = old_map;
        }


        // A series dragged out of a plot and dropped outside any of them is taken out of that plot
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            const ImGuiPayload* payload = ImGui::GetDragDropPayload();
            if (payload && payload->IsDataType(TIMELINE_SERIES_DND) && !ImGui::IsDragDropPayloadBeingAccepted()) {
                const SeriesDragPayload* dnd = (const SeriesDragPayload*)(payload->Data);
                if (dnd && dnd->src_subplot >= 0 && dnd->src_subplot < PLOT_MAX_SUBPLOTS) {
                    PlotSubplot& src = data->timeline.subplots[dnd->src_subplot];
                    plot_remove_series(src, plot_find_series(src, dnd->key));
                }
            }
        }
    }
    ImGui::End();
}

// #distribution_window

// What a getter is handed: the histogram, which member of it, and from which bin. Point i is bin u0 + i, which for periodic
// values may lie outside [0, num_bins[: the bin a whole number of periods away, at its place there (series_histogram_point_at).
struct HistogramGetterPayload {
    const SeriesHistogramView* view;
    int k;
    long long u0 = 0;
};

static ImPlotPoint histogram_getter(int i, void* payload) {
    const HistogramGetterPayload* p = (const HistogramGetterPayload*)payload;
    return series_histogram_point_at(*p->view, p->k, p->u0 + i);
}

static ImPlotPoint histogram_getter_zero(int i, void* payload) {
    const HistogramGetterPayload* p = (const HistogramGetterPayload*)payload;
    return ImPlotPoint(series_histogram_point_at(*p->view, p->k, p->u0 + i).x, 0.0);
}

// The bins of a histogram to draw, as the first (u0) and how many. All of them, unless the values are periodic: then those
// over [x_min, x_max] and one beyond each end, so that a curve runs on to the edges of the plot; or, for x_max <= x_min, the
// bins around the period from where it is cut at x_min, and one more, which closes a curve onto its beginning.
static void histogram_draw_range(long long* u0, int* count, const SeriesHistogramView& v, double x_min, double x_max) {
    const double w = series_histogram_bin_width(v);
    if (!series_periodic_valid(v.periodic) || w <= 0.0) {
        *u0 = 0;
        *count = v.num_bins;
        return;
    }
    if (x_max > x_min) {
        const long long a = (long long)floor((x_min - v.x_min) / w - 0.5) - 1;
        const long long b = (long long)ceil ((x_max - v.x_min) / w - 0.5) + 1;
        *u0 = a;
        *count = (int)MIN(b - a + 1, (long long)v.num_bins * 4);
    } else {
        *u0 = (long long)ceil((x_min - v.x_min) / w - 0.5);
        *count = v.num_bins + 1;
    }
}

// A polyline drawn dashed, carrying the dashes on across its corners and from one call to the next through phase
static void dashed_polyline(ImDrawList* dl, const ImVec2* pts, int count, ImU32 col, float* phase) {
    const float dash = 6.0f, gap = 4.0f;
    for (int i = 1; i < count; ++i) {
        const ImVec2 prev = pts[i - 1];
        const ImVec2 d = pts[i] - prev;
        const float len = sqrtf(d.x * d.x + d.y * d.y);
        if (len <= 0.0f) continue;
        const ImVec2 u = d / len;
        float t = 0.0f;
        while (t < len) {
            const float in_cycle = fmodf(*phase + t, dash + gap);
            const float step = in_cycle < dash ? MIN(dash - in_cycle, len - t) : MIN(dash + gap - in_cycle, len - t);
            if (in_cycle < dash) dl->AddLine(prev + u * t, prev + u * (t + step), col, 1.5f);
            t += MAX(step, 0.01f);
        }
        *phase = fmodf(*phase + len, dash + gap);
    }
}

// A histogram member over the plot as a preview, not a curve of the plot's own: drawn as the series is (a line, an area
// or bars), with a dashed outline and, for an area or bars, a faint fill
static void plot_histogram_preview(const SeriesHistogramView& v, int k, ImVec4 color, const PlotSeries& s, long long u0, int count) {
    if (v.num_bins <= 0 || count <= 0) return;
    ImDrawList* dl = ImPlot::GetPlotDrawList();
    ImPlot::PushPlotClipRect();
    const ImU32 line_col = ImGui::ColorConvertFloat4ToU32(color);
    const ImU32 fill_col = ImGui::ColorConvertFloat4ToU32(ImVec4(color.x, color.y, color.z, color.w * 0.25f));
    // Fills which meet edge to edge, without the seams anti-aliasing would leave between them
    const ImDrawListFlags flags = dl->Flags;
    float phase = 0.0f;

    if (s.plot_type == PlotType_Bars) {
        const double half = series_histogram_bin_width(v) * s.bar_width * 0.5;
        const bool joined = s.bar_width >= 0.999f;
        md_temp_scope_t temp = md_temp_begin();
        // Joined, the outline is one step along the tops of the bars; apart, each bar has its own
        ImVec2* outline = md_temp_alloc_array(temp, ImVec2, (size_t)count * 2 + 2);
        int n = 0;
        dl->Flags &= ~ImDrawListFlags_AntiAliasedFill;
        for (int i = 0; i < count; ++i) {
            const ImPlotPoint p = series_histogram_point_at(v, k, u0 + i);
            const ImVec2 a = ImPlot::PlotToPixels(p.x - half, 0.0);
            const ImVec2 b = ImPlot::PlotToPixels(p.x + half, p.y);
            dl->AddRectFilled(ImVec2(MIN(a.x, b.x), MIN(a.y, b.y)), ImVec2(MAX(a.x, b.x), MAX(a.y, b.y)), fill_col);
            if (joined) {
                if (n == 0) outline[n++] = a;
                outline[n++] = ImVec2(a.x, b.y);
                outline[n++] = b;
            } else {
                const ImVec2 bar[4] = { a, ImVec2(a.x, b.y), b, ImVec2(b.x, a.y) };
                dashed_polyline(dl, bar, 4, line_col, &phase);
            }
        }
        dl->Flags = flags;
        if (joined && n > 0) {
            outline[n] = ImVec2(outline[n - 1].x, outline[0].y);
            dashed_polyline(dl, outline, n + 1, line_col, &phase);
        }
        md_temp_end(temp);
    } else {
        md_temp_scope_t temp = md_temp_begin();
        ImVec2* pts = md_temp_alloc_array(temp, ImVec2, (size_t)count);
        for (int i = 0; i < count; ++i) {
            pts[i] = ImPlot::PlotToPixels(series_histogram_point_at(v, k, u0 + i));
        }
        if (s.plot_type == PlotType_Area) {
            // Down to zero, one quad per stretch between two bins
            dl->Flags &= ~ImDrawListFlags_AntiAliasedFill;
            for (int i = 1; i < count; ++i) {
                const ImVec2 base0 = ImPlot::PlotToPixels(series_histogram_point_at(v, k, u0 + i - 1).x, 0.0);
                const ImVec2 base1 = ImPlot::PlotToPixels(series_histogram_point_at(v, k, u0 + i).x, 0.0);
                dl->AddQuadFilled(pts[i - 1], pts[i], base1, base0, fill_col);
            }
            dl->Flags = flags;
        }
        dashed_polyline(dl, pts, count, line_col, &phase);
        md_temp_end(temp);
    }
    ImPlot::PopPlotClipRect();
}

// Over periods: every period's value at the bin, beside its frames, the hovered one first among equals
static void plot_periods_tooltip(const SeriesHistogramView& v, const PlotSeries& s, int bin, int hovered) {
    ImGui::BeginTooltip();
    const ImPlotPoint at = series_histogram_point(v, bin, 0);
    ImGui::Text("%s at %.4g %s", v.label, at.x, v.x_unit_str);
    if (ImGui::BeginTable("##periods", 3, ImGuiTableFlags_SizingFixedFit)) {
        for (int k = 0; k < v.num_periods; ++k) {
            const SeriesPeriodInfo& p = v.periods[k];
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            const ImVec4 color = p.own_color ? p.color : s.color;
            ImGui::ColorButton("##c", ImVec4(color.x, color.y, color.z, 1.0f), ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder,
                ImVec2(ImGui::GetTextLineHeight(), ImGui::GetTextLineHeight()));
            ImGui::SameLine();
            if (k == hovered) ImGui::TextUnformatted(p.name);
            else              ImGui::TextDisabled("%s", p.name);
            ImGui::TableNextColumn();
            if (p.valid) ImGui::Text("%.4g", series_histogram_point(v, bin, k).y);
            else         ImGui::TextDisabled("-");
            ImGui::TableNextColumn();
            if (p.num_frames < p.num_wanted) ImGui::TextDisabled("%u of %u frames", p.num_frames, p.num_wanted);
            else                             ImGui::TextDisabled("%u frames", p.num_frames);
        }
        ImGui::EndTable();
    }
    ImGui::EndTooltip();
}

// The part of a script series' legend popup which says over which periods it is shown, one curve each
static void plot_periods_popup(ApplicationState* data, PlotSeries& s, bool* refit) {
    ImGui::SeparatorText("Periods");
    const TimelinePeriods& rs = data->timeline.periods;
    auto item = [&](uint32_t id, const char* name, const ImVec4* color, const char* frames) {
        bool on = plot_series_has_period(s, id);
        const bool last = on && s.num_periods == 1;
        const bool full = !on && s.num_periods >= PLOT_MAX_PERIODS;
        ImGui::PushID((int)id);
        if (last || full) ImGui::BeginDisabled();
        if (ImGui::Checkbox("##on", &on)) {
            if (on) plot_series_add_period(s, id);
            else    plot_series_remove_period(s, id);
            *refit = true;
        }
        if (last || full) ImGui::EndDisabled();
        ImGui::SameLine();
        if (color) {
            ImGui::ColorButton("##c", *color, ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder,
                ImVec2(ImGui::GetTextLineHeight(), ImGui::GetTextLineHeight()));
            ImGui::SameLine();
        }
        ImGui::TextUnformatted(name);
        if (frames[0]) {
            ImGui::SameLine();
            ImGui::TextDisabled("%s", frames);
        }
        const uint32_t period = id != 0 ? id : s.key.period;
        if (period != 0 && ImGui::IsItemHovered()) timeline_periods_highlight(period);
        ImGui::PopID();
    };

    // The key's own frames first, in the series' colour: every frame, or its period's
    const int num_frames = (int)run_num_frames(data);
    auto frames_of = [&](char* buf, size_t cap, const TimelinePeriod& r) {
        int beg, end;
        if (timeline_period_frames(&beg, &end, r, data->animation.frame, num_frames)) snprintf(buf, cap, "%d to %d", beg, end);
    };
    {
        char frames[48] = "";
        const TimelinePeriod* own = timeline_periods_find(&rs, s.key.period);
        if (own) frames_of(frames, sizeof(frames), *own);
        item(0, s.key.period != 0 ? series_period_name(data, s.key.period) : "All frames", nullptr, frames);
    }
    for (int i = 0; i < rs.count; ++i) {
        const TimelinePeriod& r = rs.items[i];
        if (r.id == s.key.period) continue;
        char frames[48] = "";
        frames_of(frames, sizeof(frames), r);
        const ImVec4 color(r.color.x, r.color.y, r.color.z, 1.0f);
        item(r.id, r.name, &color, frames);
    }
    if (rs.count == 0) {
        ImGui::TextDisabled("Mark periods in the timeline to compare it over them");
    } else {
        ImGui::TextDisabled("Or pull a period out of the timeline onto the plot");
    }
}

static void draw_distribution_properties_menu(ApplicationState* data) {
    PlotSubplot& first = data->distributions.subplots[0];
    const uint32_t kinds = MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION;

    // A temporal property is binned over its frames, a distribution is shown as evaluated. Over the
    // frames of a timeline period too: "r (Production)".
    auto script_section = [&](uint32_t period) -> int {
        int num_listed = 0;
        series_for_each_script_property(data, SeriesSource_Script, kinds, [&](const SeriesKey& over_all) {
            SeriesKey key = over_all;
            key.period = period;
            plot_series_list_item(data, DISTRIBUTION_SERIES_DND, first, key);
            num_listed += 1;
            const bool temporal = md_script_ir_property_flags(data->script.eval_ir, series_script_ident(key)) & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL;
            if (temporal && script_property_population(data, key) > 1) {
                SeriesKey agg = key;
                agg.variant = SeriesVariant_Aggregate;
                ImGui::Indent();
                plot_series_list_item(data, DISTRIBUTION_SERIES_DND, first, agg);
                ImGui::Unindent();
            }
        });
        return num_listed;
    };

    ImGui::SeparatorText("Script");
    if (script_section(0) == 0) {
        ImGui::TextDisabled("No temporal or distribution properties, define and evaluate them in the script editor");
    } else if (data->timeline.periods.count > 0) {
        ImGui::SeparatorText("Script, over a period");
        const TimelinePeriods& rs = data->timeline.periods;
        for (int i = 0; i < rs.count; ++i) {
            const TimelinePeriod& r = rs.items[i];
            char label[TIMELINE_PERIOD_NAME_CAP + 16];
            snprintf(label, sizeof(label), "%s###period%u", r.name, r.id);
            const bool open = ImGui::BeginMenu(label);
            if (ImGui::IsItemHovered()) timeline_periods_highlight(r.id);
            if (open) {
                script_section(r.id);
                ImGui::EndMenu();
            }
        }
    }

    ImGui::SeparatorText("Loaded along the run");
    plot_system_series_menu(data, DISTRIBUTION_SERIES_DND, first);
}

static void draw_distribution_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Distributions", &data->distributions.show_window, ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_MenuBar)) {
        int& num_subplots = data->distributions.num_subplots;
        num_subplots = CLAMP(num_subplots, 1, PLOT_MAX_SUBPLOTS);

        // Refit a subplot's axes to what it holds, after it gained a series or its binning changed
        static bool fit_subplot[PLOT_MAX_SUBPLOTS] = {};

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Properties")) {
                draw_distribution_properties_menu(data);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Subplots")) {
                ImGui::SliderInt("Num Subplots", &num_subplots, 1, PLOT_MAX_SUBPLOTS);
                if (ImGui::Button("Add Subplot")) {
                    num_subplots = CLAMP(num_subplots + 1, 1, PLOT_MAX_SUBPLOTS);
                }
                if (ImGui::Button("Clear All")) {
                    plot_clear(data->distributions.subplots, PLOT_MAX_SUBPLOTS);
                }
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        // The labels of the axes are drawn within the plots (DrawInsetAxes): the plots take all of the width and height
        const ImPlotAxisFlags axis_flags   = ImPlotAxisFlags_NoSideSwitch | ImPlotAxisFlags_NoHighlight | ImPlot::InsetAxisFlags;
        const ImPlotAxisFlags axis_flags_x = axis_flags;
        const ImPlotAxisFlags axis_flags_y = axis_flags | ImPlotAxisFlags_AutoFit;
        const ImPlotFlags     plot_flags   = ImPlotFlags_NoBoxSelect | ImPlotFlags_NoFrame;

        if (ImPlot::BeginSubplots("##distribution_plots", num_subplots, 1, ImVec2(-1,-1))) {
            for (int i = 0; i < num_subplots; ++i) {
                PlotSubplot& sp = data->distributions.subplots[i];

                SeriesHistogramView views[PLOT_MAX_SERIES_PER_SUBPLOT];
                bool resolved[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                for (int j = 0; j < sp.count; ++j) {
                    PlotSeries& s = sp.series[j];
                    // A period removed takes its curve along
                    for (int p = s.num_periods - 1; p >= 0; --p) {
                        if (s.periods[p] != 0 && !timeline_periods_find(&data->timeline.periods, s.periods[p])) {
                            if (s.num_periods > 1) plot_series_remove_period(s, s.periods[p]);
                            else s.periods[0] = 0;
                        }
                    }
                    resolved[j] = series_resolve_histogram_periods(&views[j], data, s.key, s.num_bins, s.periods, s.num_periods, s.smoothing);
                }

                // Periodic values are drawn around their period from where their view cuts it (plot_periodic.h). When every
                // series here is of one period, the x axis is an axis of it: at most a period, turned by panning it, and linked
                // to every other axis of the period, in the Timelines as well.
                PeriodicView* pviews[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                PeriodicView* axis_view = nullptr;
                {
                    bool all_one = true;
                    int num = 0;
                    for (int j = 0; j < sp.count; ++j) {
                        if (!resolved[j] || views[j].num_bins <= 0) continue;
                        num += 1;
                        pviews[j] = periodic_view(views[j].periodic, views[j].x_unit_str);
                        if (!pviews[j] || (axis_view && pviews[j] != axis_view)) all_one = false;
                        if (!axis_view) axis_view = pviews[j];
                    }
                    if (!all_one || num == 0) axis_view = nullptr;
                }

                // More to show than before (a workspace loaded before anything was evaluated, an evaluation which
                // came through): fitted to it
                static int num_resolved[PLOT_MAX_SUBPLOTS] = {};
                int resolved_now = 0;
                for (int j = 0; j < sp.count; ++j) resolved_now += resolved[j] ? 1 : 0;
                if (resolved_now > num_resolved[i]) fit_subplot[i] = true;
                num_resolved[i] = resolved_now;

                // A period under the pointer in the timeline: each series over it too, as a preview
                const uint32_t preview_id = timeline_periods_hovered();
                const TimelinePeriod* preview_range = timeline_periods_find(&data->timeline.periods, preview_id);

                int remove_idx = -1;

                if (fit_subplot[i]) {
                    // An axis of a period shows it as it is cut: only the densities are fitted
                    if (axis_view) ImPlot::SetNextAxisToFit(ImAxis_Y1);
                    else           ImPlot::SetNextAxesToFit();
                    fit_subplot[i] = false;
                }

                if (ImPlot::BeginPlot("", ImVec2(-1,0), plot_flags)) {
                    // Labelled with what the series are shown in; empty when there is no unit, or several
                    char x_label[64] = "";
                    char y_label[64] = "";
                    bool first_label = true;
                    for (int j = 0; j < sp.count; ++j) {
                        if (!resolved[j]) continue;
                        if (first_label) {
                            str_copy_to_char_buf(x_label, sizeof(x_label), str_from_cstr(views[j].x_unit_str));
                            str_copy_to_char_buf(y_label, sizeof(y_label), str_from_cstr(views[j].y_unit_str));
                            first_label = false;
                        } else {
                            if (strcmp(x_label, views[j].x_unit_str) != 0) x_label[0] = '\0';
                            if (strcmp(y_label, views[j].y_unit_str) != 0) y_label[0] = '\0';
                        }
                    }

                    if (axis_view) {
                        periodic_axis_setup(ImAxis_X1, x_label, axis_view, axis_flags_x);
                    } else {
                        ImPlot::SetupAxis(ImAxis_X1, x_label, axis_flags_x);
                    }
                    ImPlot::SetupAxis(ImAxis_Y1, y_label, axis_flags_y);
                    ImPlot::SetupFinish();

                    if (axis_view && periodic_axis_popup(ImAxis_X1, axis_view)) {
                        periodic_cut_least_populated(axis_view, views, sp.count);
                    }

                    // Which bins each series draws: of periodic values, around the period from where it is cut, and on an
                    // axis of the period, over what it shows
                    const ImPlotRect limits = ImPlot::GetPlotLimits();
                    long long draw_u0[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                    int draw_count[PLOT_MAX_SERIES_PER_SUBPLOT] = {};
                    auto draw_range = [&](long long* u0, int* count, const SeriesHistogramView& hv, const PeriodicView* pv) {
                        if (axis_view)  histogram_draw_range(u0, count, hv, limits.X.Min, limits.X.Max);
                        else if (pv)    histogram_draw_range(u0, count, hv, pv->period.beg, pv->period.beg);
                        else            histogram_draw_range(u0, count, hv, 0.0, 0.0);
                    };
                    for (int j = 0; j < sp.count; ++j) {
                        if (resolved[j]) draw_range(&draw_u0[j], &draw_count[j], views[j], pviews[j]);
                    }

                    int  hovered_idx       = -1;
                    int  hovered_pop_idx   = -1;
                    int  hovered_bin       = -1;
                    char hovered_label[192] = "";

                    if (ImPlot::IsPlotHovered()) {
                        script_set_hovered_property(data, STR_LIT(""));

                        const ImPlotPoint mouse_pos = ImPlot::GetPlotMousePos();
                        const ImVec2 mouse_coord = ImPlot::PlotToPixels(mouse_pos);

                        const double max_rad = 20; // 20 pixels
                        const double area_dist = max_rad * 0.2;
                        double min_dist = max_rad;

                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j]) continue;
                            const PlotSeries& s = sp.series[j];
                            const SeriesHistogramView& v = views[j];
                            if (v.num_bins <= 0 || v.x_max <= v.x_min) continue;

                            ImPlotItem* item = ImPlot::GetItem(v.plot_id);
                            if (!item || !item->Show) continue;

                            // Among the points it draws: from bin draw_u0, which for periodic values may lie a period or more off
                            const int nb = draw_count[j];
                            if (nb <= 0) continue;
                            const double x = (mouse_pos.x - v.x_min) / series_histogram_bin_width(v) - 0.5 - (double)draw_u0[j];
                            const int xn = CLAMP((int)floor(x + 0.5), 0, nb - 1);
                            const int xf = (int)floor(x);
                            const int xi[4] = {
                                CLAMP(xf - 1, 0, nb - 1),
                                CLAMP(xf,     0, nb - 1),
                                CLAMP(xf + 1, 0, nb - 1),
                                CLAMP(xf + 2, 0, nb - 1),
                            };

                            const int dim = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            for (int k = 0; k < dim; ++k) {
                                if (v.num_periods > 0 ? !v.periods[k].valid : (dim > 1 && !s.population_mask.test(k))) continue;
                                HistogramGetterPayload payload = { &v, k, draw_u0[j] };
                                // Later entries lie on top: prefer them a little when two are equally near
                                const double layer = j + k / (double)(PLOT_MAX_POPULATION - 1);
                                const double layer_dist = ((sp.count - layer) / sp.count) * (max_rad * 0.1);
                                double d = DBL_MAX;

                                if (s.plot_type == PlotType_Bars) {
                                    const double half_width = 0.5 * series_histogram_bin_width(v) * s.bar_width;
                                    for (int l = 0; l < 3; ++l) {
                                        const ImPlotPoint top = histogram_getter(xi[l], &payload);
                                        const ImVec2 p0 = ImPlot::PlotToPixels(ImPlotPoint(top.x - half_width, 0.0));
                                        const ImVec2 p1 = ImPlot::PlotToPixels(ImPlotPoint(top.x + half_width, top.y));
                                        const ImVec2 ll = ImMin(p0, p1);
                                        const ImVec2 ur = ImMax(p0, p1);
                                        d = MIN(d, sqrt(ImLengthSqr(mouse_coord - ImClamp(mouse_coord, ll, ur))));
                                    }
                                    d += area_dist + layer_dist;
                                } else if (s.plot_type == PlotType_Area) {
                                    for (int l = 0; l < 2; ++l) {
                                        const ImVec2 a = ImPlot::PlotToPixels(histogram_getter(xi[l], &payload));
                                        const ImVec2 b = ImPlot::PlotToPixels(histogram_getter(xi[l + 1], &payload));
                                        const ImVec2 z = ImPlot::PlotToPixels(histogram_getter_zero(xi[l], &payload));
                                        const float x_min = MIN(a.x, b.x);
                                        const float x_max = MAX(a.x, b.x);
                                        const float t = (x_max > x_min) ? CLAMP((mouse_coord.x - x_min) / (x_max - x_min), 0.0f, 1.0f) : 0.0f;
                                        const float y_top = lerp(a.y, b.y, t);
                                        const ImVec2 p = {
                                            CLAMP(mouse_coord.x, x_min, x_max),
                                            CLAMP(mouse_coord.y, MIN(y_top, z.y), MAX(y_top, z.y))
                                        };
                                        d = MIN(d, sqrt(ImLengthSqr(mouse_coord - p)));
                                    }
                                    d += area_dist + layer_dist;
                                } else {
                                    ImVec2 p[4] = {
                                        ImPlot::PlotToPixels(histogram_getter(xi[0], &payload)),
                                        ImPlot::PlotToPixels(histogram_getter(xi[1], &payload)),
                                        ImPlot::PlotToPixels(histogram_getter(xi[2], &payload)),
                                        ImPlot::PlotToPixels(histogram_getter(xi[3], &payload)),
                                    };
                                    d = distance_to_linesegment(p[0], p[1], mouse_coord);
                                    d = MIN(distance_to_linesegment(p[1], p[2], mouse_coord), d);
                                    d = MIN(distance_to_linesegment(p[2], p[3], mouse_coord), d);
                                    d += layer_dist;
                                }

                                if (d < min_dist) {
                                    min_dist = d;
                                    ImPlotPoint pt = histogram_getter(xn, &payload);
                                    if (series_periodic_valid(v.periodic)) {
                                        // As the values are given, wherever the period is cut
                                        pt.x = periodic_wrap(pt.x, v.periodic.beg, series_periodic_length(v.periodic));
                                    }
                                    hovered_idx = j;
                                    hovered_pop_idx = (dim > 1 || v.num_periods > 0) ? k : -1;
                                    hovered_bin = series_histogram_bin(v, draw_u0[j] + xn);
                                    if (v.num_periods > 0) {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s, %s: %.4g at %.4g %s", v.label, v.periods[k].name, pt.y, pt.x, v.x_unit_str);
                                    } else if (dim > 1) {
                                        const str_t member = series_member_label(data, s.key, k);
                                        if (!str_empty(member)) {
                                            snprintf(hovered_label, sizeof(hovered_label), "%s[%i] (" STR_FMT "): %.4g at %.4g %s", v.label, k + 1, STR_ARG(member), pt.y, pt.x, v.x_unit_str);
                                        } else {
                                            snprintf(hovered_label, sizeof(hovered_label), "%s[%i]: %.4g at %.4g %s", v.label, k + 1, pt.y, pt.x, v.x_unit_str);
                                        }
                                    } else {
                                        snprintf(hovered_label, sizeof(hovered_label), "%s: %.4g at %.4g %s", v.label, pt.y, pt.x, v.x_unit_str);
                                    }
                                }
                            }
                        }

                        if (hovered_idx != -1) {
                            const SeriesHistogramView& v = views[hovered_idx];
                            // A period is not a member of a population: the property as a whole is what is shown
                            const int member = v.num_periods > 0 ? -1 : hovered_pop_idx;
                            if (v.script_ident[0] != '\0') {
                                script_set_hovered_property(data, str_from_cstr(v.script_ident), member);
                                script_visualize_ref(data, v.vis_ref, member, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                            }
                            if (v.num_periods > 0 && hovered_bin >= 0) {
                                if (hovered_pop_idx >= 0 && v.periods[hovered_pop_idx].id != 0) {
                                    timeline_periods_highlight(v.periods[hovered_pop_idx].id);
                                }
                                plot_periods_tooltip(v, sp.series[hovered_idx], hovered_bin, hovered_pop_idx);
                            } else if (hovered_label[0] != '\0') {
                                ImGui::SetTooltip("%s", hovered_label);
                            }
                        }
                    } else if (data->hovered_property_label[0] != '\0') {
                        for (int j = 0; j < sp.count; ++j) {
                            if (!resolved[j]) continue;
                            if (strcmp(data->hovered_property_label, views[j].script_ident) == 0) {
                                hovered_idx = j;
                                // A member of a population elsewhere is not a period here
                                hovered_pop_idx = views[j].num_periods > 0 ? -1 : data->hovered_property_pop_idx;
                                break;
                            }
                        }
                    }

                    for (int j = 0; j < sp.count; ++j) {
                        PlotSeries& s = sp.series[j];
                        const SeriesHistogramView& v = views[j];

                        if (ImPlot::IsLegendEntryHovered(v.plot_id)) {
                            if (v.script_ident[0] != '\0') {
                                script_visualize_ref(data, v.vis_ref, -1, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_GEOMETRY);
                                script_set_hovered_property(data, str_from_cstr(v.script_ident));
                            }
                            if (!resolved[j]) {
                                ImGui::SetTooltip("Nothing to show for '%s' right now", s.key.path);
                            } else if (v.state.status != SeriesStatus_Current) {
                                char text[160];
                                series_state_describe(text, sizeof(text), v.state);
                                ImGui::SetTooltip("%s", text);
                            }
                            hovered_idx = j;
                            hovered_pop_idx = -1;
                        }

                        if (ImPlot::BeginLegendPopup(v.plot_id)) {
                            if (ImGui::DeleteButton("Remove")) {
                                remove_idx = j;
                                ImGui::CloseCurrentPopup();
                            }

                            const char* plot_type_names[] = {"Line", "Area", "Bars"};
                            const PlotType plot_types[]   = {PlotType_Line, PlotType_Area, PlotType_Bars};
                            int plot_type = 0;
                            for (int k = 0; k < IM_ARRAYSIZE(plot_types); ++k) {
                                if (s.plot_type == plot_types[k]) plot_type = k;
                            }
                            if (ImGui::Combo("Plot Type", &plot_type, plot_type_names, IM_ARRAYSIZE(plot_type_names))) {
                                s.plot_type = plot_types[plot_type];
                            }
                            if (s.plot_type == PlotType_Bars) {
                                ImGui::SliderFloat("Bar Width", &s.bar_width, 0.01f, 1.0f, "%.3f");
                            }

                            // Powers of two: a script distribution only divides into those
                            int bins_exp = (int)log2((double)MAX(s.num_bins, 2));
                            char bins_fmt[32];
                            snprintf(bins_fmt, sizeof(bins_fmt), "%d", 1 << bins_exp);
                            if (ImGui::SliderInt("Num Bins", &bins_exp, 5, 10, bins_fmt)) {
                                s.num_bins = 1 << bins_exp;
                                fit_subplot[i] = true;
                            }
                            if (resolved[j] && v.num_bins != s.num_bins) {
                                ImGui::TextDisabled("Evaluated at %d bins", v.num_bins);
                            }

                            // Values binned, rather than a distribution as the script evaluated it: they may be smoothed
                            const bool binned = s.key.source != SeriesSource_Script ||
                                (md_script_ir_property_flags(data->script.eval_ir, series_script_ident(s.key)) & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL);
                            if (binned) {
                                bool smooth = s.smoothing > 0.0f;
                                if (ImGui::Checkbox("Smooth", &smooth)) {
                                    s.smoothing = smooth ? 1.0f : 0.0f;
                                }
                                ImGui::SetItemTooltip("A kernel density estimate rather than the histogram: a Gaussian about every value, as\n"
                                                      "wide as Silverman's rule of thumb makes it%s",
                                                      series_periodic_valid(v.periodic) ? ", wrapped around the period" : "");
                                if (smooth) {
                                    ImGui::SameLine();
                                    ImGui::SetNextItemWidth(ImGui::GetFontSize() * 8);
                                    ImGui::SliderFloat("Bandwidth", &s.smoothing, 0.25f, 4.0f, "%.2fx", ImGuiSliderFlags_Logarithmic);
                                    ImGui::SetItemTooltip("Of Silverman's rule of thumb, which oversmooths what has several peaks:\nless than 1 to see them apart");
                                }
                            }

                            int popup_hovered_pop = -1;
                            plot_series_style_popup(data, s, resolved[j] && v.num_periods == 0 ? v.dim : 1, v.script_ident, v.vis_ref, &popup_hovered_pop);
                            if (popup_hovered_pop != -1) {
                                hovered_idx = j;
                                hovered_pop_idx = popup_hovered_pop;
                            }
                            if (v.script_ident[0] != '\0') {
                                plot_periods_popup(data, s, &fit_subplot[i]);
                            }
                            ImPlot::EndLegendPopup();
                        }

                        if (!resolved[j] || v.num_bins <= 0) {
                            ImPlot::PlotDummy(v.plot_id);
                        } else {
                            const int population_size = CLAMP(v.dim, 1, PLOT_MAX_POPULATION);
                            auto plot = [&](int k) {
                                // A period in its own colour, so that it reads with the timeline
                                ImVec4 color = v.num_periods > 0 && v.periods[k].own_color ? v.periods[k].color : plot_series_member_color(s, k, population_size);
                                if (v.state.status >= SeriesStatus_Outdated) {
                                    // Not what the script says any more: shown, faintly, until it is evaluated
                                    color.w *= 0.35f;
                                }
                                float  fill_alpha = 1.0f;
                                float  weight = 1.0f;
                                if (hovered_idx == j) {
                                    if (hovered_pop_idx == -1 || hovered_pop_idx == k) {
                                        color = plot_highlight(color);
                                        fill_alpha = 1.25f;
                                    }
                                    if (hovered_pop_idx == k) {
                                        weight = 2.0f;
                                    }
                                }

                                HistogramGetterPayload payload = { &v, k, draw_u0[j] };
                                // Around a period from where it is cut, bars without the one which closes a curve: it would
                                // stand outside the period, beside the others
                                const bool rotated = pviews[j] && !axis_view;
                                const int count = draw_count[j] - (rotated && s.plot_type == PlotType_Bars ? 1 : 0);
                                switch (s.plot_type) {
                                case PlotType_Area:
                                    ImPlot::SetNextFillStyle(color, fill_alpha);
                                    ImPlot::PlotShadedG(v.plot_id, histogram_getter_zero, &payload, histogram_getter, &payload, count);
                                    break;
                                case PlotType_Bars:
                                    ImPlot::SetNextFillStyle(color, fill_alpha);
                                    ImPlot::PlotBarsG(v.plot_id, histogram_getter, &payload, count, series_histogram_bin_width(v) * s.bar_width);
                                    break;
                                default:
                                    ImPlot::SetNextLineStyle(color, weight);
                                    ImPlot::PlotLineG(v.plot_id, histogram_getter, &payload, count);
                                    break;
                                }
                            };

                            for (int k = 0; k < population_size; ++k) {
                                if (v.num_periods > 0 ? !v.periods[k].valid : (population_size > 1 && !s.population_mask.test(k))) continue;
                                if (hovered_idx == j && hovered_pop_idx == k) continue;
                                plot(k);
                            }
                            if (hovered_idx == j && hovered_pop_idx != -1 && hovered_pop_idx < population_size) {
                                plot(hovered_pop_idx);
                            }
                        }

                        // The period under the pointer in the timeline, drawn as the series is but dashed, unless the series is
                        // shown over it already
                        if (preview_range && v.script_ident[0] != '\0' && !plot_series_has_period(s, preview_id)) {
                            SeriesHistogramView pv;
                            if (series_resolve_histogram_periods(&pv, data, s.key, s.num_bins, &preview_id, 1, s.smoothing) && pv.num_bins > 0 && pv.periods[0].valid) {
                                const vec4_t rc = preview_range->color;
                                long long pu0 = 0;
                                int pcount = 0;
                                draw_range(&pu0, &pcount, pv, pviews[j]);
                                if (pviews[j] && !axis_view && s.plot_type == PlotType_Bars) pcount -= 1;
                                plot_histogram_preview(pv, 0, ImVec4(rc.x, rc.y, rc.z, 1.0f), s, pu0, pcount);
                            }
                        }

                        if (ImPlot::BeginDragDropSourceItem(v.plot_id)) {
                            series_set_drag_payload(DISTRIBUTION_SERIES_DND, s.key, i, v.label, s.color);
                            ImPlot::EndDragDropSource();
                        }
                    }

                    // A period pulled out of the timeline, dropped on a legend entry: that series over it too
                    if (ImPlot::BeginDragDropTargetLegend()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                            ASSERT(payload->DataSize == sizeof(uint32_t));
                            const uint32_t id = *(const uint32_t*)payload->Data;
                            for (int j = 0; j < sp.count; ++j) {
                                if (ImPlot::IsLegendEntryHovered(views[j].plot_id) && views[j].script_ident[0] != '\0') {
                                    plot_series_add_period(sp.series[j], id);
                                }
                            }
                        }
                        ImPlot::EndDragDropTarget();
                    }

                    if (ImPlot::BeginDragDropTargetPlot()) {
                        if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(DISTRIBUTION_SERIES_DND)) {
                            ASSERT(payload->DataSize == sizeof(SeriesDragPayload));
                            const SeriesDragPayload* dnd = (const SeriesDragPayload*)(payload->Data);
                            if (dnd->src_subplot >= 0 && dnd->src_subplot < PLOT_MAX_SUBPLOTS) {
                                plot_move_series(data, data->distributions.subplots[dnd->src_subplot], sp, dnd->key);
                            } else {
                                plot_add_series(data, sp, dnd->key);
                            }
                            fit_subplot[i] = true;
                        } else if (const ImGuiPayload* range_payload = ImGui::AcceptDragDropPayload(TIMELINE_PERIOD_DND)) {
                            // A period pulled out of the timeline: every script series here over it too
                            ASSERT(range_payload->DataSize == sizeof(uint32_t));
                            const uint32_t id = *(const uint32_t*)range_payload->Data;
                            for (int j = 0; j < sp.count; ++j) {
                                if (views[j].script_ident[0] != '\0') {
                                    plot_series_add_period(sp.series[j], id);
                                }
                            }
                        } else if (const ImGuiPayload* tl_payload = ImGui::AcceptDragDropPayload(TIMELINE_SERIES_DND)) {
                            // Something plotted over time, dropped here: the distribution of its
                            // values, whichever view of them (a mean, a spread) was dragged
                            ASSERT(tl_payload->DataSize == sizeof(SeriesDragPayload));
                            SeriesKey key = ((const SeriesDragPayload*)tl_payload->Data)->key;
                            key.variant = SeriesVariant_Values;
                            plot_add_series(data, sp, key);
                            fit_subplot[i] = true;
                        }
                        ImPlot::EndDragDropTarget();
                    }

                    if (ImPlot::IsPlotHovered()) {
                        ImPlotPoint plot_pos = ImPlot::GetPlotMousePos();
                        ImVec2 screen_pos = ImPlot::PlotToPixels(plot_pos);
                        ImVec2 p0 = {screen_pos.x, ImPlot::GetPlotPos().y};
                        ImVec2 p1 = {screen_pos.x, ImPlot::GetPlotPos().y + ImPlot::GetPlotSize().y};
                        ImPlot::PushPlotClipRect();
                        ImPlot::GetPlotDrawList()->AddLine(p0, p1, IM_COL32(255, 255, 255, 120));
                        ImPlot::PopPlotClipRect();
                    }

                    // Over everything else
                    ImPlot::DrawInsetAxes();

                    ImPlot::EndPlot();
                }

                if (remove_idx != -1) {
                    plot_remove_series(sp, remove_idx);
                }
            }
            ImPlot::EndSubplots();
        }

        // A series dragged out of a plot and dropped outside any of them is taken out of that plot
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
            const ImGuiPayload* payload = ImGui::GetDragDropPayload();
            if (payload && payload->IsDataType(DISTRIBUTION_SERIES_DND) && !ImGui::IsDragDropPayloadBeingAccepted()) {
                const SeriesDragPayload* dnd = (const SeriesDragPayload*)(payload->Data);
                if (dnd && dnd->src_subplot >= 0 && dnd->src_subplot < PLOT_MAX_SUBPLOTS) {
                    PlotSubplot& src = data->distributions.subplots[dnd->src_subplot];
                    plot_remove_series(src, plot_find_series(src, dnd->key));
                }
            }
        }
    }
    ImGui::End();
}

// What the profiler's times were measured at, a line each, for the head of an export (profiler::write_csv)
static void profiler_export_info(md_str_writer_t* sw, const ApplicationState* state) {
    md_str_writer_fmt(sw, "framebuffer: %u x %u\n", (unsigned)state->app.framebuffer.width, (unsigned)state->app.framebuffer.height);
    md_str_writer_fmt(sw, "viewport: %u x %u shown, render scale %.3f, %u x %u rendered\n",
                      (unsigned)state->gbuffer.out_width, (unsigned)state->gbuffer.out_height, state->settings.render_scale,
                      (unsigned)state->gbuffer.width, (unsigned)state->gbuffer.height);
    if (iso_mesh_target.width) {
        md_str_writer_fmt(sw, "isosurface meshes: %u x %u, %llu triangles, ", iso_mesh_target.width,
                          iso_mesh_target.height, (unsigned long long)iso_mesh_target.triangles);
        if (state->profiler_keep_isomeshes) {
            md_str_writer_fmt(sw, "drawn every frame, %u frames since the view or the scene changed\n", iso_mesh_target.frame);
        } else {
            md_str_writer_fmt(sw, "%u of %u frames in the mean\n", iso_mesh_target.frame, isomesh::MAX_FRAMES);
        }
    }
}

// The profiler's times (profiler::write_csv) as a CSV file
static void profiler_export_csv(const ApplicationState* state, str_t path) {
    md_file_t file = {};
    if (!md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        VIAMD_LOG_ERROR("Failed to open '" STR_FMT "' to write the profiler's times", STR_ARG(path));
        return;
    }
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    md_str_writer_t info = md_str_writer_create(md_temp_allocator(temp));
    profiler_export_info(&info, state);

    md_file_writer_t fw = md_file_writer_create(file, NULL, 0);
    profiler::write_csv(&fw.w, md_str_writer_to_cstr(&info));
    const bool ok = md_file_writer_close(&fw);
    md_file_close(&file);
    if (ok) {
        VIAMD_LOG_INFO("Exported the profiler's times to '" STR_FMT "'", STR_ARG(path));
    } else {
        VIAMD_LOG_ERROR("Failed to write the profiler's times to '" STR_FMT "'", STR_ARG(path));
    }
}

// The frame's sections and their times (gfx/profiler.h)
static void draw_profiler_window(ApplicationState* data) {
    ImGui::SetNextWindowSize(ImVec2(ImGui::GetFontSize() * 36.0f, ImGui::GetFontSize() * 32.0f), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Profiler", &data->show_profiler_window)) {
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("Per frame, averaged over 30 frames. GPU: the work a section issued; CPU: from its start to its "
                            "end. The wait for vsync is outside the frame.");
        const int made = (int)popcnt64(data->render.history);
        ImGui::TextDisabled("The scene was made in %d of the last 64 frames%s.", made,
                            data->settings.lazy_rendering ? "; shown again as it was in the others" : " (lazy rendering is off)");
        ImGui::PopTextWrapPos();

        // The table as it is shown, with what it was measured at
        ImGui::BeginDisabled(profiler::frame_cpu_ms() <= 0.0f);
        if (ImGui::Button("Export...")) {
            char path_buf[2048] = {};
            if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("csv"))) {
                profiler_export_csv(data, str_from_cstr(path_buf));
            }
        }
        ImGui::SetItemTooltip("The times as CSV, with the renderer and the sizes they were measured at");
        ImGui::SameLine();
        if (ImGui::Button("Copy")) {
            md_temp_scope_t temp = md_temp_begin();
            md_str_writer_t info = md_str_writer_create(md_temp_allocator(temp));
            profiler_export_info(&info, data);
            md_str_writer_t csv = md_str_writer_create(md_temp_allocator(temp));
            profiler::write_csv(&csv.w, md_str_writer_to_cstr(&info));
            ImGui::SetClipboardText(md_str_writer_to_cstr(&csv));
            md_temp_end(temp);
        }
        ImGui::SetItemTooltip("The same, to the clipboard");
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::Checkbox("Isosurface meshes every frame", &data->profiler_keep_isomeshes);
        ImGui::SetItemTooltip("Their passes run every frame, also once the image has settled: times to compare with the "
                              "view kept still");

        profiler::draw_table();
    }
    ImGui::End();
}

static void draw_debug_window(ApplicationState* data) {
    ASSERT(data);

    ImGui::SetNextWindowSize(ImVec2(400, 400), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Debug", &data->show_debug_window)) {       
        task_system::ID tasks[256]; 
        size_t num_tasks = task_system::pool_running_tasks(tasks, ARRAY_SIZE(tasks));
        if (num_tasks > 0) {
            ImGui::Text("Running Pool Tasks:");
            for (size_t i = 0; i < num_tasks; ++i) {
                str_t lbl = task_system::task_label(tasks[i]);
                ImGui::Text("[%i]: %.*s", (int)i, (int)lbl.len, lbl.ptr);
            }
        }

        ImGuiID active = ImGui::GetActiveID();
        ImGuiID hover  = ImGui::GetHoveredID();
        ImGui::Text("Active ID: %u, Hover ID: %u", active, hover);

        ImGui::Text("Mouse Pos: (%.3f, %.3f)", ImGui::GetMousePos().x, ImGui::GetMousePos().y);
        ImGui::Text("Camera Position: (%g, %g, %g)", data->view.camera.position.x, data->view.camera.position.y, data->view.camera.position.z);
        ImGui::Text("Camera Orientation: (%g, %g, %g, %g)", data->view.camera.orientation.x, data->view.camera.orientation.y, data->view.camera.orientation.z, data->view.camera.orientation.w);

        mat4_t P = data->view.param.matrix.curr.proj;
        ImGui::Text("proj_matrix:");
        ImGui::Text("[%g %g %g %g]", P.col[0].x, P.col[0].y, P.col[0].z, P.col[0].w);
        ImGui::Text("[%g %g %g %g]", P.col[1].x, P.col[1].y, P.col[1].z, P.col[1].w);
        ImGui::Text("[%g %g %g %g]", P.col[2].x, P.col[2].y, P.col[2].z, P.col[2].w);
        ImGui::Text("[%g %g %g %g]", P.col[3].x, P.col[3].y, P.col[3].z, P.col[3].w);
    }
    ImGui::End();
}

// Opens the script reference at topic (a procedure name, alias or heading anchor). Anything else is searched for,
// and an empty topic puts the focus in the search box (which always takes the focus). Without take_focus the keyboard
// focus stays where it is, so lookups from the script editor keep its caret.
static void open_script_reference(ApplicationState* state, str_t topic, bool take_focus) {
    ASSERT(state);
    if (str_empty(topic)) {
        script_reference::focus_search();
        take_focus = true;
    } else if (!script_reference::show(topic)) {
        script_reference::search(topic);
    }
    state->show_script_reference_window = true;
    script_reference::reveal(take_focus);
}

static void draw_script_reference_window(ApplicationState* state) {
    ASSERT(state);
    const script_reference::Action action = script_reference::draw_window(&state->show_script_reference_window);
    if (!str_empty(action.insert_code)) {
        // Examples go in as whole lines at the cursor of the script editor
        script_editor::insert_lines_at_cursor(state->editor, action.insert_code);
        state->show_script_window = true;
    }
}

// What a payload dragged over the script editor is written as in the script, false if it does not go there. An attribute
// of the system dragged as a series (from the Attributes or System window, a plot's legend, ...) goes in as the
// expression that reads it, attr("edr/potential"), text for the script as it is.
static bool script_drop_text(char* buf, size_t cap, const ApplicationState* state, const ImGuiPayload* payload) {
    ASSERT(buf && cap > 0);
    ASSERT(state);
    ASSERT(payload);
    buf[0] = '\0';

    if (payload->IsDataType(SCRIPT_TEXT_DND)) {
        if (!payload->Data || payload->DataSize <= 0) return false;
        const size_t len = strnlen((const char*)payload->Data, (size_t)payload->DataSize);
        if (len == 0 || len >= cap) return false;
        memcpy(buf, payload->Data, len);
        buf[len] = '\0';
        return true;
    }

    if ((payload->IsDataType(TIMELINE_SERIES_DND) || payload->IsDataType(DISTRIBUTION_SERIES_DND)) && payload->DataSize == (int)sizeof(SeriesDragPayload)) {
        // The values themselves: a summary of them (mean, extent) is not what attr() reads
        const SeriesKey& key = ((const SeriesDragPayload*)payload->Data)->key;
        if (key.source != SeriesSource_System || key.variant != SeriesVariant_Values) return false;
        const str_t path = str_from_cstr(key.path);
        const md_attributes_t* table = series_table(state, SeriesSource_System, path);
        const md_attribute_t* attr = table ? md_attributes_find(table, path) : nullptr;
        return attr && script_editor::attr_expression(buf, cap, table, attr);
    }

    return false;
}

static void draw_script_editor_window(ApplicationState* state) {
    ASSERT(state);

    ImGui::SetNextWindowSize({300,200}, ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Script Editor", &state->show_script_window, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
        ImGui::SetWindowSize(ImVec2(800, 600), ImGuiCond_FirstUseEver);
        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("File")) {
                char path_buf[1024] = "";
                if (ImGui::MenuItem("Load")) {
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Open, STR_LIT("txt"))) {
                        str_t txt = load_textfile(str_from_cstr(path_buf), frame_alloc);
                        std::string str(txt.ptr, txt.len);
                        state->editor.SetText(str);
                    }
                }
                if (ImGui::MenuItem("Save")) {
                    auto textToSave = state->editor.GetText();
                    if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, STR_LIT("txt"))) {
                        str_t path = str_t{path_buf, strnlen(path_buf, sizeof(path_buf))};
                        md_file_t file = {0};
                        if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                            md_file_write(file, textToSave.c_str(), textToSave.length());
                            md_file_close(&file);
                        } else {
                            VIAMD_LOG_ERROR("Failed to open file '%s' for saving script", path_buf);
                        }
                    }
                }
                if (ImGui::MenuItem("Export")) {
                    state->show_property_export_window = true;
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Edit")) {
                bool ro = state->editor.IsReadOnlyEnabled();
                if (ImGui::MenuItem("Read-only mode", nullptr, &ro))
                    state->editor.SetReadOnlyEnabled(ro);
                ImGui::Separator();

                const bool sel = state->editor.AnyCursorHasSelection();
                if (ImGui::MenuItem("Undo", "Ctrl-Z", nullptr, !ro && state->editor.CanUndo()))
                    state->editor.Undo();
                if (ImGui::MenuItem("Redo", "Ctrl-Y", nullptr, !ro && state->editor.CanRedo()))
                    state->editor.Redo();

                ImGui::Separator();

                if (ImGui::MenuItem("Copy", "Ctrl-C", nullptr, sel))
                    state->editor.Copy();
                if (ImGui::MenuItem("Cut", "Ctrl-X", nullptr, !ro && sel))
                    state->editor.Cut();
                if (ImGui::MenuItem("Delete", "Del", nullptr, !ro && sel))
                    state->editor.ReplaceTextInAllCursors("");
                if (ImGui::MenuItem("Paste", "Ctrl-V", nullptr, !ro && ImGui::GetClipboardText() != nullptr))
                    state->editor.Paste();

                ImGui::Separator();

                if (ImGui::MenuItem("Select all", "Ctrl-A"))
                    state->editor.SelectAll();
                if (ImGui::MenuItem("Find / Replace", "Ctrl-F"))
                    state->editor.OpenFindReplaceWindow();

                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Settings")) {
                if (ImGui::MenuItem("Dark palette"))
                    state->editor.SetPalette(TextEditor::GetDarkPalette());
                if (ImGui::MenuItem("Light palette"))
                    state->editor.SetPalette(TextEditor::GetLightPalette());
                if (ImGui::MenuItem("Retro blue palette"))
                    state->editor.SetPalette(script_editor::retro_blue_palette());
                ImGui::Separator();
                ImGui::ColorEdit4("Point Color",    state->script.point_color.elem);
                ImGui::ColorEdit4("Line Color",     state->script.line_color.elem);
                ImGui::ColorEdit4("Triangle Color", state->script.triangle_color.elem);
                ImGui::ColorEdit4("Text Color",     state->script.text_color.elem);
                ImGui::ColorEdit4("Text Bg Color",  state->script.text_bg_color.elem);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("Script reference", "F1")) {
                    open_script_reference(state, {});
                }
                const std::string word = script_editor::word_at_cursor(state->editor);
                char label[128];
                if (word.empty()) snprintf(label, sizeof(label), "Look up word under cursor");
                else snprintf(label, sizeof(label), "Look up '%s'", word.c_str());
                if (ImGui::MenuItem(label, "F1", nullptr, !word.empty())) {
                    open_script_reference(state, {word.data(), word.size()}, false);
                }
                ImGui::EndMenu();
            }

            ImGui::EndMenuBar();
        }

        // What Evaluate would do. It says so, and does nothing when there is nothing to do.
        const ScriptPlan plan = script_plan(*state, frame_alloc);
        char btn_text[96];
        bool can_eval = false;
        if (!plan.valid) {
            snprintf(btn_text, sizeof(btn_text), "Evaluate###evaluate");
        } else if (plan.running && plan.same) {
            snprintf(btn_text, sizeof(btn_text), "Evaluating, %u frames left###evaluate", plan.num_pending);
        } else if (plan.num_pending > 0) {
            snprintf(btn_text, sizeof(btn_text), "Evaluate %u frames###evaluate", plan.num_pending);
            can_eval = true;
        } else if (plan.changes) {
            snprintf(btn_text, sizeof(btn_text), "Evaluate###evaluate");
            can_eval = true;
        } else {
            snprintf(btn_text, sizeof(btn_text), "Up to date###evaluate");
        }

        // How each property stands, beside the line it is named on
        md_array(ScriptGutterMark) marks = nullptr;
        if (plan.valid) {
            const uint32_t* line_offsets = state->editor_markers.line_offsets;
            const size_t num_lines = md_array_size(line_offsets);
            for (size_t i = 0; i < plan.num_props; ++i) {
                const md_script_range_marker_t range = md_script_ir_property_range(state->script.ir, plan.names[i]);
                if (range.end <= range.beg || num_lines == 0) continue;
                // The last line which begins at or before it
                size_t lo = 0, hi = num_lines;
                while (hi - lo > 1) {
                    const size_t mid = (lo + hi) / 2;
                    if (line_offsets[mid] <= (uint32_t)range.beg) lo = mid; else hi = mid;
                }
                const md_script_eval_plan_t& p = plan.props[i];
                ScriptGutterMark mark = { (int)lo, p.status, 0.0f, "", "", 0 };
                str_copy_to_char_buf(mark.ident, sizeof(mark.ident), plan.names[i]);
                mark.kinds = (uint32_t)md_script_ir_property_flags(state->script.ir, plan.names[i]);
                const uint32_t needed = p.num_frames + p.num_pending;
                mark.done = p.num_pending == 0 ? 1.0f : needed ? (float)p.num_frames / (float)needed : 0.0f;
                script_plan_describe(mark.text, sizeof(mark.text), plan, i);
                md_array_push(marks, mark, frame_alloc);
            }
        }
        const size_t num_marks = md_array_size(marks);
        state->editor.SetLineDecorator(1, [state, marks, num_marks](TextEditor::Decorator& d) {
            for (size_t i = 0; i < num_marks; ++i) {
                if ((size_t)marks[i].line == d.line) {
                    draw_script_gutter_mark(state, d, marks[i]);
                    break;
                }
            }
        });

        const ImVec2 content_size = ImGui::GetContentRegionAvail();
        const ImVec2 label_size = ImGui::CalcTextSize(btn_text, NULL, true) * ImVec2(1.2, 1.0);
        const ImVec2 btn_size = ImGui::CalcItemSize(ImVec2(0,0), label_size.x + ImGui::GetStyle().FramePadding.x * 2.0f, label_size.y + ImGui::GetStyle().FramePadding.y * 2.0f);
        const ImVec2 text_size(content_size - ImVec2(0, btn_size.y + ImGui::GetStyle().ItemSpacing.y));

        // Shift+Enter evaluates the script. The editor binds it to "insert line above", so claim it before the editor
        // sees it (shortcut routes are resolved from the previous frame, hence the focus from the last draw).
        bool eval = false;
        if (ImGui::Shortcut(KEY_SCRIPT_EVALUATE, ImGuiInputFlags_RouteGlobal | ImGuiInputFlags_RouteOverFocused)) {
            eval = true;
        }

        // While a visualization is hovered, the mouse wheel steps through its elements instead of scrolling
        const script_editor::Marker* prev_hovered = script_editor::markers_hovered(&state->editor_markers);
        ImGuiWindowFlags editor_flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_HorizontalScrollbar;
        if (prev_hovered && prev_hovered->type == script_editor::MarkerType_Visualization) {
            editor_flags |= ImGuiWindowFlags_NoScrollWithMouse;
        }

        if (state->editor.Render("TextEditor", text_size, ImGuiChildFlags_None, editor_flags)) {
            state->script.compile_ir = true;
            state->script.time_since_last_change = 0;
        }
        const bool editor_hovered = ImGui::IsItemHovered();
        state->editor_focused = script_editor::has_focus_after_render();

        // Attributes dropped in go where a click would put the cursor, as the expression that reads them
        if (ImGui::BeginDragDropTarget()) {
            const ImGuiPayload* payload = ImGui::GetDragDropPayload();
            char text[512];
            if (payload && script_drop_text(text, sizeof(text), state, payload)) {
                script_editor::accept_drop(state->editor, payload->DataType, str_from_cstr(text));
            }
            ImGui::EndDragDropTarget();
        }

        const script_editor::Marker* hovered_marker = script_editor::markers_update(&state->editor_markers, &state->editor, editor_hovered);

        // F1 looks up the identifier under the mouse in the script reference, or else the one at the text cursor.
        // Hovering does not require focus: opening the reference takes the focus, and a second F1 over another
        // identifier would otherwise be ignored until the editor is clicked again.
        if ((editor_hovered || state->editor_focused) && ImGui::IsKeyPressed(ImGuiKey_F1, false)) {
            std::string word;
            if (editor_hovered) {
                word = script_editor::word_at_mouse(state->editor, ImGui::GetMousePos());
            }
            if (word.empty() && state->editor_focused) {
                word = script_editor::word_at_cursor(state->editor);
            }
            open_script_reference(state, {word.data(), word.size()}, false);
        }

        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + content_size.x - btn_size.x);
        ImGui::SetCursorPosY(ImGui::GetCursorPosY() + ImGui::GetStyle().ItemSpacing.y);

        if (!can_eval) ImGui::PushDisabled();
        if (ImGui::Button(btn_text, btn_size)) {
            eval = true;
        }
        if (!can_eval) ImGui::PopDisabled();
        if (plan.valid && ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled | ImGuiHoveredFlags_DelayShort)) {
            ImGui::BeginTooltip();
            char line[160];
            for (size_t i = 0; i < plan.num_props; ++i) {
                script_plan_describe(line, sizeof(line), plan, i);
                ImGui::TextUnformatted(line);
            }
            // What is shown which the script no longer has
            if (state->script.eval_ir && !plan.same) {
                const str_t* evaluated = md_script_ir_property_names(state->script.eval_ir);
                for (size_t j = 0; j < md_script_ir_property_count(state->script.eval_ir); ++j) {
                    bool kept = false;
                    for (size_t i = 0; i < plan.num_props && !kept; ++i) {
                        kept = str_eq(plan.names[i], evaluated[j]) || (plan.props[i].status == MD_SCRIPT_EVAL_STATUS_HELD && str_eq(plan.props[i].from, evaluated[j]));
                    }
                    if (!kept) ImGui::TextDisabled(STR_FMT ": no longer in the script", STR_ARG(evaluated[j]));
                }
            }
            ImGui::EndTooltip();
        }

        if (eval && can_eval) {
            // The script as it was evaluated: what remains of it. Otherwise a new evaluation, which takes over what
            // it can from this one.
            if (plan.same) {
                state->script.evaluate = true;
            } else {
                state->script.eval_init = true;
            }
        }

        if (editor_hovered) {
            md_bitmask_clear(&state->selection.highlight_mask);
        }

        if (hovered_marker) {
            if (hovered_marker->atoms || !md_script_vis_ref_empty(hovered_marker->vis_ref)) {
                if (hovered_marker->atoms) {
                    // Errors and warnings: the atoms the message is about
                    md_bitmask_assign(&state->selection.highlight_mask, *hovered_marker->atoms);
                }
                else if (hovered_marker->type == script_editor::MarkerType_Visualization) {
                    // Clear hovered property
                    script_set_hovered_property(state, STR_LIT(""));
                    // The markers come from the IR last compiled. One which did not compile is gone, and its
                    // markers resolve in nothing
                    const md_script_vis_ref_t ref = hovered_marker->vis_ref;
                    if (const md_script_ir_t* ir = script_ir_of(state, ref)) {
                        str_t ref_ident = md_script_vis_ref_ident(ir, ref);
                        int ref_dim  = md_script_vis_ref_dim(ir, ref);

                        if (ref_dim > 1) {
                            int delta = (int)ImGui::GetIO().MouseWheel;
                            if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
                                delta *= 10;
                            }
                            state->script.sub_idx += delta;
                            state->script.sub_idx = CLAMP(state->script.sub_idx, -1, (int)ref_dim - 1);
                        }

                        script_visualize_ref(state, ref, state->script.sub_idx, 0);
                        script_set_hovered_property(state, ref_ident, state->script.sub_idx);
                    }
                }

                bool lm_click = ImGui::IsMouseClicked(ImGuiMouseButton_Left);
                bool rm_click = ImGui::IsMouseClicked(ImGuiMouseButton_Right);
                if (ImGui::IsKeyDown(ImGuiMod_Shift) && (lm_click || rm_click)) {
                    SelectionOperator op = lm_click ? SelectionOperator::Or : SelectionOperator::AndNot;
                    modify_selection(state, md_bitmask_view(&state->selection.highlight_mask), op);
                }
            }

            // After the visualization is set: the value is computed with it
            if (!str_empty(hovered_marker->text)) {
                ImGui::BeginTooltip();
                ImGui::Text(STR_FMT, STR_ARG(hovered_marker->text));
                if (hovered_marker->type == script_editor::MarkerType_Visualization) {
                    const str_t value = script_vis_value_text(state);
                    // A value known since the compilation is in the text already, unless it was too long to keep there
                    const bool dynamic = str_find_str(NULL, hovered_marker->text, STR_LIT("[d]"));
                    const bool has_value = str_find_char(NULL, hovered_marker->text, '\n');
                    if (!str_empty(value) && (dynamic || !has_value || state->script.sub_idx != -1)) {
                        ImGui::Separator();
                        if (dynamic) ImGui::TextDisabled("At the frame shown");
                        ImGui::Text(STR_FMT, STR_ARG(value));
                    }
                }
                if (state->script.sub_idx != -1) {
                    ImGui::TextDisabled("Element %i, the mouse wheel steps through them", state->script.sub_idx + 1);
                }
                ImGui::EndTooltip();
            }
        } else {
            state->script.sub_idx = -1;
        }
    }
    ImGui::End();
}

static bool export_xvg(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename) {
    ASSERT(column_data);
    ASSERT(column_labels);

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        VIAMD_LOG_ERROR("Failed to open file '" STR_FMT "' to write data.", STR_ARG(filename));
        return false;
    }    

    time_t t;
    struct tm* info;
    time(&t);
    info = localtime(&t);

    // Print Header
    md_file_printf(file, "# This file was created %s", asctime(info));
    md_file_printf(file, "# Created by:\n");
    md_file_printf(file, "# VIAMD \n");

    // Print Legend Meta
    md_file_printf(file, "@    title \"VIAMD Properties\"\n");
    md_file_printf(file, "@    xaxis  label \"Time\"\n");
    md_file_printf(file, "@ TYPE xy\n");
    md_file_printf(file, "@ view 0.15, 0.15, 0.75, 0.85\n");
    md_file_printf(file, "@ legend on\n");
    md_file_printf(file, "@ legend box on\n");
    md_file_printf(file, "@ legend loctype view\n");
    md_file_printf(file, "@ legend 0.78, 0.8\n");
    md_file_printf(file, "@ legend length %i\n", num_columns);

    for (size_t j = 0; j < num_columns; ++j) {
        md_file_printf(file, "@ s%zu legend \"%s\"\n", j, column_labels[j]);
    }

    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_columns; ++j) {
            md_file_printf(file, "%12.6f ", column_data[j][i]);
        }
        md_file_printf(file, "\n");
    }

    md_file_close(&file);
    VIAMD_LOG_SUCCESS("Successfully exported XVG file to '%.*s'", (int)filename.len, filename.ptr);
    return true;
}

static bool export_csv(const float* column_data[], const char* column_labels[], size_t num_columns, size_t num_rows, str_t filename) {
    ASSERT(column_data);
    ASSERT(column_labels);

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        VIAMD_LOG_ERROR("Failed to open file '%.*s' to write data.", (int)filename.len, filename.ptr);
        return false;
    }
    
    for (size_t i = 0; i < num_columns; ++i) {
        md_file_printf(file, "%s,", column_labels[i]);
    }
    md_file_printf(file, "\n");

    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < num_columns; ++j) {
            md_file_printf(file, "%.6g,", column_data[j][i]);
        }
        md_file_printf(file, "\n");
    }

    md_file_close(&file);
    VIAMD_LOG_SUCCESS("Successfully exported CSV file to '%.*s'", (int)filename.len, filename.ptr);
    return true;
}

static bool export_cube(const ApplicationState& data, const SeriesVolumeView& volume, str_t filename) {
    // @NOTE: First we need to extract some meta data for the cube format, we need the atom indices/bits for any SDF
    // And the origin + extent of the volume in spatial coordinates (Ångström)

    const float* volume_values = volume.data;
    const md_script_vis_ref_t vis_ref = volume.vis_ref;
    if (!volume_values) {
        VIAMD_LOG_ERROR("Export Cube: The property to be exported did not exist");
        return false;
    }

    if (md_script_vis_ref_empty(vis_ref)) {
        VIAMD_LOG_ERROR("Export Cube: Missing input visualization data");
        return false;
    }

	md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
	defer { md_temp_end(temp); };

    const md_system_t& sys = data.mold.sys;
    md_system_state_t state = { .alloc = temp.arena };
    if (!extract_frame(&data, 0, &state)) {
        return false;
    }

    bool result = false;
    md_script_vis_t vis = { 0 };
    md_script_vis_init(&vis, temp.arena);
        
    if (md_script_ir_valid(data.script.eval_ir)) {
        md_script_vis_ctx_t ctx = {
            .ir = data.script.eval_ir,
            .sys = &data.mold.sys,
            .state = &data.mold.state,
        };
        // The structures of the volume exported: of an array, those of its member
        result = md_script_vis_eval_ref(&vis, vis_ref, (int)volume.member, &ctx, MD_SCRIPT_VISUALIZE_ATOMS | MD_SCRIPT_VISUALIZE_SDF);
    }

    if (result == true) {
        md_file_t file = { 0 };
        if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
            VIAMD_LOG_ERROR("Failed to open file '%.*s' in order to write to it.", (int)filename.len, filename.ptr);
            return false;
        }

        // Two comment lines
        md_file_printf(file, "EXPORTED DENSITY VOLUME FROM VIAMD, UNITS IN BOHR\n");
        md_file_printf(file, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");

        if (md_array_size(vis.sdf.structures) > 0) {
            const float angstrom_to_bohr = (float)(1.0 / 0.529177210903);

            // transformation matrix from world to volume
            mat4_t M = vis.sdf.matrices[0];
            const md_bitset_t* bf = &vis.sdf.structures[0];
            const int num_atoms = (int)md_bitset_count(*bf);
            const int vol_dim[3] = {(int)volume.dim[0], (int)volume.dim[1], (int)volume.dim[2]};
            const float* values = volume_values;
            const double extent = vis.sdf.extent * 2.0 * angstrom_to_bohr;
            const double voxel_ext[3] = {
                (double)extent / (double)vol_dim[0],
                (double)extent / (double)vol_dim[1],
                (double)extent / (double)vol_dim[2],
            };

            const double half_ext = extent * 0.5;

            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", -num_atoms, -half_ext, -half_ext, -half_ext);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[0], voxel_ext[0], 0.0, 0.0);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[1], 0.0, voxel_ext[1], 0.0);
            md_file_printf(file, "%5i %12.6f %12.6f %12.6f\n", vol_dim[2], 0.0, 0.0, voxel_ext[2]);

            const float scl = angstrom_to_bohr;
            M = mat4_mul(mat4_scale(scl, scl, scl), M);

            md_bitset_iter_t it = md_bitset_iter(*bf);
            while (md_bitset_iter_next(&it)) {
                size_t i = it.idx;
                vec3_t coord = md_state_coord(&state, i);
                coord = mat4_mul_vec3(M, coord, 1.0f);
                int anum     = md_atom_atomic_number(&sys.atom, i);
                float charge = (float)anum;
                md_file_printf(file, "%5i %12.6f %12.6f %12.6f %12.6f\n", anum, charge, coord.x, coord.y, coord.z);
            }

            // This entry somehow relates to the number of densities
            md_file_printf(file, "%5i %5i\n", 1, 1);

            // Write density data
            int count = 0;
            for (int x = 0; x < vol_dim[0]; ++x) {
                for (int y = 0; y < vol_dim[1]; ++y) {
                    for (int z = 0; z < vol_dim[2]; ++z) {
                        int idx = z * vol_dim[0] * vol_dim[1] + y * vol_dim[0] + x;
                        float val = values[idx];
                        md_file_printf(file, " %12.6E", val);
                        if (++count % 6 == 0) md_file_printf(file, "\n");
                    }
                }
            }
        }

        md_file_close(&file);
    } else {
        VIAMD_LOG_ERROR("Failed to visualize volume for export.");
        return false;
    }

    return true;
}

#define APPEND_BUF(buf, len, fmt, ...) (len += snprintf(buf + len, MAX(0, (int)sizeof(buf) - len), fmt, ##__VA_ARGS__) + 1)

static void draw_property_export_window(ApplicationState* data) {
    ASSERT(data);

    struct ExportFormat {
        str_t lbl;
        str_t ext;
    };

    ExportFormat table_formats[] {
        {STR_INIT("XVG"), STR_INIT("xvg")},
        {STR_INIT("CSV"), STR_INIT("csv")}
    };

    ExportFormat volume_formats[] {
        {STR_INIT("Gaussian Cube"), STR_INIT("cube")},
    };

    enum ExportKind {
        ExportKind_Temporal = 0,
        ExportKind_Distribution,
        ExportKind_Volume,
    };

    if (ImGui::Begin("Property Export", &data->show_property_export_window)) {
        static int kind = ExportKind_Temporal;
        static SeriesKey selected = {};
        static int table_format  = 0;
        static int volume_format = 0;
        static int volume_member = 0;   // of an array of volumes; -1 for each of them, in a file of its own
        static int bins_exp      = 7;   // 128 bins

        if (task_system::task_is_running(data->tasks.evaluate)) {
            ImGui::Text("The properties are currently being evaluated, please wait...");
            ImGui::End();
            return;
        }

        md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
        defer { md_temp_end(temp); };
        md_allocator_i* alloc = md_temp_allocator(temp);

        // What there is to export of the chosen kind: the same series the plots offer
        md_array(SeriesKey) candidates = 0;
        {
            const uint32_t script_kinds =
                kind == ExportKind_Temporal     ? (uint32_t)MD_SCRIPT_PROPERTY_FLAG_TEMPORAL :
                kind == ExportKind_Distribution ? (uint32_t)(MD_SCRIPT_PROPERTY_FLAG_TEMPORAL | MD_SCRIPT_PROPERTY_FLAG_DISTRIBUTION) :
                                                  (uint32_t)MD_SCRIPT_PROPERTY_FLAG_VOLUME;
            auto push_script = [&](uint32_t period) {
                series_for_each_script_property(data, SeriesSource_Script, script_kinds, [&](const SeriesKey& over_all) {
                    SeriesKey key = over_all;
                    key.period = period;
                    md_array_push(candidates, key, alloc);
                    const bool temporal = md_script_ir_property_flags(data->script.eval_ir, series_script_ident(key)) & MD_SCRIPT_PROPERTY_FLAG_TEMPORAL;
                    if (kind == ExportKind_Distribution && temporal && script_property_population(data, key) > 1) {
                        SeriesKey agg = key;
                        agg.variant = SeriesVariant_Aggregate;
                        md_array_push(candidates, agg, alloc);
                    }
                });
            };
            // Over each period too, "r (Production)": what is over time is the same over a period
            push_script(0);
            if (kind != ExportKind_Temporal) {
                for (int i = 0; i < data->timeline.periods.count; ++i) {
                    push_script(data->timeline.periods.items[i].id);
                }
            }
            if (kind != ExportKind_Volume) {
                str_t groups[64];
                const size_t num_groups = MIN(system_series_groups(groups, ARRAY_SIZE(groups), data), ARRAY_SIZE(groups));
                for (size_t g = 0; g < num_groups; ++g) {
                    const size_t num = system_series_members(nullptr, 0, data, groups[g]);
                    str_t* paths = md_temp_alloc_array(temp, str_t, num + 1);
                    system_series_members(paths, num, data, groups[g]);
                    for (size_t i = 0; i < num; ++i) {
                        md_array_push(candidates, series_key(SeriesSource_System, paths[i]), alloc);
                    }
                }
            }
        }

        bool selected_valid = false;
        for (size_t i = 0; i < md_array_size(candidates); ++i) {
            if (series_key_equal(candidates[i], selected)) { selected_valid = true; break; }
        }
        if (!selected_valid) {
            selected = md_array_size(candidates) ? candidates[0] : SeriesKey{};
            selected_valid = md_array_size(candidates) > 0;
        }

        ImGui::PushItemWidth(200);
        ImGui::Combo("Data Type", &kind, "Temporal\0Distribution\0Density Volume\0");
        ImGui::Separator();

        char selected_label[96] = "";
        if (selected_valid) {
            series_label(selected_label, sizeof(selected_label), data, selected);
        }
        if (ImGui::BeginCombo("Property", selected_label)) {
            for (size_t i = 0; i < md_array_size(candidates); ++i) {
                char label[96];
                series_label(label, sizeof(label), data, candidates[i]);
                ImGui::PushID((int)i);
                if (ImGui::Selectable(label, series_key_equal(candidates[i], selected))) {
                    selected = candidates[i];
                }
                if (ImGui::IsItemHovered() && candidates[i].source == SeriesSource_System) {
                    ImGui::SetTooltip("%s", candidates[i].path);
                }
                ImGui::PopID();
            }
            ImGui::EndCombo();
        }
        if (md_array_size(candidates) == 0) {
            ImGui::TextDisabled("Nothing of this kind, evaluate the script or load an energy file");
        }

        if (kind == ExportKind_Distribution) {
            char bins_fmt[32];
            snprintf(bins_fmt, sizeof(bins_fmt), "%d", 1 << bins_exp);
            ImGui::SliderInt("Num Bins", &bins_exp, 5, 10, bins_fmt);
        }

        str_t file_extension = {};
        if (kind == ExportKind_Volume) {
            if (ImGui::BeginCombo("File Format", volume_formats[volume_format].lbl.ptr)) {
                for (int i = 0; i < (int)ARRAY_SIZE(volume_formats); ++i) {
                    if (ImGui::Selectable(volume_formats[i].lbl.ptr, volume_format == i)) {
                        volume_format = i;
                    }
                }
                ImGui::EndCombo();
            }
            file_extension = volume_formats[volume_format].ext;

            // Of an array of volumes ('sdf(...) in residue(1:n)'), one as the script names it, or each to a file of its
            // own: each volume is about its own reference structures, which go along in the file
            const uint32_t count = selected_valid ? series_volume_count(data, selected) : 0;
            if (count > 1) {
                const str_t ident = series_script_ident(selected);
                volume_member = CLAMP(volume_member, -1, (int)count - 1);
                char preview[96];
                if (volume_member < 0) {
                    snprintf(preview, sizeof(preview), "All, one file each");
                } else {
                    snprintf(preview, sizeof(preview), STR_FMT "[%d]", STR_ARG(ident), volume_member + 1);
                }
                if (ImGui::BeginCombo("Volume", preview)) {
                    if (ImGui::Selectable("All, one file each", volume_member < 0)) {
                        volume_member = -1;
                    }
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip("'name.cube' is written as 'name_1.cube' to 'name_%u.cube'", count);
                    }
                    for (uint32_t m = 0; m < count; ++m) {
                        char item[96];
                        snprintf(item, sizeof(item), STR_FMT "[%u]", STR_ARG(ident), m + 1);
                        if (ImGui::Selectable(item, volume_member == (int)m)) {
                            volume_member = (int)m;
                        }
                        if (ImGui::IsItemHovered()) {
                            script_visualize_ref(data, md_script_ir_property_vis_ref(data->script.eval_ir, ident), (int)m, MD_SCRIPT_VISUALIZE_DEFAULT);
                        }
                    }
                    ImGui::EndCombo();
                }
            }
        } else {
            if (ImGui::BeginCombo("File Format", table_formats[table_format].lbl.ptr)) {
                for (int i = 0; i < (int)ARRAY_SIZE(table_formats); ++i) {
                    if (ImGui::Selectable(table_formats[i].lbl.ptr, table_format == i)) {
                        table_format = i;
                    }
                }
                ImGui::EndCombo();
            }
            file_extension = table_formats[table_format].ext;
        }

        if (!selected_valid) ImGui::PushDisabled();
        const bool export_clicked = ImGui::Button("Export");
        if (!selected_valid) ImGui::PopDisabled();
        ImGui::PopItemWidth();

        char path_buf[1024];
        if (export_clicked && application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, file_extension)) {
            const str_t path = {path_buf, strnlen(path_buf, sizeof(path_buf))};
            char label[96];
            series_label(label, sizeof(label), data, selected);

            if (kind == ExportKind_Volume) {
                // One volume, or of an array each to a file of its own: 'name.cube' as 'name_1.cube', 'name_2.cube', ...
                const uint32_t count = series_volume_count(data, selected);
                const bool each = count > 1 && volume_member < 0;
                const uint32_t beg = each ? 0 : (uint32_t)MAX(volume_member, 0);
                const uint32_t end = each ? count : beg + 1;
                str_t stem = path;
                str_t ext  = file_extension;
                size_t sep;
                if (each && extract_ext(&ext, path) && !str_find_char(&sep, ext, '/') && !str_find_char(&sep, ext, '\\')) {
                    extract_file_path_without_ext(&stem, path);
                } else {
                    ext = file_extension;
                }
                for (uint32_t m = beg; m < end; ++m) {
                    SeriesVolumeView view;
                    const str_t file = each ? str_printf(alloc, STR_FMT "_%u." STR_FMT, STR_ARG(stem), m + 1, STR_ARG(ext)) : path;
                    if (!series_resolve_volume(&view, data, selected, m)) {
                        VIAMD_LOG_ERROR("'%s' has no volume to export", label);
                        break;
                    }
                    if (!str_eq(file_extension, STR_LIT("cube")) || !export_cube(*data, view, file)) {
                        break;
                    }
                    VIAMD_LOG_SUCCESS("Successfully exported property '%s' to '" STR_FMT "'", view.label, STR_ARG(file));
                }
            } else {
                // Columns in the units the plots show them in, so a column and its axis label agree
                md_array(const float*) column_data   = 0;
                md_array(str_t)        column_labels = 0;
                md_array(str_t)        legends       = 0;
                str_t x_label = {};
                str_t y_label = str_from_cstr(label);
                size_t num_rows = 0;
                bool ok = false;

                auto push_columns = [&](int dim, auto value) {
                    for (int k = 0; k < dim; ++k) {
                        float* col = md_temp_alloc_array(temp, float, num_rows);
                        for (size_t r = 0; r < num_rows; ++r) {
                            col[r] = value(r, k);
                        }
                        md_array_push(column_data, col, alloc);
                        if (dim > 1) {
                            const str_t member = series_member_label(data, selected, k);
                            const str_t legend = str_empty(member) ? str_printf(alloc, "%s[%i]", label, k + 1)
                                                                   : str_printf(alloc, "%s[%i] " STR_FMT, label, k + 1, STR_ARG(member));
                            md_array_push(legends, legend, alloc);
                            md_array_push(column_labels, legend, alloc);
                        } else {
                            md_array_push(column_labels, y_label, alloc);
                        }
                    }
                };

                if (kind == ExportKind_Temporal) {
                    SeriesTemporalView view;
                    if (series_resolve_temporal(&view, data, selected)) {
                        num_rows = (size_t)view.num_samples;
                        char time_unit_str[32] = "";
                        md_unit_print(time_unit_str, sizeof(time_unit_str), data->timeline.time_unit);
                        x_label = md_unit_is_none(data->timeline.time_unit) ? STR_LIT("Frame") : str_printf(alloc, "Time (%s)", time_unit_str);
                        if (view.unit_str[0]) y_label = str_printf(alloc, "%s (%s)", label, view.unit_str);

                        float* x = md_temp_alloc_array(temp, float, num_rows);
                        MEMCPY(x, view.x, num_rows * sizeof(float));
                        md_array_push(column_data, x, alloc);
                        md_array_push(column_labels, x_label, alloc);
                        push_columns(view.dim, [&](size_t r, int k) { return (float)series_temporal_point(view, (int)r, k).y; });
                        ok = true;
                    }
                } else {
                    SeriesHistogramView view;
                    if (series_resolve_histogram(&view, data, selected, 1 << bins_exp)) {
                        num_rows = (size_t)view.num_bins;
                        x_label = str_copy(str_from_cstr(view.x_unit_str), alloc);
                        if (view.y_unit_str[0]) y_label = str_printf(alloc, "%s (%s)", label, view.y_unit_str);

                        // The centre of each bin, where the plot draws it
                        float* x = md_temp_alloc_array(temp, float, num_rows);
                        for (size_t r = 0; r < num_rows; ++r) {
                            x[r] = (float)series_histogram_point(view, (int)r, 0).x;
                        }
                        md_array_push(column_data, x, alloc);
                        md_array_push(column_labels, x_label, alloc);
                        push_columns(view.dim, [&](size_t r, int k) { return (float)series_histogram_point(view, (int)r, k).y; });
                        ok = true;
                    }
                }

                if (!ok) {
                    VIAMD_LOG_ERROR("'%s' has nothing to export right now", label);
                } else {
                    str_t out_str = {};
                    if (str_eq(file_extension, STR_LIT("xvg"))) {
                        str_t header = md_xvg_format_header(str_from_cstr(label), x_label, y_label, md_array_size(legends), legends, alloc);
                        out_str = md_xvg_format(header, md_array_size(column_data), num_rows, column_data, alloc);
                    } else if (str_eq(file_extension, STR_LIT("csv"))) {
                        out_str = md_csv_write_to_str(column_data, column_labels, md_array_size(column_data), num_rows, alloc);
                    }
                    md_file_t file = {0};
                    if (!str_empty(out_str) && md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                        md_file_write(file, out_str.ptr, out_str.len);
                        md_file_close(&file);
                        VIAMD_LOG_SUCCESS("Successfully exported property '%s' to '" STR_FMT "'", label, STR_ARG(path));
                    } else {
                        VIAMD_LOG_ERROR("Failed to write '" STR_FMT "'", STR_ARG(path));
                    }
                }
            }
        }
    }
    ImGui::End();
}

md_array(int) extract_bitfield_indices(const md_bitset_t* bf, md_allocator_i* alloc) {
    size_t count = md_bitset_count(*bf);
    md_array(int) indices = md_array_create(int, count, alloc);
    md_bitset_iter_extract_indices(indices, count, md_bitset_iter(*bf));
    return indices;
}

static void xyz_write_frame(md_file_t file, const int* atomic_nr, const md_system_state_t* state, const int32_t* atom_indices, size_t num_atoms) {
    md_file_printf(file, "%zu\n", num_atoms);
    md_file_printf(file, "XYZ format exported from VIAMD\n");
    for (size_t i = 0; i < num_atoms; ++i) {
        int idx = atom_indices ? atom_indices[i] : (int)i;
        md_file_printf(file, "%-2d %12.6f %12.6f %12.6f\n", atomic_nr[idx], state->xyz[idx].x, state->xyz[idx].y, state->xyz[idx].z);
    }
}

void draw_structure_export_window(ApplicationState* data) {
    ASSERT(data);

    if (ImGui::Begin("Structure Export", &data->structure_export.show_window)) {
        md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
        defer { md_temp_end(temp); };

        const md_system_t* sys = &data->mold.sys;
        const md_system_state_t* state = &data->mold.state;
        const bool traj = run_num_frames(data) > 0;
        auto& struct_exp = data->structure_export;

        static const char* atom_mask_options[] = {
            "Active Selection",
            "Visible Atoms",
            "Query",
            "All Atoms",
        };

        static const char* frame_mask_options[] = {
            "All Frames",
            "Current Frame",
        };

        static const char* file_formats[] = {
            "xyz",
            "pdb"
        };

        //ImGui::PushItemWidth(200);

        struct_exp.selected_atom_filter = CLAMP(struct_exp.selected_atom_filter, 0, (int)ARRAY_SIZE(atom_mask_options) - 1);
        ImGui::Combo("Atoms to Export", &struct_exp.selected_atom_filter, atom_mask_options, (int)ARRAY_SIZE(atom_mask_options));

        if (struct_exp.selected_atom_filter == 2) { // Query
            auto& query = struct_exp.query;
            query.requires_evaluation |= filter_input("Query", &query.filter, filter_context(data));
            bool preview = ImGui::IsItemFocused() || ImGui::IsItemHovered();

            if (query.requires_evaluation) {
                filter_evaluate(&query.filter, &query.mask, filter_context(data), persistent_alloc);
                query.requires_evaluation = false;
            }

            if (ImGui::IsWindowHovered()) {
                // Clear highlight mask
                md_bitmask_clear(&data->selection.highlight_mask);
            }

            if (query.filter.result.success && preview) {
                md_bitmask_assign(&data->selection.highlight_mask, query.mask);
            }
        }

        // Every frame, the current one, or those of a timeline period
        const TimelinePeriod* traj_period = timeline_periods_find(&data->timeline.periods, struct_exp.selected_traj_period);
        struct_exp.selected_traj_filter = CLAMP(struct_exp.selected_traj_filter, 0, 2);
        if (struct_exp.selected_traj_filter == 2 && !traj_period) struct_exp.selected_traj_filter = 0;
        if (traj) {
            const char* preview = struct_exp.selected_traj_filter == 2 ? traj_period->name : frame_mask_options[struct_exp.selected_traj_filter];
            if (ImGui::BeginCombo("Frames to Export", preview)) {
                for (int i = 0; i < (int)ARRAY_SIZE(frame_mask_options); ++i) {
                    if (ImGui::Selectable(frame_mask_options[i], struct_exp.selected_traj_filter == i)) struct_exp.selected_traj_filter = i;
                }
                const TimelinePeriods& rs = data->timeline.periods;
                if (rs.count > 0) ImGui::SeparatorText("Period");
                for (int i = 0; i < rs.count; ++i) {
                    const TimelinePeriod& r = rs.items[i];
                    char label[TIMELINE_PERIOD_NAME_CAP + 16];
                    snprintf(label, sizeof(label), "%s###period%u", r.name, r.id);
                    if (ImGui::Selectable(label, struct_exp.selected_traj_filter == 2 && r.id == struct_exp.selected_traj_period)) {
                        struct_exp.selected_traj_filter = 2;
                        struct_exp.selected_traj_period = r.id;
                    }
                    if (ImGui::IsItemHovered()) timeline_periods_highlight(r.id);
                }
                ImGui::EndCombo();
            }
        }

        struct_exp.selected_file_format = CLAMP(struct_exp.selected_file_format, 0, (int)ARRAY_SIZE(file_formats) - 1);
        ImGui::Combo("File Format", &struct_exp.selected_file_format, file_formats, (int)ARRAY_SIZE(file_formats));

        //ImGui::PopItemWidth();

        auto extract_atom_indices = [data](int* out_atom_indices, size_t atom_index_cap, const md_system_state_t* state) -> size_t {
            switch (data->structure_export.selected_atom_filter) {
                case 0: { // Active Selection
                    return md_bitset_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitmask_iter(&data->selection.selection_mask));
                } break;
                case 1: { // Visible Atoms
                    return md_bitset_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitmask_iter(&data->representation.visibility_mask));
                } break;
                case 2: { // Query
                    auto& query = data->structure_export.query;
                    if (query.filter.result.is_dynamic) {
                        // What it picks depends on the coordinates: picked again in each frame exported (md_filter
                        // replaces what the mask held, so nothing of the previous frame stays)
                        md_filter_context_t ctx = filter_context(data);
                        ctx.state = state;
                        filter_evaluate(&query.filter, &query.mask, ctx, persistent_alloc);
                    }
                    if (!query.filter.result.success) {
                        VIAMD_LOG_ERROR("Cannot export structure, the query expression is invalid: %s", query.filter.result.error);
                        return 0;
                    } else {
                        return md_bitset_iter_extract_indices(out_atom_indices, atom_index_cap, md_bitset_iter(query.mask));
                    }
                } break;
                case 3: { // All Atoms
                    for (size_t i = 0; i < data->mold.sys.atom.count && i < atom_index_cap; ++i) {
                        out_atom_indices[i] = (int)i;
                    }
                    return data->mold.sys.atom.count;
                } break;
                default: {
                    VIAMD_LOG_ERROR("Invalid atom filter selection.");
                } break;
            }
            return 0;
        };

        if (ImGui::Button("Export")) {
            md_array(int32_t) frame_indices = 0;
            if (traj) {
                switch (struct_exp.selected_traj_filter) {
                case 0: { // All Frames
                    size_t num_frames = run_num_frames(data);
                    for (size_t i = 0; i < num_frames; ++i) {
                        md_array_push(frame_indices, (int32_t)i, frame_alloc);
                    }
                } break;
                case 1: { // Current Frame
                    int32_t current_frame = (int)(data->animation.frame + 0.5); // Round to nearest
                    md_array_push(frame_indices, current_frame, frame_alloc);
                } break;
                case 2: { // A timeline period
                    int beg = 0, end = -1;
                    if (!traj_period || !timeline_period_frames(&beg, &end, *traj_period, data->animation.frame, (int)run_num_frames(data))) {
                        beg = 0;
                        end = -1;
                    }
                    for (int i = beg; i <= end; ++i) {
                        md_array_push(frame_indices, (int32_t)i, frame_alloc);
                    }
                } break;
                default: {
                    VIAMD_LOG_ERROR("Invalid trajectory filter selection.");
                } break;
                }
            }

            char path_buf[1024] = { 0 };
            str_t ext = str_from_cstr(file_formats[struct_exp.selected_file_format]);
            if (application::file_dialog(path_buf, sizeof(path_buf), application::FileDialogFlag_Save, ext)) {
                str_t path = str_from_cstrn(path_buf, sizeof(path_buf));
                md_file_t file = {0};
                if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
                    // Write header for file format
                    switch (struct_exp.selected_file_format) {
                        case 0: { // XYZ
                            // No header needed
                        } break;
                        case 1: { // PDB
                            md_file_printf(file, "HEADER    EXPORTED FROM VIAMD\n");
                            md_file_printf(file, "REMARK    Exported from VIAMD\n");
                            
                            // Write CRYST1 record with unit cell if available
                            if (traj) {
                                // Coordinates are not needed, only the cell of the first exported
                                // frame, which is delivered on the state.
                                md_system_state_t frame_state = {};
                                int first_frame = md_array_size(frame_indices) > 0 ? frame_indices[0] : 0;
                                if (extract_frame(data, first_frame, &frame_state)) {
                                    if (frame_state.unitcell.flags != 0) {
                                        double a, b, c, alpha, beta, gamma;
                                        md_unitcell_extract_extent_angles(&a, &b, &c, &alpha, &beta, &gamma, &frame_state.unitcell);
                                        md_file_printf(file, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", a, b, c, alpha, beta, gamma);
                                    }
                                }
                            } else if (state->unitcell.flags != 0) {
                                double a, b, c, alpha, beta, gamma;
                                md_unitcell_extract_extent_angles(&a, &b, &c, &alpha, &beta, &gamma, &state->unitcell);
                                md_file_printf(file, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", a, b, c, alpha, beta, gamma);
                            }
                        } break;
                        default: {
                            VIAMD_LOG_ERROR("Unsupported export format: '" STR_FMT "'", STR_ARG(ext));
                        } break;
                    }

                    // Write body
                    md_array(int32_t) atom_indices = 0;
                    md_array_ensure(atom_indices, sys->atom.count, frame_alloc);

                    md_array(int) atomic_numbers = md_array_create(int, sys->atom.count, frame_alloc);
                    for (size_t i = 0; i < sys->atom.count; ++i) {
                        atomic_numbers[i] = md_atom_atomic_number(&sys->atom, i);
                    }

                    size_t num_frames = md_array_size(frame_indices);
                    if (num_frames > 0) {
						md_system_state_t frame_state = { .alloc = frame_alloc };
                        md_system_state_init(&frame_state, sys->atom.count);

                        // Many frames: one context for all of them, so the run's files stay open.
                        const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
                        md_system_extract_t* ex = md_system_extract_begin(sys, str_from_cstr(data->mold.run), paths, ARRAY_SIZE(paths), md_get_heap_allocator());
                        defer { md_system_extract_end(ex); };

                        for (size_t f = 0; f < num_frames; ++f) {
                            int frame_idx = frame_indices[f];
                            if (!ex || !md_system_extract_frame(ex, frame_idx, &frame_state)) {
                                VIAMD_LOG_ERROR("Failed to load frame %d from trajectory for structure export.", frame_idx);
                                continue;
                            }

                            size_t num_atoms = extract_atom_indices(atom_indices, md_array_capacity(atom_indices), &frame_state);
                            
                            if (struct_exp.selected_file_format == 0) { // XYZ
                                xyz_write_frame(file, atomic_numbers, &frame_state, atom_indices, num_atoms);
                            } else if (struct_exp.selected_file_format == 1) { // PDB
                                int model_num = (num_frames > 1) ? (int)(f + 1) : 0;
                                md_pdb_system_write_state_to_file(file, sys, &frame_state, atom_indices, num_atoms, model_num);
                            }
                        }
                    }
                    else {
                        // Single frame from system
                        size_t num_atoms = extract_atom_indices(atom_indices, md_array_capacity(atom_indices), &data->mold.state);
                        
                        if (struct_exp.selected_file_format == 0) { // XYZ
                            xyz_write_frame(file, atomic_numbers, state, atom_indices, num_atoms);
                        } else if (struct_exp.selected_file_format == 1) { // PDB
                            md_pdb_system_write_state_to_file(file, sys, state, atom_indices, num_atoms, 0);
                        }
                    }

                    VIAMD_LOG_INFO("Successfully exported structure to: '" STR_FMT "'", STR_ARG(path));
                    md_file_close(&file);
                } else {
                    VIAMD_LOG_ERROR("Failed to open file '" STR_FMT "' for writing.", STR_ARG(path));
                }
            }
        }
    }
    ImGui::End();
}

static void update_md_buffers(ApplicationState* data) {
    ASSERT(data);
    const md_system_t& sys = data->mold.sys;
	const md_system_state_t& state = data->mold.state;

    if (sys.atom.count == 0) return;

    if (data->mold.dirty_gpu_buffers) {
        data->mold.gpu_buffers_version += 1;
    }
    if (data->mold.dirty_gpu_buffers & ~MolBit_DirtyFlags) {
        data->mold.gpu_scene_version += 1;
    }

    if (data->mold.dirty_gpu_buffers & (MolBit_DirtyPosition | MolBit_ClearVelocity)) {
        // Derived from the positions and velocities of the atoms below, when a cartoon shows them
        data->mold.cartoon_glyphs.positions_valid = false;
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyPosition) {
        vec3_t pbc_ext = { 0 };
        md_unitcell_diag_extract_float(pbc_ext.elem, &state.unitcell);
        md_gl_mol_set_atom_position(data->mold.gl_mol, 0, (uint32_t)state.num_atoms, state.xyz);
        if (!(data->mold.dirty_gpu_buffers & MolBit_ClearVelocity)) {
            md_gl_mol_compute_velocity(data->mold.gl_mol, pbc_ext.elem);
        }
#if EXPERIMENTAL_GFX_API
        vec3_t aabb_min, aabb_max;
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, state.xyz, nullptr, nullptr, state.num_atoms);
        md_gfx_structure_set_atom_position(data->mold.gfx_structure, state.xyz, (uint32_t)state.num_atoms, 0);
        md_gfx_structure_set_aabb(data->mold.gfx_structure, &aabb_min, &aabb_max);
#endif
    }

    if (data->mold.dirty_gpu_buffers & MolBit_ClearVelocity) {
        md_gl_mol_zero_velocity(data->mold.gl_mol);
    }

    if (data->mold.dirty_gpu_buffers & MolBit_ResetBackboneHistory) {
        md_gl_mol_reset_backbone_history(data->mold.gl_mol);
        if (cartoon_glyphs_valid(data->mold.cartoon_glyphs)) {
            md_gl_mol_reset_backbone_history(data->mold.cartoon_glyphs.gl_mol);
        }
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyRadius) {
        md_temp_scope_t tmp = md_temp_begin_in(frame_alloc);
        defer { md_temp_end(tmp); };

        float* radii = (float*)md_vm_arena_push(frame_alloc, sys.atom.count * sizeof(float));
        md_atom_extract_radii(radii, 0, sys.atom.count, &sys.atom);

        md_gl_mol_set_atom_radius(data->mold.gl_mol, 0, (uint32_t)sys.atom.count, radii, 0);
#if EXPERIMENTAL_GFX_API
        md_gfx_structure_set_atom_radius(data->mold.gfx_structure, 0, (uint32_t)sys.atom.count, sys.atom.radius, 0);
#endif
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyFlags) {
        md_temp_scope_t tmp = md_temp_begin_in(frame_alloc);
        defer { md_temp_end(tmp); };

        uint8_t* flags = (uint8_t*)md_vm_arena_push(frame_alloc, sys.atom.count * sizeof(uint8_t));
        MEMSET(flags, 0, sys.atom.count * sizeof(uint8_t));

        {
            md_bitset_iter_t it = md_bitmask_iter(&data->selection.highlight_mask);
            while (md_bitset_iter_next(&it)) {
                uint64_t idx = it.idx;
                flags[idx] |= AtomBit_Highlighted | AtomBit_Marked;
            }
        }
        {
            // The colour of the group of the highlight each atom is in: of the last group it is in
            const auto& groups = data->selection.highlight_groups;
            for (size_t g = 0; g < groups.count; ++g) {
                const uint8_t color = (uint8_t)((g % HIGHLIGHT_GROUP_COLORS) << AtomBit_GroupShift);
                md_bitset_iter_t it = md_bitset_iter(groups.sets[g]);
                while (md_bitset_iter_next(&it)) {
                    if (it.idx < sys.atom.count) flags[it.idx] = (uint8_t)((flags[it.idx] & ~AtomBit_GroupMask) | color);
                }
            }
        }
        {
            md_bitset_iter_t it = md_bitmask_iter(&data->selection.selection_mask);
            while (md_bitset_iter_next(&it)) {
                uint64_t idx = it.idx;
                flags[idx] |= AtomBit_Selected | AtomBit_Marked;
            }
        }
        if (!md_bitmask_empty(&data->selection.selection_mask)) {
            for (size_t i = 0; i < sys.atom.count; ++i) {
                if (!(flags[i] & AtomBit_Selected)) flags[i] |= AtomBit_Dimmed;
            }
        }
        {
            md_bitset_iter_t it = md_bitmask_iter(&data->representation.visibility_mask);
            while (md_bitset_iter_next(&it)) {
                uint64_t idx = it.idx;
                flags[idx] |= AtomBit_Visible;
            }
        }
        md_gl_mol_set_atom_flags(data->mold.gl_mol, 0, (uint32_t)sys.atom.count, flags, 0);
        cartoon_glyphs_update_flags(data->mold.cartoon_glyphs, flags, sys.atom.count);
    }

    if (data->mold.dirty_gpu_buffers & MolBit_DirtyBonds) {
        // The bonds of the graph, then those formed in the frame shown while tracking them. A broken one is set to one
        // atom twice, which is not drawn, so that the bonds keep their indices: the orders and the picking go by them.
        const md_bond_delta_t* delta = bond_tracking_delta(data);
        if (delta && (delta->num_formed > 0 || !md_bitset_empty(delta->broken))) {
            const size_t count = sys.bond.count + delta->num_formed;
            md_atom_pair_t* pairs = (md_atom_pair_t*)md_vm_arena_push(data->allocator.frame, sizeof(md_atom_pair_t) * count);
            if (sys.bond.count) MEMCPY(pairs, sys.bond.pairs, sizeof(md_atom_pair_t) * sys.bond.count);
            md_bitset_iter_t it = md_bitset_iter(delta->broken);
            while (md_bitset_iter_next(&it)) {
                if (it.idx < sys.bond.count) pairs[it.idx].idx[1] = pairs[it.idx].idx[0];
            }
            if (delta->num_formed) MEMCPY(pairs + sys.bond.count, delta->formed_pair, sizeof(md_atom_pair_t) * delta->num_formed);
            md_gl_mol_set_bonds(data->mold.gl_mol, 0, (uint32_t)count, pairs, sizeof(md_atom_pair_t));
        } else {
            md_gl_mol_set_bonds(data->mold.gl_mol, 0, (uint32_t)sys.bond.count, sys.bond.pairs, sizeof(md_atom_pair_t));
        }
        // The orders of the frame shown where the bonds tracked changed them (bond_tracking_chem)
        md_gl_mol_set_bond_orders(data->mold.gl_mol, &sys, delta, bond_tracking_chem(data));

        // The non-covalent bonds of the system are drawn dashed by the licorice and ball and stick representations
        for (size_t i = 0; i < md_array_size(data->representation.reps); ++i) {
            representation_set_weak_bonds(data, &data->representation.reps[i]);
        }
    }

    // The glyphs of the cartoon of the lipids and nucleic acids are derived on the GPU from the positions and velocities
    // of the atoms, set above: only while a cartoon shows them
    CartoonGlyphs& glyphs = data->mold.cartoon_glyphs;
    if (cartoon_glyphs_valid(glyphs) && !glyphs.positions_valid) {
        bool shown = false;
        for (size_t i = 0; i < md_array_size(data->representation.reps) && !shown; ++i) {
            const Representation& rep = data->representation.reps[i];
            shown = representation_shown(*data, rep) && rep.type == RepresentationType::Cartoon && rep.glyph_rep.id && rep.glyphs_shown;
        }
        if (shown) {
            cartoon_glyphs_update_positions(&glyphs, state);
        }
    }

    data->mold.dirty_gpu_buffers = 0;
}

void create_screenshot(str_t path) {
    md_temp_scope_t tmp = md_temp_begin_in(frame_alloc);
    defer { md_temp_end(tmp); };

    int viewport[4] = {};
    int fbo = 0;
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &fbo);

    int width  = viewport[2];
    int height = viewport[3];

    size_t bytes = width * height * sizeof(uint32_t);
    uint32_t* rgba = (uint32_t*)md_vm_arena_push(frame_alloc, bytes);
    defer { md_vm_arena_pop(frame_alloc, bytes); };

    GLenum read_buffer = fbo ? GL_COLOR_ATTACHMENT0 : GL_BACK;

    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
    glReadBuffer(read_buffer);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, rgba);

    {
        // @NOTE: Swap Rows to flip image with respect to y-axis
        const uint32_t row_byte_size = width * sizeof(uint32_t);
        uint32_t* row_t = (uint32_t*)md_alloc(frame_alloc, row_byte_size);
        defer { md_free(frame_alloc, row_t, row_byte_size); };
        for (uint32_t i = 0; i < (uint32_t)height / 2; ++i) {
            uint32_t* row_a = rgba + i * width;
            uint32_t* row_b = rgba + (height - 1 - i) * width;
            if (row_a != row_b) {
                MEMCPY(row_t, row_a, row_byte_size);  // tmp = a;
                MEMCPY(row_a, row_b, row_byte_size);  // a = b;
                MEMCPY(row_b, row_t, row_byte_size);  // b = tmp;
            }
        }
    }

    str_t ext = {};
    extract_ext(&ext, path);

    bool written = false;
    switch (image_format_from_ext(ext)) {
    case ImageFormat_JPG: {
        const int quality = 95;
        written = image_write_jpg(path, rgba, width, height, quality);
        break;
    }
    case ImageFormat_PNG:
        written = image_write_png(path, rgba, width, height);
        break;
    case ImageFormat_BMP:
        written = image_write_bmp(path, rgba, width, height);
        break;
    default:
        VIAMD_LOG_ERROR("Non supported file-extension '" STR_FMT "' when saving screenshot", STR_ARG(ext));
        return;
    }

    if (!written) {
        VIAMD_LOG_ERROR("Failed to save screenshot to: '" STR_FMT "'", STR_ARG(path));
        return;
    }

    VIAMD_LOG_SUCCESS("Screenshot saved to: '" STR_FMT "'", STR_ARG(path));
}

static void draw_coordinate_system_widget_window(ViewTransform* target, const ViewTransform& current) {
    const int def_size = 150;
    const int min_size = 100;
    const int max_size = 500;

    ImGui::SetNextWindowSize(ImVec2(def_size, def_size), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(min_size, min_size), ImVec2(max_size, max_size), [](ImGuiSizeCallbackData* data) {
        // Enforce a square aspect ratio
        data->DesiredSize = ImVec2(ImMax(data->DesiredSize.x, data->DesiredSize.y), ImMax(data->DesiredSize.x, data->DesiredSize.y));
    });

    // Scale down backround color to only show a faint outline
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImGui::GetStyleColorVec4(ImGuiCol_WindowBg) * ImVec4(1.0f, 1.0f ,1.0f ,0.1f));
    defer { ImGui::PopStyleColor(); };

    const bool editable = ImGui::IsKeyDown(ImGuiMod_Alt);

    int window_flags = ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoCollapse;
    if (!editable) {
        window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoNavFocus;
    }

    if (ImGui::Begin("Coordinate Widget", nullptr, window_flags)) {
        quat_t out_orientation = target->orientation;
        if (ImGui::CoordinateSystemWidget(&out_orientation, current.orientation, ImGui::GetContentRegionAvail())) {
            const vec3_t look_at = camera_get_look_at(*target);
            target->orientation = quat_normalize(out_orientation);
            target->position = camera_position_from_look_at(look_at, target->orientation, target->distance);
        }
    }
    ImGui::End();
}

// Input that may change what the main view draws, though not in its key: a click or a drag anywhere (a button, a
// slider), the wheel, a key that is not typed into a text field (which edits a script or a filter, whose results are in
// the key). Hovering is not: what it changes is in the keys, the highlight in the marks'.
static bool input_activity() {
    const ImGuiIO& io = ImGui::GetIO();
    for (int i = 0; i < ImGuiMouseButton_COUNT; ++i) {
        if (io.MouseDown[i] || ImGui::IsMouseReleased(i)) return true;
    }
    if (io.MouseWheel != 0.0f || io.MouseWheelH != 0.0f) return true;
    if (!io.WantTextInput) {
        for (int k = ImGuiKey_NamedKey_BEGIN; k < ImGuiKey_NamedKey_END; ++k) {
            if (ImGui::IsKeyDown((ImGuiKey)k) || ImGui::IsKeyReleased((ImGuiKey)k)) return true;
        }
    }
    return false;
}

// The frames the image takes to settle after a change: the TAA's history to within a percent of the new image (by its
// feedback on a still view) and twice through the jitter, the shadows and the ambient occlusion likewise when they are
// filtered over time; else the one frame
static uint32_t settle_frames(const ApplicationState* state) {
    uint32_t n = 1;
    if (taa_active(state)) {
        const double feedback = CLAMP((double)state->visuals.temporal_aa.feedback_max, 0.5, 0.99);
        n = (uint32_t)ceil(log(0.01) / log(feedback));
        n = MAX(n, 2 * jitter_count(state));
    }
    if (occlusion_temporal_active(state)) {
        n = MAX(n, OCCLUSION_SETTLE_FRAMES);
    }
    return CLAMP(n, 1u, 2u * JITTER_SEQUENCE_MAX);
}

// What the image of the scene is made from (render_scene), but for the marks, and for what input changes (input_activity):
// while it is the same and the image has settled, the image is shown again rather than made again
static uint64_t scene_key(const ApplicationState* state) {
    uint64_t h = 0;
    const GBuffer& gbuf = state->gbuffer;
    const uint32_t size[4] = {gbuf.width, gbuf.height, gbuf.out_width, gbuf.out_height};
    h = md_hash64(size, sizeof(size), h);

    // The view
    const Camera& c = state->view.camera;
    const float camera[11] = {c.orientation.x, c.orientation.y, c.orientation.z, c.orientation.w, c.position.x, c.position.y, c.position.z,
                              c.distance, c.near_plane, c.far_plane, c.fov_y};
    h = md_hash64(camera, sizeof(camera), h);
    h = md_hash64(&state->view.mode, sizeof(state->view.mode), h);
    const vec2_t shift = view_shift(state);
    h = md_hash64(&shift, sizeof(shift), h);
    h = md_hash64(&state->mold.unitcell_transform, sizeof(state->mold.unitcell_transform), h);

    // How it looks
    h = md_hash64(&state->visuals, sizeof(state->visuals), h);
    h = md_hash64(&state->settings, sizeof(state->settings), h);

    // What it shows: the representations and the atoms (their visibility, the dimmed by the selection: not the
    // highlighted), the volumes, what the script shows, what the components draw
    h = md_hash64_combine(h, opaque_scene_hash(state));
    h = md_hash64(&state->render.selection_hash, sizeof(state->render.selection_hash), h);
    const uint64_t volumes = volume::data_version_any();
    h = md_hash64(&volumes, sizeof(volumes), h);
    const uint64_t vis = script_vis_shown_key(state);
    h = md_hash64(&vis, sizeof(vis), h);
    const vec4_t vis_colors[3] = {state->script.point_color, state->script.line_color, state->script.triangle_color};
    h = md_hash64(vis_colors, sizeof(vis_colors), h);
    h = md_hash64(&state->render.version, sizeof(state->render.version), h);
    h = md_hash64(&use_gfx, sizeof(use_gfx), h);

    RenderKey rk = {state, h};
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderKey, viamd::EventPayloadType_RenderKey, &rk);
    return rk.key ? rk.key : 1;
}

// The image of the scene: the G-buffer, made into the image kept (gbuffer.tex.image) by post-processing, which
// render() presents
static void render_scene(ApplicationState* state, bool do_screenshot) {
    update_view_param(state);

    gbuffer_clear(&state->gbuffer);

    const GLenum draw_buffers_opaque[] = {GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
        GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // Enable all draw buffers
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, state->gbuffer.fbo);
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers_opaque), draw_buffers_opaque);

    PUSH_GPU_SECTION("G-Buffer fill")

    {
        immediate::Scope vis_scope(state->gfx.overlay, "visualization");
        immediate::Scope vis_scope_depth(state->gfx.world, "visualization with depth");

        static const md_script_vis_t nothing = {};
        const md_script_vis_t& vis = state->script.vis_shown ? *state->script.vis_shown : nothing;

        if (vis.points) {
            immediate::points(vis_scope, (immediate::Vertex*)vis.points, md_array_size(vis.points), state->script.point_color);
        }

        if (vis.triangles) {
            immediate::triangles(vis_scope, (immediate::Vertex*)vis.triangles, md_array_size(vis.triangles), state->script.triangle_color);
        }

        if (vis.lines) {
            immediate::lines(vis_scope, (immediate::Vertex*)vis.lines, md_array_size(vis.lines), state->script.line_color);
        }

        const size_t num_matrices = md_array_size(vis.sdf.matrices);
        const size_t num_structures = md_array_size(vis.sdf.structures);
        const size_t num_sdf_items = MIN(MIN(num_matrices, num_structures), 100);

        if (num_matrices != num_structures) {
            VIAMD_LOG_DEBUG("SDF visualization returned mismatched matrix and structure counts (%zu, %zu)", num_matrices, num_structures);
        }
        if (num_sdf_items > 0) {
            const vec4_t col_x = { 1, 0, 0, 0.7f };
            const vec4_t col_y = { 0, 1, 0, 0.7f };
            const vec4_t col_z = { 0, 0, 1, 0.7f };
            const float ext = vis.sdf.extent * 0.25f;
            const vec3_t box_ext = vec3_set1(vis.sdf.extent);
            for (size_t i = 0; i < num_sdf_items; ++i) {
                const mat4_t model_matrix = mat4_inverse(vis.sdf.matrices[i]);
                immediate::basis(vis_scope, model_matrix, ext, col_x, col_y, col_z);
                immediate::box_wireframe(vis_scope_depth, -box_ext, box_ext, model_matrix);
            }
        }
    }

    if (!use_gfx) {
        // DRAW VELOCITY OF STATIC OBJECTS
        PUSH_GPU_SECTION("Blit Static Velocity")
        glDrawBuffer(GL_COLOR_ATTACHMENT_VELOCITY);
        glDepthMask(0);
        postprocessing::blit_static_velocity(state->gbuffer.tex.depth, state->view.param);
        glDepthMask(1);
        POP_GPU_SECTION()
    }
    glDepthMask(1);
    glColorMask(1, 1, 1, 1);

    // DRAW REPRESENTATIONS
    PUSH_GPU_SECTION("Draw Opaque")
    glDrawBuffers((int)ARRAY_SIZE(draw_buffers_opaque), draw_buffers_opaque);
    draw_representations_opaque(state);
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderOpaque, viamd::EventPayloadType_ApplicationState, state);
    POP_GPU_SECTION()

    immediate::RenderParams params = {};
    params.view = state->view.param.matrix.curr.view;
    params.proj = state->view.param.matrix.curr.proj;

    // standard alpha blending
    PUSH_GPU_SECTION("Immediate world")
    //glEnable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    immediate::render(state->gfx.world, params);
    //glDisable(GL_BLEND);
    POP_GPU_SECTION()

    // The transparent representations test the depth of the opaque ones, but do not write it
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(0);
    glColorMask(1, 1, 1, 1);
    glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);

    const bool transparency_hdr_written = draw_representations_transparent(state);
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRenderTransparent, viamd::EventPayloadType_ApplicationState, state);

    const GLenum draw_buffers_transparent[] = { GL_COLOR_ATTACHMENT_TRANSPARENCY, 0, GL_COLOR_ATTACHMENT_VELOCITY, GL_COLOR_ATTACHMENT_PICKING };
	glDrawBuffers(ARRAY_SIZE(draw_buffers_transparent), draw_buffers_transparent);

    glDisable(GL_CULL_FACE);

    // Blend colour only. glEnable(GL_BLEND) applies to every attachment, and
    // this pass also writes velocity and picking: alpha-blending those means a
    // translucent overlay fragment scales the velocity underneath it by
    // (1 - alpha) and mixes the picking index with whatever was already there,
    // neither of which is a meaningful operation on that data.
    glEnablei(GL_BLEND, 0);
    // The transparency buffer holds premultiplied colour: straight alpha sources, accumulated premultiplied
    glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDisablei(GL_BLEND, 1);
    glDisablei(GL_BLEND, 2);
    glDisablei(GL_BLEND, 3);

    PUSH_GPU_SECTION("Immediate overlay")
    glDisable(GL_DEPTH_TEST);
    immediate::render(state->gfx.overlay, params);
    POP_GPU_SECTION()

    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);

    POP_GPU_SECTION()  // G-buffer

    PUSH_GPU_SECTION("Postprocessing")
    postprocess_pipeline::Settings settings = {};
    postprocess_pipeline::Inputs inputs = {};

    settings.background_color = state->visuals.background.color;

    settings.shading.mode = state->visuals.shading.mode;
    settings.shading.ambient = state->visuals.shading.ambient;
    settings.shading.specular = state->visuals.shading.specular;
    settings.shading.colors = state->visuals.shading.colors;
    settings.shading.light = state->visuals.shading.light;
    settings.shading.environment = state->visuals.shading.environment;

    settings.shadows.enabled = state->visuals.shadows.enabled;
    settings.shadows.softness = state->visuals.shadows.softness;
    settings.shadows.temporal = state->visuals.occlusion_temporal;

    settings.ssao.enabled = state->visuals.ssao.enabled;
    settings.ssao.intensity = state->visuals.ssao.intensity;
    settings.ssao.temporal = state->visuals.occlusion_temporal;
    settings.ssao.mode = state->visuals.ssao.mode;

    settings.tonemap.enabled = state->visuals.tonemapping.enabled;
    settings.tonemap.mode = state->visuals.tonemapping.tonemapper;
    settings.tonemap.exposure = state->visuals.tonemapping.exposure;
    settings.tonemap.gamma = state->visuals.tonemapping.gamma;

    settings.dof.enabled = state->visuals.dof.enabled;
    settings.dof.focus_depth = state->visuals.dof.focus_depth;
    settings.dof.aperture = state->visuals.dof.aperture * 0.01f;

    settings.fxaa.enabled = fxaa_active(state);

    constexpr float MOTION_BLUR_REFERENCE_DT = 1.0f / 60.0f;
    const float dt_compensation = MOTION_BLUR_REFERENCE_DT / (float)state->app.timing.delta_s;
    const float motion_scale = state->visuals.temporal_aa.motion_blur.motion_scale * dt_compensation;
    settings.taa.enabled = taa_active(state);
    settings.taa.feedback_min = state->visuals.temporal_aa.feedback_min;
    settings.taa.feedback_max = state->visuals.temporal_aa.feedback_max;
    settings.taa.motion_blur.enabled = state->visuals.temporal_aa.motion_blur.enabled;
    settings.taa.motion_blur.motion_scale = motion_scale;

    settings.sharpen.enabled = taa_active(state) && state->visuals.sharpen.enabled;
    settings.sharpen.weight = state->visuals.sharpen.weight;

    // Upscaled from the G-buffer's size to the viewport's (a screenshot is made at its size)
    if (!do_screenshot) {
        inputs.width  = (int)state->gbuffer.width;
        inputs.height = (int)state->gbuffer.height;
    }
    inputs.depth = state->gbuffer.tex.depth;
    inputs.color = state->gbuffer.tex.color;
    inputs.normal = state->gbuffer.tex.normal;
    inputs.velocity = state->gbuffer.tex.velocity;
    inputs.transparency = state->gbuffer.tex.transparency;
    inputs.transparency_hdr = transparency_hdr_written ? state->gbuffer.tex.transparency_hdr : 0;
    inputs.history = settings.taa.enabled ? state->gbuffer.tex.history : 0;
    inputs.history_prev = settings.taa.enabled ? state->gbuffer.tex.history_prev : 0;

    // Kept, and presented after (render())
    settings.present = false;
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, state->gbuffer.fbo_image);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, state->gbuffer.out_width, state->gbuffer.out_height);

    postprocess_pipeline::execute(inputs, settings, state->view.param);
    POP_GPU_SECTION()

    if (settings.taa.enabled) {
        // Ping-pong: what was written this frame is read as the previous history next frame (no copy needed)
        uint32_t tmp = state->gbuffer.tex.history;
        state->gbuffer.tex.history = state->gbuffer.tex.history_prev;
        state->gbuffer.tex.history_prev = tmp;
    }
}

// Saves the screenshot from the framebuffer bound (create_screenshot), and it is no longer pending
static void screenshot_save(ApplicationState* state) {
    create_screenshot(state->screenshot.path_to_file);
    str_free(state->screenshot.path_to_file, state->allocator.persistent);
    state->screenshot.path_to_file = {};
}

static void render(ApplicationState* state) {
    // A screenshot pending (screenshot_pending): saved below once the image has settled, without the GUI at its own size
    const bool screenshot = screenshot_pending(state);
    const bool offscreen  = screenshot_offscreen(state);

    // The image shown, and the G-buffer it is made from: the render scale of it, upscaled in post-processing.
    // A screenshot without the GUI is at its own resolution, in full.
    uint32_t out_width  = state->app.framebuffer.width;
    uint32_t out_height = state->app.framebuffer.height;
    uint32_t gbuffer_target_width  = out_width;
    uint32_t gbuffer_target_height = out_height;
    if (offscreen) {
        out_width  = gbuffer_target_width  = state->screenshot.res_x;
        out_height = gbuffer_target_height = state->screenshot.res_y;
    } else {
        render_size(state, out_width, out_height, &gbuffer_target_width, &gbuffer_target_height);
    }

    // Resize Framebuffer
    if ((state->gbuffer.width != gbuffer_target_width || state->gbuffer.height != gbuffer_target_height ||
         state->gbuffer.out_width != out_width || state->gbuffer.out_height != out_height) &&
        (gbuffer_target_width != 0 && gbuffer_target_height != 0 && out_width != 0 && out_height != 0)) {
        gbuffer_init(&state->gbuffer, gbuffer_target_width, gbuffer_target_height, out_width, out_height);
        postprocess_pipeline::initialize(state->gbuffer.width, state->gbuffer.height, out_width, out_height);
        state->render.image_valid = false;
    }

    // The scene is made when something in it changed, until it settles (and the isosurface meshes' mean has converged);
    // else the image kept is shown again. Always without lazy rendering. Input may change what the key does not show, so
    // it starts the settling again, but not while a screenshot waits for it.
    auto& lazy = state->render;
    const uint64_t key = scene_key(state);
    if (key != lazy.scene_key || !lazy.image_valid || (!screenshot && input_activity())) {
        lazy.scene_key = key;
        lazy.settle_frames = settle_frames(state);
    }
    const bool scene = !state->settings.lazy_rendering || lazy.settle_frames > 0 || !lazy.image_valid || !iso_mesh_converged;
    lazy.history = (lazy.history << 1) | (scene ? 1 : 0);
    if (scene) {
        PUSH_GPU_SECTION("Scene")
        render_scene(state, offscreen);
        POP_GPU_SECTION()
        lazy.image_valid = true;
        if (lazy.settle_frames > 0) lazy.settle_frames -= 1;
    }
    // Settled: as made as it gets. The meshes' mean that the profiler keeps drawing never converges: then its frames.
    const bool meshes_settled = iso_mesh_converged || iso_mesh_target.frame >= isomesh::MAX_FRAMES;
    const bool settled = lazy.image_valid && lazy.settle_frames == 0 && meshes_settled;

    // The marks of the selection and the highlight: the selected and highlighted atoms drawn once, both, into a target
    // of the size of the image made, unjittered, each fragment marked in view or behind the scene (its depth, from the
    // G-buffer); put over the image as it is presented. They are not part of the scene: changing them changes nothing in
    // it, nor in its history. Drawn again when they change, or the scene was made.
    const bool selection_empty = md_bitmask_empty(&state->selection.selection_mask);
    const bool highlight_empty = md_bitmask_empty(&state->selection.highlight_mask);
    const bool draw_marks = !use_gfx && (!selection_empty || !highlight_empty);
    if (draw_marks) {
        uint64_t marks_key = md_hash64(&lazy.selection_hash, sizeof(lazy.selection_hash), lazy.scene_key);
        marks_key = md_hash64(&lazy.highlight_hash, sizeof(lazy.highlight_hash), marks_key);
        marks_key = md_hash64(&state->selection.highlight_groups.key, sizeof(state->selection.highlight_groups.key), marks_key);
        marks_key = marks_key ? marks_key : 1;
        if (scene || marks_key != lazy.marks_key) {
            PUSH_GPU_SECTION("Marks")
            gbuffer_marks_begin(&state->gbuffer);
            draw_representations_marks(state);
            gbuffer_marks_end(&state->gbuffer);
            POP_GPU_SECTION()
            lazy.marks_key = marks_key;
        }
    } else {
        lazy.marks_key = 0;
    }

    // The image kept, the marks over it, into the framebuffer bound
    auto present = [&]() {
        PUSH_GPU_SECTION("Present")
        postprocess_pipeline::Settings settings = {};
        if (draw_marks) {
            // The highlight fades in, the outline is as wide on a display of more pixels per point
            constexpr double HIGHLIGHT_FADE_IN_S = 0.12;
            const float t = (float)CLAMP((ImGui::GetTime() - state->selection.highlight_time) / HIGHLIGHT_FADE_IN_S, 0.0, 1.0);
            const int outline_step = MAX(1, (int)roundf(ImGui::GetIO().DisplayFramebufferScale.x));
            marks_settings(&settings, *state, t * t * (3.0f - 2.0f * t), outline_step);
        }
        postprocess_pipeline::present(state->gbuffer.tex.image, draw_marks ? state->gbuffer.tex.marks : 0, draw_marks ? state->gbuffer.tex.marks_group : 0, settings);
        POP_GPU_SECTION()
    };

    // Without the GUI: once settled, the image presented into the G-buffer's target (its size) and saved from there
    if (offscreen && settled) {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, state->gbuffer.fbo);
        glViewport(0, 0, state->gbuffer.width, state->gbuffer.height);
        glDrawBuffer(GL_COLOR_ATTACHMENT0);
        present();
        screenshot_save(state);
    }

    // The window: the image, or while a screenshot without the GUI is made (at a size of its own), how far it has come
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glViewport(0, 0, state->app.framebuffer.width, state->app.framebuffer.height);
    glDrawBuffer(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT);
    if (!offscreen) {
        present();
    } else {
        const uint32_t total = MAX(settle_frames(state), 1u);
        float fraction = screenshot_pending(state) ? 1.0f - (float)lazy.settle_frames / (float)total : 1.0f;
        if (!meshes_settled) fraction = MIN(fraction, (float)iso_mesh_target.frame / (float)isomesh::MAX_FRAMES);
        if (screenshot_pending(state)) ImGui::OpenPopup("Capturing Frame");
        ImVec2 size = ImGui::CalcTextSize("Capturing Frame");
        ImGui::SetNextWindowSize(size * ImVec2(2, 5), ImGuiCond_Always);
        if (ImGui::BeginPopupModal("Capturing Frame", 0, ImGuiWindowFlags_NoResize)) {
            ImGui::ProgressBar(CLAMP(fraction, 0.0f, 1.0f));
            if (!screenshot_pending(state)) ImGui::CloseCurrentPopup();     // saved this frame
            ImGui::EndPopup();
        }
    }

    PUSH_GPU_SECTION("Imgui render")
    application::render_imgui(&state->app);
    POP_GPU_SECTION()

    // With the GUI: the window as shown, once settled
    if (screenshot && !offscreen && settled) {
        screenshot_save(state);
    }
}

// The longest bond md_gl draws, covalent or weak (md_gl_draw_args_t max_bond_length): a longer one is taken to be across
// a periodic boundary, from an atom to the far side of the cell, and is left out.
//   In a cell: half its smallest width (the minimum image), at least 3 Å.
//   Without a cell the boundary is not known, but the coordinates may still be wrapped (a periodic trajectory without
//   its cell). Then: BOND_LENGTH_MARGIN times the longest bond of the system in the frame where its topology last
//   changed (where its bonds were made, read or edited), at least BOND_LENGTH_MIN for the weak bonds, which are not of
//   the topology. A bond stays near the length it had, one wrapped across a cell since is about as long as the cell.
//   A bond already wrapped in that frame cannot be told from a long one, and lifts the limit with it.
static constexpr float BOND_LENGTH_MARGIN = 1.5f;
static constexpr float BOND_LENGTH_MIN    = 5.0f;

static float bond_length_limit(ApplicationState* state) {
    const md_unitcell_t& cell = state->mold.state.unitcell;
    if (md_unitcell_flags(&cell) != 0) {
        vec3_t ext = {0};
        md_unitcell_diag_extract_float(ext.elem, &cell);
        return MAX(3.0f, vec3_reduce_min(ext) * 0.5f);
    }

    const md_system_t& sys = state->mold.sys;
    if (state->mold.bond_length_limit_version != sys.topology_version || !(state->mold.bond_length_limit > 0.0f)) {
        const md_system_state_t& frame = state->mold.state;
        float max_d2 = 0.0f;
        if (frame.xyz && frame.num_atoms == sys.atom.count) {
            for (size_t i = 0; i < sys.bond.count; ++i) {
                const md_atom_pair_t pair = sys.bond.pairs[i];
                if ((size_t)pair.idx[0] < frame.num_atoms && (size_t)pair.idx[1] < frame.num_atoms) {
                    max_d2 = MAX(max_d2, vec3_distance_squared(frame.xyz[pair.idx[0]], frame.xyz[pair.idx[1]]));
                }
            }
        }
        state->mold.bond_length_limit = MAX(BOND_LENGTH_MIN, BOND_LENGTH_MARGIN * sqrtf(max_d2));
        state->mold.bond_length_limit_version = sys.topology_version;
    }
    return state->mold.bond_length_limit;
}

// The copies of the scene the periodic images show (view.periodic_images): its own, as it is, then the translations of
// the scene into the neighboring cells, through the turn the coordinates carry if the orientation is kept (as the box).
// None without a cell or images.
static void compute_periodic_transforms(md_array(mat4_t)* transforms, const ApplicationState& state, md_allocator_i* alloc) {
    md_array_shrink(*transforms, 0);
    const md_unitcell_t& cell = state.mold.state.unitcell;
    const int* n = state.view.periodic_images;
    if (md_unitcell_flags(&cell) == 0 || (n[0] | n[1] | n[2]) == 0 || state.view.assembly >= 0) return;

    mat3_t A = {0};
    md_unitcell_A_extract_float(A.elem, &cell);
    const mat4_t M = state.operations.state_rotation * mat4_from_mat3(A);
    md_array_push(*transforms, mat4_ident(), alloc);
    for (int k = -n[2]; k <= n[2]; ++k) {
        for (int j = -n[1]; j <= n[1]; ++j) {
            for (int i = -n[0]; i <= n[0]; ++i) {
                if (i == 0 && j == 0 && k == 0) continue;
                const vec3_t t = vec3_from_vec4(mat4_mul_vec4(M, vec4_set((float)i, (float)j, (float)k, 0.0f)));
                md_array_push(*transforms, mat4_translate(t.x, t.y, t.z), alloc);
            }
        }
    }
}

// The copies a representation is drawn in without an assembly: the periodic images if it has them, by instancing.
// md_gl culls those out of view on the CPU. Read from the view state, which outlives the draw.
static md_gl_copies_t rep_copies(const ApplicationState& state, const Representation& rep) {
    md_gl_copies_t copies = {};
    const size_t n = md_array_size(state.view.periodic_transforms);
    if (rep.periodic_images && n > 1) {
        copies.transforms = &state.view.periodic_transforms[0].elem[0][0];
        copies.prev_transforms = md_array_size(state.view.periodic_transforms_prev) == n ? &state.view.periodic_transforms_prev[0].elem[0][0] : NULL;
        copies.count = (uint32_t)n;
    }
    return copies;
}

// The ops from beg on, of a representation, in the copies of the scene. With an assembly shown, each op once per part
// of it: of the instances of the part (its subset, of the molecule or of the glyphs of the cartoon, assembly_subsets_update)
// in the copies its operators make. Otherwise once, in the periodic images of the representation. Read from the system
// and the view state, which outlive the draw.
static void place_ops(md_array(md_gl_draw_op_t)* ops, size_t beg, const ApplicationState& state, const Representation& rep, bool glyphs, md_allocator_i* alloc) {
    const size_t end = md_array_size(*ops);
    const md_assembly_data_t& a = state.mold.sys.assembly;
    const bool shown = state.view.assembly >= 0 && (size_t)state.view.assembly < a.count;
    const size_t part_beg = shown ? a.group_offset[state.view.assembly] : 0;
    const size_t part_end = shown ? a.group_offset[state.view.assembly + 1] : 0;
    const size_t num_parts = md_array_size(state.view.assembly_copy_offset) > 0 ? md_array_size(state.view.assembly_copy_offset) - 1 : 0;
    const md_gl_subset_t* subsets = glyphs ? state.view.assembly_glyph_subsets : state.view.assembly_subsets;
    // Without an assembly, or its copies or subsets not of it (yet)
    if (num_parts == 0 || num_parts != part_end - part_beg || state.view.assembly_subsets_of != state.view.assembly || md_array_size(subsets) != num_parts) {
        const md_gl_copies_t copies = rep_copies(state, rep);
        for (size_t i = beg; i < end; ++i) (*ops)[i].copies = copies;
        return;
    }
    for (size_t p = 0; p < num_parts; ++p) {
        const uint32_t c0 = state.view.assembly_copy_offset[p];
        const uint32_t c1 = state.view.assembly_copy_offset[p + 1];
        for (size_t i = beg; i < end; ++i) {
            md_gl_draw_op_t op = (*ops)[i];
            op.copies = {&state.view.assembly_transforms[c0].elem[0][0], NULL, c1 - c0};
            op.subset = subsets[p];
            if (p == 0) {
                (*ops)[i] = op;
            } else {
                md_array_push(*ops, op, alloc);
            }
        }
    }
}

// The edges of the unit cell of the frame shown (a SimulationBox representation), through the turn the coordinates
// carry if the orientation is kept (the cell itself cannot hold one). Capsules of the thickness, shaded as the rest of
// the scene (the world queue writes the G-buffer), never thinner than a pixel where they are, so a box seen from afar
// does not break up; or lines of a pixel. With axes, the edges along a, b and c from the origin red, green and blue.
// The cells out of view are left out.
static void draw_simulation_box(ApplicationState* state, const Representation& rep) {
    const md_unitcell_t& cell = state->mold.state.unitcell;
    if (md_unitcell_flags(&cell) == 0) return;

    mat3_t A = {0};
    md_unitcell_A_extract_float(A.elem, &cell);
    const mat4_t M = state->operations.state_rotation * mat4_from_mat3(A);
    const uint32_t color = convert_color(vec4_set(rep.box.color.x, rep.box.color.y, rep.box.color.z, 1.0f));
    const uint32_t axis_color[3] = { convert_color(vec4_set(0.85f, 0.15f, 0.15f, 1.0f)), convert_color(vec4_set(0.15f, 0.70f, 0.15f, 1.0f)), convert_color(vec4_set(0.15f, 0.30f, 0.90f, 1.0f)) };

    immediate::Scope scope(state->gfx.world, "simulation box");

    // The corners: bit k of i is the coordinate along cell vector k
    vec3_t corner[8];
    for (int i = 0; i < 8; ++i) {
        corner[i] = vec3_from_vec4(mat4_mul_vec4(M, vec4_set((float)(i & 1), (float)((i >> 1) & 1), (float)((i >> 2) & 1), 1.0f)));
    }

    // A pixel at the depth d of view space: 2 d / (P11 h) in a perspective, 2 / (P11 h) in an orthographic projection
    const mat4_t& V = state->view.param.matrix.curr.view;
    const mat4_t& P = state->view.param.matrix.curr.proj;
    const float   h = (float)MAX(state->gbuffer.out_height, 1U);    // the pixels shown, whatever the render scale
    const bool    ortho = P.elem[3][3] == 1.0f;
    const float   px_per_depth = 2.0f / (P.elem[1][1] * h);

    // The planes of the view volume, inward, and the sphere of a cell about its center: a cell wholly outside of one is
    // out of view. The capsules are at most a few pixels thick, a margin of a hundredth of the cell covers them.
    const mat4_t PV = P * V;
    vec4_t plane[6];
    for (int p = 0; p < 6; ++p) {
        const int r = p / 2;
        const float s = (p & 1) ? -1.0f : 1.0f;
        const vec4_t row_w = {PV.elem[0][3], PV.elem[1][3], PV.elem[2][3], PV.elem[3][3]};
        const vec4_t row_r = {PV.elem[0][r], PV.elem[1][r], PV.elem[2][r], PV.elem[3][r]};
        plane[p] = vec4_add(row_w, vec4_mul1(row_r, s));
        const float len = vec3_length(vec3_from_vec4(plane[p]));
        plane[p] = len > 0.0f ? vec4_mul1(plane[p], 1.0f / len) : vec4_t{0, 0, 0, 1};
    }
    const vec3_t cell_center = (corner[0] + corner[7]) * 0.5f;
    float cell_radius = 0.0f;
    for (int i = 0; i < 8; ++i) cell_radius = MAX(cell_radius, vec3_length(corner[i] - cell_center));
    cell_radius = cell_radius * 1.01f + 0.5f * rep.box.thickness;

    // Its own cell, then those of the periodic images (whose shared edges are drawn twice); the axes in its own only
    const md_array(mat4_t) transforms = state->view.periodic_transforms;
    const size_t num_cells = rep.periodic_images ? MAX(md_array_size(transforms), (size_t)1) : 1;
    for (size_t c = 0; c < num_cells; ++c) {
        const vec3_t t = c == 0 ? vec3_t{0, 0, 0} : vec3_from_vec4(transforms[c].col[3]);
        const vec4_t center = vec4_from_vec3(cell_center + t, 1.0f);
        bool out_of_view = false;
        for (int p = 0; p < 6 && !out_of_view; ++p) {
            out_of_view = vec4_dot(plane[p], center) < -cell_radius;
        }
        if (out_of_view) continue;
        for (int i = 0; i < 8; ++i) {
            for (int k = 0; k < 3; ++k) {
                const int bit = 1 << k;
                if (i & bit) continue;
                const int j = i | bit;
                const uint32_t col = (rep.box.axes && c == 0 && i == 0) ? axis_color[k] : color;
                const vec3_t p0 = corner[i] + t;
                const vec3_t p1 = corner[j] + t;
                if (rep.box.thickness > 0.0f) {
                    const vec3_t mid  = (p0 + p1) * 0.5f;
                    const float depth = ortho ? 1.0f : MAX(-vec3_from_vec4(mat4_mul_vec4(V, vec4_from_vec3(mid, 1.0f))).z, 0.0f);
                    const float radius = MAX(0.5f * rep.box.thickness, 0.6f * px_per_depth * depth);
                    immediate::capsule(scope, p0, p1, radius, col);
                } else {
                    immediate::line(scope, p0, p1, col);
                }
            }
        }
    }
}

static void draw_representations_opaque(ApplicationState* state) {
    ASSERT(state);

    if (state->mold.sys.atom.count == 0) {
        return;
    }

#if EXPERIMENTAL_GFX_API
    if (use_gfx) {
        const uint32_t instance_count = 10000;
        static mat4_t* transforms = 0;

        if (transforms == 0) {
            auto rnd = []() -> float {
                return (float)rand() / RAND_MAX;
            };
            for (uint32_t i = 0; i < instance_count; ++i) {
                vec3_t axis = {rnd(), rnd(), rnd()};
                quat_t ori = quat_angle_axis(rnd() * TWO_PI, vec3_normalize(axis));
                mat4_t R = mat4_from_quat(ori);
                mat4_t T = mat4_translate(rnd() * 4000, rnd() * 4000, rnd() * 4000);
                mat4_t M = T * R;
                md_urange_t range = {0, (int32_t)state->mold.sys.atom.count};
                md_array_push(transforms, M, persistent_alloc);
            }
        }

        md_gfx_draw_op_t* draw_ops = 0;
        for (int64_t i = 0; i < md_array_size(state->representation.reps); ++i) {
            if (representation_shown(*state, state->representation.reps[i])) {
                md_gfx_draw_op_t op;
                op.structure = state->mold.gfx_structure;
                op.representation = state->representation.reps[i].gfx_rep;
                op.model_mat = NULL;
                md_array_push(draw_ops, op, frame_alloc);
                
                for (uint32_t j = 0; j < instance_count; ++j) {
                    md_gfx_draw_op_t op;
                    op.structure = state->mold.gfx_structure;
                    op.representation = state->representation.reps[i].gfx_rep;
                    op.model_mat = &transforms[j];
                    md_array_push(draw_ops, op, frame_alloc);
                }
                
            }
        }

        md_gfx_draw((uint32_t)md_array_size(draw_ops), draw_ops, &state->view.param.matrix.curr.proj, &state->view.param.matrix.curr.view, &state->view.param.matrix.inv.proj, &state->view.param.matrix.inv.view);
    } else {
#endif
        const size_t num_representations = md_array_size(state->representation.reps);
        if (num_representations == 0) return;

		glEnable(GL_DEPTH_TEST);

        md_array(md_gl_draw_op_t) draw_ops = 0;
        // The glyphs of the cartoon of the lipids and nucleic acids, a molecule of their own with picking of their own
        md_array(md_gl_draw_op_t) glyph_ops = 0;
        for (size_t i = 0; i < num_representations; ++i) {
            const Representation& rep = state->representation.reps[i];

            if (rep.type == RepresentationType::Cartoon && representation_shown(*state, rep) && rep.type_is_valid) {
                const size_t beg = md_array_size(glyph_ops);
                representation_glyph_draw_ops(&glyph_ops, *state, rep, frame_alloc);
                place_ops(&glyph_ops, beg, *state, rep, true, frame_alloc);
                // The cartoon of the proteins needs a protein backbone to draw
                if (state->mold.sys.protein_backbone.range.count == 0) continue;
            }

            if (RepresentationType::SpaceFill <= rep.type && rep.type <= RepresentationType::Cartoon) {
                if (representation_shown(*state, rep) && rep.type_is_valid) {
                    const md_gl_draw_op_t op = representation_draw_op(rep);
                    const size_t beg = md_array_size(draw_ops);
                    md_array_push(draw_ops, op, frame_alloc);
                    // The weak bonds it draws over its sticks (the non-covalent bonds of the system), in its copies
                    representation_weak_bond_draw_ops(&draw_ops, *state, rep, frame_alloc);
                    place_ops(&draw_ops, beg, *state, rep, false, frame_alloc);
                }
            } else if (rep.type == RepresentationType::Bonds) {
                if (representation_shown(*state, rep) && rep.type_is_valid) {
                    const size_t beg = md_array_size(draw_ops);
                    representation_weak_bond_draw_ops(&draw_ops, *state, rep, frame_alloc);
                    place_ops(&draw_ops, beg, *state, rep, false, frame_alloc);
                }
            } else if (rep.type == RepresentationType::DipoleMoment) {
                if (representation_shown(*state, rep)) {
                    // immediate draw of dipole moment as arrow
                    vec3_t dipole_vec = { 0, 0, 0 };
                    vec3_t dipole_org = { 0, 0, 0 };
                    if (dipole_moment_read(&dipole_vec, &dipole_org, state->mold.sys, rep.dipole.dipole_key, rep.dipole.dipole_index)) {
                        // The representation's own (key, index) IS the picking address: the group's
                        // range was reserved under that key this frame, and the shader adds the base to
                        // the per primitive index, so the element index is what goes on the primitive.
                        // A group that did not fit in the picking space draws with INVALID_PICKING_IDX,
                        // which the shader passes through untouched - visible, just not pickable.
                        const PickingSpace* space = picking_handler_current_space(&state->picking_handler);
                        const PickingRange* range = space ? picking_space_find_range(*space, PickingDomain_Dipole, rep.dipole.dipole_key) : nullptr;

                        immediate::Scope scope(state->gfx.world, "debug_dipole_moment");
                        immediate::set_picking_base_idx(scope, range ? range->beg : 0);

                        const vec3_t org = dipole_org;
                        const vec3_t vec = dipole_vec * (float)rep.dipole.scale;

                        // cylinder body
                        const float body_scale = 0.8f;

                        const float body_radius = rep.dipole.radius;
                        const float head_radius = body_radius * 1.5f;

                        vec3_t cyl_beg = rep.dipole.offset + org;
                        vec3_t cyl_end = rep.dipole.offset + org + vec * body_scale;
                        vec3_t arrow_end = rep.dipole.offset + org + vec;

                        uint32_t color_u32 = convert_color(rep.dipole.color);

                        uint32_t picking_idx = range ? rep.dipole.dipole_index : INVALID_PICKING_IDX;

                        immediate::cylinder(scope, cyl_beg, cyl_end, body_radius, color_u32, picking_idx);
                        immediate::cone(scope, cyl_end, arrow_end, head_radius, color_u32, picking_idx);
                    }
                }
            } else if (rep.type == RepresentationType::SimulationBox) {
                if (representation_shown(*state, rep) && rep.type_is_valid) {
                    draw_simulation_box(state, rep);
                }
            }
        }

        const float max_bond_length = bond_length_limit(state);

        md_gl_draw_args_t args = {
            .shaders = state->gl.shaders,
            .draw_operations = {
                .count = (uint32_t)md_array_size(draw_ops),
                .ops = draw_ops,
            },
            .view_transform = {
                .view_matrix = &state->view.param.matrix.curr.view.elem[0][0],
                .proj_matrix = &state->view.param.matrix.curr.proj.elem[0][0],
                // These two are for temporal anti-aliasing reprojection (optional)
                .prev_view_matrix = &state->view.param.matrix.prev.view.elem[0][0],
                .prev_proj_matrix = &state->view.param.matrix.prev.proj.elem[0][0],
            },
            .picking_offset = gl_picking_offset(*state),
            .max_bond_length = max_bond_length,
        };

        draw_ops_timed(args, draw_ops, md_array_size(draw_ops), false, frame_alloc);

        if (md_array_size(glyph_ops) > 0) {
            args.picking_offset = gl_picking_offset_glyphs(*state);
            draw_ops_timed(args, glyph_ops, md_array_size(glyph_ops), true, frame_alloc);
        }
#if EXPERIMENTAL_GFX_API
    }
#endif
}

// What the opaque image under the isosurfaces depends on, as far as the mesh mode's accumulation is concerned:
// the depth of what the representations draw (draw_representations_opaque): not what is highlighted or selected
static uint64_t opaque_scene_hash(const ApplicationState* state) {
    uint64_t h = md_hash64(&state->mold.gpu_scene_version, sizeof(state->mold.gpu_scene_version), state->representation.gl_version);
    h = md_hash64(&state->representation.visibility_mask_hash, sizeof(state->representation.visibility_mask_hash), h);
    h = md_hash64(&state->animation.frame, sizeof(state->animation.frame), h);
    h = md_hash64(state->view.periodic_images, sizeof(state->view.periodic_images), h);
    h = md_hash64(&state->view.assembly, sizeof(state->view.assembly), h);
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        const Representation& rep = state->representation.reps[i];
        const bool shown = representation_shown(*state, rep) && rep.type_is_valid;
        h = md_hash64(&shown, sizeof(shown), h);
        if (!shown) continue;
        h = md_hash64(&rep.type, sizeof(rep.type), h);
        h = md_hash64(&rep.md_rep, sizeof(rep.md_rep), h);
        h = md_hash64(&rep.scale, sizeof(rep.scale), h);
        h = md_hash64(&rep.bond_orders, sizeof(rep.bond_orders), h);
        h = md_hash64(&rep.chain_gaps, sizeof(rep.chain_gaps), h);
        h = md_hash64(&rep.elastic_bonds, sizeof(rep.elastic_bonds), h);
        h = md_hash64(&rep.noncovalent_bonds, sizeof(rep.noncovalent_bonds), h);
        if (rep.type == RepresentationType::Bonds) {
            h = md_hash64(&rep.bonds_dashed, sizeof(rep.bonds_dashed), h);
        }
        h = md_hash64(&rep.periodic_images, sizeof(rep.periodic_images), h);
        if (rep.type == RepresentationType::Cartoon) {
            h = md_hash64(&rep.two_head_spheres, sizeof(rep.two_head_spheres), h);
            h = md_hash64(&rep.nucleic_band_radius, sizeof(rep.nucleic_band_radius), h);
            h = md_hash64(&rep.nucleic_ring_radius, sizeof(rep.nucleic_ring_radius), h);
        }
        if (rep.type == RepresentationType::SimulationBox) {
            h = md_hash64(&rep.box.thickness, sizeof(rep.box.thickness), h);
            h = md_hash64(&rep.box.axes, sizeof(rep.box.axes), h);
        }
        if (rep.type == RepresentationType::DipoleMoment) {
            h = md_hash64(&rep.dipole.dipole_key, sizeof(rep.dipole.dipole_key), h);
            h = md_hash64(&rep.dipole.dipole_index, sizeof(rep.dipole.dipole_index), h);
            h = md_hash64(&rep.dipole.offset, sizeof(rep.dipole.offset), h);
            h = md_hash64(&rep.dipole.scale, sizeof(rep.dipole.scale), h);
            h = md_hash64(&rep.dipole.radius, sizeof(rep.dipole.radius), h);
        }
    }
    return h;
}

// Pixels of the render target to a point of the UI, for what is sized in points over the view (the width of an
// isoline): the framebuffer scale, and in a screenshot at another resolution as much larger as it is
static float isoline_px_per_point(const ApplicationState* state) {
    const float h = ImGui::GetIO().DisplaySize.y;
    return h > 0.0f ? (float)state->gbuffer.height / h : 1.0f;
}

// The isosurfaces of every electronic structure representation as meshes, in one composite
static bool draw_representations_transparent_mesh(ApplicationState* state) {
    static constexpr size_t MAX_SURFACES = 256;
    isomesh::Surface surfaces[MAX_SURFACES];
    size_t num_surfaces = 0;
    uint64_t scene_hash = opaque_scene_hash(state);

    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        Representation& rep = state->representation.reps[i];
        if (rep.type != RepresentationType::ElectronicStructure || !representation_shown(*state, rep) || !rep.type_is_valid) continue;
        ElectronicStructureRepresentation& es = rep.electronic_structure;

        IsoDesc iso;
        electronic_structure_iso_desc_init(&iso, es);
        const size_t n = electronic_structure_mesh_surfaces(surfaces + num_surfaces, MAX_SURFACES - num_surfaces, state, &rep, iso);

        // A field is evaluated at the vertices (Coloring_VertexField, set by electronic_structure_mesh_surfaces);
        // the colour volume spans the box of the density volume
        const mat4_t world_to_volume = mat4_inverse(es.density_vol.texture_to_world);
        const bool use_field  = es.coloring == SurfaceColoring::Field;
        const bool use_colors = es.coloring == SurfaceColoring::AtomColors && es.color_vol.tex_id;
        const uint32_t colormap = use_field ? surface_field_colormap_texture(&es.field_vol, es.field_map.colormap) : 0;
        const IsolineDesc isolines = color_scale_isoline_desc(es.field_isolines, isoline_px_per_point(state));
        for (size_t j = num_surfaces; j < num_surfaces + n; ++j) {
            isomesh::Surface& s = surfaces[j];
            s.world_to_volume = world_to_volume;
            if (s.coloring == isomesh::Coloring_VertexField) {
                s.field_colormap = colormap;
                s.field_range_beg = es.field_map.range_beg;
                s.field_range_end = es.field_map.range_end;
                s.isolines = isolines;
            } else if (use_colors) {
                s.coloring = isomesh::Coloring_ColorVolume;
                s.color_volume = es.color_vol.tex_id;
            }
        }
        // The contents of the colour volume (the vertex values are in the mesh's version)
        if (use_colors) scene_hash = md_hash64(&es.col_hash, sizeof(es.col_hash), scene_hash);
        if (use_field) scene_hash = md_hash64(&es.field_map.colormap, sizeof(es.field_map.colormap), scene_hash);
        
        num_surfaces += n;
    }

    const isomesh::RenderDesc desc = {
        .render_target = {
            .depth = state->gbuffer.tex.depth,
            .color = state->gbuffer.tex.transparency_hdr,
            .width = state->gbuffer.width,
            .height = state->gbuffer.height,
            .clear_color = true,
        },
        .matrix = {
            .view = state->view.param.matrix.curr.view,
            .proj = state->view.param.matrix.curr.proj,
            .proj_no_jitter = state->view.param.matrix.curr.proj_no_jitter,
        },
        .surfaces = surfaces,
        .num_surfaces = num_surfaces,
        .shading = {
            .roughness = 0.3f,
            .ior = 1.5f,
            .mode = state->visuals.shading.mode,
            .ambient = state->visuals.shading.ambient,
            .specular = state->visuals.shading.specular,
            .colors = state->visuals.shading.colors,
            .light = state->visuals.shading.light,
            .environment = state->visuals.shading.environment,
        },
        .scene_hash = scene_hash,
        .keep_drawing = state->profiler_keep_isomeshes,
    };
    if (num_surfaces == 0) {
        isomesh::target_free(&iso_mesh_target);
        return false;
    }
    return isomesh::render(&iso_mesh_target, desc, &iso_mesh_converged);
}

static bool draw_representations_transparent(ApplicationState* state) {
    ASSERT(state);
    if (state->mold.sys.atom.count == 0) return false;

    const size_t num_representations = md_array_size(state->representation.reps);
    if (num_representations == 0) return false;

    iso_mesh_converged = true;
    if (state->settings.isosurface_rendering == IsosurfaceRendering_Mesh) {
        return draw_representations_transparent_mesh(state);
    }
    isomesh::target_free(&iso_mesh_target);

    bool written = false;

    for (size_t i = 0; i < num_representations; ++i) {
        Representation& rep = state->representation.reps[i];
        if (!representation_shown(*state, rep)) continue;
        if (rep.type == RepresentationType::ElectronicStructure) {
            const bool use_field = rep.electronic_structure.coloring == SurfaceColoring::Field && rep.electronic_structure.field_vol.tex_id;
            IsoDesc iso;
            electronic_structure_iso_desc_init(&iso, rep.electronic_structure);

#if VIAMD_RECOMPUTE_ORBITAL_PER_FRAME
            flag_representation_as_dirty(&state->representation.reps[i]);
#endif

            volume::IsoRenderDesc desc = {
                .render_target = {
                    .depth = state->gbuffer.tex.depth,
                    .color = state->gbuffer.tex.transparency_hdr,
                    .width = state->gbuffer.width,
                    .height = state->gbuffer.height,
                    .clear_color = !written,
                },
                .texture = {
                    .density_volume = rep.electronic_structure.density_vol.tex_id,
                    .color_volume = rep.electronic_structure.color_vol.tex_id,
                    .field_volume = use_field ? rep.electronic_structure.field_vol.tex_id : 0,
                    .field_colormap = use_field ? surface_field_colormap_texture(&rep.electronic_structure.field_vol, rep.electronic_structure.field_map.colormap) : 0,
                },
                .matrix = {
                    .model = rep.electronic_structure.density_vol.texture_to_world,
                    .view = state->view.param.matrix.curr.view,
                    .proj = state->view.param.matrix.curr.proj,
                    .inv_proj = state->view.param.matrix.inv.proj,
                },
                .clip_volume = {
                    .min = {0,0,0},
                    .max = {1,1,1},
                },
                .iso = {
                    .count = iso.count,
                    .values = iso.values,
                    .colors = iso.colors,
                    .optical_densities = iso.optical_densities,
                    .use_color_volume = rep.electronic_structure.coloring == SurfaceColoring::AtomColors,
                    .use_field = use_field,
                    .exact = state->settings.exact_isosurfaces,
                },
                .field = {
                    .range_beg = rep.electronic_structure.field_map.range_beg,
                    .range_end = rep.electronic_structure.field_map.range_end,
                },
                .isolines = use_field ? color_scale_isoline_desc(rep.electronic_structure.field_isolines, isoline_px_per_point(state)) : IsolineDesc{},
                .shading = {
                    .roughness = 0.3f,
                    .ior = 1.5f,
                    .mode = state->visuals.shading.mode,
                    .ambient = state->visuals.shading.ambient,
                    .specular = state->visuals.shading.specular,
                    .colors = state->visuals.shading.colors,
                    .light = state->visuals.shading.light,
                    .environment = state->visuals.shading.environment,
                },
            };

            if (volume::render_isosurfaces(desc)) {
                written = true;
            }

#if DEBUG
            {
				immediate::Scope scope(state->gfx.world, "debug_electronic_structure");
                immediate::box_wireframe(scope, { 0,0,0 }, { 1,1,1 }, rep.electronic_structure.density_vol.texture_to_world, immediate::COLOR_BLACK);
            }
#endif
        }
    }
    return written;
}

static void draw_representations_marks(ApplicationState* data) {
    md_gl_draw_op_t* draw_ops = 0;
    md_gl_draw_op_t* glyph_ops = 0;
    // In the copies as the opaque pass draws them (the periodic images, an assembly), so the selected and highlighted
    // atoms show in them too
    for (size_t i = 0; i < md_array_size(data->representation.reps); ++i) {
        const Representation& rep = data->representation.reps[i];

        if (rep.type > RepresentationType::Cartoon && rep.type != RepresentationType::Bonds) continue;
        if (rep.type == RepresentationType::Cartoon && representation_shown(*data, rep) && rep.type_is_valid) {
            const size_t beg = md_array_size(glyph_ops);
            representation_glyph_draw_ops(&glyph_ops, *data, rep, frame_alloc);
            place_ops(&glyph_ops, beg, *data, rep, true, frame_alloc);
            if (data->mold.sys.protein_backbone.range.count == 0) continue;
        }

        if (representation_shown(*data, rep) && rep.type_is_valid) {
            // The op of the opaque pass: its depth tells the marks in view from those behind it, so the geometry must
            // be the same, the bond orders too
            const size_t beg = md_array_size(draw_ops);
            if (rep.type != RepresentationType::Bonds) {
                md_array_push(draw_ops, representation_draw_op(rep), frame_alloc);
            }
            representation_weak_bond_draw_ops(&draw_ops, *data, rep, frame_alloc);
            place_ops(&draw_ops, beg, *data, rep, false, frame_alloc);
        }
    }

    const float max_bond_length = bond_length_limit(data);

    // Unjittered: the marks are put over the image made from the jittered frames, which shows the scene still
    md_gl_draw_args_t args = {
        .shaders = data->gl.shaders_marks,
        .draw_operations = {
            .count = md_array_size(draw_ops),
            .ops = draw_ops,
        },
        .view_transform = {
            .view_matrix = &data->view.param.matrix.curr.view.elem[0][0],
            .proj_matrix = &data->view.param.matrix.curr.proj_no_jitter.elem[0][0],
        },
        .picking_offset = gl_picking_offset(*data),
        .atom_mask = AtomBit_Marked | AtomBit_Visible,
		.max_bond_length = max_bond_length,
    };

    md_gl_draw(&args);

    if (md_array_size(glyph_ops) > 0) {
        args.draw_operations = { .count = md_array_size(glyph_ops), .ops = glyph_ops };
        args.picking_offset = gl_picking_offset_glyphs(*data);
        md_gl_draw(&args);
    }
}
