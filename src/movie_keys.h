#pragma once

#include <gfx/camera.h>
#include <movie_overlay.h>

#include <stddef.h>
#include <stdint.h>
#include <array>
#include <string>
#include <vector>

struct md_system_t;
struct md_system_state_t;
struct md_bitfield_t;
struct md_allocator_i;

// Centre of a selected atom/group, made whole across periodic boundaries and transformed to camera space.
// Returns false for an empty selection or indices outside the current system. alloc must be a VM arena.
bool movie_target_center(vec3_t* out, const md_system_t& system, const md_system_state_t& state,
                         const md_bitfield_t& mask, const mat4_t& world_transform, md_allocator_i* alloc);

// ## What a key looks at and what is sharp
//
// A key has a pose (where the eye is), a look-at and a focus, and the look-at and the focus are set on the key without touching
// the rest. The atoms that a key tracks are kept in a table of the movie, by id, so that a key stays a plain value.

struct MovieAtomSet {
    uint32_t id = 0;
    std::vector<uint32_t> atoms;   // Sorted, without duplicates
    std::string expr;              // A selection in the language of the filters (protein, resname("AIN")): then 'atoms' is not used
};

// The id of a set with these atoms: one the table has already, or a new one. 0 for no atoms.
uint32_t movie_atoms_add(std::vector<MovieAtomSet>* sets, std::vector<uint32_t> atoms);
// The id of a set that is this selection expression (one the table has already, or a new one). 0 for an empty expression.
uint32_t movie_atoms_add_expr(std::vector<MovieAtomSet>* sets, std::string expr);
const MovieAtomSet* movie_atoms_find(const std::vector<MovieAtomSet>& sets, uint32_t id);
// Takes out the sets that no key uses
void movie_atoms_prune(std::vector<MovieAtomSet>* sets, const CameraKeyframe* keys, size_t count);

// The focus of a key as numbers, for the workspace: target, point (3), distance, blur, transition, ease, set. The decoder is
// false for values that no key could have made (not finite, an unknown target or ease).
std::array<float, 9> movie_focus_encode(const CameraKeyframe& key);
bool movie_focus_decode(CameraKeyframe* key, const float (&values)[9]);

// The look-at at 'time': through the points where the keys aim ('looks', aligned with 'keys', which are sorted by time and have
// different times; the caller resolves those that track atoms), with the easing of the key each stretch leads to.
vec3_t movie_look_evaluate(const CameraKeyframe* keys, size_t count, const vec3_t* looks, double time);

// The eye at 'time': through the keys' positions, with their easing.
vec3_t movie_position_evaluate(const CameraKeyframe* keys, size_t count, double time, bool loop);

// The camera at 'time': the eye from the positions of the keys, the aim from the look-at points, the field of view, the roll
// and the spin from the keys. A spin turns the eye about the look-at. With 'upright' the camera is kept level about it.
void movie_keys_pose(ViewTransform* out, float* fov_y, const CameraKeyframe* keys, size_t count, double time,
                       const vec3_t* looks, const vec3_t* upright, vec3_t up_axis);

// Turns a camera to look at 'look' from the eye of 'unspun' turned about 'look' by 'spin' (or by how the camera is turned from
// 'unspun' when none is given), kept level about 'up' and tilted by 'roll' when there is one.
void movie_camera_independent(ViewTransform* camera, const ViewTransform& unspun, vec3_t look, const vec3_t* up, float roll, const quat_t* spin = nullptr);

// ## Workspaces from before a key had its own look-at and focus

// A key of the separate look-at or focus track of an earlier version: its own time, a target and, for the focus, a blur.
struct MovieLegacyTarget {
    double time = 0.0;
    FocusTarget target = FocusTarget::Point;
    vec3_t point = {};
    float distance = 10.0f;
    float blur = 0.0f;
    float transition = 1.0f;
    KeyEase ease = KeyEase::EaseInOut;
    std::vector<uint32_t> atoms;   // Selection
};
bool movie_legacy_target_decode(MovieLegacyTarget* key, const float (&values)[10], bool focus);

// Puts the look-at and focus tracks of an earlier version on the camera keys: a camera key is added (on the path, so that the
// movie goes the same way) at every time a track has a key, and each key then looks at what the look-at track gives at its time
// (a set of atoms where the track key is a selection) and sets the focus of the focus key that is at its time. 'keys' must be
// sorted and not empty.
void movie_keys_from_tracks(std::vector<CameraKeyframe>* keys, std::vector<MovieAtomSet>* sets, const std::vector<MovieLegacyTarget>& look,
                            const std::vector<MovieLegacyTarget>& focus, const vec3_t* upright, vec3_t up_axis);

// Keys that followed the movie's follow target (or an atom of their own) of an earlier version look at the centre of those atoms
void movie_keys_from_follow(std::vector<CameraKeyframe>* keys, std::vector<MovieAtomSet>* sets, const std::vector<uint32_t>& target_atoms);

// Whether any key sets the focus
bool movie_focus_keys_exist(const CameraKeyframe* keys, size_t count);

// What is sharp at 'time'. 'points' are aligned with 'keys': the live Selection centres.
// A key without an explicit target returns focus to look-at. Explicit transitions start at the key's time
// and finish no later than the next camera key.
// How blurred the rest is belongs to the look parameter of the blur.
vec3_t movie_focus_point(const CameraKeyframe* keys, size_t count, const vec3_t* points, double time, const ViewTransform& camera);

// The depth in front of the camera of that point
float movie_focus_depth(const CameraKeyframe* keys, size_t count, const vec3_t* points, double time, const ViewTransform& camera);

// What is keyed on a movie's timeline besides the camera: look parameters (background, depth of field,
// clipping ...), each with its own keys. The application has a table of the parameters, a key refers to
// one by its index there.
struct ParamKey {
    int     param = 0;
    double  time = 0.0;
    float   value[3] = {0, 0, 0};   // A colour uses all three, a scalar the first
    KeyEase ease = KeyEase::Smooth; // Shapes the stretch leading to this key
};

// The id of the look parameter for the blur of the depth of field (it is saved, and also the id of the application's table)
constexpr int MOVIE_PARAM_DOF_APERTURE = 9;

// Keys of earlier versions carried the blur of the focus (focus_blur of 0 or more; -1 is no blur of its own). This moves it to
// keys of the blur's look parameter, with the same transitions, and sets focus_blur to -1.
std::vector<ParamKey> movie_blur_keys_from_focus(CameraKeyframe* keys, size_t count);

// The value of a parameter at 'time' from the keys that are for it, which need not be sorted. Returns
// false if there are none. Passes through every key and holds outside the first and last.
bool param_keys_evaluate(float* out, int comps, const ParamKey* keys, size_t count, int param, double time);

// What can be keyed on a representation. The numbers are saved in workspaces: never renumber, add at the end.
// Scale0..2 are the components of the representation's scale (what each is depends on its type).
enum class RepProp : int {
    Visible = 0,
    Scale0,
    Scale1,
    Scale2,
    TintScale,
    Saturation,
    BaseColor,   // Shown where the colour mapping is Uniform
    TintColor,   // Shown where the tint scale is above zero
    Count,
};

// How many numbers a property has: a colour has three
inline int rep_prop_comps(int prop) { return prop == (int)RepProp::BaseColor || prop == (int)RepProp::TintColor ? 3 : 1; }

// One key of one property of one representation, which is named by its id (not its place in the list)
struct RepKey {
    uint32_t rep = 0;
    int      prop = 0;               // A RepProp
    double   time = 0.0;
    float    value[3] = {0, 0, 0};   // A colour uses all three, the others the first; Visible is 0 or 1
    KeyEase  ease = KeyEase::Smooth; // Shapes the stretch leading to this key
};

// The value of a property of a representation at 'time' from the keys for it, which need not be sorted.
// Returns false if there are none. Visible is held: it changes at its keys, it is not blended between them.
// 'out' has rep_prop_comps(prop) numbers.
bool rep_keys_evaluate(float* out, const RepKey* keys, size_t count, uint32_t rep, int prop, double time);

// How much of a representation is there at 'time', 0 (gone) to 1 (all of it), from its Visible keys. At a key the
// factor starts to go to the key's value (1 shown, 0 hidden), taking 'transition' seconds, smoothly; from wherever
// it was if the previous change has not finished. The first key holds from the start. With a transition of 0 it
// is held and changes at its keys, as rep_keys_evaluate does. Returns false if there are no keys.
bool rep_visible_factor(float* out, const RepKey* keys, size_t count, uint32_t rep, double time, double transition);

// ## Representations over time as stretches
//
// What the Visible keys of a representation say, as the times it is shown: from the key that shows it to the key that hides
// it, or to the end of the movie. (The transition of Visible is not part of this: it comes after the key.)

struct RepInterval {
    double begin = 0.0;
    double end = 0.0;
    int    begin_key = -1;            // The key (an index into the keys it was made from) that shows it
    int    end_key = -1;              // The key that hides it again, -1 when nothing does: it lasts to the end of the movie
    bool   begin_is_start = false;    // Shown from the start by its first key, which can be later: before it the first key holds
};

std::vector<RepInterval> rep_shown_intervals(const std::vector<RepKey>& keys, uint32_t rep, double duration);
// Without visibility keys, an enabled representation covers the whole movie.
std::vector<RepInterval> rep_effective_intervals(const std::vector<RepKey>& keys, uint32_t rep, double duration, bool enabled);

// Moves the ends of a stretch to new times, with the keys it is made of (and a hidden key at the start of the movie or at the
// new end where there was none, so that nothing else changes). The new times are kept inside the neighbouring stretches and the
// movie; a stretch that was only moved keeps its length. Not sorted afterwards.
void rep_move_interval(std::vector<RepKey>* keys, uint32_t rep, const RepInterval& iv, double new_begin, double new_end, double duration);

// Adds a stretch from 'begin' to 'end' where the representation is hidden, ending before the next one. False if it is shown at
// 'begin' or there is not room (less than 'min_length').
bool rep_add_interval(std::vector<RepKey>* keys, uint32_t rep, double begin, double end, double duration);

// Takes out the keys that make a stretch: the representation is hidden there afterwards
void rep_remove_interval(std::vector<RepKey>* keys, const RepInterval& iv);

// Transfers one interval, merging target visibility overlaps while retaining other property keys.
void rep_transfer_interval(std::vector<RepKey>* keys, uint32_t from, uint32_t to, const RepInterval& iv, double duration);

// At time t, 'from' goes and 'to' comes: a hidden key for one and a shown key for the other. A representation without keys is
// taken to have been shown ('from') or hidden ('to') up to then.
void rep_swap_at(std::vector<RepKey>* keys, uint32_t from, uint32_t to, double t);

// The union of stretches (a group's: shown whenever any of its members is)
std::vector<RepInterval> rep_union_intervals(std::vector<RepInterval> intervals);

// "protein-cpk": the group 'protein' and the member 'cpk'. A name without a hyphen in it (or one that starts or ends with it)
// is a group of its own: the member is empty and false is returned.
bool rep_name_split(const char* name, std::string* group, std::string* member);

// A row of the representation overview: one system, named by the part of the representation names before the first hyphen
struct RepRow {
    bool        header = false;   // The row of a group
    int         rep = -1;         // The place in the list of a representation row; for a header the first member
    std::string label;
    std::string group;
    int         members = 1;
    bool        indented = false;
};
std::vector<RepRow> rep_system_rows(const std::vector<std::string>& names);

struct RepBlock {
    int rep = -1;
    int interval_index = 0;
    RepInterval interval;
    int slot = 0;
};
int rep_pack_blocks(std::vector<RepBlock>* blocks, double transition);

// ## The camera keys as a lane
//
// What the camera does between its keys, as spans on the timeline. Keys are sorted by time.

enum class CameraBandKind : int {
    LookAtSet,      // The key looks at the centre of a set of atoms, tracked through the trajectory
    Spin,           // Extra turns around what is looked at, in the stretch leading to a key
};

struct CameraBand {
    CameraBandKind kind = CameraBandKind::LookAtSet;
    int      first = 0;    // The keys it starts and ends at (the same one for a single key that tracks)
    int      last = 0;
    double   begin = 0.0;
    double   end = 0.0;
    uint32_t set = 0;      // LookAtSet
    int      turns = 0;    // Spin
};

// One band for each run of keys that look at the same set, one spin band for each key with turns
std::vector<CameraBand> camera_bands(const std::vector<CameraKeyframe>& keys);

// "3", or "3 intro" for a key with a name. 'index' counts from 0.
std::string camera_key_label(const CameraKeyframe& key, int index);

// A key at 'time' with the pose 'pose' (and field of view), that the path of 'keys' (not empty, sorted) has there, so that adding it
// does not move the camera. It looks at what the keys on both sides look at, if that is the same set of atoms; it has no frame and
// no spin, and leaves the focus as it is. With upright (the movie keeps the camera level about it), it has the roll of the pose.
CameraKeyframe camera_key_on_path(const ViewTransform& pose, float fov_y, const std::vector<CameraKeyframe>& keys, double time, const vec3_t* upright = nullptr);

// ## The camera path in the viewport

// Where the eye and what it looks at are along the movie, sampled in time (sorted, at least one sample)
struct CameraPathSamples {
    std::vector<double> time;
    std::vector<vec3_t> eye, look;
};

// The eye and the look-at point at 'time', between the samples (held before the first and after the last)
void camera_path_at(const CameraPathSamples& samples, double time, vec3_t* eye, vec3_t* look);

// A round step (0.1, 0.2, 0.5, 1, 2, 5, 10 ... seconds) for ticks along 'span' seconds that gives at most 'max_ticks' of them
double camera_tick_step(double span, int max_ticks);

// Cuts a line between two points in clip space at the camera (where w is 0), keeping the part in front of it. False if none is.
bool clip_segment_near(vec4_t* a, vec4_t* b);

// The ring that the eye of key 'to' goes round when the camera spins in the stretch from key 'from' to it: its centre, two
// perpendicular directions in its plane (going from u to v is a positive turn) and its radius. False if the eye is on the axis.
bool camera_spin_ring(const CameraKeyframe& from, const CameraKeyframe& to, vec3_t* center, vec3_t* u, vec3_t* v, float* radius);

// The point where a ray meets a plane. False if they are parallel or the plane is behind the ray.
bool ray_plane_hit(vec3_t origin, vec3_t dir, vec3_t plane_point, vec3_t plane_normal, vec3_t* out);

// The point of a polyline (points on the screen; 'ok' tells which of them can be used, and a segment needs both of its ends)
// that is nearest to 'q': the distance, the segment and how far along it. FLT_MAX when there is none.
float polyline_nearest(const std::vector<vec2_t>& points, const std::vector<char>& ok, vec2_t q, int* segment, float* along);

// Editing a key by its handles in the viewport: the eye moves and what it looks at stays, the other way round, or both move
bool camera_key_set_eye(CameraKeyframe* key, vec3_t eye);
bool camera_key_set_look(CameraKeyframe* key, vec3_t look);
void camera_key_translate(CameraKeyframe* key, vec3_t delta);

// Everything on the timeline that the user edits, so that it can be undone as one
struct MovieKeys {
    std::vector<CameraKeyframe> camera;
    std::vector<MovieAtomSet> sets;     // The atoms that the keys look at and focus on
    std::vector<ParamKey> params;
    std::vector<RepKey> reps;
    std::vector<MovieOverlay> overlays;
    std::vector<MovieMarker> markers;
    bool loop = false;
    bool keep_upright = false;
    int  up_axis = 1;
    // The timing, which is scaled together with the keys
    float  duration = 0.0f;
    float  traj_begin = 0.0f;
    float  traj_end = 0.0f;
    double start_frame = 0.0;
    double end_frame = 0.0;
};

bool movie_keys_equal(const MovieKeys& a, const MovieKeys& b);
void movie_keys_scale_time(MovieKeys* keys, double scale);

// ## Selecting keys and moving them together
//
// A key is picked out by what it belongs to and its time (a time is unique among the keys of one thing), so a selection stays
// valid when the keys are re-sorted. A camera key is one key however many lanes it is drawn in.

// Overlay: a bar of the overlay lane (subject: its place in the list, time: when it starts, end: when it stops).
// Block: a stretch of the representation overview where a representation is shown (subject: its id, time and end as the stretch).
enum class KeyKind : int { Camera, Param, Rep, Overlay, Block };

// What a key of a parameter or of a property of a representation belongs to
inline int64_t rep_key_subject(uint32_t rep, int prop) { return ((int64_t)rep << 8) | (int64_t)prop; }

struct KeyId {
    KeyKind kind = KeyKind::Camera;
    int64_t subject = 0;       // The parameter, rep_key_subject() of a representation, 0 for the camera
    double  time = 0.0;
    double  end = 0.0;         // Overlay and Block: when it stops
};

struct KeySelection {
    std::vector<KeyId> ids;

    bool   contains(KeyKind kind, int64_t subject, double time) const;
    void   add(KeyKind kind, int64_t subject, double time, double end = 0.0);
    void   toggle(KeyKind kind, int64_t subject, double time, double end = 0.0);
    void   set(KeyKind kind, int64_t subject, double time, double end = 0.0);
    void   clear() { ids.clear(); }
    bool   empty() const { return ids.empty(); }
    size_t size() const { return ids.size(); }
};

// Takes out what no key is there for any more (after an undo or a removal)
void key_selection_prune(KeySelection* sel, const MovieKeys& keys);

// All the camera keys, the keys of one parameter (-1: none) and the keys of one representation property (-1: none)
void key_selection_all(KeySelection* sel, const MovieKeys& keys, int64_t param, int64_t rep_subject);

// The value part of a move: which quantity of the keys changes, in the units the lane shows
enum class KeyLane : int { None, Frame, Distance, Fov, Param, Rep };

struct KeyShift {
    KeyLane lane = KeyLane::None;
    int64_t subject = 0;        // Param: the parameter, Rep: rep_key_subject()
    double  dy = 0.0;           // Added to the value (with 'ratio': multiplies it)
    bool    ratio = false;
    double  lo = 0.0, hi = 0.0; // Limits of Frame, Param and Rep
    double  unit = 1.0;         // Distance: the display unit per internal unit
    double  step = 1.0;         // How far an arrow key moves the value (not used by movie_keys_shift)
    bool    step_ratio = false; // An arrow key multiplies the value by 1.05 instead of adding 'step'
};

// Moves the keys, overlay bars and representation blocks of 'start' that are in 'start_sel' by 'dt' seconds and, in the lane, the values
// of the keys by shift.dy: the result is put in 'keys' (the camera, parameter and representation keys and the overlays; the order is
// kept, nothing is sorted) and the selection in 'sel'. A bar takes what is timed inside it along; a block stops at the other blocks
// of its representation (so it may move less than the rest).
// The selection is kept inside 0 .. duration as a whole: the dt that was used is returned. A distance changes with the
// camera looking at the same point, and every value stays within its limits.
double movie_keys_shift(MovieKeys* keys, KeySelection* sel, const MovieKeys& start, const KeySelection& start_sel, double dt, const KeyShift& shift, double duration);

// The first and last time of what is selected (keys, bars, blocks). False if nothing is.
bool key_selection_extent(const MovieKeys& keys, const KeySelection& sel, double* t0, double* t1);

// Stretches the selected items in time about 'anchor': a time t goes to anchor + (t - anchor) * factor (the length of a bar or a
// block too, and what is timed inside a bar). The factor is kept so that nothing leaves the movie, and the one used is returned.
// Made like movie_keys_shift: from 'start', into 'keys' (camera, parameter, representation keys and overlays) and 'sel', not sorted.
double movie_keys_scale(MovieKeys* keys, KeySelection* sel, const MovieKeys& start, const KeySelection& start_sel, double anchor, double factor, double duration);

// The easing of the picked keys that have one that can be chosen (the camera keys but the first, the parameter keys and the
// properties of representations but Visible, which holds). Returns how many there are; 'mixed' is whether they differ, and
// 'common' the easing of the first.
int key_selection_ease(const CameraKeyframe* camera, size_t num_camera, const std::vector<ParamKey>& params, const std::vector<RepKey>& reps,
    const KeySelection& sel, KeyEase* common, bool* mixed);

// Sets it on those keys
void movie_keys_set_ease(MovieKeys* keys, const KeySelection& sel, KeyEase ease);

// Sorts the keys. Where a key lies on another (within a millisecond) one of them goes: the selected one stays.
void movie_keys_resolve(MovieKeys* keys, const KeySelection& sel);

// Removes the selected keys and representation blocks (an overlay is not removed this way)
void movie_keys_delete(MovieKeys* keys, KeySelection* sel);

// Selected keys remembered to be put in somewhere else
struct KeyClip {
    std::vector<CameraKeyframe> camera;
    std::vector<ParamKey> params;
    std::vector<RepKey> reps;
    double begin = 0.0;         // The time of the first
    bool empty() const { return camera.empty() && params.empty() && reps.empty(); }
};

KeyClip movie_keys_copy(const MovieKeys& keys, const KeySelection& sel);

// Puts a clip in with its first key at 'time' (earlier if it would go past the end of the movie), a key on another replacing it.
// The keys that came in are the selection afterwards; the keys are sorted.
void movie_keys_paste(MovieKeys* keys, KeySelection* sel, const KeyClip& clip, double time, double duration);

// Undo and redo of edits to the keys. It is not told what changed: it is shown the keys every frame, and
// whenever they differ from the last committed state and no edit is under way (a slider being dragged, a
// point being moved) that is one step. A drag is therefore one undo, not hundreds.
class MovieHistory {
public:
    // Forgets everything, e.g. when another workspace is loaded
    void clear(const MovieKeys& current);

    void update(const MovieKeys& current, bool editing);

    bool can_undo() const { return !undo_.empty(); }
    bool can_redo() const { return !redo_.empty(); }

    // 'state' is the present keys. On true it holds the ones to go back to (or forward to).
    bool undo(MovieKeys* state);
    bool redo(MovieKeys* state);

private:
    static constexpr size_t MAX_STEPS = 200;

    std::vector<MovieKeys> undo_;
    std::vector<MovieKeys> redo_;
    MovieKeys committed_;
    bool has_committed_ = false;
};

// The frames of a movie that a recording makes, first to last (inclusive), out of 'num_frames' at 'fps'. With a
// range, the frames between begin and end seconds, which always include at least one.
void movie_frame_range(int num_frames, double fps, bool use_range, double begin, double end, int* first, int* last);

// The size of a frame at a percentage of the full size, kept even for the video encoder
void movie_scaled_size(int* w, int* h, int percent);

// How long a recording has left from its pace so far: 'done' of 'total' frames in 'active_seconds'. False
// until there is enough to go on.
bool movie_time_left(int done, int total, double active_seconds, double* seconds);

// A time moved to the nearest frame of the movie (at 'fps'), kept within 0 .. duration
double movie_snap_to_frame(double time, double fps, double duration);
