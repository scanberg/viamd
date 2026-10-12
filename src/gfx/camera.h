#pragma once

#include <core/md_vec_math.h>

struct ViewTransform {
    quat_t orientation = {0, 0, 0, 1};
    vec3_t position = {0, 0, 0};
	float  distance = 10.0f;
};

struct Camera : ViewTransform {
    float near_plane = 1.0f;
    float far_plane = 10000.0f;
    float fov_y = (3.1415926534f / 4.0f);

    Camera& operator = (const ViewTransform& t) {
        this->orientation = t.orientation;
        this->position = t.position;
        this->distance = t.distance;
        return *this;
	}
};

// How the camera, and the trajectory frame, move in the stretch leading to a keyframe
enum class KeyEase : int {
    Smooth,     // Through the key without stopping
    EaseInOut,  // Starts and ends slowly, so it stops at the keys
    Linear,     // Constant speed in a straight line
    Hold,       // Stays as it was until the key, where it jumps
    Count,
};

// What a keyframe's extra turns go around
enum class SpinAxis : int {
    ViewUp,   // The up direction of the camera at the start of the segment
    WorldY,
    WorldX,
    WorldZ,
    Count,
};

// What is sharp from a key on
enum class FocusTarget : int {
    Point,       // A fixed point
    Selection,   // The centre of a set of atoms, tracked through the trajectory
    LookAt,      // What the camera looks at
    Distance,    // A distance from the camera
    Count,
};

// A camera pose that a movie passes through at 'time' seconds on the movie timeline.
struct CameraKeyframe {
    ViewTransform transform = {};
    float  fov_y = (3.1415926534f / 4.0f);
    double time = 0.0;

    KeyEase ease = KeyEase::Smooth;
    // The trajectory frame shown at this time. The keys that have one decide how the trajectory plays in the movie.
    bool   use_frame = false;
    double frame = 0.0;

    // Extra whole turns of the camera around what it looks at, made over the segment that ends at this
    // key (positive is counter-clockwise seen from the tip of the axis). A whole number so that the camera
    // still arrives at this key's pose.
    int      spin_turns = 0;
    SpinAxis spin_axis = SpinAxis::ViewUp;
    bool     spin_constant_speed = false;  // Otherwise it eases in and out

    // The look-at point moves with a target (the movie's follow target) rather than staying put. follow_center
    // is where the target was when the key was made, so the key's look-at is kept relative to it.
    bool   follow = false;
    vec3_t follow_center = {0, 0, 0};
    // With follow set: the atom this key's look-at tracks. Negative means the movie's follow target.
    int32_t follow_atom = -1;

    // What the user calls this key ("intro", "close-up"), shown in the lane and the table. Empty: only its number.
    char name[24] = {};

    // The tilt of the view about the direction it looks in, in radians: positive leans the camera's up to the left,
    // so the image turns clockwise. Used when the movie keeps the camera upright: then the key's stored orientation
    // only gives where it looks.
    float roll = 0.0f;

    // Where the key looks at is the point its pose aims at, unless 'look_set' names a set of atoms (an id in the movie's table
    // of atom sets): then it looks at the centre of the set wherever it is in the trajectory, and the point the pose aims at
    // is only where that was when the key was made.
    uint32_t look_set = 0;

    // Without an explicit target, focus follows look-at. Explicit changes take focus_transition seconds from this key.
    bool        focus_on = false;
    FocusTarget focus_target = FocusTarget::LookAt;
    uint32_t    focus_set = 0;                 // Selection: the atoms
    vec3_t      focus_point = {0, 0, 0};       // Point, and where the set was when Selection was made
    float       focus_distance = 10.0f;        // Distance
    float       focus_blur = -1.0f;            // Only keys of earlier versions have one (0 or more); the blur is a look parameter now
    float       focus_transition = 1.0f;
    KeyEase     focus_ease = KeyEase::EaseInOut;
};
