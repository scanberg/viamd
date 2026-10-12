#include "movie_keys.h"

#include <gfx/camera_utils.h>
#include <md_system.h>
#include <md_util.h>
#include <core/md_bitfield.h>
#include <core/md_allocator.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>

bool movie_target_center(vec3_t* out, const md_system_t& system, const md_system_state_t& state,
                         const md_bitfield_t& mask, const mat4_t& world_transform, md_allocator_i* alloc) {
    const size_t count = md_bitfield_popcount(&mask);
    uint64_t first = 0, last = 0;
    if (!out || count == 0 || !state.xyz || state.num_atoms != system.atom.count ||
        !md_bitfield_get_range(&first, &last, &mask) || last >= system.atom.count) return false;

    md_temp_scope_t temp = md_temp_begin_in(alloc);
    defer { md_temp_end(temp); };
    vec4_t* xyzw = md_temp_alloc_array(temp, vec4_t, count);
    md_util_system_extract_xyzw_from_mask(xyzw, &mask, &system, &state);
    vec3_t center = vec3_zero();
    md_util_deperiodize_self_vec4(xyzw, count, &state.unitcell, &center);
    *out = mat4_mul_vec3(world_transform, center, 1.0f);
    return true;
}

uint32_t movie_atoms_add(std::vector<MovieAtomSet>* sets, std::vector<uint32_t> atoms) {
    std::sort(atoms.begin(), atoms.end());
    atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
    if (atoms.empty()) return 0;
    uint32_t next = 1;
    for (const MovieAtomSet& s : *sets) {
        if (s.expr.empty() && s.atoms == atoms) return s.id;
        next = std::max(next, s.id + 1);
    }
    sets->push_back({next, std::move(atoms)});
    return next;
}

uint32_t movie_atoms_add_expr(std::vector<MovieAtomSet>* sets, std::string expr) {
    const size_t first = expr.find_first_not_of(" \t\r\n"), last = expr.find_last_not_of(" \t\r\n");
    if (first == std::string::npos) return 0;
    expr = expr.substr(first, last - first + 1);
    uint32_t next = 1;
    for (const MovieAtomSet& s : *sets) {
        if (s.expr == expr) return s.id;
        next = std::max(next, s.id + 1);
    }
    MovieAtomSet set;
    set.id = next;
    set.expr = std::move(expr);
    sets->push_back(std::move(set));
    return next;
}

const MovieAtomSet* movie_atoms_find(const std::vector<MovieAtomSet>& sets, uint32_t id) {
    if (id == 0) return nullptr;
    for (const MovieAtomSet& s : sets) if (s.id == id) return &s;
    return nullptr;
}

void movie_atoms_prune(std::vector<MovieAtomSet>* sets, const CameraKeyframe* keys, size_t count) {
    sets->erase(std::remove_if(sets->begin(), sets->end(), [&](const MovieAtomSet& s) {
        for (size_t i = 0; i < count; ++i) {
            if (keys[i].look_set == s.id || (keys[i].focus_on && keys[i].focus_target == FocusTarget::Selection && keys[i].focus_set == s.id)) return false;
        }
        return true;
    }), sets->end());
}

static bool finite_bits(float v) {
    // Floating-point predicates can be optimized away by the application's fast-math flags.
    uint32_t bits;
    memcpy(&bits, &v, sizeof(bits));
    return (bits & 0x7f800000u) != 0x7f800000u;
}

std::array<float, 9> movie_focus_encode(const CameraKeyframe& k) {
    return {(float)(int)k.focus_target, k.focus_point.x, k.focus_point.y, k.focus_point.z, k.focus_distance, k.focus_blur,
        k.focus_transition, (float)(int)k.focus_ease, (float)k.focus_set};
}

bool movie_focus_decode(CameraKeyframe* key, const float (&v)[9]) {
    for (float value : v) if (!finite_bits(value)) return false;
    if (v[0] < 0 || v[0] >= (float)FocusTarget::Count || v[0] != floorf(v[0]) ||
        v[7] < 0 || v[7] >= (float)KeyEase::Count || v[7] != floorf(v[7]) || v[8] < 0 || v[8] != floorf(v[8])) return false;
    key->focus_on = true;
    key->focus_target = (FocusTarget)(int)v[0];
    key->focus_point = vec3_set(v[1], v[2], v[3]);
    key->focus_distance = MAX(v[4], 0.001f);
    key->focus_blur = v[5] >= 0.0f ? CLAMP(v[5], 0.0f, 4.0f) : -1.0f;
    key->focus_transition = MAX(v[6], 0.0f);
    key->focus_ease = (KeyEase)(int)v[7];
    key->focus_set = (uint32_t)v[8];
    return true;
}

vec3_t movie_look_evaluate(const CameraKeyframe* keys, size_t count, const vec3_t* looks, double time) {
    ASSERT(count > 0);
    std::vector<double> times(count), values(count);
    std::vector<KeyEase> eases(count);
    for (size_t i = 0; i < count; ++i) { times[i] = keys[i].time; eases[i] = keys[i].ease; }
    vec3_t result = {};
    for (int c = 0; c < 3; ++c) {
        for (size_t i = 0; i < count; ++i) values[i] = looks[i].elem[c];
        result.elem[c] = (float)keyed_curve_evaluate(times.data(), values.data(), eases.data(), count, time);
    }
    return result;
}

void movie_camera_independent(ViewTransform* camera, const ViewTransform& unspun, vec3_t look, const vec3_t* up, float roll, const quat_t* spin) {
    const quat_t orbit = spin ? *spin : quat_normalize(camera->orientation * quat_conj(unspun.orientation));
    camera->position = look + orbit * (unspun.position - look);
    if (camera_aim_at(camera, look) && up) camera_level(camera, *up, roll);
}

vec3_t movie_position_evaluate(const CameraKeyframe* keys, size_t count, double time, bool loop) {
    if (!count) return {};
    std::vector<CameraKeyframe> ordered(keys, keys + count);
    std::stable_sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) { return a.time < b.time; });
    ordered.erase(std::unique(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) { return a.time == b.time; }), ordered.end());
    const double start = ordered.front().time, span = ordered.back().time - start;
    if (loop && span > 0.0 && time > ordered.back().time) time = start + fmod(time - start, span);
    std::vector<double> times, values;
    std::vector<KeyEase> eases;
    for (const auto& key : ordered) { times.push_back(key.time); eases.push_back(key.ease); }
    vec3_t point = {};
    for (int axis = 0; axis < 3; ++axis) {
        values.clear();
        for (const auto& key : ordered) values.push_back(key.transform.position.elem[axis]);
        point.elem[axis] = (float)keyed_curve_evaluate(times.data(), values.data(), eases.data(), times.size(), time);
    }
    return point;
}

void movie_keys_pose(ViewTransform* out, float* fov_y, const CameraKeyframe* keys, size_t count, double time,
                       const vec3_t* looks, const vec3_t* upright, vec3_t up_axis) {
    ASSERT(count > 0);
    // The field of view, the tilt and the turns of the spin come from the keys' own interpolation; where the eye and the aim are
    // is decided separately, below
    std::vector<CameraKeyframe> plain(keys, keys + count);
    for (auto& key : plain) key.follow = false;
    quat_t spin;
    camera_keyframes_evaluate(out, fov_y, plain.data(), count, time, false, nullptr, nullptr, upright, &spin);
    for (auto& key : plain) key.spin_turns = 0;
    ViewTransform unspun;
    float unused;
    camera_keyframes_evaluate(&unspun, &unused, plain.data(), count, time, false, nullptr, nullptr, upright);
    const float roll = camera_roll(unspun, up_axis);
    unspun.position = movie_position_evaluate(keys, count, time, false);
    movie_camera_independent(out, unspun, movie_look_evaluate(keys, count, looks, time), upright, roll, &spin);
}

bool movie_legacy_target_decode(MovieLegacyTarget* key, const float (&v)[10], bool focus) {
    for (float value : v) if (!finite_bits(value)) return false;
    if (v[9] != 1 || v[0] < 0 || v[1] < 0 || v[1] >= (float)FocusTarget::Count || v[1] != floorf(v[1]) ||
        (!focus && v[1] > (float)FocusTarget::Selection) || v[8] < 0 || v[8] >= (float)KeyEase::Count || v[8] != floorf(v[8])) return false;
    MovieLegacyTarget result;
    result.time = v[0];
    result.target = (FocusTarget)(int)v[1];
    result.point = vec3_set(v[2], v[3], v[4]);
    result.distance = MAX(v[5], 0.001f);
    result.blur = CLAMP(v[6], 0.0f, 4.0f);
    result.transition = MAX(v[7], 0.0f);
    result.ease = (KeyEase)(int)v[8];
    *key = std::move(result);
    return true;
}

static vec3_t legacy_look_evaluate(const std::vector<MovieLegacyTarget>& look, double time) {
    std::vector<size_t> order(look.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return look[a].time < look[b].time; });
    std::vector<double> times, values;
    std::vector<KeyEase> eases;
    std::vector<size_t> used;
    for (size_t i : order) {
        if (!times.empty() && look[i].time <= times.back()) continue;
        times.push_back(look[i].time);
        eases.push_back(look[i].ease);
        used.push_back(i);
    }
    vec3_t result = {};
    for (int c = 0; c < 3; ++c) {
        values.clear();
        for (size_t i : used) values.push_back(look[i].point.elem[c]);
        result.elem[c] = (float)keyed_curve_evaluate(times.data(), values.data(), eases.data(), times.size(), time);
    }
    return result;
}

void movie_keys_from_tracks(std::vector<CameraKeyframe>* keys, std::vector<MovieAtomSet>* sets, const std::vector<MovieLegacyTarget>& look,
                            const std::vector<MovieLegacyTarget>& focus, const vec3_t* upright, vec3_t up_axis) {
    if (keys->empty()) return;
    auto at_time = [](const std::vector<MovieLegacyTarget>& track, double t) -> const MovieLegacyTarget* {
        for (const auto& k : track) if (fabs(k.time - t) < 1.0e-6) return &k;
        return nullptr;
    };
    auto look_at = [&](double t, const CameraKeyframe& key) {
        return look.empty() ? camera_get_look_at(key.transform) : legacy_look_evaluate(look, t);
    };

    // What the camera did: aimed where the look-at track says. The keys' own aim is not used.
    const std::vector<CameraKeyframe> old = *keys;
    auto old_looks = [&]() {
        std::vector<vec3_t> looks;
        for (const auto& k : old) looks.push_back(look_at(k.time, k));
        return looks;
    };
    const std::vector<vec3_t> looks = old_looks();

    std::vector<double> times;
    for (const auto& k : old) times.push_back(k.time);
    for (const auto* track : {&look, &focus}) for (const auto& k : *track) {
        bool have = false;
        for (double t : times) have |= fabs(t - k.time) < 1.0e-6;
        if (!have) times.push_back(k.time);
    }
    std::sort(times.begin(), times.end());

    std::vector<CameraKeyframe> result;
    for (double t : times) {
        const CameraKeyframe* existing = nullptr;
        for (const auto& k : old) if (fabs(k.time - t) < 1.0e-6) existing = &k;
        CameraKeyframe key;
        if (existing) key = *existing;
        else {
            ViewTransform pose;
            float fov;
            movie_keys_pose(&pose, &fov, old.data(), old.size(), t, looks.data(), upright, up_axis);
            key = camera_key_on_path(pose, fov, old, t, upright);
            key.ease = KeyEase::Smooth;
        }
        // Aimed where the look-at track puts it at this time
        key.follow = false;
        key.follow_atom = -1;
        key.look_set = 0;
        const vec3_t aim = look_at(t, key);
        if (!camera_key_set_look(&key, aim)) key.transform.distance = MAX(key.transform.distance, 0.001f);
        if (const MovieLegacyTarget* l = at_time(look, t)) {
            if (l->target == FocusTarget::Selection) key.look_set = movie_atoms_add(sets, l->atoms);
        }
        if (const MovieLegacyTarget* f = at_time(focus, t)) {
            key.focus_on = true;
            key.focus_target = f->target;
            key.focus_point = f->point;
            key.focus_distance = f->distance;
            key.focus_blur = f->blur;
            key.focus_transition = f->transition;
            key.focus_ease = f->ease;
            key.focus_set = f->target == FocusTarget::Selection ? movie_atoms_add(sets, f->atoms) : 0;
        }
        result.push_back(key);
    }
    *keys = std::move(result);
}

void movie_keys_from_follow(std::vector<CameraKeyframe>* keys, std::vector<MovieAtomSet>* sets, const std::vector<uint32_t>& target_atoms) {
    for (CameraKeyframe& k : *keys) {
        if (!k.follow) continue;
        k.look_set = k.follow_atom >= 0 ? movie_atoms_add(sets, {(uint32_t)k.follow_atom}) : movie_atoms_add(sets, target_atoms);
        k.follow = false;
        k.follow_atom = -1;
    }
}

bool movie_focus_keys_exist(const CameraKeyframe* keys, size_t count) {
    for (size_t i = 0; i < count; ++i) if (keys[i].focus_on) return true;
    return false;
}

vec3_t movie_focus_point(const CameraKeyframe* keys, size_t count, const vec3_t* points, double time, const ViewTransform& camera) {
    struct Focus {
        double time;
        FocusTarget target;
        vec3_t point;
        float distance, transition;
        KeyEase ease;
    };
    // What is sharp before any key says: what the camera looks at
    std::vector<Focus> f;
    f.push_back({-DBL_MAX, FocusTarget::LookAt, {}, 10.0f, 0.0f, KeyEase::EaseInOut});
    for (size_t i = 0; i < count; ++i) {
        f.push_back({keys[i].time, keys[i].focus_on ? keys[i].focus_target : FocusTarget::LookAt,
            points ? points[i] : keys[i].focus_point, keys[i].focus_distance,
            keys[i].focus_on ? keys[i].focus_transition : 0.0f, keys[i].focus_ease});
    }
    std::stable_sort(f.begin(), f.end(), [](const Focus& a, const Focus& b) { return a.time < b.time; });
    size_t current = 0, previous = 0;
    for (size_t i = 0; i < f.size(); ++i) {
        if (f[i].time > time) break;
        previous = current;
        current = i;
    }
    auto point = [&](const Focus& k) {
        switch (k.target) {
        case FocusTarget::LookAt: return camera_get_look_at(camera);
        case FocusTarget::Distance: return camera.position + camera.orientation * vec3_t{0, 0, -k.distance};
        default: return k.point;
        }
    };
    float u = 1.0f;
    if (previous != current && f[current].transition > 0.0f) {
        double duration = f[current].transition;
        if (current + 1 < f.size()) duration = std::min(duration, f[current + 1].time - f[current].time);
        u = CLAMP((float)((time - f[current].time) / duration), 0.0f, 1.0f);
        if (f[current].ease == KeyEase::Hold) u = u >= 1.0f ? 1.0f : 0.0f;
        else if (f[current].ease != KeyEase::Linear) u = u * u * (3.0f - 2.0f * u);
    }
    return point(f[previous]) * (1.0f - u) + point(f[current]) * u;
}

float movie_focus_depth(const CameraKeyframe* keys, size_t count, const vec3_t* points, double time, const ViewTransform& camera) {
    return MAX(camera_depth_of_point(camera, movie_focus_point(keys, count, points, time, camera)), 1.0e-3f);
}

std::vector<ParamKey> movie_blur_keys_from_focus(CameraKeyframe* keys, size_t count) {
    std::vector<size_t> order;
    for (size_t i = 0; i < count; ++i) if (keys[i].focus_on) order.push_back(i);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return keys[a].time < keys[b].time; });
    std::vector<ParamKey> result;
    float previous = -1.0f;
    double last_time = -DBL_MAX;
    for (size_t n = 0; n < order.size(); ++n) {
        CameraKeyframe& k = keys[order[n]];
        if (k.focus_blur < 0.0f) continue;
        auto add = [&](double time, float value, KeyEase ease) {
            ParamKey p;
            p.param = MOVIE_PARAM_DOF_APERTURE;
            p.time = time;
            p.value[0] = value;
            p.ease = ease;
            result.push_back(p);
            last_time = time;
        };
        const bool ramp = previous >= 0.0f && k.focus_transition > 0.0f && previous != k.focus_blur;
        if (ramp) {
            double end = k.time + k.focus_transition;
            for (size_t m = n + 1; m < order.size(); ++m) {
                if (keys[order[m]].time > k.time) { end = std::min(end, keys[order[m]].time); break; }
            }
            if (last_time < k.time - 1.0e-9) add(k.time, previous, KeyEase::Linear);
            add(std::max(end, k.time + 1.0e-6), k.focus_blur, k.focus_ease);
        } else {
            add(k.time, k.focus_blur, KeyEase::Hold);
        }
        previous = k.focus_blur;
        k.focus_blur = -1.0f;
    }
    return result;
}

bool param_keys_evaluate(float* out, int comps, const ParamKey* keys, size_t count, int param, double time) {
    std::vector<const ParamKey*> mine;
    for (size_t i = 0; i < count; ++i) {
        if (keys[i].param == param) mine.push_back(&keys[i]);
    }
    if (mine.empty()) return false;
    std::stable_sort(mine.begin(), mine.end(), [](const ParamKey* a, const ParamKey* b) { return a->time < b->time; });

    // Keys on the same time are one key
    std::vector<double> times;
    std::vector<KeyEase> eases;
    std::vector<const ParamKey*> unique;
    for (const ParamKey* k : mine) {
        if (!times.empty() && k->time <= times.back()) continue;
        times.push_back(k->time);
        eases.push_back(k->ease);
        unique.push_back(k);
    }

    std::vector<double> values(unique.size());
    for (int c = 0; c < comps && c < 3; ++c) {
        for (size_t i = 0; i < unique.size(); ++i) values[i] = unique[i]->value[c];
        out[c] = (float)keyed_curve_evaluate(times.data(), values.data(), eases.data(), times.size(), time);
    }
    return true;
}

static bool equal(const CameraKeyframe& a, const CameraKeyframe& b) {
    return a.transform.orientation.x == b.transform.orientation.x && a.transform.orientation.y == b.transform.orientation.y &&
           a.transform.orientation.z == b.transform.orientation.z && a.transform.orientation.w == b.transform.orientation.w &&
           a.transform.position.x == b.transform.position.x && a.transform.position.y == b.transform.position.y &&
           a.transform.position.z == b.transform.position.z && a.transform.distance == b.transform.distance &&
           a.fov_y == b.fov_y && a.time == b.time && a.ease == b.ease &&
           a.use_frame == b.use_frame && a.frame == b.frame &&
           a.spin_turns == b.spin_turns && a.spin_axis == b.spin_axis && a.spin_constant_speed == b.spin_constant_speed &&
           a.follow == b.follow && a.follow_center.x == b.follow_center.x && a.follow_center.y == b.follow_center.y &&
           a.follow_center.z == b.follow_center.z && a.follow_atom == b.follow_atom && strcmp(a.name, b.name) == 0 && a.roll == b.roll &&
           a.look_set == b.look_set && a.focus_on == b.focus_on && a.focus_target == b.focus_target && a.focus_set == b.focus_set &&
           a.focus_point.x == b.focus_point.x && a.focus_point.y == b.focus_point.y && a.focus_point.z == b.focus_point.z &&
           a.focus_distance == b.focus_distance && a.focus_blur == b.focus_blur && a.focus_transition == b.focus_transition &&
           a.focus_ease == b.focus_ease;
}

static bool equal(const ParamKey& a, const ParamKey& b) {
    return a.param == b.param && a.time == b.time && a.ease == b.ease &&
           a.value[0] == b.value[0] && a.value[1] == b.value[1] && a.value[2] == b.value[2];
}

bool rep_keys_evaluate(float* out, const RepKey* keys, size_t count, uint32_t rep, int prop, double time) {
    std::vector<const RepKey*> mine;
    for (size_t i = 0; i < count; ++i) {
        if (keys[i].rep == rep && keys[i].prop == prop) mine.push_back(&keys[i]);
    }
    if (mine.empty()) return false;
    std::stable_sort(mine.begin(), mine.end(), [](const RepKey* a, const RepKey* b) { return a->time < b->time; });

    std::vector<double> times;
    std::vector<KeyEase> eases;
    std::vector<const RepKey*> unique;
    for (const RepKey* k : mine) {
        if (!times.empty() && k->time <= times.back()) continue;
        times.push_back(k->time);
        eases.push_back(prop == (int)RepProp::Visible ? KeyEase::Hold : k->ease);
        unique.push_back(k);
    }

    std::vector<double> values(unique.size());
    for (int c = 0; c < rep_prop_comps(prop); ++c) {
        for (size_t i = 0; i < unique.size(); ++i) values[i] = unique[i]->value[c];
        out[c] = (float)keyed_curve_evaluate(times.data(), values.data(), eases.data(), times.size(), time);
    }
    return true;
}

bool rep_visible_factor(float* out, const RepKey* keys, size_t count, uint32_t rep, double time, double transition) {
    std::vector<const RepKey*> mine;
    for (size_t i = 0; i < count; ++i) {
        if (keys[i].rep == rep && keys[i].prop == (int)RepProp::Visible) mine.push_back(&keys[i]);
    }
    if (mine.empty()) return false;
    std::stable_sort(mine.begin(), mine.end(), [](const RepKey* a, const RepKey* b) { return a->time < b->time; });

    // Keys on the same time are one key
    std::vector<const RepKey*> unique;
    for (const RepKey* k : mine) {
        if (unique.empty() || k->time > unique.back()->time) unique.push_back(k);
    }

    auto smooth = [](double u) { u = std::clamp(u, 0.0, 1.0); return u * u * (3.0 - 2.0 * u); };
    auto target = [](const RepKey* k) { return k->value[0] >= 0.5f ? 1.0 : 0.0; };

    if (time < unique[0]->time) {
        *out = (float)target(unique[0]);
        return true;
    }
    size_t last = 0;
    while (last + 1 < unique.size() && unique[last + 1]->time <= time) ++last;

    // Each change starts from where the one before it had got to by then
    auto ramp = [&](double since) { return transition > 0.0 ? smooth(since / transition) : 1.0; };
    double start = target(unique[0]);
    for (size_t i = 1; i <= last; ++i) {
        start += (target(unique[i - 1]) - start) * ramp(unique[i]->time - unique[i - 1]->time);
    }
    const double end = target(unique[last]);
    *out = (float)(start + (end - start) * ramp(time - unique[last]->time));
    return true;
}

// ## Representations over time as stretches

namespace {
constexpr double REP_EPS = 1.0e-6;
constexpr double REP_MIN_LENGTH = 0.05;

// The indices of the Visible keys of a representation by time, keys on the same time being one (the first)
std::vector<int> visible_key_order(const std::vector<RepKey>& keys, uint32_t rep) {
    std::vector<int> idx;
    for (int i = 0; i < (int)keys.size(); ++i) {
        if (keys[i].rep == rep && keys[i].prop == (int)RepProp::Visible) idx.push_back(i);
    }
    std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return keys[a].time < keys[b].time; });
    std::vector<int> unique;
    for (int i : idx) {
        if (unique.empty() || keys[i].time > keys[unique.back()].time) unique.push_back(i);
    }
    return unique;
}

RepKey hidden_key(uint32_t rep, double time) {
    RepKey k;
    k.rep = rep;
    k.prop = (int)RepProp::Visible;
    k.time = time;
    k.value[0] = 0.0f;
    k.ease = KeyEase::Hold;
    return k;
}

RepKey shown_key(uint32_t rep, double time) {
    RepKey k = hidden_key(rep, time);
    k.value[0] = 1.0f;
    return k;
}
}

std::vector<RepInterval> rep_shown_intervals(const std::vector<RepKey>& keys, uint32_t rep, double duration) {
    std::vector<RepInterval> out;
    const std::vector<int> order = visible_key_order(keys, rep);
    if (order.empty()) return out;
    bool shown = keys[order[0]].value[0] >= 0.5f;
    RepInterval cur;
    if (shown) {
        cur.begin = 0.0;
        cur.begin_key = order[0];
        cur.begin_is_start = true;
    }
    for (size_t i = 1; i < order.size(); ++i) {
        const RepKey& k = keys[order[i]];
        const bool s = k.value[0] >= 0.5f;
        if (s && !shown) {
            cur = RepInterval();
            cur.begin = k.time;
            cur.begin_key = order[i];
            shown = true;
        } else if (!s && shown) {
            cur.end = k.time;
            cur.end_key = order[i];
            out.push_back(cur);
            shown = false;
        }
    }
    if (shown) {
        cur.end = duration;
        cur.end_key = -1;
        out.push_back(cur);
    }
    return out;
}

void rep_move_interval(std::vector<RepKey>* keys, uint32_t rep, const RepInterval& iv, double new_begin, double new_end, double duration) {
    // The room it has between its neighbours and the movie
    double lo = 0.0, hi = duration;
    const std::vector<RepInterval> all = rep_shown_intervals(*keys, rep, duration);
    for (size_t i = 0; i < all.size(); ++i) {
        if (all[i].begin_key != iv.begin_key || all[i].end_key != iv.end_key) continue;
        // A gap is left between stretches: keys on the same time are one key, so they would melt together
        if (i > 0) lo = all[i - 1].end + REP_MIN_LENGTH;
        if (i + 1 < all.size()) hi = all[i + 1].begin - REP_MIN_LENGTH;
    }
    const double length = iv.end - iv.begin;
    const bool only_moved = fabs((new_end - new_begin) - length) < 1.0e-6;
    if (only_moved) {
        if (new_begin < lo) { new_end += lo - new_begin; new_begin = lo; }
        if (new_end > hi) { new_begin -= new_end - hi; new_end = hi; }
    }
    new_begin = std::clamp(new_begin, lo, hi - REP_MIN_LENGTH);
    new_end = std::clamp(new_end, new_begin + REP_MIN_LENGTH, hi);

    if (iv.begin_key >= 0) {
        if (iv.begin_is_start) {
            // Shown from the start by its first key: to start later it needs a hidden key at the start
            if (new_begin > REP_EPS) {
                (*keys)[iv.begin_key].time = new_begin;
                keys->push_back(hidden_key(rep, 0.0));
            }
        } else {
            (*keys)[iv.begin_key].time = new_begin;
        }
    }
    if (iv.end_key >= 0) {
        (*keys)[iv.end_key].time = new_end;
    } else if (new_end < duration - REP_EPS) {
        keys->push_back(hidden_key(rep, new_end));
    }
}

std::vector<RepInterval> rep_effective_intervals(const std::vector<RepKey>& keys, uint32_t rep, double duration, bool enabled) {
    if (!visible_key_order(keys, rep).empty()) return rep_shown_intervals(keys, rep, duration);
    if (!enabled || duration <= 0.0) return {};
    RepInterval interval;
    interval.end = duration;
    interval.begin_is_start = true;
    return {interval};
}

bool rep_add_interval(std::vector<RepKey>* keys, uint32_t rep, double begin, double end, double duration) {
    const std::vector<RepInterval> all = rep_shown_intervals(*keys, rep, duration);
    double next = duration;
    for (const RepInterval& iv : all) {
        if (begin >= iv.begin - REP_MIN_LENGTH && begin < iv.end + REP_MIN_LENGTH) return false;
        if (iv.begin > begin) next = std::min(next, iv.begin - REP_MIN_LENGTH);
    }
    end = std::min(end, next);
    if (end - begin < REP_MIN_LENGTH) return false;
    if (all.empty() && begin > REP_EPS) keys->push_back(hidden_key(rep, 0.0));
    keys->push_back(shown_key(rep, begin));
    if (end < duration - REP_EPS) keys->push_back(hidden_key(rep, end));
    return true;
}

void rep_remove_interval(std::vector<RepKey>* keys, const RepInterval& iv) {
    const uint32_t rep = iv.begin_key >= 0 && iv.begin_key < (int)keys->size() ? (*keys)[iv.begin_key].rep : 0;
    int a = iv.begin_key, b = iv.end_key;
    if (a < b) std::swap(a, b);
    if (a >= 0 && a < (int)keys->size()) keys->erase(keys->begin() + a);
    if (b >= 0 && b < (int)keys->size()) keys->erase(keys->begin() + b);
    if (rep) {
        const auto order = visible_key_order(*keys, rep);
        if (order.empty() || ((*keys)[order.front()].time > REP_EPS && (*keys)[order.front()].value[0] >= 0.5f)) {
            keys->push_back(hidden_key(rep, 0.0));
        }
    }
}

void rep_transfer_interval(std::vector<RepKey>* keys, uint32_t from, uint32_t to, const RepInterval& iv, double duration) {
    if (from == to) return;
    std::vector<RepInterval> intervals = rep_shown_intervals(*keys, to, duration);
    intervals.push_back(iv);
    intervals = rep_union_intervals(intervals);
    rep_remove_interval(keys, iv);
    keys->erase(std::remove_if(keys->begin(), keys->end(), [to](const RepKey& key) {
        return key.rep == to && key.prop == (int)RepProp::Visible;
    }), keys->end());
    if (intervals.front().begin > REP_EPS) keys->push_back(hidden_key(to, 0.0));
    for (const RepInterval& span : intervals) {
        keys->push_back(shown_key(to, span.begin));
        if (span.end < duration - REP_EPS) keys->push_back(hidden_key(to, span.end));
    }
}

int rep_pack_blocks(std::vector<RepBlock>* blocks, double transition) {
    std::vector<double> ends;
    std::vector<size_t> order(blocks->size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [blocks](size_t a, size_t b) {
        return (*blocks)[a].interval.begin < (*blocks)[b].interval.begin;
    });
    for (size_t i : order) {
        RepBlock& block = (*blocks)[i];
        size_t slot = 0;
        while (slot < ends.size() && ends[slot] > block.interval.begin + REP_EPS) ++slot;
        const double end = block.interval.end + std::max(transition, 0.0);
        if (slot == ends.size()) ends.push_back(end);
        else ends[slot] = end;
        block.slot = (int)slot;
    }
    return std::max((int)ends.size(), 1);
}

void rep_swap_at(std::vector<RepKey>* keys, uint32_t from, uint32_t to, double t) {
    auto put = [&](uint32_t rep, bool shown_before, bool shown_after) {
        bool has_keys = false;
        for (RepKey& k : *keys) {
            if (k.rep != rep || k.prop != (int)RepProp::Visible) continue;
            has_keys = true;
            if (fabs(k.time - t) < REP_EPS) {
                k.value[0] = shown_after ? 1.0f : 0.0f;
                k.ease = KeyEase::Hold;
                return;
            }
        }
        if (!has_keys && t > REP_EPS) keys->push_back(shown_before ? shown_key(rep, 0.0) : hidden_key(rep, 0.0));
        keys->push_back(shown_after ? shown_key(rep, t) : hidden_key(rep, t));
    };
    put(from, true, false);
    put(to, false, true);
}

std::vector<RepInterval> rep_union_intervals(std::vector<RepInterval> intervals) {
    std::sort(intervals.begin(), intervals.end(), [](const RepInterval& a, const RepInterval& b) { return a.begin < b.begin; });
    std::vector<RepInterval> out;
    for (const RepInterval& iv : intervals) {
        if (!out.empty() && iv.begin <= out.back().end + REP_EPS) {
            out.back().end = std::max(out.back().end, iv.end);
        } else {
            RepInterval u;
            u.begin = iv.begin;
            u.end = iv.end;
            out.push_back(u);
        }
    }
    return out;
}

bool rep_name_split(const char* name, std::string* group, std::string* member) {
    const char* h = strchr(name, '-');
    if (h && h != name && h[1] != '\0') {
        group->assign(name, (size_t)(h - name));
        member->assign(h + 1);
        return true;
    }
    group->assign(name);
    member->clear();
    return false;
}

std::vector<RepRow> rep_system_rows(const std::vector<std::string>& names) {
    std::vector<RepRow> rows;
    for (size_t i = 0; i < names.size(); ++i) {
        std::string group, member;
        rep_name_split(names[i].c_str(), &group, &member);
        auto it = std::find_if(rows.begin(), rows.end(), [&group](const RepRow& row) { return row.group == group; });
        if (it == rows.end()) {
            RepRow row;
            row.rep = (int)i;
            row.label = group;
            row.group = group;
            rows.push_back(row);
        } else {
            ++it->members;
        }
    }
    return rows;
}

static bool equal(const RepKey& a, const RepKey& b) {
    return a.rep == b.rep && a.prop == b.prop && a.time == b.time && a.ease == b.ease &&
           a.value[0] == b.value[0] && a.value[1] == b.value[1] && a.value[2] == b.value[2];
}

bool movie_keys_equal(const MovieKeys& a, const MovieKeys& b) {
    if (a.sets.size() != b.sets.size()) return false;
    for (size_t i = 0; i < a.sets.size(); ++i) if (a.sets[i].id != b.sets[i].id || a.sets[i].atoms != b.sets[i].atoms || a.sets[i].expr != b.sets[i].expr) return false;
    if (a.loop != b.loop || a.keep_upright != b.keep_upright || a.up_axis != b.up_axis || a.duration != b.duration || a.traj_begin != b.traj_begin || a.traj_end != b.traj_end ||
        a.start_frame != b.start_frame || a.end_frame != b.end_frame || a.camera.size() != b.camera.size() ||
        a.params.size() != b.params.size() || a.reps.size() != b.reps.size() || a.overlays.size() != b.overlays.size()) return false;
    for (size_t i = 0; i < a.camera.size(); ++i) {
        if (!equal(a.camera[i], b.camera[i])) return false;
    }
    for (size_t i = 0; i < a.params.size(); ++i) {
        if (!equal(a.params[i], b.params[i])) return false;
    }
    for (size_t i = 0; i < a.reps.size(); ++i) {
        if (!equal(a.reps[i], b.reps[i])) return false;
    }
    for (size_t i = 0; i < a.overlays.size(); ++i) {
        const MovieOverlay& x = a.overlays[i];
        const MovieOverlay& y = b.overlays[i];
        if (x.type != y.type || x.enabled != y.enabled || x.begin != y.begin || x.end != y.end ||
            x.fade_in != y.fade_in || x.fade_out != y.fade_out || x.anchor != y.anchor ||
            x.num_bins != y.num_bins ||
            x.size != y.size || x.size_unit != y.size_unit || x.width != y.width || x.show_elapsed != y.show_elapsed || x.show_speed != y.show_speed || x.legacy_subplot_mask != y.legacy_subplot_mask || x.plot_axis != y.plot_axis || x.font_points != y.font_points || x.line_points != y.line_points || x.palette != y.palette || x.show_markers != y.show_markers || x.show_titles != y.show_titles || x.panels.size() != y.panels.size() || x.reveal != y.reveal || x.show_value != y.show_value || x.length != y.length || strcmp(x.text, y.text) != 0 || strcmp(x.path, y.path) != 0) return false;
        for (int c = 0; c < 4; ++c) {
            if (x.color[c] != y.color[c] || x.background[c] != y.background[c]) return false;
        }
        for (size_t p = 0; p < x.panels.size(); ++p) {
            if (x.panels[p].view != y.panels[p].view || x.panels[p].subplot != y.panels[p].subplot ||
                x.panels[p].begin != y.panels[p].begin || x.panels[p].end != y.panels[p].end || strcmp(x.panels[p].title, y.panels[p].title) != 0) return false;
        }
    }
    if (a.markers.size() != b.markers.size()) return false;
    for (size_t i = 0; i < a.markers.size(); ++i) {
        if (a.markers[i].time != b.markers[i].time || a.markers[i].subplot != b.markers[i].subplot ||
            strcmp(a.markers[i].label, b.markers[i].label) != 0 || memcmp(a.markers[i].color, b.markers[i].color, sizeof(a.markers[i].color)) != 0) return false;
    }
    return true;
}

void movie_keys_scale_time(MovieKeys* keys, double scale) {
    for (CameraKeyframe& k : keys->camera) {
        k.time *= scale;
        k.focus_transition *= (float)scale;
    }
    for (ParamKey& k : keys->params) k.time *= scale;
    for (RepKey& k : keys->reps) k.time *= scale;
    for (MovieOverlay& o : keys->overlays) {
        o.begin *= scale;
        o.end *= scale;
        o.fade_in *= (float)scale;
        o.fade_out *= (float)scale;
        for (MoviePlotPanel& p : o.panels) {
            p.begin *= scale;
            p.end *= scale;
        }
    }
    for (MovieMarker& k : keys->markers) k.time *= scale;
    keys->duration = (float)(keys->duration * scale);
    keys->traj_begin = (float)(keys->traj_begin * scale);
    keys->traj_end = (float)(keys->traj_end * scale);
}

void MovieHistory::clear(const MovieKeys& current) {
    undo_.clear();
    redo_.clear();
    committed_ = current;
    has_committed_ = true;
}

void MovieHistory::update(const MovieKeys& current, bool editing) {
    if (!has_committed_) {
        committed_ = current;
        has_committed_ = true;
        return;
    }
    if (editing || movie_keys_equal(committed_, current)) return;

    undo_.push_back(committed_);
    if (undo_.size() > MAX_STEPS) undo_.erase(undo_.begin());
    redo_.clear();
    committed_ = current;
}

bool MovieHistory::undo(MovieKeys* state) {
    // An edit that is still under way counts as done before it is undone
    update(*state, false);
    if (undo_.empty()) return false;
    redo_.push_back(*state);
    *state = undo_.back();
    undo_.pop_back();
    committed_ = *state;
    return true;
}

bool MovieHistory::redo(MovieKeys* state) {
    update(*state, false);
    if (redo_.empty()) return false;
    undo_.push_back(*state);
    *state = redo_.back();
    redo_.pop_back();
    committed_ = *state;
    return true;
}

void movie_frame_range(int num_frames, double fps, bool use_range, double begin, double end, int* first, int* last) {
    const int n = num_frames < 1 ? 1 : num_frames;
    *first = 0;
    *last = n - 1;
    if (!use_range) return;
    *first = std::clamp((int)std::ceil(begin * fps - 1.0e-6), 0, n - 1);
    *last  = std::clamp((int)std::floor(end * fps + 1.0e-6), *first, n - 1);
}

void movie_scaled_size(int* w, int* h, int percent) {
    const int s = std::clamp(percent, 10, 100);
    if (s >= 100) return;
    *w = std::max(2, (*w * s / 100) & ~1);
    *h = std::max(2, (*h * s / 100) & ~1);
}

bool movie_time_left(int done, int total, double active_seconds, double* seconds) {
    if (done < 2 || active_seconds < 1.0) return false;
    *seconds = active_seconds / (double)done * (double)std::max(total - done, 0);
    return true;
}

double movie_snap_to_frame(double time, double fps, double duration) {
    if (fps > 0.0) time = std::round(time * fps) / fps;
    return std::clamp(time, 0.0, std::max(duration, 0.0));
}

std::vector<CameraBand> camera_bands(const std::vector<CameraKeyframe>& keys) {
    std::vector<CameraBand> bands;
    const int n = (int)keys.size();
    for (int i = 0; i < n;) {
        if (keys[i].look_set == 0) { ++i; continue; }
        int j = i;
        while (j + 1 < n && keys[j + 1].look_set == keys[i].look_set) ++j;
        CameraBand b;
        b.kind = CameraBandKind::LookAtSet;
        b.first = i;
        b.last = j;
        b.begin = keys[i].time;
        b.end = keys[j].time;
        b.set = keys[i].look_set;
        bands.push_back(b);
        i = j + 1;
    }
    for (int i = 1; i < n; ++i) {
        if (keys[i].spin_turns == 0) continue;
        CameraBand b;
        b.kind = CameraBandKind::Spin;
        b.first = i - 1;
        b.last = i;
        b.begin = keys[i - 1].time;
        b.end = keys[i].time;
        b.turns = keys[i].spin_turns;
        bands.push_back(b);
    }
    return bands;
}

std::string camera_key_label(const CameraKeyframe& key, int index) {
    std::string s = std::to_string(index + 1);
    if (key.name[0] != '\0') s += std::string(" ") + key.name;
    return s;
}

CameraKeyframe camera_key_on_path(const ViewTransform& pose, float fov_y, const std::vector<CameraKeyframe>& keys, double time, const vec3_t* upright) {
    CameraKeyframe key = {};
    key.transform = pose;
    if (upright) key.roll = camera_roll(pose, *upright);
    key.fov_y = fov_y;
    key.time = time;

    const CameraKeyframe* prev = nullptr;
    const CameraKeyframe* next = nullptr;
    for (const CameraKeyframe& k : keys) {
        if (k.time <= time) prev = &k;
        else if (!next) next = &k;
    }
    const CameraKeyframe* ref = prev ? prev : next;
    if (ref && ref->look_set != 0 && (!prev || !next || prev->look_set == next->look_set)) key.look_set = ref->look_set;
    return key;
}

void camera_path_at(const CameraPathSamples& s, double time, vec3_t* eye, vec3_t* look) {
    *eye = vec3_t{0, 0, 0};
    *look = vec3_t{0, 0, 0};
    const size_t n = s.time.size();
    if (n == 0 || s.eye.size() < n || s.look.size() < n) return;
    if (n == 1 || time <= s.time[0]) { *eye = s.eye[0]; *look = s.look[0]; return; }
    if (time >= s.time[n - 1]) { *eye = s.eye[n - 1]; *look = s.look[n - 1]; return; }
    const size_t j = (size_t)(std::upper_bound(s.time.begin(), s.time.end(), time) - s.time.begin());
    const size_t i = j - 1;
    const float u = (float)((time - s.time[i]) / std::max(s.time[j] - s.time[i], 1.0e-12));
    *eye = s.eye[i] + (s.eye[j] - s.eye[i]) * u;
    *look = s.look[i] + (s.look[j] - s.look[i]) * u;
}

double camera_tick_step(double span, int max_ticks) {
    if (span <= 0.0 || max_ticks < 1) return 1.0;
    // In hundredths of a second, so that the steps are exact
    int64_t decade = 1;
    for (int k = 0; k < 10; ++k, decade *= 10) {
        for (int64_t m : {1, 2, 5}) {
            const double step = (double)(m * decade) / 100.0;
            if (span / step <= (double)max_ticks) return step;
        }
    }
    return span;
}

bool clip_segment_near(vec4_t* a, vec4_t* b) {
    const float eps = 1.0e-3f;
    const bool in_a = a->w > eps, in_b = b->w > eps;
    if (!in_a && !in_b) return false;
    if (in_a && in_b) return true;
    vec4_t* out = in_a ? b : a;
    const vec4_t in = in_a ? *a : *b;
    const float t = (in.w - eps) / (in.w - out->w);
    out->x = in.x + (out->x - in.x) * t;
    out->y = in.y + (out->y - in.y) * t;
    out->z = in.z + (out->z - in.z) * t;
    out->w = eps;
    return true;
}

bool camera_spin_ring(const CameraKeyframe& from, const CameraKeyframe& to, vec3_t* center, vec3_t* u, vec3_t* v, float* radius) {
    vec3_t axis = {0, 1, 0};
    switch (to.spin_axis) {
    case SpinAxis::ViewUp: axis = from.transform.orientation * vec3_t{0, 1, 0}; break;
    case SpinAxis::WorldX: axis = {1, 0, 0}; break;
    case SpinAxis::WorldZ: axis = {0, 0, 1}; break;
    default: break;
    }
    axis = vec3_normalize(axis);
    const vec3_t look = camera_get_look_at(to.transform);
    const vec3_t d = to.transform.position - look;
    const float h = vec3_dot(d, axis);
    const vec3_t w = d - axis * h;
    const float r = vec3_length(w);
    if (r < 1.0e-4f) return false;
    *center = look + axis * h;
    *u = w / r;
    *v = vec3_cross(axis, *u);
    *radius = r;
    return true;
}

bool ray_plane_hit(vec3_t origin, vec3_t dir, vec3_t plane_point, vec3_t plane_normal, vec3_t* out) {
    const float denom = vec3_dot(dir, plane_normal);
    if (fabsf(denom) < 1.0e-6f) return false;
    const float t = vec3_dot(plane_point - origin, plane_normal) / denom;
    if (t < 0.0f) return false;
    *out = origin + dir * t;
    return true;
}

float polyline_nearest(const std::vector<vec2_t>& pts, const std::vector<char>& ok, vec2_t q, int* segment, float* along) {
    float best = FLT_MAX;
    for (size_t i = 1; i < pts.size() && i < ok.size(); ++i) {
        if (!ok[i - 1] || !ok[i]) continue;
        const float ax = pts[i - 1].x, ay = pts[i - 1].y;
        const float dx = pts[i].x - ax, dy = pts[i].y - ay;
        const float len2 = dx * dx + dy * dy;
        float u = len2 > 1.0e-12f ? ((q.x - ax) * dx + (q.y - ay) * dy) / len2 : 0.0f;
        u = std::min(std::max(u, 0.0f), 1.0f);
        const float ex = q.x - (ax + dx * u), ey = q.y - (ay + dy * u);
        const float dist = sqrtf(ex * ex + ey * ey);
        if (dist < best) {
            best = dist;
            if (segment) *segment = (int)i - 1;
            if (along) *along = u;
        }
    }
    return best;
}

bool camera_key_set_eye(CameraKeyframe* key, vec3_t eye) {
    const vec3_t look = camera_get_look_at(key->transform);
    ViewTransform t = key->transform;
    t.position = eye;
    if (!camera_aim_at(&t, look)) return false;
    key->transform = t;
    return true;
}

bool camera_key_set_look(CameraKeyframe* key, vec3_t look) {
    ViewTransform t = key->transform;
    if (!camera_aim_at(&t, look)) return false;
    key->transform = t;
    return true;
}

void camera_key_translate(CameraKeyframe* key, vec3_t delta) {
    key->transform.position = key->transform.position + delta;
    if (key->follow) key->follow_center = key->follow_center + delta;
    if (key->focus_on && key->focus_target == FocusTarget::Point) key->focus_point = key->focus_point + delta;
}

static bool same_key(const KeyId& k, KeyKind kind, int64_t subject, double time) {
    return k.kind == kind && k.subject == subject && fabs(k.time - time) < 1.0e-9;
}

bool KeySelection::contains(KeyKind kind, int64_t subject, double time) const {
    for (const KeyId& k : ids) {
        if (same_key(k, kind, subject, time)) return true;
    }
    return false;
}

void KeySelection::add(KeyKind kind, int64_t subject, double time, double end) {
    if (!contains(kind, subject, time)) ids.push_back({kind, subject, time, end});
}

void KeySelection::toggle(KeyKind kind, int64_t subject, double time, double end) {
    for (size_t i = 0; i < ids.size(); ++i) {
        if (same_key(ids[i], kind, subject, time)) {
            ids.erase(ids.begin() + (ptrdiff_t)i);
            return;
        }
    }
    ids.push_back({kind, subject, time, end});
}

void KeySelection::set(KeyKind kind, int64_t subject, double time, double end) {
    ids.clear();
    ids.push_back({kind, subject, time, end});
}

static bool has_key(const MovieKeys& keys, const KeyId& id) {
    switch (id.kind) {
    case KeyKind::Camera:
        for (const CameraKeyframe& k : keys.camera) if (fabs(k.time - id.time) < 1.0e-9) return true;
        return false;
    case KeyKind::Param:
        for (const ParamKey& k : keys.params) if ((int64_t)k.param == id.subject && fabs(k.time - id.time) < 1.0e-9) return true;
        return false;
    case KeyKind::Rep:
        for (const RepKey& k : keys.reps) if (rep_key_subject(k.rep, k.prop) == id.subject && fabs(k.time - id.time) < 1.0e-9) return true;
        return false;
    case KeyKind::Overlay:
        return id.subject >= 0 && (size_t)id.subject < keys.overlays.size() && fabs(keys.overlays[(size_t)id.subject].begin - id.time) < 1.0e-9;
    case KeyKind::Block:
        for (const RepInterval& iv : rep_shown_intervals(keys.reps, (uint32_t)id.subject, (double)keys.duration)) {
            if (iv.begin_key >= 0 && fabs(iv.begin - id.time) < 1.0e-9) return true;
        }
        return false;
    }
    return false;
}

void key_selection_prune(KeySelection* sel, const MovieKeys& keys) {
    sel->ids.erase(std::remove_if(sel->ids.begin(), sel->ids.end(), [&](const KeyId& id) { return !has_key(keys, id); }), sel->ids.end());
}

void key_selection_all(KeySelection* sel, const MovieKeys& keys, int64_t param, int64_t rep_subject) {
    for (const CameraKeyframe& k : keys.camera) sel->add(KeyKind::Camera, 0, k.time);
    if (param >= 0) {
        for (const ParamKey& k : keys.params) if ((int64_t)k.param == param) sel->add(KeyKind::Param, param, k.time);
    }
    if (rep_subject >= 0) {
        for (const RepKey& k : keys.reps) if (rep_key_subject(k.rep, k.prop) == rep_subject) sel->add(KeyKind::Rep, rep_subject, k.time);
    }
}

static double shifted_value(double v, const KeyShift& s) {
    return s.ratio ? v * s.dy : v + s.dy;
}

struct BlockMove { size_t id; uint32_t rep; RepInterval iv; };

// The first and the last time of what is selected (nothing selected: t_min > t_max), and the stretches of the selected blocks
static void selection_extent(const MovieKeys& keys, const KeySelection& sel, double duration, double* t_min, double* t_max, std::vector<BlockMove>* blocks) {
    *t_min = DBL_MAX;
    *t_max = -DBL_MAX;
    auto take = [&](double a, double b) { *t_min = std::min(*t_min, a); *t_max = std::max(*t_max, b); };
    for (const CameraKeyframe& k : keys.camera) if (sel.contains(KeyKind::Camera, 0, k.time)) take(k.time, k.time);
    for (const ParamKey& k : keys.params) if (sel.contains(KeyKind::Param, (int64_t)k.param, k.time)) take(k.time, k.time);
    for (const RepKey& k : keys.reps) if (sel.contains(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time)) take(k.time, k.time);
    for (size_t i = 0; i < sel.ids.size(); ++i) {
        const KeyId& id = sel.ids[i];
        if (id.kind == KeyKind::Overlay && id.subject >= 0 && (size_t)id.subject < keys.overlays.size()) {
            const MovieOverlay& o = keys.overlays[(size_t)id.subject];
            if (fabs(o.begin - id.time) < 1.0e-9) take(o.begin, o.end);
        } else if (id.kind == KeyKind::Block) {
            for (const RepInterval& iv : rep_shown_intervals(keys.reps, (uint32_t)id.subject, duration)) {
                if (iv.begin_key < 0 || fabs(iv.begin - id.time) >= 1.0e-9) continue;
                if (blocks) blocks->push_back({i, (uint32_t)id.subject, iv});
                take(iv.begin, iv.end);
            }
        }
    }
}

bool key_selection_extent(const MovieKeys& keys, const KeySelection& sel, double* t0, double* t1) {
    selection_extent(keys, sel, (double)keys.duration, t0, t1, nullptr);
    return *t0 <= *t1;
}

double movie_keys_shift(MovieKeys* keys, KeySelection* sel, const MovieKeys& start, const KeySelection& start_sel, double dt, const KeyShift& s, double duration) {
    std::vector<CameraKeyframe> camera = start.camera;
    std::vector<ParamKey> params = start.params;
    std::vector<RepKey> reps = start.reps;
    std::vector<MovieOverlay> overlays = start.overlays;

    auto cam_sel = [&](const CameraKeyframe& k) { return start_sel.contains(KeyKind::Camera, 0, k.time); };
    auto par_sel = [&](const ParamKey& k) { return start_sel.contains(KeyKind::Param, (int64_t)k.param, k.time); };
    auto rep_sel = [&](const RepKey& k) { return start_sel.contains(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time); };

    double t_min = 0.0, t_max = 0.0;
    std::vector<BlockMove> blocks;
    selection_extent(start, start_sel, duration, &t_min, &t_max, &blocks);

    double d = 0.0;
    if (t_min <= t_max) {
        const double lo_d = -t_min;
        const double hi_d = std::max(duration - t_max, lo_d);
        d = std::min(std::max(dt, lo_d), hi_d);
    }

    for (CameraKeyframe& k : camera) {
        if (!cam_sel(k)) continue;
        k.time += d;
        if (s.lane == KeyLane::Frame && k.use_frame) {
            k.frame = std::min(std::max(shifted_value(k.frame, s), s.lo), s.hi);
        } else if (s.lane == KeyLane::Distance) {
            const double unit = s.unit > 0.0 ? s.unit : 1.0;
            const double shown = shifted_value((double)k.transform.distance * unit, s);
            const float dist = (float)std::max(shown / unit, 0.01);
            const vec3_t look = camera_get_look_at(k.transform);
            k.transform.distance = dist;
            k.transform.position = camera_position_from_look_at(look, k.transform.orientation, dist);
        } else if (s.lane == KeyLane::Fov) {
            const double deg = shifted_value((double)k.fov_y * 57.29577951308232, s);
            k.fov_y = (float)(std::min(std::max(deg, 1.0), 170.0) * 0.017453292519943295);
        }
    }
    for (ParamKey& k : params) {
        if (!par_sel(k)) continue;
        k.time += d;
        if (s.lane == KeyLane::Param && s.subject == (int64_t)k.param) k.value[0] = (float)std::min(std::max(shifted_value((double)k.value[0], s), s.lo), s.hi);
    }
    for (RepKey& k : reps) {
        if (!rep_sel(k)) continue;
        k.time += d;
        if (s.lane == KeyLane::Rep && s.subject == rep_key_subject(k.rep, k.prop)) k.value[0] = (float)std::min(std::max(shifted_value((double)k.value[0], s), s.lo), s.hi);
    }

    for (size_t i = 0; i < overlays.size(); ++i) {
        MovieOverlay& o = overlays[i];
        if (!start_sel.contains(KeyKind::Overlay, (int64_t)i, o.begin)) continue;
        o.begin += d;
        o.end += d;
        // What is timed inside a bar moves with it
        for (MoviePlotPanel& panel : o.panels) {
            panel.begin = std::max(panel.begin + d, 0.0);
            if (panel.end > 0.0) panel.end = std::max(panel.end + d, 0.0);
        }
    }
    // The blocks one after the other, the ones ahead of the move first so that they make room
    std::sort(blocks.begin(), blocks.end(), [&](const BlockMove& a, const BlockMove& b) { return d > 0.0 ? a.iv.begin > b.iv.begin : a.iv.begin < b.iv.begin; });
    for (const BlockMove& b : blocks) rep_move_interval(&reps, b.rep, b.iv, b.iv.begin + d, b.iv.end + d, duration);

    std::vector<KeyId> ids = start_sel.ids;
    for (KeyId& id : ids) {
        if (id.kind == KeyKind::Block) continue;
        id.time += d;
        if (id.kind == KeyKind::Overlay) id.end += d;
    }
    for (const BlockMove& b : blocks) {
        for (const RepInterval& iv : rep_shown_intervals(reps, b.rep, duration)) {
            if (iv.begin_key == b.iv.begin_key) { ids[b.id].time = iv.begin; ids[b.id].end = iv.end; }
        }
    }
    keys->overlays = std::move(overlays);
    keys->camera = std::move(camera);
    keys->params = std::move(params);
    keys->reps = std::move(reps);
    sel->ids = std::move(ids);
    return d;
}

// Of keys within a millisecond of each other (same thing, sorted) one is kept: the selected one
template <typename T, typename Same, typename IsSel>
static void keep_one_per_time(std::vector<T>& v, Same same, IsSel is_sel) {
    std::vector<T> out;
    size_t i = 0;
    while (i < v.size()) {
        size_t keep = i, j = i + 1;
        bool keep_sel = is_sel(v[i]);
        while (j < v.size() && same(v[i], v[j]) && v[j].time - v[j - 1].time < 1.0e-3) {
            if (!keep_sel && is_sel(v[j])) { keep = j; keep_sel = true; }
            ++j;
        }
        out.push_back(v[keep]);
        i = j;
    }
    v = std::move(out);
}

void movie_keys_resolve(MovieKeys* keys, const KeySelection& sel) {
    std::stable_sort(keys->camera.begin(), keys->camera.end(), [](const CameraKeyframe& a, const CameraKeyframe& b) { return a.time < b.time; });
    keep_one_per_time(keys->camera, [](const CameraKeyframe&, const CameraKeyframe&) { return true; },
        [&](const CameraKeyframe& k) { return sel.contains(KeyKind::Camera, 0, k.time); });

    std::stable_sort(keys->params.begin(), keys->params.end(), [](const ParamKey& a, const ParamKey& b) {
        return a.param != b.param ? a.param < b.param : a.time < b.time;
    });
    keep_one_per_time(keys->params, [](const ParamKey& a, const ParamKey& b) { return a.param == b.param; },
        [&](const ParamKey& k) { return sel.contains(KeyKind::Param, (int64_t)k.param, k.time); });

    std::stable_sort(keys->reps.begin(), keys->reps.end(), [](const RepKey& a, const RepKey& b) {
        if (a.rep != b.rep) return a.rep < b.rep;
        return a.prop != b.prop ? a.prop < b.prop : a.time < b.time;
    });
    keep_one_per_time(keys->reps, [](const RepKey& a, const RepKey& b) { return a.rep == b.rep && a.prop == b.prop; },
        [&](const RepKey& k) { return sel.contains(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time); });
}

void movie_keys_delete(MovieKeys* keys, KeySelection* sel) {
    for (const KeyId& id : sel->ids) {
        if (id.kind != KeyKind::Block) continue;
        for (const RepInterval& iv : rep_shown_intervals(keys->reps, (uint32_t)id.subject, (double)keys->duration)) {
            if (iv.begin_key >= 0 && fabs(iv.begin - id.time) < 1.0e-9) {
                rep_remove_interval(&keys->reps, iv);
                break;
            }
        }
    }
    keys->camera.erase(std::remove_if(keys->camera.begin(), keys->camera.end(),
        [&](const CameraKeyframe& k) { return sel->contains(KeyKind::Camera, 0, k.time); }), keys->camera.end());
    keys->params.erase(std::remove_if(keys->params.begin(), keys->params.end(),
        [&](const ParamKey& k) { return sel->contains(KeyKind::Param, (int64_t)k.param, k.time); }), keys->params.end());
    keys->reps.erase(std::remove_if(keys->reps.begin(), keys->reps.end(),
        [&](const RepKey& k) { return sel->contains(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time); }), keys->reps.end());
    sel->clear();
}

KeyClip movie_keys_copy(const MovieKeys& keys, const KeySelection& sel) {
    KeyClip clip;
    double begin = DBL_MAX;
    for (const CameraKeyframe& k : keys.camera) if (sel.contains(KeyKind::Camera, 0, k.time)) { clip.camera.push_back(k); begin = std::min(begin, k.time); }
    for (const ParamKey& k : keys.params) if (sel.contains(KeyKind::Param, (int64_t)k.param, k.time)) { clip.params.push_back(k); begin = std::min(begin, k.time); }
    for (const RepKey& k : keys.reps) if (sel.contains(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time)) { clip.reps.push_back(k); begin = std::min(begin, k.time); }
    clip.begin = clip.empty() ? 0.0 : begin;
    return clip;
}

void movie_keys_paste(MovieKeys* keys, KeySelection* sel, const KeyClip& clip, double time, double duration) {
    if (clip.empty()) return;
    double last = clip.begin;
    for (const CameraKeyframe& k : clip.camera) last = std::max(last, k.time);
    for (const ParamKey& k : clip.params) last = std::max(last, k.time);
    for (const RepKey& k : clip.reps) last = std::max(last, k.time);
    const double offset = std::max(std::min(time, duration - (last - clip.begin)), 0.0) - clip.begin;

    sel->clear();
    for (CameraKeyframe k : clip.camera) { k.time += offset; keys->camera.push_back(k); sel->add(KeyKind::Camera, 0, k.time); }
    for (ParamKey k : clip.params) { k.time += offset; keys->params.push_back(k); sel->add(KeyKind::Param, (int64_t)k.param, k.time); }
    for (RepKey k : clip.reps) { k.time += offset; keys->reps.push_back(k); sel->add(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time); }
    movie_keys_resolve(keys, *sel);
}

double movie_keys_scale(MovieKeys* keys, KeySelection* sel, const MovieKeys& start, const KeySelection& start_sel, double anchor, double factor, double duration) {
    std::vector<CameraKeyframe> camera = start.camera;
    std::vector<ParamKey> params = start.params;
    std::vector<RepKey> reps = start.reps;
    std::vector<MovieOverlay> overlays = start.overlays;

    double t_min = 0.0, t_max = 0.0;
    std::vector<BlockMove> blocks;
    selection_extent(start, start_sel, duration, &t_min, &t_max, &blocks);

    double f = 1.0;
    if (t_min <= t_max) {
        f = std::max(factor, 0.01);
        if (t_max > anchor) f = std::min(f, (duration - anchor) / (t_max - anchor));
        if (t_min < anchor) f = std::min(f, anchor / (anchor - t_min));
        f = std::max(f, 0.0);
    }
    auto at = [&](double t) { return anchor + (t - anchor) * f; };

    for (CameraKeyframe& k : camera) if (start_sel.contains(KeyKind::Camera, 0, k.time)) { k.time = at(k.time); k.focus_transition *= (float)f; }
    for (ParamKey& k : params) if (start_sel.contains(KeyKind::Param, (int64_t)k.param, k.time)) k.time = at(k.time);
    for (RepKey& k : reps) if (start_sel.contains(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time)) k.time = at(k.time);
    for (size_t i = 0; i < overlays.size(); ++i) {
        MovieOverlay& o = overlays[i];
        if (!start_sel.contains(KeyKind::Overlay, (int64_t)i, o.begin)) continue;
        o.begin = at(o.begin);
        o.end = at(o.end);
        for (MoviePlotPanel& panel : o.panels) {
            panel.begin = std::max(at(panel.begin), 0.0);
            if (panel.end > 0.0) panel.end = std::max(at(panel.end), 0.0);
        }
    }
    // The blocks that move away from the anchor first, so that they make room
    std::sort(blocks.begin(), blocks.end(), [&](const BlockMove& a, const BlockMove& b) { return f >= 1.0 ? a.iv.begin > b.iv.begin : a.iv.begin < b.iv.begin; });
    for (const BlockMove& b : blocks) rep_move_interval(&reps, b.rep, b.iv, at(b.iv.begin), at(b.iv.end), duration);

    std::vector<KeyId> ids = start_sel.ids;
    for (KeyId& id : ids) {
        if (id.kind == KeyKind::Block) continue;
        id.time = at(id.time);
        if (id.kind == KeyKind::Overlay) id.end = at(id.end);
    }
    for (const BlockMove& b : blocks) {
        for (const RepInterval& iv : rep_shown_intervals(reps, b.rep, duration)) {
            if (iv.begin_key == b.iv.begin_key) { ids[b.id].time = iv.begin; ids[b.id].end = iv.end; }
        }
    }
    keys->overlays = std::move(overlays);
    keys->camera = std::move(camera);
    keys->params = std::move(params);
    keys->reps = std::move(reps);
    sel->ids = std::move(ids);
    return f;
}

int key_selection_ease(const CameraKeyframe* camera, size_t num_camera, const std::vector<ParamKey>& params, const std::vector<RepKey>& reps,
    const KeySelection& sel, KeyEase* common, bool* mixed) {
    double first = DBL_MAX;
    for (size_t i = 0; i < num_camera; ++i) first = std::min(first, camera[i].time);
    int count = 0;
    *mixed = false;
    auto take = [&](KeyEase e) {
        if (count == 0) *common = e;
        else if (e != *common) *mixed = true;
        ++count;
    };
    for (size_t i = 0; i < num_camera; ++i) if (camera[i].time > first && sel.contains(KeyKind::Camera, 0, camera[i].time)) take(camera[i].ease);
    for (const ParamKey& k : params) if (sel.contains(KeyKind::Param, (int64_t)k.param, k.time)) take(k.ease);
    for (const RepKey& k : reps) if (k.prop != (int)RepProp::Visible && sel.contains(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time)) take(k.ease);
    return count;
}

void movie_keys_set_ease(MovieKeys* keys, const KeySelection& sel, KeyEase ease) {
    double first = DBL_MAX;
    for (const CameraKeyframe& k : keys->camera) first = std::min(first, k.time);
    for (CameraKeyframe& k : keys->camera) if (k.time > first && sel.contains(KeyKind::Camera, 0, k.time)) k.ease = ease;
    for (ParamKey& k : keys->params) if (sel.contains(KeyKind::Param, (int64_t)k.param, k.time)) k.ease = ease;
    for (RepKey& k : keys->reps) if (k.prop != (int)RepProp::Visible && sel.contains(KeyKind::Rep, rep_key_subject(k.rep, k.prop), k.time)) k.ease = ease;
}
