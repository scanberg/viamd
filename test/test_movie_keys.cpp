#include "utest.h"

#include <gfx/camera_utils.h>
#include <movie_keys.h>
#include <md_system.h>
#include <md_unitcell.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitfield.h>
#include <md_xyz.h>
#include <serialization_utils.h>

#include <float.h>
#include <math.h>
#include <filesystem>
#include <fstream>
#include <iterator>

UTEST(viamd_movie_keys, independent_position_does_not_depend_on_legacy_aim) {
    CameraKeyframe keys[2];
    keys[0].time = 0; keys[1].time = 10;
    keys[0].transform.position = vec3_set(0, 0, 40);
    keys[1].transform.position = vec3_set(10, 0, 40);
    keys[1].ease = KeyEase::Linear;
    keys[1].transform.orientation = quat_axis_angle(vec3_set(0, 1, 0), 1.5f);
    keys[1].transform.distance = 200;
    vec3_t eye = movie_position_evaluate(keys, 2, 5, false);
    EXPECT_NEAR(5.0f, eye.x, 1.0e-5f);
    EXPECT_NEAR(40.0f, eye.z, 1.0e-5f);
    keys[0].transform.distance = 999;
    keys[1].transform.orientation = quat_axis_angle(vec3_set(1, 0, 0), 1.0f);
    vec3_t unchanged = movie_position_evaluate(keys, 2, 5, false);
    EXPECT_NEAR(0.0f, vec3_length(eye - unchanged), 1.0e-5f);
    EXPECT_NEAR(10.0f, movie_position_evaluate(keys, 2, 10, false).x, 1.0e-5f);
}

static CameraKeyframe plain_key(double time, vec3_t eye, float distance) {
    CameraKeyframe key;
    key.time = time;
    key.transform.orientation = quat_t{0, 0, 0, 1};
    key.transform.position = eye;
    key.transform.distance = distance;
    return key;
}

UTEST(viamd_movie_keys, atom_sets_are_shared_sorted_and_pruned_when_no_key_uses_them) {
    std::vector<MovieAtomSet> sets;
    const uint32_t a = movie_atoms_add(&sets, {4, 1, 4});
    EXPECT_EQ(a, movie_atoms_add(&sets, {1, 4}));
    ASSERT_EQ((size_t)1, sets.size());
    ASSERT_EQ((size_t)2, sets[0].atoms.size());
    EXPECT_EQ((uint32_t)1, sets[0].atoms[0]);
    EXPECT_EQ((uint32_t)4, sets[0].atoms[1]);
    const uint32_t b = movie_atoms_add(&sets, {7});
    EXPECT_NE(a, b);
    EXPECT_EQ((uint32_t)0, movie_atoms_add(&sets, {}));
    EXPECT_TRUE(movie_atoms_find(sets, 0) == nullptr);
    ASSERT_TRUE(movie_atoms_find(sets, b) != nullptr);

    CameraKeyframe keys[2];
    keys[0].look_set = a;
    keys[1].focus_on = true;
    keys[1].focus_target = FocusTarget::Selection;
    keys[1].focus_set = b;
    movie_atoms_prune(&sets, keys, 2);
    EXPECT_EQ((size_t)2, sets.size());
    keys[1].focus_target = FocusTarget::LookAt;   // The set is no longer the focus
    movie_atoms_prune(&sets, keys, 2);
    ASSERT_EQ((size_t)1, sets.size());
    EXPECT_EQ(a, sets[0].id);
}

UTEST(viamd_movie_keys, look_at_is_keyed_through_the_points_the_keys_aim_at) {
    CameraKeyframe keys[2];
    keys[0].time = 2; keys[1].time = 6;
    keys[1].ease = KeyEase::Linear;
    vec3_t looks[2] = {vec3_set(0, 0, 0), vec3_set(8, 0, 0)};
    EXPECT_NEAR(0.0f, movie_look_evaluate(keys, 2, looks, 0).x, 1.0e-5f);
    EXPECT_NEAR(4.0f, movie_look_evaluate(keys, 2, looks, 4).x, 1.0e-5f);
    looks[1].x = 12;   // The atoms the second key tracks have moved
    EXPECT_NEAR(6.0f, movie_look_evaluate(keys, 2, looks, 4).x, 1.0e-5f);
    EXPECT_NEAR(12.0f, movie_look_evaluate(keys, 2, looks, 20).x, 1.0e-5f);
}

UTEST(viamd_movie_keys, the_eye_goes_through_the_keys_and_the_aim_through_their_look_at) {
    CameraKeyframe keys[2] = {plain_key(0, vec3_set(0, 0, 40), 32), plain_key(10, vec3_set(10, 0, 40), 32)};
    keys[1].ease = KeyEase::Linear;
    vec3_t looks[2] = {vec3_set(0, 0, 8), vec3_set(10, 0, 8)};
    ViewTransform view;
    float fov;
    movie_keys_pose(&view, &fov, keys, 2, 5, looks, nullptr, vec3_set(0, 1, 0));
    EXPECT_NEAR(5.0f, view.position.x, 1.0e-4f);
    EXPECT_NEAR(40.0f, view.position.z, 1.0e-4f);
    EXPECT_NEAR(5.0f, camera_get_look_at(view).x, 1.0e-3f);
    EXPECT_NEAR(8.0f, camera_get_look_at(view).z, 1.0e-3f);
    // Where the second key looks decides the aim, not the eye
    looks[1] = vec3_set(10, 0, 0);
    movie_keys_pose(&view, &fov, keys, 2, 10, looks, nullptr, vec3_set(0, 1, 0));
    EXPECT_NEAR(10.0f, view.position.x, 1.0e-4f);
    EXPECT_NEAR(40.0f, view.position.z, 1.0e-4f);
    EXPECT_NEAR(0.0f, camera_get_look_at(view).z, 1.0e-3f);
}

static ViewTransform front_camera() {
    ViewTransform camera = {};
    camera.orientation = quat_t{0, 0, 0, 1};
    camera.position = vec3_set(0, 0, 40);
    camera.distance = 32;
    return camera;
}

UTEST(viamd_movie_keys, the_focus_is_what_the_camera_looks_at_until_a_key_sets_it) {
    const ViewTransform camera = front_camera();
    CameraKeyframe keys[2];
    EXPECT_FALSE(movie_focus_keys_exist(keys, 2));
    EXPECT_NEAR(32.0f, movie_focus_depth(keys, 2, nullptr, 3, camera), 1.0e-5f);

    keys[1].time = 4;
    keys[1].focus_on = true;
    keys[1].focus_target = FocusTarget::Distance;
    keys[1].focus_distance = 48;
    keys[1].focus_transition = 2;
    keys[1].focus_ease = KeyEase::Linear;
    EXPECT_TRUE(movie_focus_keys_exist(keys, 2));
    EXPECT_NEAR(32.0f, movie_focus_depth(keys, 2, nullptr, 4, camera), 1.0e-5f);
    EXPECT_NEAR(40.0f, movie_focus_depth(keys, 2, nullptr, 5, camera), 1.0e-5f);
    EXPECT_NEAR(48.0f, movie_focus_depth(keys, 2, nullptr, 9, camera), 1.0e-5f);
}

UTEST(viamd_movie_keys, default_focus_returns_to_look_at_after_an_explicit_target) {
    ViewTransform camera = front_camera();
    CameraKeyframe keys[3];
    keys[0].time = 0;
    keys[0].focus_on = true;
    keys[0].focus_target = FocusTarget::Selection;
    keys[0].focus_transition = 0;
    keys[1].time = 5;
    keys[1].focus_on = false;
    keys[2].time = 10;
    keys[2].focus_on = true;
    keys[2].focus_target = FocusTarget::Distance;
    keys[2].focus_distance = 48;
    keys[2].focus_transition = 2;
    keys[2].focus_ease = KeyEase::Linear;
    vec3_t points[3] = {vec3_set(0, 0, -8), {}, {}};
    EXPECT_NEAR(48.0f, movie_focus_depth(keys, 3, points, 4, camera), 1.0e-5f);
    EXPECT_NEAR(32.0f, movie_focus_depth(keys, 3, points, 5, camera), 1.0e-5f);
    camera.distance = 24;
    EXPECT_NEAR(24.0f, movie_focus_depth(keys, 3, points, 7, camera), 1.0e-5f);
    EXPECT_NEAR(36.0f, movie_focus_depth(keys, 3, points, 11, camera), 1.0e-5f);
    EXPECT_NEAR(48.0f, movie_focus_depth(keys, 3, points, 12, camera), 1.0e-5f);
    EXPECT_NEAR(48.0f, movie_focus_depth(keys, 3, points, 4, camera), 1.0e-5f);
    EXPECT_NEAR(24.0f, movie_focus_depth(keys, 3, points, 5, camera), 1.0e-5f);
}

UTEST(viamd_movie_keys, default_focus_interrupts_a_transition_and_holds_look_at) {
    ViewTransform camera = front_camera();
    CameraKeyframe keys[3];
    keys[0].time = 0;
    keys[1].time = 10;
    keys[1].focus_on = true;
    keys[1].focus_target = FocusTarget::Distance;
    keys[1].focus_distance = 48;
    keys[1].focus_transition = 8;
    keys[1].focus_ease = KeyEase::Linear;
    keys[2].time = 12;
    EXPECT_NEAR(40.0f, movie_focus_depth(keys, 3, nullptr, 11, camera), 1.0e-5f);
    EXPECT_NEAR(32.0f, movie_focus_depth(keys, 3, nullptr, 12, camera), 1.0e-5f);
    EXPECT_NEAR(32.0f, movie_focus_depth(keys, 3, nullptr, 40, camera), 1.0e-5f);
    camera.distance = 20;
    EXPECT_NEAR(20.0f, movie_focus_depth(keys, 3, nullptr, 40, camera), 1.0e-5f);
    EXPECT_NEAR(20.0f, movie_focus_depth(nullptr, 0, nullptr, 0, camera), 1.0e-5f);
}

UTEST(viamd_movie_keys, the_focus_switches_between_targets_and_follows_atoms) {
    const ViewTransform camera = front_camera();
    CameraKeyframe keys[3];
    keys[0].focus_on = true;
    keys[0].focus_target = FocusTarget::LookAt;
    keys[0].focus_transition = 0;
    keys[1].time = 10;
    keys[1].focus_on = true;
    keys[1].focus_target = FocusTarget::Selection;
    keys[1].focus_transition = 2;
    keys[1].focus_ease = KeyEase::Linear;
    keys[2].time = 15;
    keys[2].focus_on = true;
    keys[2].focus_target = FocusTarget::Distance;
    keys[2].focus_distance = 48;
    keys[2].focus_transition = 0;
    vec3_t points[3] = {{}, vec3_set(4, 0, -8), {}};   // Where the atoms of the second key are
    EXPECT_NEAR(32.0f, movie_focus_depth(keys, 3, points, 0, camera), 1.0e-5f);
    EXPECT_NEAR(32.0f, movie_focus_depth(keys, 3, points, 10, camera), 1.0e-5f);
    EXPECT_NEAR(40.0f, movie_focus_depth(keys, 3, points, 11, camera), 1.0e-5f);
    EXPECT_NEAR(48.0f, movie_focus_depth(keys, 3, points, 12, camera), 1.0e-5f);
    points[1].z = 0;   // The atoms moved
    EXPECT_NEAR(40.0f, movie_focus_depth(keys, 3, points, 14, camera), 1.0e-5f);
    EXPECT_NEAR(48.0f, movie_focus_depth(keys, 3, points, 15, camera), 1.0e-5f);
    EXPECT_NEAR(40.0f, camera.position.z, 1.0e-5f);   // The camera itself is not touched
    EXPECT_NEAR(32.0f, camera.distance, 1.0e-5f);
}

UTEST(viamd_movie_keys, a_focus_transition_finishes_before_the_next_key_that_sets_the_focus) {
    const ViewTransform camera = front_camera();
    CameraKeyframe keys[3];
    keys[0].focus_on = true;
    keys[0].focus_target = FocusTarget::LookAt;
    keys[0].focus_transition = 0;
    keys[1].time = 10; keys[1].focus_on = true; keys[1].focus_target = FocusTarget::Point;
    keys[1].focus_transition = 8; keys[1].focus_ease = KeyEase::Linear;
    keys[2].time = 12; keys[2].focus_on = true; keys[2].focus_target = FocusTarget::Point;
    keys[2].focus_transition = 2; keys[2].focus_ease = KeyEase::Linear;
    vec3_t points[3] = {vec3_set(0, 0, 8), vec3_set(0, 0, 0), vec3_set(0, 0, 8)};
    // The first key's change would take until 18, but the next key that sets the focus is at 12: it is done by then
    EXPECT_NEAR(36.0f, movie_focus_depth(keys, 3, points, 11, camera), 1.0e-5f);
    EXPECT_NEAR(40.0f, movie_focus_depth(keys, 3, points, 12, camera), 1.0e-5f);
    EXPECT_NEAR(36.0f, movie_focus_depth(keys, 3, points, 13, camera), 1.0e-5f);
}

UTEST(viamd_movie_keys, a_set_can_be_a_selection_expression_that_is_shared_and_trimmed) {
    std::vector<MovieAtomSet> sets;
    const uint32_t a = movie_atoms_add_expr(&sets, "  resname(\"AIN\") ");
    EXPECT_EQ(a, movie_atoms_add_expr(&sets, "resname(\"AIN\")"));
    EXPECT_EQ((uint32_t)0, movie_atoms_add_expr(&sets, "   "));
    const uint32_t list = movie_atoms_add(&sets, {1, 2});
    EXPECT_NE(a, list);
    EXPECT_EQ(list, movie_atoms_add(&sets, {2, 1}));
    ASSERT_EQ((size_t)2, sets.size());
    ASSERT_TRUE(movie_atoms_find(sets, a) != nullptr);
    EXPECT_STREQ("resname(\"AIN\")", movie_atoms_find(sets, a)->expr.c_str());
    MovieKeys x, y;
    x.sets = sets; y.sets = sets;
    EXPECT_TRUE(movie_keys_equal(x, y));
    y.sets[0].expr = "protein";
    EXPECT_FALSE(movie_keys_equal(x, y));
}

UTEST(viamd_movie_keys, the_blur_of_focus_keys_becomes_keys_of_the_blur_parameter_with_the_same_transition) {
    std::vector<CameraKeyframe> keys(4);
    for (int i = 0; i < 4; ++i) keys[i].time = 5.0 * i;
    keys[0].focus_on = true; keys[0].focus_blur = 0.0f; keys[0].focus_transition = 0;
    keys[1].focus_on = true; keys[1].focus_blur = 3.0f; keys[1].focus_transition = 2; keys[1].focus_ease = KeyEase::Linear;
    keys[2].focus_on = true;   // Sets the target only: no blur of its own
    keys[3].focus_on = true; keys[3].focus_blur = 3.0f; keys[3].focus_transition = 0;
    const std::vector<ParamKey> blur = movie_blur_keys_from_focus(keys.data(), keys.size());
    float v = -1.0f;
    ASSERT_TRUE(param_keys_evaluate(&v, 1, blur.data(), blur.size(), MOVIE_PARAM_DOF_APERTURE, 0));
    EXPECT_NEAR(0.0f, v, 1.0e-5f);
    param_keys_evaluate(&v, 1, blur.data(), blur.size(), MOVIE_PARAM_DOF_APERTURE, 5);
    EXPECT_NEAR(0.0f, v, 1.0e-5f);
    param_keys_evaluate(&v, 1, blur.data(), blur.size(), MOVIE_PARAM_DOF_APERTURE, 6);
    EXPECT_NEAR(1.5f, v, 1.0e-4f);
    param_keys_evaluate(&v, 1, blur.data(), blur.size(), MOVIE_PARAM_DOF_APERTURE, 12);
    EXPECT_NEAR(3.0f, v, 1.0e-5f);
    for (const CameraKeyframe& k : keys) EXPECT_TRUE(k.focus_blur < 0.0f);
    EXPECT_TRUE(movie_blur_keys_from_focus(keys.data(), keys.size()).empty());
}

UTEST(viamd_movie_keys, camera_bands_follow_the_sets_a_run_of_keys_look_at_and_their_spins) {
    std::vector<CameraKeyframe> keys(4);
    for (int i = 0; i < 4; ++i) keys[i].time = 2.0 * i;
    keys[0].look_set = keys[1].look_set = 3;
    keys[3].spin_turns = 2;
    int look = 0, spin = 0;
    for (const CameraBand& b : camera_bands(keys)) {
        if (b.kind == CameraBandKind::LookAtSet) {
            ++look;
            EXPECT_EQ((uint32_t)3, b.set);
            EXPECT_EQ(0, b.first);
            EXPECT_EQ(1, b.last);
            EXPECT_NEAR(0.0, b.begin, 1.0e-9);
            EXPECT_NEAR(2.0, b.end, 1.0e-9);
        } else {
            ++spin;
            EXPECT_EQ(2, b.turns);
            EXPECT_EQ(3, b.last);
        }
    }
    EXPECT_EQ(1, look);
    EXPECT_EQ(1, spin);
}

UTEST(viamd_movie_keys, keys_that_followed_a_target_of_an_older_workspace_look_at_its_atoms) {
    std::vector<CameraKeyframe> keys(3);
    keys[0].follow = true;
    keys[1].follow = true; keys[1].follow_atom = 5;
    std::vector<MovieAtomSet> sets;
    movie_keys_from_follow(&keys, &sets, {3, 2});
    EXPECT_FALSE(keys[0].follow);
    EXPECT_FALSE(keys[1].follow);
    ASSERT_TRUE(keys[0].look_set != 0 && keys[1].look_set != 0);
    EXPECT_NE(keys[0].look_set, keys[1].look_set);
    EXPECT_EQ((uint32_t)0, keys[2].look_set);
    const MovieAtomSet* target = movie_atoms_find(sets, keys[0].look_set);
    ASSERT_TRUE(target != nullptr);
    ASSERT_EQ((size_t)2, target->atoms.size());
    EXPECT_EQ((uint32_t)2, target->atoms[0]);
    EXPECT_EQ((uint32_t)5, movie_atoms_find(sets, keys[1].look_set)->atoms[0]);
}

UTEST(viamd_movie_keys, the_focus_of_a_key_is_saved_and_values_no_key_could_have_are_refused) {
    md_allocator_i* alloc = md_vm_arena_create(1024 * 1024);
    viamd::serialization_state_t writer = {};
    writer.sb = md_strb_create(alloc);
    CameraKeyframe key;
    key.focus_on = true;
    key.focus_target = FocusTarget::Selection;
    key.focus_point = vec3_set(1, 2, 3);
    key.focus_distance = 7.5f;
    key.focus_blur = 3.5f;
    key.focus_transition = 2;
    key.focus_ease = KeyEase::Linear;
    key.focus_set = 4;
    viamd::write_section_header(writer, STR_LIT("Movie"));
    const auto encoded = movie_focus_encode(key);
    viamd::write_flt_vec(writer, STR_LIT("KeyframeFocus"), encoded.data(), encoded.size());
    viamd::deserialization_state_t reader = {};
    reader.text = md_strb_to_str(writer.sb);
    str_t section, ident, arg;
    ASSERT_TRUE(viamd::next_section_header(section, reader));
    ASSERT_TRUE(viamd::next_entry(ident, arg, reader));
    float values[9];
    ASSERT_TRUE(viamd::extract_flt_vec(values, 9, arg));
    CameraKeyframe restored;
    ASSERT_TRUE(movie_focus_decode(&restored, values));
    EXPECT_TRUE(restored.focus_on);
    EXPECT_EQ((int)key.focus_target, (int)restored.focus_target);
    EXPECT_NEAR(key.focus_point.z, restored.focus_point.z, 1.0e-6f);
    EXPECT_NEAR(key.focus_distance, restored.focus_distance, 1.0e-6f);
    EXPECT_NEAR(key.focus_blur, restored.focus_blur, 1.0e-6f);
    EXPECT_NEAR(key.focus_transition, restored.focus_transition, 1.0e-6f);
    EXPECT_EQ((int)key.focus_ease, (int)restored.focus_ease);
    EXPECT_EQ((uint32_t)4, restored.focus_set);
    values[0] = 99;
    EXPECT_FALSE(movie_focus_decode(&restored, values));
    values[0] = (float)FocusTarget::Point; values[7] = 99;
    EXPECT_FALSE(movie_focus_decode(&restored, values));
    values[7] = 0; values[1] = NAN;
    EXPECT_FALSE(movie_focus_decode(&restored, values));
    md_vm_arena_destroy(alloc);
}

UTEST(viamd_movie_keys, look_and_focus_are_part_of_the_history_and_scale_with_the_length) {
    MovieKeys keys;
    CameraKeyframe key;
    key.time = 3; key.look_set = 1; key.focus_on = true; key.focus_target = FocusTarget::Selection;
    key.focus_set = 2; key.focus_transition = 2;
    keys.camera.push_back(key);
    keys.sets.push_back({1, {1, 4}});
    keys.sets.push_back({2, {6}});
    MovieHistory history;
    history.clear(keys);
    movie_keys_scale_time(&keys, 2);
    EXPECT_NEAR(6.0, keys.camera[0].time, 1.0e-6);
    EXPECT_NEAR(4.0f, keys.camera[0].focus_transition, 1.0e-6f);
    history.update(keys, false);
    ASSERT_TRUE(history.undo(&keys));
    EXPECT_NEAR(3.0, keys.camera[0].time, 1.0e-6);
    EXPECT_NEAR(2.0f, keys.camera[0].focus_transition, 1.0e-6f);
    ASSERT_TRUE(history.redo(&keys));
    EXPECT_NEAR(4.0f, keys.camera[0].focus_transition, 1.0e-6f);
    EXPECT_EQ((uint32_t)4, keys.sets[0].atoms[1]);
    MovieKeys changed = keys;
    changed.sets[0].atoms[1] = 5;
    EXPECT_FALSE(movie_keys_equal(changed, keys));
    changed = keys;
    changed.camera[0].focus_blur += 1.0f;
    EXPECT_FALSE(movie_keys_equal(changed, keys));
    changed = keys;
    changed.camera[0].look_set = 2;
    EXPECT_FALSE(movie_keys_equal(changed, keys));
}

UTEST(viamd_movie_keys, a_key_keeps_its_look_and_focus_when_moved_copied_and_deleted) {
    MovieKeys start;
    start.duration = 20;
    CameraKeyframe key;
    key.time = 2; key.look_set = 1; key.focus_on = true; key.focus_target = FocusTarget::Selection;
    key.focus_set = 2; key.focus_transition = 2; key.focus_blur = 3;
    start.camera.push_back(key);
    key.time = 8; key.look_set = 0; key.focus_on = false;
    start.camera.push_back(key);
    KeySelection selection;
    selection.add(KeyKind::Camera, 0, 2);
    MovieKeys shifted = start;
    KeySelection moved = selection;
    movie_keys_shift(&shifted, &moved, start, selection, 3, KeyShift{}, 20);
    EXPECT_NEAR(5.0, shifted.camera[0].time, 1.0e-6);
    EXPECT_EQ((uint32_t)1, shifted.camera[0].look_set);
    EXPECT_EQ((uint32_t)2, shifted.camera[0].focus_set);
    EXPECT_TRUE(moved.contains(KeyKind::Camera, 0, 5));
    KeyClip copied = movie_keys_copy(shifted, moved);
    ASSERT_FALSE(copied.empty());
    movie_keys_paste(&shifted, &moved, copied, 12, 20);
    ASSERT_EQ((size_t)3, shifted.camera.size());
    int pasted = -1;
    for (size_t i = 0; i < shifted.camera.size(); ++i) if (fabs(shifted.camera[i].time - 12.0) < 1.0e-6) pasted = (int)i;
    ASSERT_TRUE(pasted >= 0);
    EXPECT_EQ((uint32_t)1, shifted.camera[pasted].look_set);
    EXPECT_TRUE(shifted.camera[pasted].focus_on);
    EXPECT_NEAR(3.0f, shifted.camera[pasted].focus_blur, 1.0e-6f);
    movie_keys_delete(&shifted, &moved);
    EXPECT_EQ((size_t)2, shifted.camera.size());
    MovieKeys scaled = start;
    movie_keys_scale(&scaled, &moved, start, selection, 0, 2, 20);
    EXPECT_NEAR(4.0, scaled.camera[0].time, 1.0e-6);
    EXPECT_NEAR(4.0f, scaled.camera[0].focus_transition, 1.0e-6f);
}

UTEST(viamd_movie_keys, independent_spin_orbits_the_new_look_target) {
    ViewTransform unspun = {};
    unspun.orientation = quat_t{0, 0, 0, 1};
    unspun.position = vec3_set(0, 0, 10);
    unspun.distance = 10;
    ViewTransform spun = unspun;
    spun.orientation = quat_axis_angle(vec3_set(0, 1, 0), 3.14159265358979f / 2);
    vec3_t look = vec3_set(3, 0, 0);
    movie_camera_independent(&spun, unspun, look, nullptr, 0);
    EXPECT_NEAR(13.0f, spun.position.x, 1.0e-4f);
    EXPECT_NEAR(3.0f, spun.position.z, 1.0e-4f);
    EXPECT_NEAR(vec3_length(unspun.position - look), vec3_length(spun.position - look), 1.0e-4f);
    EXPECT_NEAR(0.0f, vec3_length(camera_get_look_at(spun) - look), 1.0e-4f);
    ViewTransform stationary = unspun;
    movie_camera_independent(&stationary, unspun, look, nullptr, 0);
    EXPECT_NEAR(0.0f, vec3_length(stationary.position - unspun.position), 1.0e-4f);
}

UTEST(viamd_movie_keys, independent_spin_exports_axis_even_with_upright_roll) {
    CameraKeyframe keys[2];
    keys[0].time = 0; keys[1].time = 10;
    for (auto& key : keys) {
        key.transform.position = vec3_set(0, 0, 10);
        key.transform.orientation = quat_t{0, 0, 0, 1};
        key.transform.distance = 10;
    }
    keys[1].spin_turns = 1;
    keys[1].spin_axis = SpinAxis::ViewUp;
    keys[1].spin_constant_speed = true;
    keys[1].roll = 0.4f;
    const vec3_t up = vec3_set(0, 0, 1);
    ViewTransform camera;
    quat_t spin;
    float fov;
    camera_keyframes_evaluate(&camera, &fov, keys, 2, 2.5, false, nullptr, nullptr, &up, &spin);
    const vec3_t rotated = spin * vec3_set(1, 0, 0);
    EXPECT_NEAR(0.0f, rotated.x, 1.0e-4f);
    EXPECT_NEAR(1.0f, rotated.y, 1.0e-4f);
    camera_keyframes_evaluate(&camera, &fov, keys, 2, 10, false, nullptr, nullptr, &up, &spin);
    EXPECT_NEAR(1.0f, spin.w, 1.0e-5f);
}

UTEST(viamd_movie_keys, optics_example_of_separate_tracks_becomes_keys_with_their_own_look_and_focus) {
    const auto dir = std::filesystem::path(__FILE__).parent_path().parent_path() / "docs/examples";
    std::ifstream input(dir / "fov_and_focus.via");
    ASSERT_TRUE(input.good());
    const std::string text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    viamd::deserialization_state_t reader = {};
    reader.text = {text.data(), text.size()};
    md_allocator_i* alloc = md_vm_arena_create(1024 * 1024);
    md_bitfield_t mask = {};
    md_bitfield_init(&mask, alloc);
    std::vector<CameraKeyframe> keys;
    std::vector<MovieLegacyTarget> look, focus;
    str_t section, ident, arg;
    int focus_mode = -1;
    while (viamd::next_section_header(section, reader)) {
        while (viamd::next_entry(ident, arg, reader)) {
            if (str_eq(section, STR_LIT("Movie")) && str_eq(ident, STR_LIT("FocusAtoms"))) {
                ASSERT_TRUE(viamd::extract_bitfield(&mask, arg));
            } else if (str_eq(section, STR_LIT("Movie")) && (str_eq(ident, STR_LIT("FocusKey")) || str_eq(ident, STR_LIT("LookKey")))) {
                const bool is_focus = str_eq(ident, STR_LIT("FocusKey"));
                float v[10];
                ASSERT_TRUE(viamd::extract_flt_vec(v, 10, arg));
                MovieLegacyTarget key;
                ASSERT_TRUE(movie_legacy_target_decode(&key, v, is_focus));
                (is_focus ? focus : look).push_back(key);
            } else if (str_eq(section, STR_LIT("RenderSettings")) && str_eq(ident, STR_LIT("DofFocusMode"))) {
                ASSERT_TRUE(viamd::extract_int(focus_mode, arg));
            } else if (str_eq(section, STR_LIT("Movie")) && str_eq(ident, STR_LIT("KeyframeV3"))) {
                float v[21];
                ASSERT_TRUE(viamd::extract_flt_vec(v, 21, arg));
                CameraKeyframe key;
                key.time = v[0];
                key.fov_y = v[1];
                key.transform.distance = v[2];
                key.transform.position = vec3_set(v[3], v[4], v[5]);
                key.transform.orientation = quat_t{v[6], v[7], v[8], v[9]};
                key.use_frame = v[10] != 0;
                key.frame = v[11];
                key.ease = (KeyEase)(int)v[15];
                EXPECT_NEAR(0.0f, v[16], 1.0e-6f);
                keys.push_back(key);
            }
        }
    }
    EXPECT_EQ(0, focus_mode);
    ASSERT_EQ((size_t)1, look.size());
    ASSERT_EQ((size_t)4, focus.size());
    EXPECT_EQ((size_t)1, md_bitfield_popcount(&mask));
    EXPECT_TRUE(md_bitfield_test_bit(&mask, 2));
    ASSERT_EQ((size_t)8, keys.size());
    for (MovieLegacyTarget& f : focus) if (f.target == FocusTarget::Selection) f.atoms = {2};

    // The tracks become keys: the camera keys, and one where each focus key is that no camera key has
    std::vector<MovieAtomSet> sets;
    movie_keys_from_tracks(&keys, &sets, look, focus, nullptr, vec3_set(0, 1, 0));
    ASSERT_EQ((size_t)10, keys.size());
    ASSERT_EQ((size_t)1, sets.size());
    int sets_focus = 0;
    for (const CameraKeyframe& k : keys) {
        EXPECT_EQ((uint32_t)0, k.look_set);
        sets_focus += k.focus_on ? 1 : 0;
    }
    EXPECT_EQ(4, sets_focus);

    // The blur of the focus keys is the keys of the blur look parameter
    const std::vector<ParamKey> blur_keys = movie_blur_keys_from_focus(keys.data(), keys.size());
    for (const CameraKeyframe& k : keys) EXPECT_TRUE(k.focus_blur < 0.0f);
    std::vector<vec3_t> looks, points;
    for (const CameraKeyframe& k : keys) {
        looks.push_back(camera_get_look_at(k.transform));
        points.push_back(k.focus_point);
    }
    for (double time : {0.0, 7.0, 13.0, 18.0, 22.5, 31.5, 36.0}) {
        ViewTransform view;
        float fov;
        movie_keys_pose(&view, &fov, keys.data(), keys.size(), time, looks.data(), nullptr, vec3_set(0, 1, 0));
        EXPECT_NEAR(40.0f, view.position.z, 1.0e-3f);
        EXPECT_NEAR(0.0f, view.position.x, 1.0e-3f);
        EXPECT_NEAR(8.0f, camera_get_look_at(view).z, 1.0e-3f);
        const float degrees = fov * 180.0f / 3.14159265358979f;
        EXPECT_NEAR(time == 7.0 ? 20.0f : time == 13.0 ? 70.0f : 45.0f, degrees, 1.0e-2f);
        float blur = -1.0f;
        ASSERT_TRUE(param_keys_evaluate(&blur, 1, blur_keys.data(), blur_keys.size(), MOVIE_PARAM_DOF_APERTURE, time));
        EXPECT_NEAR(time <= 18.0 || time >= 35.0 ? 0.0f : 3.5f, blur, 1.0e-4f);
        (void)movie_focus_point(keys.data(), keys.size(), points.data(), time, view);
    }
    md_xyz_data_t xyz = {};
    const std::string path = (dir / "fov_and_focus.xyz").string();
    ASSERT_TRUE(md_xyz_data_parse_file(&xyz, {path.data(), path.size()}, alloc));
    ASSERT_EQ((size_t)9, xyz.num_models);
    for (size_t i = 0; i < xyz.num_models; ++i) {
        const md_xyz_model_t& model = xyz.models[i];
        ASSERT_EQ((uint32_t)5, model.end_coord_index - model.beg_coord_index);
        const md_xyz_coordinate_t* atoms = xyz.coordinates + model.beg_coord_index;
        EXPECT_NEAR(8.0f, atoms[0].z, 1.0e-5f);
        EXPECT_NEAR(-8.0f, atoms[4].z, 1.0e-5f);
        EXPECT_EQ(6, atoms[2].atomic_number);
    }
    EXPECT_NEAR(8.0f, xyz.coordinates[xyz.models[2].beg_coord_index + 2].z, 1.0e-5f);
    EXPECT_NEAR(-8.0f, xyz.coordinates[xyz.models[6].beg_coord_index + 2].z, 1.0e-5f);
    md_vm_arena_destroy(alloc);
}

UTEST(viamd_movie_keys, independent_focus_tracks_target_without_changing_camera) {
    md_allocator_i* alloc = md_vm_arena_create(1024 * 1024);
    md_bitfield_t mask = {};
    md_bitfield_init(&mask, alloc);
    md_bitfield_set_bit(&mask, 0);
    md_atom_type_idx_t types[] = {0, 0};
    float mass[] = {1.0f};
    md_system_t system = {};
    system.atom.count = 2;
    system.atom.type_idx = types;
    system.atom.type.count = 1;
    system.atom.type.mass = mass;
    vec3_t positions[] = {vec3_set(3, 0, -7), vec3_set(-3, 0, -12)};
    md_system_state_t state = {};
    state.num_atoms = 2;
    state.xyz = positions;
    ViewTransform camera = {};
    camera.orientation = quat_ident();
    camera.distance = 20.0f;
    const vec3_t look_at = camera_get_look_at(camera);
    vec3_t center;
    ASSERT_TRUE(movie_target_center(&center, system, state, mask, mat4_ident(), alloc));
    EXPECT_NEAR(7.0f, camera_depth_of_point(camera, center), 1.0e-5f);
    positions[0] = vec3_set(6, 0, -9);
    ASSERT_TRUE(movie_target_center(&center, system, state, mask, mat4_ident(), alloc));
    EXPECT_NEAR(9.0f, camera_depth_of_point(camera, center), 1.0e-5f);
    EXPECT_NEAR(look_at.z, camera_get_look_at(camera).z, 1.0e-5f);
    EXPECT_NEAR(20.0f, camera.distance, 1.0e-5f);
    EXPECT_NEAR(0.0f, camera.position.z, 1.0e-5f);
    md_vm_arena_destroy(alloc);
}

UTEST(viamd_movie_keys, focus_group_is_mass_weighted_periodic_and_transformed) {
    md_allocator_i* alloc = md_vm_arena_create(1024 * 1024);
    md_bitfield_t mask = {};
    md_bitfield_init(&mask, alloc);
    md_bitfield_set_range(&mask, 0, 2);
    md_atom_type_idx_t types[] = {0, 1};
    float masses[] = {1.0f, 3.0f};
    md_system_t system = {};
    system.atom.count = 2;
    system.atom.type_idx = types;
    system.atom.type.count = 2;
    system.atom.type.mass = masses;
    vec3_t positions[] = {vec3_set(9, 0, 0), vec3_set(1, 0, 0)};
    md_system_state_t state = {};
    state.num_atoms = 2;
    state.xyz = positions;
    state.unitcell = md_unitcell_from_extent(10, 10, 10);
    vec3_t center;
    ASSERT_TRUE(movie_target_center(&center, system, state, mask, mat4_translate(0, 0, -8), alloc));
    const float image = center.x - 10.0f * floorf(center.x / 10.0f);
    EXPECT_NEAR(0.5f, image, 1.0e-4f);
    EXPECT_NEAR(-8.0f, center.z, 1.0e-4f);
    EXPECT_NEAR(9.0f, positions[0].x, 1.0e-5f);
    EXPECT_NEAR(1.0f, positions[1].x, 1.0e-5f);
    md_vm_arena_destroy(alloc);
}

UTEST(viamd_movie_keys, focus_target_rejects_empty_and_invalid_selections) {
    md_allocator_i* alloc = md_vm_arena_create(1024 * 1024);
    md_bitfield_t mask = {};
    md_bitfield_init(&mask, alloc);
    md_system_t system = {};
    system.atom.count = 1;
    vec3_t position = vec3_set(0, 0, -5);
    md_system_state_t state = {};
    state.num_atoms = 1;
    state.xyz = &position;
    vec3_t center = vec3_set1(42);
    EXPECT_FALSE(movie_target_center(&center, system, state, mask, mat4_ident(), alloc));
    md_bitfield_set_bit(&mask, 1);
    EXPECT_FALSE(movie_target_center(&center, system, state, mask, mat4_ident(), alloc));
    md_bitfield_clear(&mask);
    md_bitfield_set_bit(&mask, 0);
    state.num_atoms = 0;
    EXPECT_FALSE(movie_target_center(&center, system, state, mask, mat4_ident(), alloc));
    EXPECT_NEAR(42.0f, center.x, 1.0e-5f);
    md_vm_arena_destroy(alloc);
}

UTEST(viamd_movie_keys, changing_length_scales_all_timing_and_can_be_undone) {
    MovieKeys keys;
    keys.duration = 20.0f;
    keys.traj_begin = 4.0f;
    keys.traj_end = 18.0f;
    keys.start_frame = 0.0;
    keys.end_frame = 600.0;
    CameraKeyframe camera;
    camera.time = 10.0;
    camera.use_frame = true;
    camera.frame = 300.0;
    camera.spin_turns = 2;
    keys.camera.push_back(camera);
    ParamKey param;
    param.time = 12.0;
    param.value[0] = 7.0f;
    keys.params.push_back(param);
    MovieOverlay overlay;
    overlay.begin = 2.0;
    overlay.end = 16.0;
    overlay.fade_in = 1.0f;
    overlay.fade_out = 2.0f;
    keys.overlays.push_back(overlay);
    const MovieKeys original = keys;
    MovieHistory history;
    history.clear(keys);

    movie_keys_scale_time(&keys, 2.0);
    EXPECT_EQ(40.0f, keys.duration);
    EXPECT_EQ(8.0f, keys.traj_begin);
    EXPECT_EQ(36.0f, keys.traj_end);
    EXPECT_EQ(20.0, keys.camera[0].time);
    EXPECT_EQ(300.0, keys.camera[0].frame);
    EXPECT_EQ(2, keys.camera[0].spin_turns);
    EXPECT_EQ(24.0, keys.params[0].time);
    EXPECT_EQ(7.0f, keys.params[0].value[0]);
    EXPECT_EQ(4.0, keys.overlays[0].begin);
    EXPECT_EQ(32.0, keys.overlays[0].end);
    EXPECT_EQ(2.0f, keys.overlays[0].fade_in);
    EXPECT_EQ(4.0f, keys.overlays[0].fade_out);
    EXPECT_EQ(600.0, keys.end_frame);
    const MovieKeys doubled = keys;
    history.update(keys, false);
    ASSERT_TRUE(history.undo(&keys));
    EXPECT_TRUE(movie_keys_equal(original, keys));
    ASSERT_TRUE(history.redo(&keys));
    EXPECT_TRUE(movie_keys_equal(doubled, keys));
    movie_keys_scale_time(&keys, 0.5);
    EXPECT_TRUE(movie_keys_equal(original, keys));
}

/* Curve through keyed values with the easing of the key a segment leads to */

UTEST(viamd_movie_keys, a_keyed_curve_passes_through_its_keys_and_holds_outside) {
    const double t[3] = {0.0, 2.0, 4.0};
    const double v[3] = {0.0, 10.0, 20.0};
    const KeyEase e[3] = {KeyEase::Smooth, KeyEase::Smooth, KeyEase::Smooth};

    for (int i = 0; i < 3; ++i) EXPECT_NEAR(v[i], keyed_curve_evaluate(t, v, e, 3, t[i]), 1.0e-9);
    EXPECT_NEAR(5.0, keyed_curve_evaluate(t, v, e, 3, 1.0), 1.0e-9);   /* even spacing of a line stays a line */
    EXPECT_NEAR(0.0, keyed_curve_evaluate(t, v, e, 3, -3.0), 1.0e-12);
    EXPECT_NEAR(20.0, keyed_curve_evaluate(t, v, e, 3, 9.0), 1.0e-12);
    EXPECT_NEAR(0.0, keyed_curve_evaluate(t, v, e, 1, 3.0), 1.0e-12);   /* a single key is its value, whenever */
    EXPECT_NEAR(0.0, keyed_curve_evaluate(t, v, e, 0, 3.0), 1.0e-12);
}

UTEST(viamd_movie_keys, the_easing_of_a_key_shapes_the_stretch_leading_to_it) {
    const double t[3] = {0.0, 2.0, 4.0};
    const double v[3] = {0.0, 10.0, 20.0};

    KeyEase e[3] = {KeyEase::Smooth, KeyEase::Hold, KeyEase::Smooth};
    EXPECT_NEAR(0.0, keyed_curve_evaluate(t, v, e, 3, 1.99), 1.0e-12);    /* stays until the key */
    EXPECT_NEAR(10.0, keyed_curve_evaluate(t, v, e, 3, 2.0), 1.0e-9);

    e[1] = KeyEase::EaseInOut;
    const double u = 0.25;
    EXPECT_NEAR(10.0 * (u * u * (3.0 - 2.0 * u)), keyed_curve_evaluate(t, v, e, 3, 0.5), 1.0e-9);

    e[1] = KeyEase::Smooth;
    e[2] = KeyEase::Linear;
    EXPECT_NEAR(15.0, keyed_curve_evaluate(t, v, e, 3, 3.0), 1.0e-9);
}

UTEST(viamd_movie_keys, a_smooth_curve_does_not_overshoot_a_key) {
    const double t[4] = {0.0, 1.0, 2.0, 3.0};
    const double v[4] = {0.0, 10.0, 0.0, 5.0};   /* the second is a peak */
    const KeyEase e[4] = {KeyEase::Smooth, KeyEase::Smooth, KeyEase::Smooth, KeyEase::Smooth};
    for (double x = 0.0; x <= 3.0; x += 0.01) {
        const double y = keyed_curve_evaluate(t, v, e, 4, x);
        EXPECT_LE(y, 10.0 + 1.0e-9);
        EXPECT_GE(y, -1.0e-9);
    }
}

/* Look parameters */

static ParamKey pk(int param, double time, float a, float b = 0.0f, float c = 0.0f, KeyEase ease = KeyEase::Smooth) {
    ParamKey k;
    k.param = param;
    k.time = time;
    k.value[0] = a; k.value[1] = b; k.value[2] = c;
    k.ease = ease;
    return k;
}

UTEST(viamd_movie_keys, a_parameter_without_keys_is_left_alone) {
    ParamKey keys[1] = { pk(3, 1.0, 5.0f) };
    float v[3] = {-1, -1, -1};
    EXPECT_FALSE(param_keys_evaluate(v, 1, keys, 1, 2, 1.0));
    EXPECT_FALSE(param_keys_evaluate(v, 1, keys, 0, 3, 1.0));
    EXPECT_EQ(-1.0f, v[0]);
}

UTEST(viamd_movie_keys, parameters_are_keyed_independently_of_each_other_and_of_order) {
    /* Out of order, and with keys for another parameter in between */
    ParamKey keys[5] = {
        pk(1, 4.0, 40.0f),
        pk(2, 1.0, 999.0f),
        pk(1, 0.0, 0.0f),
        pk(2, 3.0, -999.0f),
        pk(1, 2.0, 20.0f),
    };
    float v[3] = {};
    ASSERT_TRUE(param_keys_evaluate(v, 1, keys, 5, 1, 1.0));
    EXPECT_NEAR(10.0f, v[0], 1.0e-4f);
    ASSERT_TRUE(param_keys_evaluate(v, 1, keys, 5, 1, 3.0));
    EXPECT_NEAR(30.0f, v[0], 1.0e-4f);
    ASSERT_TRUE(param_keys_evaluate(v, 1, keys, 5, 2, 2.0));
    EXPECT_NEAR(0.0f, v[0], 1.0f);                    /* between its own two keys, the other's do not matter */
    ASSERT_TRUE(param_keys_evaluate(v, 1, keys, 5, 2, 10.0));
    EXPECT_NEAR(-999.0f, v[0], 1.0e-3f);
}

UTEST(viamd_movie_keys, a_colour_is_keyed_in_all_three_channels) {
    ParamKey keys[2] = { pk(0, 0.0, 0.0f, 1.0f, 0.5f, KeyEase::Smooth), pk(0, 2.0, 1.0f, 0.0f, 0.5f, KeyEase::Linear) };
    float v[3] = {};
    ASSERT_TRUE(param_keys_evaluate(v, 3, keys, 2, 0, 1.0));
    EXPECT_NEAR(0.5f, v[0], 1.0e-5f);
    EXPECT_NEAR(0.5f, v[1], 1.0e-5f);
    EXPECT_NEAR(0.5f, v[2], 1.0e-5f);
}

/* Undo */

static MovieKeys keys_with(int n, bool loop = false) {
    MovieKeys k;
    for (int i = 0; i < n; ++i) {
        CameraKeyframe c;
        c.time = (double)i;
        k.camera.push_back(c);
    }
    k.loop = loop;
    return k;
}

UTEST(viamd_movie_keys, equality_looks_at_everything_that_is_edited) {
    MovieKeys a = keys_with(2), b = keys_with(2);
    EXPECT_TRUE(movie_keys_equal(a, b));
    b.camera[1].spin_turns = 1;
    EXPECT_FALSE(movie_keys_equal(a, b));
    b = a;
    b.camera[0].ease = KeyEase::Hold;
    EXPECT_FALSE(movie_keys_equal(a, b));
    b = a;
    b.loop = true;
    EXPECT_FALSE(movie_keys_equal(a, b));
    b = a;
    b.params.push_back(pk(0, 0.0, 1.0f));
    EXPECT_FALSE(movie_keys_equal(a, b));
}

UTEST(viamd_movie_keys, an_edit_is_one_undo_step_and_redo_brings_it_back) {
    MovieHistory h;
    MovieKeys cur = keys_with(1);
    h.clear(cur);
    EXPECT_FALSE(h.can_undo());

    cur = keys_with(2);
    h.update(cur, false);
    EXPECT_TRUE(h.can_undo());
    EXPECT_FALSE(h.can_redo());

    ASSERT_TRUE(h.undo(&cur));
    EXPECT_EQ(1u, (uint32_t)cur.camera.size());
    EXPECT_TRUE(h.can_redo());
    EXPECT_FALSE(h.can_undo());

    h.update(cur, false);   /* being shown the restored keys is not an edit */
    EXPECT_TRUE(h.can_redo());

    ASSERT_TRUE(h.redo(&cur));
    EXPECT_EQ(2u, (uint32_t)cur.camera.size());
    EXPECT_FALSE(h.redo(&cur));
}

UTEST(viamd_movie_keys, a_drag_is_one_step_not_many) {
    MovieHistory h;
    MovieKeys cur = keys_with(2);
    h.clear(cur);

    for (int i = 1; i <= 50; ++i) {
        cur.camera[1].time = 1.0 + 0.01 * i;
        h.update(cur, true);          /* the mouse is down */
    }
    EXPECT_FALSE(h.can_undo());
    h.update(cur, false);             /* released */
    ASSERT_TRUE(h.can_undo());

    ASSERT_TRUE(h.undo(&cur));
    EXPECT_NEAR(1.0, cur.camera[1].time, 1.0e-12);
    EXPECT_FALSE(h.can_undo());       /* it was the only step */
}

UTEST(viamd_movie_keys, undoing_in_the_middle_of_an_edit_counts_the_edit_first) {
    MovieHistory h;
    MovieKeys cur = keys_with(2);
    h.clear(cur);

    cur.camera[0].time = 0.5;
    h.update(cur, true);              /* not committed yet */
    ASSERT_TRUE(h.undo(&cur));
    EXPECT_NEAR(0.0, cur.camera[0].time, 1.0e-12);
}

UTEST(viamd_movie_keys, a_new_edit_after_an_undo_drops_what_could_be_redone) {
    MovieHistory h;
    MovieKeys cur = keys_with(1);
    h.clear(cur);
    cur = keys_with(2); h.update(cur, false);
    cur = keys_with(3); h.update(cur, false);

    ASSERT_TRUE(h.undo(&cur));
    EXPECT_TRUE(h.can_redo());
    cur = keys_with(5, true);
    h.update(cur, false);
    EXPECT_FALSE(h.can_redo());
    ASSERT_TRUE(h.undo(&cur));
    EXPECT_EQ(2u, (uint32_t)cur.camera.size());
}

UTEST(viamd_movie_keys, clearing_forgets_the_past) {
    MovieHistory h;
    MovieKeys cur = keys_with(1);
    h.clear(cur);
    cur = keys_with(2); h.update(cur, false);
    EXPECT_TRUE(h.can_undo());

    h.clear(cur);
    EXPECT_FALSE(h.can_undo());
    EXPECT_FALSE(h.can_redo());
    h.update(cur, false);
    EXPECT_FALSE(h.can_undo());       /* what was cleared to is the present */
}

/* Which frames a recording makes, how big they are, how long it has left */

UTEST(viamd_movie_keys, without_a_range_every_frame_is_rendered) {
    int a = -1, b = -1;
    movie_frame_range(121, 24.0, false, 2.0, 3.0, &a, &b);
    EXPECT_EQ(0, a);
    EXPECT_EQ(120, b);
}

UTEST(viamd_movie_keys, a_range_picks_the_frames_inside_it) {
    int a = -1, b = -1;
    movie_frame_range(121, 24.0, true, 1.0, 2.0, &a, &b);   /* 24 .. 48 */
    EXPECT_EQ(24, a);
    EXPECT_EQ(48, b);

    movie_frame_range(121, 24.0, true, 1.01, 1.99, &a, &b); /* between frames: the ones inside */
    EXPECT_EQ(25, a);
    EXPECT_EQ(47, b);
}

UTEST(viamd_movie_keys, a_range_always_holds_a_frame_and_stays_inside_the_movie) {
    int a = -1, b = -1;
    movie_frame_range(121, 24.0, true, 3.0, 3.0, &a, &b);
    EXPECT_EQ(72, a);
    EXPECT_EQ(72, b);

    movie_frame_range(121, 24.0, true, 1.01, 1.02, &a, &b); /* nothing exactly inside: the first after begin */
    EXPECT_EQ(25, a);
    EXPECT_GE(b, a);

    movie_frame_range(121, 24.0, true, 4.0, 90.0, &a, &b);
    EXPECT_EQ(96, a);
    EXPECT_EQ(120, b);

    movie_frame_range(121, 24.0, true, 90.0, 95.0, &a, &b);
    EXPECT_EQ(120, a);
    EXPECT_EQ(120, b);
}

UTEST(viamd_movie_keys, a_scaled_size_is_even_and_never_larger) {
    int w = 1920, h = 1080;
    movie_scaled_size(&w, &h, 100);
    EXPECT_EQ(1920, w);
    EXPECT_EQ(1080, h);

    movie_scaled_size(&w, &h, 50);
    EXPECT_EQ(960, w);
    EXPECT_EQ(540, h);

    w = 1001; h = 667;
    movie_scaled_size(&w, &h, 75);
    EXPECT_EQ(0, w % 2);
    EXPECT_EQ(0, h % 2);
    EXPECT_LE(w, 751);

    w = 20; h = 10;
    movie_scaled_size(&w, &h, 10);
    EXPECT_GE(w, 2);
    EXPECT_GE(h, 2);
}

UTEST(viamd_movie_keys, time_left_follows_the_pace_so_far) {
    double s = 0.0;
    EXPECT_FALSE(movie_time_left(1, 100, 5.0, &s));    /* too early to tell */
    EXPECT_FALSE(movie_time_left(50, 100, 0.2, &s));
    ASSERT_TRUE(movie_time_left(25, 100, 50.0, &s));   /* 2 s a frame, 75 to go */
    EXPECT_NEAR(150.0, s, 1.0e-9);
    ASSERT_TRUE(movie_time_left(100, 100, 50.0, &s));
    EXPECT_NEAR(0.0, s, 1.0e-12);
}

UTEST(viamd_movie_keys, a_time_snaps_to_the_nearest_frame_inside_the_movie) {
    EXPECT_NEAR(24.0 / 24.0, movie_snap_to_frame(1.0 + 0.4 / 24.0, 24.0, 5.0), 1.0e-9);
    EXPECT_NEAR(25.0 / 24.0, movie_snap_to_frame(1.0 + 0.6 / 24.0, 24.0, 5.0), 1.0e-9);
    EXPECT_NEAR(0.0, movie_snap_to_frame(-3.0, 24.0, 5.0), 1.0e-12);
    EXPECT_NEAR(5.0, movie_snap_to_frame(9.0, 24.0, 5.0), 1.0e-12);
    /* A length that is not a whole number of frames: the end is the length, not a frame past it */
    EXPECT_NEAR(5.01, movie_snap_to_frame(5.2, 10.0, 5.01), 1.0e-12);
}

/* Keys of properties of representations */

static RepKey rk(uint32_t rep, RepProp prop, double time, float value, KeyEase ease = KeyEase::Smooth) {
    RepKey k;
    k.rep = rep;
    k.prop = (int)prop;
    k.time = time;
    k.value[0] = value;
    k.ease = ease;
    return k;
}

UTEST(viamd_movie_keys, rep_keys_follow_their_own_representation_and_property) {
    const RepKey keys[] = {
        rk(1, RepProp::Scale0, 0.0, 1.0f, KeyEase::Linear), rk(1, RepProp::Scale0, 4.0, 3.0f, KeyEase::Linear),
        rk(2, RepProp::Scale0, 0.0, 10.0f), rk(1, RepProp::Saturation, 0.0, 0.5f),
    };
    float v = 0.0f;
    ASSERT_TRUE(rep_keys_evaluate(&v, keys, 4, 1, (int)RepProp::Scale0, 2.0));
    EXPECT_NEAR(2.0f, v, 1.0e-6f);
    ASSERT_TRUE(rep_keys_evaluate(&v, keys, 4, 2, (int)RepProp::Scale0, 2.0));
    EXPECT_NEAR(10.0f, v, 1.0e-6f);
    ASSERT_TRUE(rep_keys_evaluate(&v, keys, 4, 1, (int)RepProp::Saturation, 9.0));
    EXPECT_NEAR(0.5f, v, 1.0e-6f);
    EXPECT_FALSE(rep_keys_evaluate(&v, keys, 4, 3, (int)RepProp::Scale0, 2.0));
    EXPECT_FALSE(rep_keys_evaluate(&v, keys, 4, 1, (int)RepProp::TintScale, 2.0));
}

UTEST(viamd_movie_keys, a_visible_key_changes_at_its_key_and_is_not_blended) {
    /* Even with a smooth ease on the keys: shown, hidden at 2 s, shown again at 4 s */
    const RepKey keys[] = {
        rk(1, RepProp::Visible, 4.0, 1.0f), rk(1, RepProp::Visible, 0.0, 1.0f), rk(1, RepProp::Visible, 2.0, 0.0f),
    };
    float v = -1.0f;
    ASSERT_TRUE(rep_keys_evaluate(&v, keys, 3, 1, (int)RepProp::Visible, 1.99));
    EXPECT_NEAR(1.0f, v, 1.0e-6f);
    ASSERT_TRUE(rep_keys_evaluate(&v, keys, 3, 1, (int)RepProp::Visible, 2.0));
    EXPECT_NEAR(0.0f, v, 1.0e-6f);
    ASSERT_TRUE(rep_keys_evaluate(&v, keys, 3, 1, (int)RepProp::Visible, 3.99));
    EXPECT_NEAR(0.0f, v, 1.0e-6f);
    ASSERT_TRUE(rep_keys_evaluate(&v, keys, 3, 1, (int)RepProp::Visible, 4.0));
    EXPECT_NEAR(1.0f, v, 1.0e-6f);
}

UTEST(viamd_movie_keys, rep_keys_are_part_of_the_undo_state_and_scale_with_time) {
    MovieKeys a, b;
    a.reps.push_back(rk(1, RepProp::Scale0, 2.0, 1.0f));
    EXPECT_FALSE(movie_keys_equal(a, b));
    b.reps.push_back(rk(1, RepProp::Scale0, 2.0, 1.0f));
    EXPECT_TRUE(movie_keys_equal(a, b));
    b.reps[0].rep = 2;
    EXPECT_FALSE(movie_keys_equal(a, b));

    movie_keys_scale_time(&a, 2.0);
    EXPECT_NEAR(4.0, a.reps[0].time, 1.0e-12);
}

UTEST(viamd_movie_keys, a_color_key_blends_all_three_components) {
    RepKey a = rk(1, RepProp::BaseColor, 0.0, 0.0f, KeyEase::Linear);
    RepKey b = rk(1, RepProp::BaseColor, 2.0, 1.0f, KeyEase::Linear);
    a.value[1] = 1.0f; a.value[2] = 0.5f;
    b.value[1] = 0.0f; b.value[2] = 0.5f;
    const RepKey keys[] = {a, b};
    EXPECT_EQ(3, rep_prop_comps((int)RepProp::BaseColor));
    EXPECT_EQ(3, rep_prop_comps((int)RepProp::TintColor));
    EXPECT_EQ(1, rep_prop_comps((int)RepProp::Saturation));

    float v[3] = {};
    ASSERT_TRUE(rep_keys_evaluate(v, keys, 2, 1, (int)RepProp::BaseColor, 1.0));
    EXPECT_NEAR(0.5f, v[0], 1.0e-6f);
    EXPECT_NEAR(0.5f, v[1], 1.0e-6f);
    EXPECT_NEAR(0.5f, v[2], 1.0e-6f);
    ASSERT_TRUE(rep_keys_evaluate(v, keys, 2, 1, (int)RepProp::BaseColor, 9.0));
    EXPECT_NEAR(1.0f, v[0], 1.0e-6f);
    EXPECT_NEAR(0.0f, v[1], 1.0e-6f);
}

/* A representation that grows in and shrinks away at its Visible keys */

static RepKey vk(double time, bool shown) {
    return rk(1, RepProp::Visible, time, shown ? 1.0f : 0.0f);
}

UTEST(viamd_movie_keys, a_visible_key_starts_a_smooth_transition) {
    const RepKey keys[] = { vk(0.0, false), vk(10.0, true), vk(30.0, false) };
    float f = -1.0f;
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 5.0, 2.0));   /* before the first change */
    EXPECT_NEAR(0.0f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 10.0, 2.0));  /* it starts at the key */
    EXPECT_NEAR(0.0f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 11.0, 2.0));  /* halfway */
    EXPECT_NEAR(0.5f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 12.0, 2.0));  /* done */
    EXPECT_NEAR(1.0f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 20.0, 2.0));
    EXPECT_NEAR(1.0f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 31.0, 2.0));  /* shrinking away */
    EXPECT_NEAR(0.5f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 99.0, 2.0));
    EXPECT_NEAR(0.0f, f, 1.0e-6f);
}

UTEST(viamd_movie_keys, the_transition_never_overshoots_and_only_goes_one_way) {
    const RepKey keys[] = { vk(0.0, false), vk(10.0, true) };
    float prev = 0.0f, f = 0.0f;
    for (double t = 10.0; t <= 13.0; t += 0.01) {
        ASSERT_TRUE(rep_visible_factor(&f, keys, 2, 1, t, 2.0));
        EXPECT_GE(f, prev - 1.0e-6f);
        EXPECT_LE(f, 1.0f);
        prev = f;
    }
}

UTEST(viamd_movie_keys, a_transition_of_zero_changes_at_the_key) {
    const RepKey keys[] = { vk(0.0, true), vk(4.0, false) };
    float f = -1.0f;
    ASSERT_TRUE(rep_visible_factor(&f, keys, 2, 1, 3.99, 0.0));
    EXPECT_NEAR(1.0f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 2, 1, 4.0, 0.0));
    EXPECT_NEAR(0.0f, f, 1.0e-6f);
}

UTEST(viamd_movie_keys, a_change_that_comes_before_the_last_is_done_turns_around_from_where_it_was) {
    /* Shown at 10 (takes 4 s), hidden again at 12: at 12 it was halfway up, and goes down from there */
    const RepKey keys[] = { vk(0.0, false), vk(10.0, true), vk(12.0, false) };
    float f = -1.0f;
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 12.0, 4.0));
    EXPECT_NEAR(0.5f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 14.0, 4.0));
    EXPECT_NEAR(0.25f, f, 1.0e-6f);
    ASSERT_TRUE(rep_visible_factor(&f, keys, 3, 1, 16.0, 4.0));
    EXPECT_NEAR(0.0f, f, 1.0e-6f);
}

UTEST(viamd_movie_keys, no_visible_keys_means_no_factor) {
    const RepKey keys[] = { rk(1, RepProp::Scale0, 0.0, 2.0f), vk(0.0, true) };
    float f = 0.0f;
    EXPECT_FALSE(rep_visible_factor(&f, keys, 2, 2, 1.0, 2.0));
    EXPECT_FALSE(rep_visible_factor(&f, keys, 1, 1, 1.0, 2.0));
}

UTEST(viamd_movie_keys, an_overlay_background_is_part_of_the_undo_state) {
    MovieKeys a, b;
    a.overlays.push_back(MovieOverlay{});
    b.overlays.push_back(MovieOverlay{});
    EXPECT_TRUE(movie_keys_equal(a, b));
    b.overlays[0].background[3] = 0.5f;
    EXPECT_FALSE(movie_keys_equal(a, b));
}

UTEST(viamd_movie_keys, an_overlay_image_file_is_part_of_the_undo_state) {
    MovieKeys a, b;
    a.overlays.push_back(MovieOverlay{});
    b.overlays.push_back(MovieOverlay{});
    snprintf(b.overlays[0].path, sizeof(b.overlays[0].path), "/pics/group.png");
    EXPECT_FALSE(movie_keys_equal(a, b));
    snprintf(a.overlays[0].path, sizeof(a.overlays[0].path), "/pics/group.png");
    EXPECT_TRUE(movie_keys_equal(a, b));
}

UTEST(viamd_movie_keys, the_panels_and_the_look_of_a_figure_are_part_of_the_undo_state) {
    MovieKeys a, b;
    a.overlays.push_back(MovieOverlay{});
    b.overlays.push_back(MovieOverlay{});
    b.overlays[0].panels.push_back({MoviePlotView::Timeline, 3});
    EXPECT_FALSE(movie_keys_equal(a, b));
    a.overlays[0].panels.push_back({MoviePlotView::Timeline, 4});
    EXPECT_FALSE(movie_keys_equal(a, b));   /* another subplot */
    a.overlays[0].panels[0].subplot = 3;
    EXPECT_TRUE(movie_keys_equal(a, b));
    a.overlays[0].panels[0].view = MoviePlotView::Distribution;
    EXPECT_FALSE(movie_keys_equal(a, b));
    a.overlays[0].panels[0].view = MoviePlotView::Timeline;
    b.overlays[0].font_points = 24.0f;
    EXPECT_FALSE(movie_keys_equal(a, b));
    b.overlays[0].font_points = 0.0f;
    b.overlays[0].palette = 2;
    EXPECT_FALSE(movie_keys_equal(a, b));
}

UTEST(viamd_movie_keys, markers_are_part_of_the_undo_state_and_scale_with_the_length) {
    MovieKeys a, b;
    a.markers.push_back({10.0, "water"});
    EXPECT_FALSE(movie_keys_equal(a, b));
    b.markers.push_back({10.0, "water"});
    EXPECT_TRUE(movie_keys_equal(a, b));
    b.markers[0].label[0] = 'W';
    EXPECT_FALSE(movie_keys_equal(a, b));
    movie_keys_scale_time(&a, 2.0);
    EXPECT_NEAR(20.0, a.markers[0].time, 1e-9);
}

UTEST(viamd_movie_keys, bin_counts_and_marker_targets_can_be_undone_and_redone) {
    MovieKeys cur;
    cur.overlays.push_back(MovieOverlay{});
    cur.overlays[0].type = MovieOverlayType::Distribution;
    cur.markers.push_back({10.0, "water"});
    MovieHistory history;
    history.clear(cur);

    cur.overlays[0].num_bins = 37;
    history.update(cur, false);
    cur.markers[0].subplot = 4;
    history.update(cur, false);

    ASSERT_TRUE(history.undo(&cur));
    EXPECT_EQ(0u, cur.markers[0].subplot);
    EXPECT_EQ(37, cur.overlays[0].num_bins);
    ASSERT_TRUE(history.undo(&cur));
    EXPECT_EQ(0, cur.overlays[0].num_bins);
    ASSERT_TRUE(history.redo(&cur));
    ASSERT_TRUE(history.redo(&cur));
    EXPECT_EQ(37, cur.overlays[0].num_bins);
    EXPECT_EQ(4u, cur.markers[0].subplot);

    movie_keys_scale_time(&cur, 2.0);
    EXPECT_NEAR(20.0, cur.markers[0].time, 1e-9);
    EXPECT_EQ(4u, cur.markers[0].subplot);
    EXPECT_EQ(37, cur.overlays[0].num_bins);
}

UTEST(viamd_movie_keys, the_time_of_a_panel_is_part_of_the_undo_state_and_scales_with_the_length) {
    MovieKeys a, b;
    a.overlays.push_back(MovieOverlay{});
    b.overlays.push_back(MovieOverlay{});
    a.overlays[0].panels.push_back({MoviePlotView::Timeline, 3, 0.0, 0.0});
    b.overlays[0].panels.push_back({MoviePlotView::Timeline, 3, 0.0, 0.0});
    EXPECT_TRUE(movie_keys_equal(a, b));
    b.overlays[0].panels[0].begin = 12.0;
    EXPECT_FALSE(movie_keys_equal(a, b));
    a.overlays[0].panels[0].begin = 12.0;
    a.overlays[0].panels[0].end = 20.0;
    EXPECT_FALSE(movie_keys_equal(a, b));
    movie_keys_scale_time(&a, 2.0);
    EXPECT_NEAR(24.0, a.overlays[0].panels[0].begin, 1e-9);
    EXPECT_NEAR(40.0, a.overlays[0].panels[0].end, 1e-9);
}

UTEST(viamd_movie_keys, the_title_of_a_panel_is_part_of_the_undo_state) {
    MovieKeys a, b;
    a.overlays.push_back(MovieOverlay{});
    b.overlays.push_back(MovieOverlay{});
    a.overlays[0].panels.push_back({MoviePlotView::Timeline, 3, 0.0, 0.0});
    b.overlays[0].panels.push_back({MoviePlotView::Timeline, 3, 0.0, 0.0});
    EXPECT_TRUE(movie_keys_equal(a, b));
    snprintf(b.overlays[0].panels[0].title, sizeof(b.overlays[0].panels[0].title), "Distance to the ligand");
    EXPECT_FALSE(movie_keys_equal(a, b));
    snprintf(a.overlays[0].panels[0].title, sizeof(a.overlays[0].panels[0].title), "Distance to the ligand");
    EXPECT_TRUE(movie_keys_equal(a, b));
}

/* Representations over time as stretches */

static RepKey vis_key(uint32_t rep, double time, bool shown) {
    RepKey k;
    k.rep = rep;
    k.prop = (int)RepProp::Visible;
    k.time = time;
    k.value[0] = shown ? 1.0f : 0.0f;
    k.ease = KeyEase::Hold;
    return k;
}

UTEST(viamd_movie_keys, the_stretches_a_representation_is_shown_are_read_from_its_keys) {
    const std::vector<RepKey> keys = {vis_key(8, 0, false), vis_key(8, 2.4, true), vis_key(8, 9.6, false), vis_key(8, 60, true), vis_key(3, 0, true)};
    const std::vector<RepInterval> iv = rep_shown_intervals(keys, 8, 80.0);
    ASSERT_EQ(2, (int)iv.size());
    EXPECT_NEAR(2.4, iv[0].begin, 1e-9);
    EXPECT_NEAR(9.6, iv[0].end, 1e-9);
    EXPECT_EQ(1, iv[0].begin_key);
    EXPECT_EQ(2, iv[0].end_key);
    EXPECT_NEAR(60.0, iv[1].begin, 1e-9);
    EXPECT_NEAR(80.0, iv[1].end, 1e-9);
    EXPECT_EQ(-1, iv[1].end_key);
}

UTEST(viamd_movie_keys, a_first_key_that_shows_it_holds_from_the_start_of_the_movie) {
    const std::vector<RepKey> keys = {vis_key(1, 3, true), vis_key(1, 20, false)};
    const std::vector<RepInterval> iv = rep_shown_intervals(keys, 1, 50.0);
    ASSERT_EQ(1, (int)iv.size());
    EXPECT_NEAR(0.0, iv[0].begin, 1e-9);
    EXPECT_TRUE(iv[0].begin_is_start);
    EXPECT_NEAR(20.0, iv[0].end, 1e-9);
}

UTEST(viamd_movie_keys, a_representation_without_keys_has_no_stretches) {
    const std::vector<RepKey> keys = {vis_key(2, 5, true)};
    EXPECT_TRUE(rep_shown_intervals(keys, 1, 50.0).empty());
}

UTEST(viamd_movie_keys, the_end_of_a_stretch_is_moved_with_its_key) {
    std::vector<RepKey> keys = {vis_key(8, 0, false), vis_key(8, 2.4, true), vis_key(8, 9.6, false), vis_key(8, 60, true)};
    const RepInterval first = rep_shown_intervals(keys, 8, 80.0)[0];
    rep_move_interval(&keys, 8, first, 2.4, 12.0, 80.0);
    const std::vector<RepInterval> iv = rep_shown_intervals(keys, 8, 80.0);
    EXPECT_NEAR(2.4, iv[0].begin, 1e-9);
    EXPECT_NEAR(12.0, iv[0].end, 1e-9);
    EXPECT_EQ(4, (int)keys.size());
}

UTEST(viamd_movie_keys, starting_a_stretch_that_began_at_the_start_later_adds_a_hidden_key_before_it) {
    std::vector<RepKey> keys = {vis_key(1, 3, true), vis_key(1, 20, false)};
    const RepInterval iv = rep_shown_intervals(keys, 1, 50.0)[0];
    rep_move_interval(&keys, 1, iv, 5.0, 20.0, 50.0);
    const std::vector<RepInterval> after = rep_shown_intervals(keys, 1, 50.0);
    ASSERT_EQ(1, (int)after.size());
    EXPECT_NEAR(5.0, after[0].begin, 1e-9);
    EXPECT_NEAR(20.0, after[0].end, 1e-9);
    EXPECT_FALSE(after[0].begin_is_start);
    EXPECT_EQ(3, (int)keys.size());
}

UTEST(viamd_movie_keys, ending_a_stretch_that_lasted_to_the_end_earlier_adds_a_hidden_key) {
    std::vector<RepKey> keys = {vis_key(8, 0, false), vis_key(8, 60, true)};
    const RepInterval iv = rep_shown_intervals(keys, 8, 80.0)[0];
    rep_move_interval(&keys, 8, iv, 60.0, 70.0, 80.0);
    const std::vector<RepInterval> after = rep_shown_intervals(keys, 8, 80.0);
    ASSERT_EQ(1, (int)after.size());
    EXPECT_NEAR(70.0, after[0].end, 1e-9);
}

UTEST(viamd_movie_keys, a_stretch_stays_between_its_neighbours) {
    std::vector<RepKey> keys = {vis_key(8, 0, false), vis_key(8, 2.4, true), vis_key(8, 9.6, false), vis_key(8, 60, true), vis_key(8, 70, false)};
    const RepInterval second = rep_shown_intervals(keys, 8, 80.0)[1];
    rep_move_interval(&keys, 8, second, 5.0, 15.0, 80.0);    /* moved only, into the first one */
    const std::vector<RepInterval> iv = rep_shown_intervals(keys, 8, 80.0);
    ASSERT_EQ(2, (int)iv.size());
    EXPECT_NEAR(9.65, iv[1].begin, 1e-9);   /* a little gap is left */
    EXPECT_NEAR(19.65, iv[1].end, 1e-9);    /* its length of 10 s is kept */
}

UTEST(viamd_movie_keys, a_stretch_keeps_a_least_length) {
    std::vector<RepKey> keys = {vis_key(1, 0, false), vis_key(1, 10, true), vis_key(1, 20, false)};
    const RepInterval iv = rep_shown_intervals(keys, 1, 50.0)[0];
    rep_move_interval(&keys, 1, iv, 10.0, 5.0, 50.0);
    const std::vector<RepInterval> after = rep_shown_intervals(keys, 1, 50.0);
    EXPECT_GT(after[0].end - after[0].begin, 0.0);
}

UTEST(viamd_movie_keys, a_stretch_is_added_where_it_is_hidden_and_ends_before_the_next) {
    std::vector<RepKey> keys = {vis_key(8, 0, false), vis_key(8, 2.4, true), vis_key(8, 9.6, false), vis_key(8, 60, true)};
    EXPECT_TRUE(rep_add_interval(&keys, 8, 20.0, 80.0, 100.0));
    std::vector<RepInterval> iv = rep_shown_intervals(keys, 8, 100.0);
    ASSERT_EQ(3, (int)iv.size());
    EXPECT_NEAR(20.0, iv[1].begin, 1e-9);
    EXPECT_NEAR(59.95, iv[1].end, 1e-9);    /* it stops a little before the next one */
    EXPECT_FALSE(rep_add_interval(&keys, 8, 3.0, 5.0, 100.0));   /* shown there already */
}

UTEST(viamd_movie_keys, a_stretch_added_to_a_representation_without_keys_is_hidden_before_it) {
    std::vector<RepKey> keys;
    EXPECT_TRUE(rep_add_interval(&keys, 4, 10.0, 20.0, 50.0));
    const std::vector<RepInterval> iv = rep_shown_intervals(keys, 4, 50.0);
    ASSERT_EQ(1, (int)iv.size());
    EXPECT_NEAR(10.0, iv[0].begin, 1e-9);
    EXPECT_NEAR(20.0, iv[0].end, 1e-9);
    float f = 1.0f;
    ASSERT_TRUE(rep_visible_factor(&f, keys.data(), keys.size(), 4, 5.0, 0.0));
    EXPECT_NEAR(0.0f, f, 1e-6);
}

UTEST(viamd_movie_keys, removing_a_stretch_leaves_the_others) {
    std::vector<RepKey> keys = {vis_key(8, 0, false), vis_key(8, 2.4, true), vis_key(8, 9.6, false), vis_key(8, 60, true)};
    rep_remove_interval(&keys, rep_shown_intervals(keys, 8, 80.0)[0]);
    const std::vector<RepInterval> iv = rep_shown_intervals(keys, 8, 80.0);
    ASSERT_EQ(1, (int)iv.size());
    EXPECT_NEAR(60.0, iv[0].begin, 1e-9);
}

UTEST(viamd_movie_keys, a_swap_hides_one_and_shows_the_other_at_the_same_time) {
    std::vector<RepKey> keys;    /* neither has keys: 'from' was shown, 'to' hidden */
    rep_swap_at(&keys, 1, 2, 10.0);
    const std::vector<RepInterval> a = rep_shown_intervals(keys, 1, 50.0);
    const std::vector<RepInterval> b = rep_shown_intervals(keys, 2, 50.0);
    ASSERT_EQ(1, (int)a.size());
    ASSERT_EQ(1, (int)b.size());
    EXPECT_NEAR(0.0, a[0].begin, 1e-9);
    EXPECT_NEAR(10.0, a[0].end, 1e-9);
    EXPECT_NEAR(10.0, b[0].begin, 1e-9);
    EXPECT_NEAR(50.0, b[0].end, 1e-9);
}

UTEST(viamd_movie_keys, a_swap_on_a_key_that_is_there_changes_it_instead_of_adding_one) {
    std::vector<RepKey> keys = {vis_key(1, 0, true), vis_key(1, 10, true)};
    rep_swap_at(&keys, 1, 2, 10.0);
    int at_ten = 0;
    for (const RepKey& k : keys) if (k.rep == 1 && fabs(k.time - 10.0) < 1e-9) { at_ten += 1; EXPECT_NEAR(0.0f, k.value[0], 1e-9); }
    EXPECT_EQ(1, at_ten);
}

UTEST(viamd_movie_keys, stretches_that_touch_or_overlap_are_one_in_the_union) {
    std::vector<RepInterval> in(3);
    in[0].begin = 10; in[0].end = 20;
    in[1].begin = 20; in[1].end = 30;
    in[2].begin = 40; in[2].end = 50;
    const std::vector<RepInterval> u = rep_union_intervals(in);
    ASSERT_EQ(2, (int)u.size());
    EXPECT_NEAR(10.0, u[0].begin, 1e-9);
    EXPECT_NEAR(30.0, u[0].end, 1e-9);
    EXPECT_NEAR(40.0, u[1].begin, 1e-9);
}

UTEST(viamd_movie_keys, a_name_is_split_at_its_first_hyphen) {
    std::string g, m;
    EXPECT_TRUE(rep_name_split("protein-cpk", &g, &m));
    EXPECT_STREQ("protein", g.c_str());
    EXPECT_STREQ("cpk", m.c_str());
    EXPECT_TRUE(rep_name_split("protein-cpk-blue", &g, &m));
    EXPECT_STREQ("protein", g.c_str());
    EXPECT_STREQ("cpk-blue", m.c_str());
    EXPECT_FALSE(rep_name_split("water", &g, &m));
    EXPECT_STREQ("water", g.c_str());
    EXPECT_TRUE(m.empty());
    EXPECT_FALSE(rep_name_split("-odd", &g, &m));
    EXPECT_STREQ("-odd", g.c_str());
    EXPECT_FALSE(rep_name_split("odd-", &g, &m));
    EXPECT_STREQ("odd-", g.c_str());
}

UTEST(viamd_movie_keys, system_rows_have_one_label_per_system_in_first_occurrence_order) {
    const auto rows = rep_system_rows({"protein-cartoon", "ligand", "protein-cpk", "water-vdw", "ligand-vdw"});
    ASSERT_EQ(3, (int)rows.size());
    EXPECT_STREQ("protein", rows[0].label.c_str());
    EXPECT_EQ(2, rows[0].members);
    EXPECT_EQ(0, rows[0].rep);
    EXPECT_STREQ("ligand", rows[1].label.c_str());
    EXPECT_EQ(2, rows[1].members);
    EXPECT_STREQ("water", rows[2].label.c_str());
    EXPECT_TRUE(rep_system_rows({}).empty());
}

UTEST(viamd_movie_keys, overlapping_blocks_stack_and_nonoverlapping_blocks_reuse_slots) {
    std::vector<RepBlock> blocks = {
        {0, 0, {0.0, 10.0}}, {1, 0, {5.0, 15.0}}, {2, 0, {10.0, 20.0}}, {0, 1, {20.0, 30.0}},
    };
    EXPECT_EQ(2, rep_pack_blocks(&blocks, 0.0));
    EXPECT_EQ(0, blocks[0].slot);
    EXPECT_EQ(1, blocks[1].slot);
    EXPECT_EQ(0, blocks[2].slot);
    EXPECT_EQ(0, blocks[3].slot);
    EXPECT_EQ(1, blocks[1].rep);
    EXPECT_EQ(1, blocks[3].interval_index);
}

UTEST(viamd_movie_keys, block_packing_accounts_for_transition_tails_and_unsorted_input) {
    std::vector<RepBlock> blocks = {{1, 0, {10.0, 20.0}}, {0, 0, {0.0, 10.0}}};
    EXPECT_EQ(2, rep_pack_blocks(&blocks, 2.0));
    EXPECT_EQ(0, blocks[1].slot);
    EXPECT_EQ(1, blocks[0].slot);
    std::vector<RepBlock> empty;
    EXPECT_EQ(1, rep_pack_blocks(&empty, 2.0));
}

UTEST(viamd_movie_keys, unkeyed_enabled_representations_cover_the_movie_but_hidden_ones_do_not) {
    const std::vector<RepKey> keys;
    const auto spans = rep_effective_intervals(keys, 1, 80.0, true);
    ASSERT_EQ(1, (int)spans.size());
    EXPECT_EQ(0.0, spans[0].begin);
    EXPECT_EQ(80.0, spans[0].end);
    EXPECT_EQ(-1, spans[0].begin_key);
    EXPECT_TRUE(rep_effective_intervals(keys, 1, 80.0, false).empty());
    const std::vector<RepKey> hidden = {vis_key(1, 0, false)};
    EXPECT_TRUE(rep_effective_intervals(hidden, 1, 80.0, true).empty());
}

UTEST(viamd_movie_keys, switching_a_block_merges_target_overlaps_and_preserves_other_representations) {
    std::vector<RepKey> keys = {
        vis_key(1, 0, false), vis_key(1, 10, true), vis_key(1, 20, false),
        vis_key(2, 0, false), vis_key(2, 15, true), vis_key(2, 25, false),
        vis_key(3, 0, true),
    };
    RepKey scale;
    scale.rep = 2;
    scale.prop = (int)RepProp::Scale0;
    scale.value[0] = 2.5f;
    keys.push_back(scale);
    rep_transfer_interval(&keys, 1, 2, rep_shown_intervals(keys, 1, 80.0)[0], 80.0);
    EXPECT_TRUE(rep_shown_intervals(keys, 1, 80.0).empty());
    const auto target = rep_shown_intervals(keys, 2, 80.0);
    ASSERT_EQ(1, (int)target.size());
    EXPECT_EQ(10.0, target[0].begin);
    EXPECT_EQ(25.0, target[0].end);
    EXPECT_EQ(80.0, rep_shown_intervals(keys, 3, 80.0)[0].end);
    float value[3] = {};
    ASSERT_TRUE(rep_keys_evaluate(value, keys.data(), keys.size(), 2, (int)RepProp::Scale0, 12.0));
    EXPECT_EQ(2.5f, value[0]);
}

UTEST(viamd_movie_keys, removing_initial_or_only_blocks_keeps_the_system_hidden_in_the_gap) {
    std::vector<RepKey> keys = {vis_key(1, 0, true), vis_key(1, 10, false), vis_key(1, 20, true)};
    rep_remove_interval(&keys, rep_shown_intervals(keys, 1, 80.0)[0]);
    const auto spans = rep_shown_intervals(keys, 1, 80.0);
    ASSERT_EQ(1, (int)spans.size());
    EXPECT_EQ(20.0, spans[0].begin);
    rep_remove_interval(&keys, spans[0]);
    EXPECT_TRUE(rep_effective_intervals(keys, 1, 80.0, true).empty());
    std::vector<RepKey> only = {vis_key(2, 0, true)};
    rep_remove_interval(&only, rep_shown_intervals(only, 2, 80.0)[0]);
    EXPECT_TRUE(rep_effective_intervals(only, 2, 80.0, true).empty());
}

UTEST(viamd_movie_keys, switching_representation_blocks_is_undoable_without_changing_other_keys) {
    MovieKeys cur;
    cur.reps = {vis_key(1, 0, false), vis_key(1, 10, true), vis_key(1, 20, false)};
    const MovieKeys before = cur;
    MovieHistory history;
    history.clear(cur);
    rep_transfer_interval(&cur.reps, 1, 2, rep_shown_intervals(cur.reps, 1, 80.0)[0], 80.0);
    const MovieKeys after = cur;
    history.update(cur, false);
    ASSERT_TRUE(history.undo(&cur));
    EXPECT_TRUE(movie_keys_equal(before, cur));
    ASSERT_TRUE(history.redo(&cur));
    EXPECT_TRUE(movie_keys_equal(after, cur));
}

static CameraKeyframe cam_key(double time, float distance, uint32_t look_set = 0) {
    CameraKeyframe k = {};
    k.time = time;
    k.transform.distance = distance;
    k.look_set = look_set;
    return k;
}

UTEST(viamd_movie_keys, camera_bands_are_runs_of_keys_that_look_at_the_same_atoms) {
    std::vector<CameraKeyframe> keys = {
        cam_key(0, 10), cam_key(1, 10, 1), cam_key(2, 10, 1), cam_key(3, 10, 2), cam_key(4, 10, 2), cam_key(5, 10), cam_key(6, 10, 1),
    };
    const std::vector<CameraBand> bands = camera_bands(keys);
    ASSERT_EQ(bands.size(), (size_t)3);
    EXPECT_EQ(bands[0].kind, CameraBandKind::LookAtSet);
    EXPECT_EQ(bands[0].set, (uint32_t)1);
    EXPECT_EQ(bands[0].first, 1);
    EXPECT_EQ(bands[0].last, 2);
    EXPECT_EQ(bands[1].kind, CameraBandKind::LookAtSet);
    EXPECT_EQ(bands[1].set, (uint32_t)2);
    EXPECT_EQ(bands[1].begin, 3.0);
    EXPECT_EQ(bands[1].end, 4.0);
    EXPECT_EQ(bands[2].set, (uint32_t)1);
    EXPECT_EQ(bands[2].first, 6);
    EXPECT_EQ(bands[2].last, 6);
}

UTEST(viamd_movie_keys, a_spin_is_a_band_over_the_stretch_leading_to_its_key) {
    std::vector<CameraKeyframe> keys = {cam_key(0, 10), cam_key(2, 10), cam_key(5, 10)};
    keys[0].spin_turns = 3;   // Nothing leads to the first key
    keys[2].spin_turns = -2;
    const std::vector<CameraBand> bands = camera_bands(keys);
    ASSERT_EQ(bands.size(), (size_t)1);
    EXPECT_EQ(bands[0].kind, CameraBandKind::Spin);
    EXPECT_EQ(bands[0].begin, 2.0);
    EXPECT_EQ(bands[0].end, 5.0);
    EXPECT_EQ(bands[0].turns, -2);
    EXPECT_TRUE(camera_bands({}).empty());
}

UTEST(viamd_movie_keys, a_key_is_labelled_by_its_number_and_its_name) {
    CameraKeyframe k = {};
    EXPECT_STREQ(camera_key_label(k, 0).c_str(), "1");
    strcpy(k.name, "close-up");
    EXPECT_STREQ(camera_key_label(k, 2).c_str(), "3 close-up");
}

UTEST(viamd_movie_keys, a_key_on_the_path_does_not_move_the_camera) {
    std::vector<CameraKeyframe> keys = {cam_key(0, 10), cam_key(4, 30), cam_key(8, 20)};
    keys[1].spin_turns = 1;
    keys[1].use_frame = true;
    keys[1].focus_on = true;
    strcpy(keys[1].name, "mid");
    ViewTransform vt;
    float fov;
    camera_keyframes_evaluate(&vt, &fov, keys.data(), keys.size(), 2.0, false);
    const CameraKeyframe k = camera_key_on_path(vt, fov, keys, 2.0);
    EXPECT_EQ(k.time, 2.0);
    EXPECT_EQ(k.transform.distance, vt.distance);
    EXPECT_EQ(k.fov_y, fov);
    EXPECT_EQ(k.spin_turns, 0);
    EXPECT_FALSE(k.use_frame);
    EXPECT_FALSE(k.focus_on);   // The focus holds from the key before: nothing to say here
    EXPECT_EQ(k.name[0], '\0');
    EXPECT_EQ(k.look_set, (uint32_t)0);
}

UTEST(viamd_movie_keys, a_key_between_keys_that_look_at_the_same_atoms_looks_at_them_too) {
    std::vector<CameraKeyframe> keys = {cam_key(0, 10, 5), cam_key(4, 10, 5), cam_key(8, 10)};
    const ViewTransform pose = {};
    EXPECT_EQ(camera_key_on_path(pose, 0.8f, keys, 1.0).look_set, (uint32_t)5);
    // Between keys that look at different things it looks at a point, and so it does when the keys on one side look at none
    EXPECT_EQ(camera_key_on_path(pose, 0.8f, keys, 6.0).look_set, (uint32_t)0);
    // Outside the keys it holds what the end key does
    keys[2].look_set = 5;
    EXPECT_EQ(camera_key_on_path(pose, 0.8f, keys, 10.0).look_set, (uint32_t)5);
}

UTEST(viamd_movie_keys, renaming_a_key_is_an_edit) {
    MovieKeys a = keys_with(2), b = keys_with(2);
    EXPECT_TRUE(movie_keys_equal(a, b));
    strcpy(b.camera[1].name, "end");
    EXPECT_FALSE(movie_keys_equal(a, b));
}

UTEST(viamd_movie_keys, the_path_is_found_between_its_samples_and_held_outside) {
    CameraPathSamples s;
    s.time = {0.0, 2.0, 6.0};
    s.eye = {{0, 0, 0}, {2, 0, 0}, {2, 4, 0}};
    s.look = {{0, 0, 1}, {0, 0, 1}, {0, 0, 5}};
    vec3_t e, l;
    camera_path_at(s, 1.0, &e, &l);
    EXPECT_NEAR(e.x, 1.0f, 1.0e-6f);
    camera_path_at(s, 4.0, &e, &l);
    EXPECT_NEAR(e.y, 2.0f, 1.0e-6f);
    EXPECT_NEAR(l.z, 3.0f, 1.0e-6f);
    camera_path_at(s, -3.0, &e, &l);
    EXPECT_EQ(e.x, 0.0f);
    camera_path_at(s, 99.0, &e, &l);
    EXPECT_EQ(e.y, 4.0f);
}

UTEST(viamd_movie_keys, ticks_have_a_round_step_that_gives_few_enough_of_them) {
    EXPECT_EQ(camera_tick_step(10.0, 30), 0.5);
    EXPECT_EQ(camera_tick_step(30.0, 30), 1.0);
    EXPECT_EQ(camera_tick_step(60.0, 30), 2.0);
    EXPECT_EQ(camera_tick_step(300.0, 30), 10.0);
    EXPECT_EQ(camera_tick_step(0.0, 30), 1.0);
    for (double span : {0.7, 3.3, 12.0, 47.0, 1000.0}) EXPECT_LE(span / camera_tick_step(span, 20), 20.0);
}

UTEST(viamd_movie_keys, a_line_is_cut_at_the_camera) {
    vec4_t a = {0, 0, 0, 2.0f}, b = {4, 0, 0, -2.0f};
    EXPECT_TRUE(clip_segment_near(&a, &b));
    EXPECT_EQ(a.w, 2.0f);
    EXPECT_NEAR(b.w, 1.0e-3f, 1.0e-6f);
    EXPECT_NEAR(b.x, 2.0f, 0.01f);
    vec4_t c = {0, 0, 0, 1.0f}, d = {1, 0, 0, 3.0f};
    EXPECT_TRUE(clip_segment_near(&c, &d));
    EXPECT_EQ(d.x, 1.0f);
    vec4_t e = {0, 0, 0, -1.0f}, f = {1, 0, 0, -3.0f};
    EXPECT_FALSE(clip_segment_near(&e, &f));
}

UTEST(viamd_movie_keys, the_spin_ring_goes_through_the_eye_and_turns_the_way_of_the_axis) {
    CameraKeyframe from = cam_key(0, 10), to = cam_key(4, 10);
    to.spin_turns = 1;
    to.spin_axis = SpinAxis::WorldY;
    to.transform.position = {3, 5, 4};   // The look-at point is 10 in front of it, along -z in the camera's own frame
    vec3_t c, u, v;
    float r;
    ASSERT_TRUE(camera_spin_ring(from, to, &c, &u, &v, &r));
    const vec3_t eye = c + u * r;
    EXPECT_NEAR(eye.x, 3.0f, 1.0e-4f);
    EXPECT_NEAR(eye.y, 5.0f, 1.0e-4f);
    EXPECT_NEAR(eye.z, 4.0f, 1.0e-4f);
    EXPECT_NEAR(vec3_dot(u, vec3_t{0, 1, 0}), 0.0f, 1.0e-5f);
    // A positive turn is counter-clockwise seen from the tip of the axis
    const vec3_t axis = vec3_cross(u, v);
    EXPECT_NEAR(axis.y, 1.0f, 1.0e-5f);

    // The eye on the axis has no ring
    to.transform.orientation = quat_axis_angle(vec3_t{1, 0, 0}, -1.5707963f);   // Looks straight down, so the eye is above what it looks at
    EXPECT_FALSE(camera_spin_ring(from, to, &c, &u, &v, &r));
}

UTEST(viamd_movie_keys, a_ray_meets_a_plane_in_front_of_it) {
    vec3_t p;
    EXPECT_TRUE(ray_plane_hit({0, 0, 0}, {0, 0, -1}, {0, 0, -5}, {0, 0, 1}, &p));
    EXPECT_NEAR(p.z, -5.0f, 1.0e-6f);
    EXPECT_FALSE(ray_plane_hit({0, 0, 0}, {1, 0, 0}, {0, 0, -5}, {0, 0, 1}, &p));
    EXPECT_FALSE(ray_plane_hit({0, 0, 0}, {0, 0, 1}, {0, 0, -5}, {0, 0, 1}, &p));
}

UTEST(viamd_movie_keys, the_nearest_point_of_a_polyline_skips_what_cannot_be_used) {
    std::vector<vec2_t> pts = {{0, 0}, {10, 0}, {10, 10}, {20, 10}};
    std::vector<char> ok = {1, 1, 1, 1};
    int seg = -1;
    float u = -1.0f;
    float d = polyline_nearest(pts, ok, {4, 3}, &seg, &u);
    EXPECT_NEAR(d, 3.0f, 1.0e-5f);
    EXPECT_EQ(seg, 0);
    EXPECT_NEAR(u, 0.4f, 1.0e-5f);
    d = polyline_nearest(pts, ok, {12, 7}, &seg, &u);
    EXPECT_EQ(seg, 1);
    ok[2] = 0;   // The middle segments lose an end
    d = polyline_nearest(pts, ok, {12, 7}, &seg, &u);
    EXPECT_EQ(seg, 0);
    EXPECT_NEAR(d, 7.0f * 1.0f + 0.0f, 8.0f);
    EXPECT_EQ(polyline_nearest({}, {}, {0, 0}, &seg, &u), FLT_MAX);
}

UTEST(viamd_movie_keys, a_key_is_edited_by_its_eye_or_by_what_it_looks_at) {
    CameraKeyframe k = cam_key(0, 10);
    k.transform.position = {0, 0, 10};   // Looks along -z at the origin
    const vec3_t look = camera_get_look_at(k.transform);

    CameraKeyframe a = k;
    ASSERT_TRUE(camera_key_set_eye(&a, {5, 0, 10}));
    const vec3_t look_a = camera_get_look_at(a.transform);
    EXPECT_NEAR(look_a.x, look.x, 1.0e-3f);
    EXPECT_NEAR(look_a.y, look.y, 1.0e-3f);
    EXPECT_NEAR(look_a.z, look.z, 1.0e-3f);
    EXPECT_NEAR(a.transform.position.x, 5.0f, 1.0e-6f);

    CameraKeyframe b = k;
    ASSERT_TRUE(camera_key_set_look(&b, {4, 0, 0}));
    EXPECT_NEAR(b.transform.position.z, 10.0f, 1.0e-6f);
    const vec3_t look_b = camera_get_look_at(b.transform);
    EXPECT_NEAR(look_b.x, 4.0f, 1.0e-3f);

    CameraKeyframe c = k;
    c.follow = true;
    c.follow_center = {1, 1, 1};
    camera_key_translate(&c, {1, 2, 3});
    const vec3_t look_c = camera_get_look_at(c.transform);
    EXPECT_NEAR(look_c.x, look.x + 1.0f, 1.0e-4f);
    EXPECT_NEAR(look_c.z, look.z + 3.0f, 1.0e-4f);
    EXPECT_NEAR(c.follow_center.y, 3.0f, 1.0e-6f);

    // The eye on the point it looks at cannot aim
    CameraKeyframe d = k;
    EXPECT_FALSE(camera_key_set_eye(&d, look));
}

UTEST(viamd_movie_keys, roll_and_keeping_upright_are_edits_that_can_be_undone) {
    MovieKeys a, b;
    a.camera = {cam_key(0, 10)};
    b = a;
    b.camera[0].roll = 0.2f;
    EXPECT_FALSE(movie_keys_equal(a, b));
    b = a;
    b.keep_upright = !a.keep_upright;
    EXPECT_FALSE(movie_keys_equal(a, b));
    b = a;
    b.up_axis = 2;
    EXPECT_FALSE(movie_keys_equal(a, b));
}

UTEST(viamd_movie_keys, a_key_on_an_upright_path_takes_the_roll_there) {
    std::vector<CameraKeyframe> keys = {cam_key(0, 10), cam_key(4, 30)};
    keys[0].roll = 0.0f;
    keys[1].roll = 0.6f;
    const vec3_t up = {0, 1, 0};
    ViewTransform vt;
    float fov;
    camera_keyframes_evaluate(&vt, &fov, keys.data(), keys.size(), 4.0, false, nullptr, nullptr, &up);
    const CameraKeyframe k = camera_key_on_path(vt, fov, keys, 4.0, &up);
    EXPECT_NEAR(k.roll, 0.6f, 1.0e-4f);
    EXPECT_EQ(camera_key_on_path(vt, fov, keys, 4.0).roll, 0.0f);
}

/* Selecting and moving keys together */

static RepKey rk(uint32_t rep, int prop, double time, float value = 1.0f) {
    RepKey k;
    k.rep = rep;
    k.prop = prop;
    k.time = time;
    k.value[0] = value;
    return k;
}

static MovieKeys select_keys() {
    MovieKeys k;
    k.duration = 20.0f;
    for (double t : {2.0, 5.0, 9.0, 14.0}) k.camera.push_back(cam_key(t, 10));
    k.params.push_back(pk(1, 3.0, 1.0f));
    k.params.push_back(pk(1, 7.0, 2.0f));
    k.params.push_back(pk(2, 3.0, 5.0f));
    k.reps.push_back(rk(4, 2, 6.0, 0.5f));
    return k;
}

UTEST(viamd_movie_keys, a_selection_picks_keys_by_what_they_belong_to_and_their_time) {
    KeySelection s;
    s.add(KeyKind::Camera, 0, 5.0);
    s.add(KeyKind::Camera, 0, 5.0);
    s.add(KeyKind::Param, 1, 5.0);
    EXPECT_EQ(s.size(), (size_t)2);
    EXPECT_TRUE(s.contains(KeyKind::Param, 1, 5.0));
    EXPECT_FALSE(s.contains(KeyKind::Param, 2, 5.0));
    s.toggle(KeyKind::Camera, 0, 5.0);
    EXPECT_FALSE(s.contains(KeyKind::Camera, 0, 5.0));
    s.toggle(KeyKind::Camera, 0, 6.0);
    EXPECT_TRUE(s.contains(KeyKind::Camera, 0, 6.0));
    s.set(KeyKind::Rep, rep_key_subject(4, 2), 1.0);
    EXPECT_EQ(s.size(), (size_t)1);
}

UTEST(viamd_movie_keys, a_selection_forgets_keys_that_are_gone_and_can_take_all_of_a_lane) {
    const MovieKeys keys = select_keys();
    KeySelection s;
    key_selection_all(&s, keys, 1, -1);
    EXPECT_EQ(s.size(), (size_t)6);   // 4 camera keys and 2 of parameter 1
    key_selection_all(&s, keys, -1, rep_key_subject(4, 2));
    EXPECT_EQ(s.size(), (size_t)7);
    s.add(KeyKind::Camera, 0, 99.0);
    key_selection_prune(&s, keys);
    EXPECT_EQ(s.size(), (size_t)7);
    EXPECT_FALSE(s.contains(KeyKind::Camera, 0, 99.0));
}

UTEST(viamd_movie_keys, selected_keys_move_together_and_the_others_stay) {
    const MovieKeys start = select_keys();
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 5.0);
    sel.add(KeyKind::Camera, 0, 9.0);
    sel.add(KeyKind::Param, 1, 7.0);
    MovieKeys out;
    KeySelection out_sel;
    const double used = movie_keys_shift(&out, &out_sel, start, sel, 1.5, KeyShift{}, 20.0);
    EXPECT_EQ(used, 1.5);
    EXPECT_EQ(out.camera[0].time, 2.0);
    EXPECT_EQ(out.camera[1].time, 6.5);
    EXPECT_EQ(out.camera[2].time, 10.5);
    EXPECT_EQ(out.camera[3].time, 14.0);
    EXPECT_EQ(out.params[0].time, 3.0);
    EXPECT_EQ(out.params[1].time, 8.5);
    EXPECT_EQ(out.params[2].time, 3.0);
    EXPECT_EQ(out.reps[0].time, 6.0);
    EXPECT_TRUE(out_sel.contains(KeyKind::Camera, 0, 6.5));
    EXPECT_TRUE(out_sel.contains(KeyKind::Param, 1, 8.5));
    EXPECT_FALSE(out_sel.contains(KeyKind::Camera, 0, 5.0));
    // What it was made from is left alone
    EXPECT_EQ(start.camera[1].time, 5.0);
}

UTEST(viamd_movie_keys, a_group_stops_at_the_ends_of_the_movie_as_a_whole) {
    const MovieKeys start = select_keys();
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 9.0);
    sel.add(KeyKind::Camera, 0, 14.0);
    MovieKeys out;
    KeySelection out_sel;
    double used = movie_keys_shift(&out, &out_sel, start, sel, 100.0, KeyShift{}, 20.0);
    EXPECT_EQ(used, 6.0);
    EXPECT_EQ(out.camera[2].time, 15.0);   // The spacing of the group is kept
    EXPECT_EQ(out.camera[3].time, 20.0);
    used = movie_keys_shift(&out, &out_sel, start, sel, -100.0, KeyShift{}, 20.0);
    EXPECT_EQ(used, -9.0);
    EXPECT_EQ(out.camera[2].time, 0.0);
    EXPECT_EQ(out.camera[3].time, 5.0);
    // Nothing selected: nothing moves
    KeySelection none;
    used = movie_keys_shift(&out, &out_sel, start, none, 3.0, KeyShift{}, 20.0);
    EXPECT_EQ(used, 0.0);
    EXPECT_EQ(out.camera[1].time, 5.0);
    EXPECT_EQ(out.params[1].time, 7.0);
}

UTEST(viamd_movie_keys, the_values_of_a_lane_move_with_the_keys_and_stay_in_their_limits) {
    MovieKeys start = select_keys();
    start.camera[0].use_frame = true;
    start.camera[0].frame = 40.0;
    start.camera[1].use_frame = true;
    start.camera[1].frame = 90.0;
    start.camera[0].fov_y = 0.5f;
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 2.0);
    sel.add(KeyKind::Camera, 0, 5.0);
    sel.add(KeyKind::Param, 1, 3.0);
    sel.add(KeyKind::Rep, rep_key_subject(4, 2), 6.0);
    MovieKeys out;
    KeySelection out_sel;

    KeyShift frame;
    frame.lane = KeyLane::Frame;
    frame.dy = 25.0;
    frame.lo = 0.0;
    frame.hi = 100.0;
    movie_keys_shift(&out, &out_sel, start, sel, 0.0, frame, 20.0);
    EXPECT_EQ(out.camera[0].frame, 65.0);
    EXPECT_EQ(out.camera[1].frame, 100.0);   // Limited

    KeyShift fov;
    fov.lane = KeyLane::Fov;
    fov.dy = 1000.0;
    movie_keys_shift(&out, &out_sel, start, sel, 0.0, fov, 20.0);
    EXPECT_NEAR(out.camera[0].fov_y, 170.0f * 3.14159265f / 180.0f, 1.0e-4f);

    KeyShift param;
    param.lane = KeyLane::Param;
    param.subject = 1;
    param.dy = 3.0;
    param.lo = 0.0;
    param.hi = 10.0;
    movie_keys_shift(&out, &out_sel, start, sel, 0.0, param, 20.0);
    EXPECT_EQ(out.params[0].value[0], 4.0f);
    EXPECT_EQ(out.params[1].value[0], 2.0f);   // Not selected
    EXPECT_EQ(out.params[2].value[0], 5.0f);   // Another parameter

    KeyShift rep;
    rep.lane = KeyLane::Rep;
    rep.subject = rep_key_subject(4, 2);
    rep.dy = 2.0;
    rep.ratio = true;
    rep.lo = 0.0;
    rep.hi = 10.0;
    movie_keys_shift(&out, &out_sel, start, sel, 0.0, rep, 20.0);
    EXPECT_EQ(out.reps[0].value[0], 1.0f);
}

UTEST(viamd_movie_keys, a_distance_moves_along_the_line_of_sight_and_the_camera_keeps_looking_at_the_same_point) {
    MovieKeys start;
    CameraKeyframe k = cam_key(1.0, 10.0f);
    k.transform.position = {2, 3, 10};
    k.look_set = 7;
    k.focus_on = true;
    k.focus_target = FocusTarget::Selection;
    k.focus_set = 8;
    k.focus_point = {4, 5, 6};
    start.camera.push_back(k);
    const vec3_t look = camera_get_look_at(k.transform);
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 1.0);
    KeyShift s;
    s.lane = KeyLane::Distance;
    s.unit = 2.0;    // Shown in units that are twice the internal ones
    s.dy = 6.0;      // 6 shown = 3 internal
    MovieKeys out;
    KeySelection out_sel;
    movie_keys_shift(&out, &out_sel, start, sel, 0.0, s, 20.0);
    EXPECT_NEAR(out.camera[0].transform.distance, 13.0f, 1.0e-5f);
    const vec3_t look_after = camera_get_look_at(out.camera[0].transform);
    EXPECT_NEAR(look_after.x, look.x, 1.0e-4f);
    EXPECT_NEAR(look_after.y, look.y, 1.0e-4f);
    EXPECT_NEAR(look_after.z, look.z, 1.0e-4f);
    EXPECT_EQ(out.camera[0].look_set, 7u);
    EXPECT_TRUE(out.camera[0].focus_on);
    EXPECT_EQ(out.camera[0].focus_target, FocusTarget::Selection);
    EXPECT_EQ(out.camera[0].focus_set, 8u);
    EXPECT_NEAR(out.camera[0].focus_point.x, 4.0f, 1.0e-5f);
    EXPECT_NEAR(out.camera[0].focus_point.y, 5.0f, 1.0e-5f);
    EXPECT_NEAR(out.camera[0].focus_point.z, 6.0f, 1.0e-5f);
    s.dy = -1000.0;
    movie_keys_shift(&out, &out_sel, start, sel, 0.0, s, 20.0);
    EXPECT_GT(out.camera[0].transform.distance, 0.0f);
}

UTEST(viamd_movie_keys, a_moved_key_that_lands_on_another_replaces_it_and_everything_is_sorted) {
    MovieKeys keys = select_keys();
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 5.0);
    sel.add(KeyKind::Param, 1, 3.0);
    keys.camera[1].time = 9.0;           // On the camera key at 9
    keys.camera[1].fov_y = 0.123f;
    sel.ids[0].time = 9.0;
    keys.params[0].time = 7.0;           // On the parameter key at 7
    keys.params[0].value[0] = 42.0f;
    sel.ids[1].time = 7.0;
    keys.camera.push_back(cam_key(1.0, 10));   // Out of order
    movie_keys_resolve(&keys, sel);
    ASSERT_EQ(keys.camera.size(), (size_t)4);
    EXPECT_EQ(keys.camera[0].time, 1.0);
    EXPECT_EQ(keys.camera[2].time, 9.0);
    EXPECT_EQ(keys.camera[2].fov_y, 0.123f);   // The selected one stayed
    ASSERT_EQ(keys.params.size(), (size_t)2);
    EXPECT_EQ(keys.params[0].value[0], 42.0f);
    EXPECT_EQ(keys.params[1].param, 2);
}

UTEST(viamd_movie_keys, selected_keys_are_deleted) {
    MovieKeys keys = select_keys();
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 2.0);
    sel.add(KeyKind::Param, 2, 3.0);
    sel.add(KeyKind::Rep, rep_key_subject(4, 2), 6.0);
    movie_keys_delete(&keys, &sel);
    EXPECT_EQ(keys.camera.size(), (size_t)3);
    EXPECT_EQ(keys.params.size(), (size_t)2);
    EXPECT_TRUE(keys.reps.empty());
    EXPECT_TRUE(sel.empty());
}

UTEST(viamd_movie_keys, a_copy_is_put_in_with_its_spacing_and_inside_the_movie) {
    MovieKeys keys = select_keys();
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 2.0);
    sel.add(KeyKind::Camera, 0, 5.0);
    sel.add(KeyKind::Param, 1, 7.0);
    const KeyClip clip = movie_keys_copy(keys, sel);
    EXPECT_EQ(clip.camera.size(), (size_t)2);
    EXPECT_EQ(clip.params.size(), (size_t)1);
    EXPECT_EQ(clip.begin, 2.0);

    KeySelection pasted;
    movie_keys_paste(&keys, &pasted, clip, 10.0, 20.0);
    EXPECT_EQ(keys.camera.size(), (size_t)6);
    EXPECT_TRUE(pasted.contains(KeyKind::Camera, 0, 10.0));
    EXPECT_TRUE(pasted.contains(KeyKind::Camera, 0, 13.0));
    EXPECT_TRUE(pasted.contains(KeyKind::Param, 1, 15.0));
    EXPECT_EQ(pasted.size(), (size_t)3);

    // Near the end it is moved earlier as a whole
    KeySelection late;
    movie_keys_paste(&keys, &late, clip, 19.0, 20.0);
    EXPECT_TRUE(late.contains(KeyKind::Camera, 0, 12.0) || late.contains(KeyKind::Camera, 0, 15.0));
    for (const KeyId& id : late.ids) EXPECT_LE(id.time, 20.0);

    // An empty clip changes nothing
    const size_t count = keys.camera.size();
    movie_keys_paste(&keys, &late, KeyClip{}, 3.0, 20.0);
    EXPECT_EQ(keys.camera.size(), count);
}

static MovieOverlay bar(double begin, double end) {
    MovieOverlay o;
    o.begin = begin;
    o.end = end;
    return o;
}

UTEST(viamd_movie_keys, overlay_bars_move_with_what_is_timed_inside_them) {
    MovieKeys start;
    start.duration = 20.0f;
    start.overlays = {bar(1.0, 4.0), bar(6.0, 9.0), bar(10.0, 12.0)};
    MoviePlotPanel panel;
    panel.begin = 7.0;
    panel.end = 8.0;
    start.overlays[1].panels.push_back(panel);
    start.camera.push_back(cam_key(5.0, 10));
    KeySelection sel;
    sel.add(KeyKind::Overlay, 1, 6.0, 9.0);
    sel.add(KeyKind::Camera, 0, 5.0);
    MovieKeys out;
    KeySelection out_sel;
    double used = movie_keys_shift(&out, &out_sel, start, sel, 2.5, KeyShift{}, 20.0);
    EXPECT_EQ(used, 2.5);
    EXPECT_EQ(out.overlays[0].begin, 1.0);
    EXPECT_EQ(out.overlays[1].begin, 8.5);
    EXPECT_EQ(out.overlays[1].end, 11.5);
    EXPECT_EQ(out.overlays[1].panels[0].begin, 9.5);
    EXPECT_EQ(out.overlays[1].panels[0].end, 10.5);
    EXPECT_EQ(out.camera[0].time, 7.5);
    EXPECT_TRUE(out_sel.contains(KeyKind::Overlay, 1, 8.5));
    EXPECT_EQ(out_sel.ids[0].end, 11.5);

    // The end of the bar is what stops it at the end of the movie
    used = movie_keys_shift(&out, &out_sel, start, sel, 100.0, KeyShift{}, 20.0);
    EXPECT_EQ(used, 11.0);
    EXPECT_EQ(out.overlays[1].end, 20.0);
    used = movie_keys_shift(&out, &out_sel, start, sel, -100.0, KeyShift{}, 20.0);
    EXPECT_EQ(used, -5.0);
    EXPECT_EQ(out.overlays[1].begin, 1.0);
}

UTEST(viamd_movie_keys, representation_blocks_move_together_and_follow_their_keys) {
    MovieKeys start;
    start.duration = 20.0f;
    start.reps = {rk(7, 0, 0.0, 0.0f), rk(7, 0, 2.0, 1.0f), rk(7, 0, 4.0, 0.0f), rk(9, 0, 0.0, 0.0f), rk(9, 0, 10.0, 1.0f), rk(9, 0, 14.0, 0.0f)};
    KeySelection sel;
    sel.add(KeyKind::Block, 7, 2.0, 4.0);
    sel.add(KeyKind::Block, 9, 10.0, 14.0);
    MovieKeys out;
    KeySelection out_sel;
    double used = movie_keys_shift(&out, &out_sel, start, sel, 3.0, KeyShift{}, 20.0);
    EXPECT_EQ(used, 3.0);
    EXPECT_EQ(out.reps[1].time, 5.0);
    EXPECT_EQ(out.reps[2].time, 7.0);
    EXPECT_EQ(out.reps[4].time, 13.0);
    EXPECT_EQ(out.reps[5].time, 17.0);
    EXPECT_TRUE(out_sel.contains(KeyKind::Block, 7, 5.0));
    EXPECT_TRUE(out_sel.contains(KeyKind::Block, 9, 13.0));
    EXPECT_EQ(out_sel.ids[1].end, 17.0);

    // The group stops at the end of the movie by the block that is the furthest
    used = movie_keys_shift(&out, &out_sel, start, sel, 100.0, KeyShift{}, 20.0);
    EXPECT_EQ(used, 6.0);
    EXPECT_EQ(out.reps[5].time, 20.0);
    EXPECT_EQ(out.reps[2].time, 10.0);

    // A block of the selection that no stretch is there for is left alone
    KeySelection stale;
    stale.add(KeyKind::Block, 7, 3.0, 4.0);
    used = movie_keys_shift(&out, &out_sel, start, stale, 3.0, KeyShift{}, 20.0);
    EXPECT_EQ(out.reps[1].time, 2.0);
}

UTEST(viamd_movie_keys, a_picked_block_is_removed_with_its_keys_and_stale_bars_and_blocks_are_forgotten) {
    MovieKeys keys;
    keys.duration = 20.0f;
    keys.reps = {rk(7, 0, 0.0, 0.0f), rk(7, 0, 2.0, 1.0f), rk(7, 0, 4.0, 0.0f)};
    keys.overlays = {bar(1.0, 4.0)};
    KeySelection sel;
    sel.add(KeyKind::Block, 7, 2.0, 4.0);
    sel.add(KeyKind::Overlay, 0, 1.0, 4.0);
    sel.add(KeyKind::Overlay, 3, 1.0, 4.0);      // No such bar
    sel.add(KeyKind::Block, 7, 8.0, 9.0);        // No such stretch
    key_selection_prune(&sel, keys);
    EXPECT_EQ(sel.size(), (size_t)2);
    movie_keys_delete(&keys, &sel);
    for (const RepKey& k : keys.reps) EXPECT_LT(k.value[0], 0.5f);   // Nothing is shown any more
    EXPECT_EQ(keys.overlays.size(), (size_t)1);                      // The bar stays
    EXPECT_TRUE(sel.empty());
}

UTEST(viamd_movie_keys, the_picked_items_are_stretched_about_an_anchor) {
    MovieKeys start = select_keys();
    start.overlays = {bar(6.0, 9.0)};
    MoviePlotPanel panel;
    panel.begin = 7.0;
    panel.end = 8.0;
    start.overlays[0].panels.push_back(panel);
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 5.0);
    sel.add(KeyKind::Camera, 0, 9.0);
    sel.add(KeyKind::Param, 1, 7.0);
    sel.add(KeyKind::Overlay, 0, 6.0, 9.0);
    double t0 = 0.0, t1 = 0.0;
    ASSERT_TRUE(key_selection_extent(start, sel, &t0, &t1));
    EXPECT_EQ(t0, 5.0);
    EXPECT_EQ(t1, 9.0);

    MovieKeys out;
    KeySelection out_sel;
    double used = movie_keys_scale(&out, &out_sel, start, sel, 5.0, 2.0, 20.0);
    EXPECT_EQ(used, 2.0);
    EXPECT_EQ(out.camera[0].time, 2.0);     // Not picked
    EXPECT_EQ(out.camera[1].time, 5.0);     // The anchor stays
    EXPECT_EQ(out.camera[2].time, 13.0);
    EXPECT_EQ(out.params[1].time, 9.0);
    EXPECT_EQ(out.overlays[0].begin, 7.0);
    EXPECT_EQ(out.overlays[0].end, 13.0);
    EXPECT_EQ(out.overlays[0].panels[0].begin, 9.0);
    EXPECT_EQ(out.overlays[0].panels[0].end, 11.0);
    EXPECT_TRUE(out_sel.contains(KeyKind::Camera, 0, 13.0));
    EXPECT_TRUE(out_sel.contains(KeyKind::Overlay, 0, 7.0));
    EXPECT_EQ(out_sel.ids[3].end, 13.0);

    // It is kept from going out of the movie, and a squeeze keeps the order
    used = movie_keys_scale(&out, &out_sel, start, sel, 5.0, 100.0, 20.0);
    EXPECT_EQ(used, 15.0 / 4.0);
    EXPECT_EQ(out.camera[2].time, 20.0);
    used = movie_keys_scale(&out, &out_sel, start, sel, 5.0, 0.5, 20.0);
    EXPECT_EQ(out.camera[2].time, 7.0);

    // About a time inside the selection: both sides are limited
    used = movie_keys_scale(&out, &out_sel, start, sel, 8.0, 100.0, 20.0);
    EXPECT_NEAR(used, 8.0 / 3.0, 1.0e-9);
    EXPECT_NEAR(out.camera[1].time, 0.0, 1.0e-9);

    // Nothing picked: the factor is 1
    KeySelection none;
    used = movie_keys_scale(&out, &out_sel, start, none, 5.0, 3.0, 20.0);
    EXPECT_EQ(used, 1.0);
    EXPECT_EQ(out.camera[1].time, 5.0);
    EXPECT_FALSE(key_selection_extent(start, none, &t0, &t1));
}

UTEST(viamd_movie_keys, blocks_are_stretched_with_their_keys) {
    MovieKeys start;
    start.duration = 20.0f;
    start.reps = {rk(7, 0, 0.0, 0.0f), rk(7, 0, 2.0, 1.0f), rk(7, 0, 4.0, 0.0f), rk(9, 0, 0.0, 0.0f), rk(9, 0, 6.0, 1.0f), rk(9, 0, 8.0, 0.0f)};
    KeySelection sel;
    sel.add(KeyKind::Block, 7, 2.0, 4.0);
    sel.add(KeyKind::Block, 9, 6.0, 8.0);
    MovieKeys out;
    KeySelection out_sel;
    const double used = movie_keys_scale(&out, &out_sel, start, sel, 2.0, 2.0, 20.0);
    EXPECT_EQ(used, 2.0);
    EXPECT_EQ(out.reps[1].time, 2.0);
    EXPECT_EQ(out.reps[2].time, 6.0);    // Twice as long
    EXPECT_EQ(out.reps[4].time, 10.0);
    EXPECT_EQ(out.reps[5].time, 14.0);
    EXPECT_TRUE(out_sel.contains(KeyKind::Block, 9, 10.0));
    EXPECT_EQ(out_sel.ids[1].end, 14.0);
}

UTEST(viamd_movie_keys, the_easing_of_the_picked_keys_is_read_and_set_together) {
    MovieKeys keys = select_keys();
    keys.camera[0].ease = KeyEase::Hold;          // The first key has no easing to choose
    keys.camera[2].ease = KeyEase::Linear;
    keys.params[0].ease = KeyEase::Linear;
    keys.reps.push_back(rk(4, 0, 8.0, 1.0f));     // Visible holds
    KeySelection sel;
    sel.add(KeyKind::Camera, 0, 2.0);
    sel.add(KeyKind::Camera, 0, 9.0);
    sel.add(KeyKind::Param, 1, 3.0);
    sel.add(KeyKind::Rep, rep_key_subject(4, 0), 8.0);
    KeyEase common = KeyEase::Smooth;
    bool mixed = true;
    int count = key_selection_ease(keys.camera.data(), keys.camera.size(), keys.params, keys.reps, sel, &common, &mixed);
    EXPECT_EQ(count, 2);
    EXPECT_FALSE(mixed);
    EXPECT_EQ(common, KeyEase::Linear);

    sel.add(KeyKind::Camera, 0, 5.0);             // Smooth
    count = key_selection_ease(keys.camera.data(), keys.camera.size(), keys.params, keys.reps, sel, &common, &mixed);
    EXPECT_EQ(count, 3);
    EXPECT_TRUE(mixed);

    movie_keys_set_ease(&keys, sel, KeyEase::EaseInOut);
    EXPECT_EQ(keys.camera[0].ease, KeyEase::Hold);
    EXPECT_EQ(keys.camera[1].ease, KeyEase::EaseInOut);
    EXPECT_EQ(keys.camera[2].ease, KeyEase::EaseInOut);
    EXPECT_EQ(keys.camera[3].ease, KeyEase::Smooth);
    EXPECT_EQ(keys.params[0].ease, KeyEase::EaseInOut);
    EXPECT_EQ(keys.reps.back().ease, KeyEase::Smooth);

    KeySelection none;
    EXPECT_EQ(key_selection_ease(keys.camera.data(), keys.camera.size(), keys.params, keys.reps, none, &common, &mixed), 0);
}
