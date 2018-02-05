#pragma once

#include <core/camera.h>

void camera_controller_fps(Camera* camera, bool key_fwd, bool key_bwd, bool key_lft, bool key_rht, vec2 mouse_vel,
                           float delta_time, float move_speed = 1.f, float rot_speed = 1.f);

void camera_controller_trackball(Camera* camera, bool lft_mouse, bool mdl_mouse, bool rht_mouse, vec2 mouse_vel);