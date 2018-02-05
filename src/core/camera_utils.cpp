#include <core/camera_utils.h>
#include <core/common.h>
#include <core/math_utils.h>

void camera_controller_fps(Camera* camera, bool key_fwd, bool key_bwd, bool key_lft, bool key_rht, vec2 mouse_vel, float delta_time, float move_speed,
                           float rot_speed) {
    ASSERT(camera);

    camera->orientation = quat(mouse_vel.x, vec3(0, 1, 0)) * camera->orientation * quat(mouse_vel.y, vec3(1, 0, 0));
    camera->orientation = normalize(camera->orientation);

    mat3 mat = glm::mat3_cast(camera->orientation);
    vec3 move_vec(0);

    if (key_fwd) {
        move_vec.z -= move_speed;
    }
    if (key_bwd) {
        move_vec.z += move_speed;
    }
    if (key_lft) {
        move_vec.x -= move_speed;
    }
    if (key_rht) {
        move_vec.x += move_speed;
    }

    camera->position += mat * move_vec;
}

void camera_controller_trackball(Camera* camera, bool lft_mouse, bool mdl_mouse, bool rht_mouse, vec2 mouse_vel) {
    
}