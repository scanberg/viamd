#include <core/camera_utils.h>
#include <core/common.h>
#include <core/math_utils.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/euler_angles.hpp>

void camera_controller_fps(Camera* camera, bool key_fwd, bool key_bwd, bool key_lft, bool key_rht, vec2 mouse_vel, float delta_time, float move_speed,
                           float rot_speed) {
    ASSERT(camera);

	mouse_vel *= 0.01f;
	quat pitch = glm::quat(glm::vec3(-mouse_vel.y, 0, 0));
	quat yaw = glm::quat(glm::vec3(0, mouse_vel.x, 0));
	//glm::clamp(euler.x, -math::PI * 0.5f, math::PI * 0.5f);
	camera->orientation = yaw * camera->orientation * pitch;
	camera->orientation = glm::normalize(camera->orientation);
    //camera->orientation = quat(mouse_vel.x, vec3(0, 1, 0)) * camera->orientation * quat(mouse_vel.y, vec3(1, 0, 0));

    vec3 move_vec(0);
	move_speed *= 10.f * delta_time;

    if (key_fwd) move_vec.z -= move_speed;
    if (key_bwd) move_vec.z += move_speed;
    if (key_lft) move_vec.x -= move_speed;
    if (key_rht) move_vec.x += move_speed;

	mat3 mat = glm::mat3_cast(camera->orientation);
	camera->position += mat * move_vec;
}

void camera_controller_trackball(Camera* camera, bool lft_mouse, bool mdl_mouse, bool rht_mouse, vec2 mouse_vel, int screen_w, int screen_h) {
    
}