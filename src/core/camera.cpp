#include <core/camera.h>
#include <core/common.h>
#include <core/math_utils.h>

mat4 compute_view_to_world_matrix(const Camera& camera) {
    auto r = glm::mat4_cast(camera.orientation);
	auto t = glm::translate(mat4(1), camera.position);
    return  t * r;
}

mat4 compute_world_to_view_matrix(const Camera& camera) {
	auto r = glm::mat4_cast(glm::conjugate(camera.orientation));
	auto t = glm::translate(mat4(1), -camera.position);
	return r * t;
}

mat4 compute_perspective_projection_matrix(const Camera& camera, int width, int height) {
    float aspect = (float)width / (float)height;
    return glm::perspective(camera.fov_y, aspect, camera.near_plane, camera.far_plane);
}

void camera_look_at(Camera* camera, vec3 look_at, vec3 look_up) {
    ASSERT(camera);
    vec3 look = camera->position - look_at;
    float len2 = dot(look, look);
    if (len2 > 0.f) {
        look /= sqrtf(len2);
        vec3 right = normalize(cross(look_up, look));
        // TODO: Make sure look and look_up dont coincide
        look_up = cross(look, right);

        camera->orientation = glm::quat_cast(mat3(right, look_up, look));
    }

    // TODO: make sure look_up is kept.
}