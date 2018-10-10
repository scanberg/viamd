#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/math_utils.h>

namespace radial_basis {
enum class FunctionType : uint8 { Wendland_3_1, Gaussian };

void compute_radial_basis(Array<float> rbf_weights, Array<const vec3> positions, Array<const float> values,
                          FunctionType type = FunctionType::Wendland_3_1, float cutoff = 10.f);

void compute_radial_basis(Array<vec2> rbf_weights, Array<const vec3> positions, Array<const vec2> values,
                          FunctionType type = FunctionType::Wendland_3_1, float cutoff = 10.f);

void compute_radial_basis(Array<vec3> rbf_weights, Array<const vec3> positions, Array<const vec3> values,
                          FunctionType type = FunctionType::Wendland_3_1, float cutoff = 10.f);

void compute_radial_basis(Array<vec4> rbf_weights, Array<const vec3> positions, Array<const vec4> values,
                          FunctionType type = FunctionType::Wendland_3_1, float cutoff = 10.f);

float evaluate_radial_basis(const vec3& point, Array<const vec3> rbf_positions, Array<const float> rbf_wights, FunctionType type,
                            float cutoff = 10.f);
vec2 evaluate_radial_basis(const vec3& point, Array<const vec3> rbf_positions, Array<const vec2> rbf_wights, FunctionType type, float cutoff = 10.f);
vec3 evaluate_radial_basis(const vec3& point, Array<const vec3> rbf_positions, Array<const vec3> rbf_wights, FunctionType type, float cutoff = 10.f);
vec4 evaluate_radial_basis(const vec3& point, Array<const vec3> rbf_positions, Array<const vec4> rbf_wights, FunctionType type, float cutoff = 10.f);

}  // namespace radial_basis
