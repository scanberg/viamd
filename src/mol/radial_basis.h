#include <core/types.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/math_utils.h>

enum class RadialBasisFunctionType : uint8 { Wendland_3_1, Gaussian };

struct RadialBasis {
    struct ControlPoint {
        vec3 position;
        vec3 weight;
    };
    DynamicArray<ControlPoint> control_points = {};

    struct {
        RadialBasisFunctionType type = RadialBasisFunctionType::Wendland_3_1;
        float param = 0.f;
    } function;
};

// RBF functions
inline float Wendland_3_1(float r, float radius = 1.f) {
    r = r / radius;
    if (r <= 0.f) return 0.f;
    if (r >= 1.f) return 0.f;

    const float x = 1.f - r;
    const float x2 = x * x;
    return x2 * x2 * (4.f * r + 1.f);
}

inline float Gaussian(float r) { return math::exp(-math::pow(0.5f * r, 2)); }

void compute_radial_basis(RadialBasis* dst, Array<const vec3> points, Array<const vec3> values,
                          RadialBasisFunctionType type = RadialBasisFunctionType::Wendland_3_1, float function_param = 10.f);
RadialBasis compute_radial_basis(Array<const vec3> points, Array<const vec3> values,
                                 RadialBasisFunctionType type = RadialBasisFunctionType::Wendland_3_1, float function_param = 10.f);

inline vec3 evaluate_radial_basis(const RadialBasis& basis, const vec3& point) {
    vec3 result{0};
    switch (basis.function.type) {
        case RadialBasisFunctionType::Wendland_3_1:
            for (const auto& cp : basis.control_points) {
                result += cp.weight * Wendland_3_1(math::length(point - cp.position), basis.function.param);
            }
            break;
        case RadialBasisFunctionType::Gaussian:
            for (const auto& cp : basis.control_points) {
                result += cp.weight * Gaussian(math::length(point - cp.position));
            }
            break;
        default:
            break;
    }

    return result;
}
