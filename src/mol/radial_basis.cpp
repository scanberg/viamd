#include "radial_basis.h"
#include <Eigen/Eigen>

namespace radial_basis {

// RBF functions
inline float Wendland_3_1(float r, float radius = 1.f) {
    r = r / radius;
    if (r <= 0.f) return 0.f;
    if (r >= 1.f) return 0.f;

    const float x = 1.f - r;
    const float x2 = x * x;
    return x2 * x2 * (4.f * r + 1.f);
}

inline float Gaussian(float r, float cutoff) {
    if (r > cutoff) return 0.f;
    return math::exp(-math::pow(0.5f * r, 2));
}

void compute_radial_basis(Array<vec3> rbf_weights, Array<const vec3> positions, Array<const vec3> values, FunctionType type, float cutoff) {
    ASSERT(rbf_weights.count == positions.count);
    ASSERT(positions.count == values.count);
    const int N = (int)rbf_weights.count;

    // Create matrix A NxN
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N, N);
    if (type == FunctionType::Wendland_3_1) {
        for (int32 i = 0; i < N; i++) {
            for (int32 j = 0; j < N; j++) {
                A(i, j) = Wendland_3_1(math::length(positions[i] - positions[j]), cutoff);
            }
        }
    } else if (type == FunctionType::Gaussian) {
        for (int32 i = 0; i < N; i++) {
            for (int32 j = 0; j < N; j++) {
                A(i, j) = Gaussian(math::length(positions[i] - positions[j]), cutoff);
            }
        }
    }

    // Create Vector b
    Eigen::MatrixXf b = Eigen::MatrixXf::Zero(N, 3);
    for (int32 i = 0; i < N; i++) {
        b.row(i) = Eigen::Vector3f(values[i].x, values[i].y, values[i].z);
    }

    Eigen::MatrixXf x = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    for (int32 i = 0; i < N; i++) {
        rbf_weights[i] = {x(i, 0), x(i, 1), x(i, 2)};
    }
}

vec3 evaluate_radial_basis(const vec3& point, Array<const vec3> rbf_positions, Array<const vec3> rbf_weights, FunctionType type, float cutoff) {
    ASSERT(rbf_positions.count == rbf_weights.count);
    const auto N = rbf_positions.count;
    vec3 result{0};

    switch (type) {
        case FunctionType::Wendland_3_1:
            for (int64 i = 0; i < N; i++) {
                result += rbf_weights[i] * Wendland_3_1(math::length(point - rbf_positions[i]), cutoff);
            }
            break;
        case FunctionType::Gaussian:
            for (int64 i = 0; i < N; i++) {
                result += rbf_weights[i] * Gaussian(math::length(point - rbf_positions[i]), cutoff);
            }
            break;
        default:
            break;
    }

    return result;
}

}  // namespace radial_basis
