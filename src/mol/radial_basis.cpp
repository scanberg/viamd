#include "radial_basis.h"
#include <Eigen/Eigen>

void compute_radial_basis(RadialBasis* basis, Array<const vec3> points, Array<const vec3> values, RadialBasisFunctionType type,
                          float function_param) {
    ASSERT(basis);
    ASSERT(points.count == values.count);
    const int N = points.count;

    // Create matrix A NxN
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N, N);
    for (int32 i = 0; i < N; i++) {
        for (int32 j = 0; j < N; j++) {
            if (type == RadialBasisFunctionType::Wendland_3_1) {
                A(i, j) = Wendland_3_1(math::length(points[i] - points[j]), function_param);
            } else if (type == RadialBasisFunctionType::Gaussian) {
                A(i, j) = Gaussian(math::length(points[i] - points[j]));
            }
        }
    }

    // Create Vector b
    Eigen::MatrixXf b = Eigen::MatrixXf::Zero(N, 3);
    for (int32 i = 0; i < N; i++) {
        b.row(i) = Eigen::Vector3f(values[i].x, values[i].y, values[i].z);
    }

    Eigen::MatrixXf x = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    basis->control_points.resize(N);
    for (int32 i = 0; i < N; i++) {
        basis->control_points[i] = {points[i], vec3(x(i, 0), x(i, 1), x(i, 2))};
    }
}

RadialBasis compute_radial_basis(Array<const vec3> points, Array<const vec3> values, RadialBasisFunctionType type, float function_param) {
    RadialBasis basis;
    compute_radial_basis(&basis, points, values);
    return basis;
}

vec3 evaluate_radial_basis(const RadialBasis& basis, const vec3& point) {
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