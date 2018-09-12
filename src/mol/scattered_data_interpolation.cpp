#include "scattered_data_interpolation.h"
#include <core/math_utils.h>
#include <iostream>
#include <Eigen/Eigen>

// RBF functions
inline float Wendland_3_1(float r, float radius = 1.f) {
    r = r / radius;
    if (0.f <= r && r <= 1.f) {
        const float x = 1.f - r;
        const float x2 = x * x;
        return x2 * x2 * (4.f * r + 1.f);
    }
    return 0.f;
}

inline float Gaussian(float r) { return math::exp(-math::pow(0.5f * r, 2)); }

void compute_displacement_basis(RadialBasis* basis, Array<const vec3> orig_points, Array<const vec3> disp_points) {
    ASSERT(basis);
    ASSERT(orig_points.count == disp_points.count);
    const int N = math::min(orig_points.count, 2048LL);

    DynamicArray<int32> indices(N);
    /*
srand(0);
for (int32 i = 0; i < N; i++) {
    indices[i] = rand() % orig_points.count;
}
    */

    for (int32 i = 0; i < N; i++) {
        indices[i] = i;
    }

    DynamicArray<vec3> disp(N);
    for (int32 i = 0; i < N; i++) {
        int32 idx = indices[i];
        disp[i] = orig_points[idx] - disp_points[idx];
    }

    // Create matrix A NxN
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N, N);
    for (int32 i = 0; i < N; i++) {
        for (int32 j = 0; j < N; j++) {
            A(i, j) = Wendland_3_1(math::length(disp_points[indices[i]] - disp_points[indices[j]]), 10.f);
            // A(i, j) = Gaussian(math::length(disp_points[indices[i]] - disp_points[indices[j]]));
        }
    }

    // std::cout << "This is A:\n" << A << "\n";
    /*
// Create matrix P     Nx4
Eigen::MatrixXf P = Eigen::MatrixXf::Zero(N, 4);
for (int32 i = 0; i < N; i++) {
    P(i, 0) = 1.f;
    P(i, 1) = disp_points[i].x;
    P(i, 2) = disp_points[i].y;
    P(i, 3) = disp_points[i].z;
}

Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N + 4, N + 4);
A.block(0, 0, N, N) = Theta;
A.block(0, N, N, 4) = P;
A.block(N, 0, 4, N) = P.transpose();
    */

    // Create Vector b
    Eigen::MatrixXf b = Eigen::MatrixXf::Zero(N, 3);
    for (int32 i = 0; i < N; i++) {
        b.row(i) = Eigen::Vector3f(disp[i].x, disp[i].y, disp[i].z);
    }

    /*
std::cout << "This is b:\n" << b << "\n";

std::cout << "This is x:\n" << x << "\n";
    */
    Eigen::MatrixXf x = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    basis->control_points.resize(N);
    basis->weights.resize(N);

    for (int32 i = 0; i < N; i++) {
        basis->control_points[i] = disp_points[indices[i]];
        basis->weights[i] = vec3(x(i, 0), x(i, 1), x(i, 2));
    }

    /*
int32 idx = 8;
vec3 exact = disp[idx];
vec3 approx = evaluate_displacement(*basis, disp_points[indices[idx]]);

std::cout << "Testing\n";

printf("Exact:  %f %f %f\n", exact.x, exact.y, exact.z);

printf("Approx: %f %f %f\n", approx.x, approx.y, approx.z);
    */
}

RadialBasis compute_displacement_basis(Array<const vec3> orig_points, Array<const vec3> disp_points) {
    RadialBasis basis;
    compute_displacement_basis(&basis, orig_points, disp_points);
    return basis;
}

vec3 evaluate_displacement(const RadialBasis& basis, const vec3& x) {
    ASSERT(basis.control_points.count == basis.weights.count);
    const int N = basis.control_points.count;

    vec3 result{0.f};
    for (int32 i = 0; i < N; i++) {
        result += basis.weights[i] * Wendland_3_1(math::length(x - basis.control_points[i]), 10.f);
        // result += basis.weights[i] * Gaussian(math::length(x - basis.control_points[i]));
    }

    return result;
}
