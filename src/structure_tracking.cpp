#include "structure_tracking.h"

#include <core/log.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#include <svd3/svd3.h>

//#pragma warning(disable : 4127)  // disable warnings about expressions which could be constexpr in Eigen
//#include <Eigen/Eigen>

namespace structure_tracking {

struct Structure {
    ID id = 0;
    int32 num_frames = 0;

    struct {
        struct {
            vec3* com = nullptr;
            struct {
                quat* absolute = nullptr;
                quat* relative = nullptr;
                quat* corrected = nullptr;
            } rotation;
        } transform;

        struct {
            mat3* vector = nullptr;
            vec3* value = nullptr;
        } eigen;

        struct {
            float* abs = nullptr;
            float* rel = nullptr;
        } determinant;
    } frame_data;
};

struct Entry {
    ID id;
    Structure* ptr;
};

struct Context {
    uint32 next_hash = 0xdeadf00d;
    DynamicArray<Entry> entries{};
};

Context* context = nullptr;

// from here https://stackoverflow.com/questions/4372224/fast-method-for-computing-3x3-symmetric-matrix-spectral-decomposition
// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0
static void Diagonalize(const float (&A)[3][3], float (&Q)[3][3], float (&D)[3][3]) {
    // A must be a symmetric matrix.
    // returns Q and D such that
    // Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
    const int maxsteps = 24;  // certainly wont need that many.
    int k0, k1, k2;
    float o[3], m[3];
    float q[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    float jr[4];
    float sqw, sqx, sqy, sqz;
    float tmp1, tmp2, mq;
    float AQ[3][3];
    float thet, sgn, t, c;
    for (int i = 0; i < maxsteps; ++i) {
        // quat to matrix
        sqx = q[0] * q[0];
        sqy = q[1] * q[1];
        sqz = q[2] * q[2];
        sqw = q[3] * q[3];
        Q[0][0] = (sqx - sqy - sqz + sqw);
        Q[1][1] = (-sqx + sqy - sqz + sqw);
        Q[2][2] = (-sqx - sqy + sqz + sqw);
        tmp1 = q[0] * q[1];
        tmp2 = q[2] * q[3];
        Q[1][0] = 2.0f * (tmp1 + tmp2);
        Q[0][1] = 2.0f * (tmp1 - tmp2);
        tmp1 = q[0] * q[2];
        tmp2 = q[1] * q[3];
        Q[2][0] = 2.0f * (tmp1 - tmp2);
        Q[0][2] = 2.0f * (tmp1 + tmp2);
        tmp1 = q[1] * q[2];
        tmp2 = q[0] * q[3];
        Q[2][1] = 2.0f * (tmp1 + tmp2);
        Q[1][2] = 2.0f * (tmp1 - tmp2);

        // AQ = A * Q
        AQ[0][0] = Q[0][0] * A[0][0] + Q[1][0] * A[0][1] + Q[2][0] * A[0][2];
        AQ[0][1] = Q[0][1] * A[0][0] + Q[1][1] * A[0][1] + Q[2][1] * A[0][2];
        AQ[0][2] = Q[0][2] * A[0][0] + Q[1][2] * A[0][1] + Q[2][2] * A[0][2];
        AQ[1][0] = Q[0][0] * A[0][1] + Q[1][0] * A[1][1] + Q[2][0] * A[1][2];
        AQ[1][1] = Q[0][1] * A[0][1] + Q[1][1] * A[1][1] + Q[2][1] * A[1][2];
        AQ[1][2] = Q[0][2] * A[0][1] + Q[1][2] * A[1][1] + Q[2][2] * A[1][2];
        AQ[2][0] = Q[0][0] * A[0][2] + Q[1][0] * A[1][2] + Q[2][0] * A[2][2];
        AQ[2][1] = Q[0][1] * A[0][2] + Q[1][1] * A[1][2] + Q[2][1] * A[2][2];
        AQ[2][2] = Q[0][2] * A[0][2] + Q[1][2] * A[1][2] + Q[2][2] * A[2][2];
        // D = Qt * AQ
        D[0][0] = AQ[0][0] * Q[0][0] + AQ[1][0] * Q[1][0] + AQ[2][0] * Q[2][0];
        D[0][1] = AQ[0][0] * Q[0][1] + AQ[1][0] * Q[1][1] + AQ[2][0] * Q[2][1];
        D[0][2] = AQ[0][0] * Q[0][2] + AQ[1][0] * Q[1][2] + AQ[2][0] * Q[2][2];
        D[1][0] = AQ[0][1] * Q[0][0] + AQ[1][1] * Q[1][0] + AQ[2][1] * Q[2][0];
        D[1][1] = AQ[0][1] * Q[0][1] + AQ[1][1] * Q[1][1] + AQ[2][1] * Q[2][1];
        D[1][2] = AQ[0][1] * Q[0][2] + AQ[1][1] * Q[1][2] + AQ[2][1] * Q[2][2];
        D[2][0] = AQ[0][2] * Q[0][0] + AQ[1][2] * Q[1][0] + AQ[2][2] * Q[2][0];
        D[2][1] = AQ[0][2] * Q[0][1] + AQ[1][2] * Q[1][1] + AQ[2][2] * Q[2][1];
        D[2][2] = AQ[0][2] * Q[0][2] + AQ[1][2] * Q[1][2] + AQ[2][2] * Q[2][2];
        o[0] = D[1][2];
        o[1] = D[0][2];
        o[2] = D[0][1];
        m[0] = fabs(o[0]);
        m[1] = fabs(o[1]);
        m[2] = fabs(o[2]);

        k0 = (m[0] > m[1] && m[0] > m[2]) ? 0 : (m[1] > m[2]) ? 1 : 2;  // index of largest element of offdiag
        k1 = (k0 + 1) % 3;
        k2 = (k0 + 2) % 3;
        if (o[k0] == 0.0f) {
            break;  // diagonal already
        }
        thet = (D[k2][k2] - D[k1][k1]) / (2.0f * o[k0]);
        sgn = (thet > 0.0f) ? 1.0f : -1.0f;
        thet *= sgn;                                                             // make it positive
        t = sgn / (thet + ((thet < 1.E6f) ? sqrtf(thet * thet + 1.0f) : thet));  // sign(T)/(|T|+sqrt(T^2+1))
        c = 1.0f / sqrtf(t * t + 1.0f);                                          //  c= 1/(t^2+1) , t=s/c
        if (c == 1.0f) {
            break;  // no room for improvement - reached machine precision.
        }
        jr[0] = jr[1] = jr[2] = jr[3] = 0.0f;
        jr[k0] = sgn * sqrtf((1.0f - c) / 2.0f);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)
        jr[k0] *= -1.0f;                          // since our quat-to-matrix convention was for v*M instead of M*v
        jr[3] = sqrtf(1.0f - jr[k0] * jr[k0]);
        if (jr[3] == 1.0f) {
            break;  // reached limits of floating point precision
        }
        q[0] = (q[3] * jr[0] + q[0] * jr[3] + q[1] * jr[2] - q[2] * jr[1]);
        q[1] = (q[3] * jr[1] - q[0] * jr[2] + q[1] * jr[3] + q[2] * jr[0]);
        q[2] = (q[3] * jr[2] + q[0] * jr[1] - q[1] * jr[0] + q[2] * jr[3]);
        q[3] = (q[3] * jr[3] - q[0] * jr[0] - q[1] * jr[1] - q[2] * jr[2]);
        mq = sqrtf(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
        q[0] /= mq;
        q[1] /= mq;
        q[2] /= mq;
        q[3] /= mq;
    }
}

static void diagonalize(const mat3& M, mat3* Q, mat3* D) {
    ASSERT(Q);
    ASSERT(D);
    Diagonalize((const float(&)[3][3])M, (float(&)[3][3]) * Q, (float(&)[3][3]) * D);
}

static void decompose(const mat3& M, mat3* R, mat3* S) {
    ASSERT(R);
    ASSERT(S);
    mat3 AtA = math::transpose(M) * M;
    mat3 Q, D;
    diagonalize(AtA, &Q, &D);
    // const float det = math::determinant(AtA); // For debugging
    D = mat3(sqrtf(D[0][0]), 0, 0, 0, sqrtf(D[1][1]), 0, 0, 0, sqrtf(D[2][2]));
    /*
        D[0][0] = sqrtf(D[0][0]);
    D[1][1] = sqrtf(D[1][1]);
    D[2][2] = sqrtf(D[2][2]);
        */
    // @NOTE: Should one zero every non diagonal element???
    *S = math::transpose(Q) * D * Q;
    *R = M * math::inverse(*S);
}

/*
// clang-format off
static mat3 compute_mass_weighted_cross_covariance_matrix(const float* x0, const float* y0, const float* z0, const float* mass0,
                                                                                                              const float* x1, const float* y1, const float* z1, const float* mass1,
                                                                                                                  int64 count, const vec3& com0, const vec3& com1)
// clang-format on
{
    mat3 A{0};
    float mass_sum = 0.0f;
    for (int64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float qx = x0[i] - com0.x;
        const float qy = y0[i] - com0.y;
        const float qz = z0[i] - com0.z;

        const float px = x1[i] - com1.x;
        const float py = y1[i] - com1.y;
        const float pz = z1[i] - com1.z;

        const float m = mass0[i] * mass1[i];
        mass_sum += m;

        A[0][0] += m * px * qx;
        A[0][1] += m * py * qx;
        A[0][2] += m * pz * qx;
        A[1][0] += m * px * qy;
        A[1][1] += m * py * qy;
        A[1][2] += m * pz * qy;
        A[2][0] += m * px * qz;
        A[2][1] += m * py * qz;
        A[2][2] += m * pz * qz;
    }

    return A;
}
*/

// clang-format off
static mat3 compute_mass_weighted_cross_covariance_matrix(const float* x0, const float* y0, const float* z0,
                                                          const float* x1, const float* y1, const float* z1,
                                                          const float* mass,
                                                          int64 count, const vec3& com0, const vec3& com1)
// clang-format on
{
    mat3 A{0};
    float mass_sum = 0.0f;
    for (int64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float qx = x0[i] - com0.x;
        const float qy = y0[i] - com0.y;
        const float qz = z0[i] - com0.z;

        const float px = x1[i] - com1.x;
        const float py = y1[i] - com1.y;
        const float pz = z1[i] - com1.z;

        const float m = mass[i];
        mass_sum += m;

        A[0][0] += m * px * qx;
        A[0][1] += m * py * qx;
        A[0][2] += m * pz * qx;
        A[1][0] += m * px * qy;
        A[1][1] += m * py * qy;
        A[1][2] += m * pz * qy;
        A[2][0] += m * px * qz;
        A[2][1] += m * py * qz;
        A[2][2] += m * pz * qz;
    }

    return A;
}

// clang-format off
static mat3 compute_mass_weighted_cross_covariance_matrix(const float* x0, const float* y0, const float* z0,
                                                          const float* x1, const float* y1, const float* z1,
                                                          const float* mass, int64 count)
// clang-format on
{
    mat3 A{0};
    float mass_sum = 0.0f;
    for (int64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float qx = x0[i];
        const float qy = y0[i];
        const float qz = z0[i];

        const float px = x1[i];
        const float py = y1[i];
        const float pz = z1[i];

        const float m = mass[i];
        mass_sum += m;

        A[0][0] += m * px * qx;
        A[0][1] += m * py * qx;
        A[0][2] += m * pz * qx;
        A[1][0] += m * px * qy;
        A[1][1] += m * py * qy;
        A[1][2] += m * pz * qy;
        A[2][0] += m * px * qz;
        A[2][1] += m * py * qz;
        A[2][2] += m * pz * qz;
    }

    return A;
}

// clang-format off
static mat3 compute_mass_weighted_covariance_matrix(const float* x, const float* y, const float* z,
                                                    const float* mass, int64 count, const vec3& com)
// clang-format on
{
    mat3 A{0};
    float mass_sum = 0.0f;
    for (int64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float qx = x[i] - com.x;
        const float qy = y[i] - com.y;
        const float qz = z[i] - com.z;
        const float m = mass[i];
        mass_sum += m;

        A[0][0] += m * qx * qx;
        A[0][1] += m * qy * qx;
        A[0][2] += m * qz * qx;
        A[1][0] += m * qx * qy;
        A[1][1] += m * qy * qy;
        A[1][2] += m * qz * qy;
        A[2][0] += m * qx * qz;
        A[2][1] += m * qy * qz;
        A[2][2] += m * qz * qz;
    }

    return A;
}

// clang-format off
static mat3 compute_mass_weighted_covariance_matrix(const float* x, const float* y, const float* z,
                                                    const float* mass, int64 count)
// clang-format on
{
    mat3 A{0};
    float mass_sum = 0.0f;
    for (int64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float qx = x[i];
        const float qy = y[i];
        const float qz = z[i];
        const float m = mass[i];
        mass_sum += m;

        A[0][0] += m * qx * qx;
        A[0][1] += m * qy * qx;
        A[0][2] += m * qz * qx;
        A[1][0] += m * qx * qy;
        A[1][1] += m * qy * qy;
        A[1][2] += m * qz * qy;
        A[2][0] += m * qx * qz;
        A[2][1] += m * qy * qz;
        A[2][2] += m * qz * qz;
    }

    return A;
}

#define ARGS(M) M[0][0], M[1][0], M[2][0], M[0][1], M[1][1], M[2][1], M[0][2], M[1][2], M[2][2]

static void compute_eigen(const mat3& M, vec3 (&vectors)[3], float (&value)[3]) {
    /*
    Eigen::Matrix3f A = Eigen::Matrix3f({
                                    {M[0][0], M[1][0], M[2][0]},
                                    {M[0][1], M[1][1], M[2][1]},
                                    {M[0][2], M[1][2], M[2][2]} });

    Eigen::EigenSolver<Eigen::Matrix3f> es(A);

    const auto swap = [](int& x, int& y) {
            int tmp = x;
            x = y;
            y = tmp;
    };

    const auto& e_vec = es.eigenvectors().real();
    const auto& e_val = es.eigenvalues().real();

    int l0 = 0, l1 = 1, l2 = 2;
    if (e_val[l0] < e_val[l1]) swap(l0, l1);
    if (e_val[l1] < e_val[l2]) swap(l1, l2);
    if (e_val[l0] < e_val[l1]) swap(l0, l1);

    vectors[0] = { e_vec(0, l0), e_vec(1, l0), e_vec(2, l0) };
    vectors[1] = { e_vec(0, l1), e_vec(1, l1), e_vec(2, l1) };
    vectors[2] = { e_vec(0, l2), e_vec(1, l2), e_vec(2, l2) };
    values[0]  = e_val[l0] / e_val[l0];
    values[1]  = e_val[l1] / e_val[l0];
    values[2]  = e_val[l2] / e_val[l0];
    */

    mat3 U, S, V;
    svd(ARGS(M), ARGS(U), ARGS(S), ARGS(V));
    float max_val = math::max(S[0][0], math::max(S[1][1], S[2][2]));
    const mat3 Ut = glm::transpose(U);

    const float e_val[] = {S[0][0] / max_val, S[1][1] / max_val, S[2][2] / max_val};
    const vec3 e_vec[] = {Ut[0], Ut[1], Ut[2]};
    int l[3] = {0, 1, 2};

    const auto swap = [](int& x, int& y) {
        int tmp = x;
        x = y;
        y = tmp;
    };

    if (e_val[l[0]] < e_val[l[1]]) swap(l[0], l[1]);
    if (e_val[l[1]] < e_val[l[2]]) swap(l[1], l[2]);
    if (e_val[l[0]] < e_val[l[1]]) swap(l[0], l[1]);

    value[0] = e_val[l[0]];
    value[1] = e_val[l[1]];
    value[2] = e_val[l[2]];

    vectors[0] = e_vec[l[0]];
    vectors[1] = e_vec[l[1]];
    vectors[2] = e_vec[l[2]];
}

static void compute_svd(const mat3& A, mat3& U, mat3& S, mat3& V) { svd(ARGS(A), ARGS(U), ARGS(S), ARGS(V)); }

#undef ARGS

vec3 compute_eigen_values(const float* RESTRICT x, const float* RESTRICT y, const float* RESTRICT z, const float* RESTRICT mass, int64 count) {
    const vec3 com = compute_com(x, y, z, mass, count);
    const mat3 A = compute_mass_weighted_covariance_matrix(x, y, z, mass, count, com);
    vec3 vecs[3];
    vec3 vals;
    compute_eigen(A, vecs, (float(&)[3])vals);
    return vals;
}

static mat3 extract_rotation(const mat3& M) {
    // mat3 R, S;
    // decompose(M, &R, &S);
    // return R;

    mat3 U, S, V;
    compute_svd(M, U, S, V);

    const mat3 Ut = math::transpose(U);
    const float d = math::determinant(V * Ut);
    const mat3 D = {1, 0, 0, 0, 1, 0, 0, 0, d};
    const mat3 R = V * D * Ut;
    return R;
}

// clang-format off
void compute_direction_vectors(float* out_x, float* out_y, float* out_z,
                               const float* in_x, const float* in_y, const float* in_z,
                               int64 count, const vec3& com)
// clang-format on
{
    for (int64 i = 0; i < count; i++) {
        const vec3 v = vec3(in_x[i], in_y[i], in_z[i]) - com;
        const vec3 dir = math::normalize(v);
        out_x[i] = dir.x;
        out_y[i] = dir.y;
        out_z[i] = dir.z;
    }
}

// clang-format off
mat3 compute_rotation(const float* x0, const float* y0, const float* z0,
                      const float* x1, const float* y1, const float* z1,
                      const float* mass, int64 count, const vec3& com0, const vec3& com1)
// clang-format on
{
#ifdef USE_DIRECTION_VECTORS
    void* tmp_mem = TMP_MALLOC(count * sizeof(float) * 6);
    float* dx0 = (float*)tmp_mem;
    float* dy0 = (float*)dx0 + count;
    float* dz0 = (float*)dy0 + count;
    float* dx1 = (float*)dz0 + count;
    float* dy1 = (float*)dx1 + count;
    float* dz1 = (float*)dy1 + count;

    compute_direction_vectors(dx0, dy0, dz0, x0, y0, z0, count, com0);
    compute_direction_vectors(dx1, dy1, dz1, x1, y1, z1, count, com1);

    const mat3 Apq = compute_mass_weighted_cross_covariance_matrix(dx0, dy0, dz0, dx1, dy1, dz1, mass, count);
    const mat3 Aqq = compute_mass_weighted_covariance_matrix(dx0, dy0, dz0, mass, count);
#else
    const mat3 Apq = compute_mass_weighted_cross_covariance_matrix(x0, y0, z0, x1, y1, z1, mass, count, com0, com1);
    const mat3 Aqq = compute_mass_weighted_covariance_matrix(x0, y0, z0, mass, count);
#endif

    const mat3 A = Apq;
    // return A; // Return complete linear transform with skewing and all
    return extract_rotation(A);  // Return rotational part
}

// clang-format off
static void compute_residual_error(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
                                   const float* RESTRICT src_x, const float* RESTRICT src_y, const float* RESTRICT src_z,
                                   const float* RESTRICT ref_x, const float* RESTRICT ref_y, const float* RESTRICT ref_z,
                                   int64 count, const mat4& matrix)
// clang-format on
{
    for (int32 i = 0; i < count; i++) {
        // @TODO: Vectorize this...

        const vec4 r = {ref_x[i], ref_y[i], ref_z[i], 1.0f};
        const vec4 v = {src_x[i], src_y[i], src_z[i], 1.0f};
        const vec4 u = matrix * v;
        const vec4 d = u - r;

        out_x[i] = d.x;
        out_y[i] = d.y;
        out_z[i] = d.z;
    }
}

#if 0
// RBF functions
inline float Wendland_3_1(float r) {
    const float x = 1.f - r;
    const float x2 = x * x;
    return x2 * x2 * (4.f * r + 1.f);
}

inline float Gaussian(float r) {
    const float a = 0.5f * r;
    return math::exp(-a*a);
}

#define RBF_FUNCTION Wendland_3_1

static void compute_rbf_weights(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
                         const float* RESTRICT in_pos_x, const float* RESTRICT in_pos_y, const float* RESTRICT in_pos_z,
                         const float* RESTRICT in_val_x, const float* RESTRICT in_val_y, const float* RESTRICT in_val_z,
                         int64 count, const float radial_cutoff = 10.f) {

    const float d2_max = radial_cutoff * radial_cutoff;
    const float scl = 1.f / radial_cutoff;
    const int32 N = (int32)count;

    // Create matrix A NxN
    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(N, N);
    for (int32 i = 0; i < N; i++) {
        for (int32 j = 0; j < N; j++) {
            const float dx = in_pos_x[i] - in_pos_x[j];
            const float dy = in_pos_y[i] - in_pos_y[j];
            const float dz = in_pos_z[i] - in_pos_z[j];
            const float d2 = dx * dx + dy * dy + dz * dz;
            if (d2 > d2_max) continue;

            const float d = math::sqrt(d2);
            A(i, j) = RBF_FUNCTION(d * scl);
        }
    }

    // Create Vector b
    Eigen::MatrixXf b = Eigen::MatrixXf::Zero(N, 3);
    for (int32 i = 0; i < N; i++) {
        b.row(i) = Eigen::Vector3f(in_val_x[i], in_val_y[i], in_val_z[i]);
    }

    Eigen::MatrixXf x = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    for (int32 i = 0; i < N; i++) {
        out_x[i] = x(i, 0);
        out_y[i] = x(i, 1);
out_z[i] = x(i, 2);
    };
}
#endif

static void free_structure_data(Structure* s) {
    ASSERT(s);
    if (s->frame_data.transform.rotation.absolute) FREE(s->frame_data.transform.rotation.absolute);
    if (s->frame_data.transform.rotation.relative) FREE(s->frame_data.transform.rotation.relative);
    if (s->frame_data.transform.rotation.corrected) FREE(s->frame_data.transform.rotation.corrected);
    if (s->frame_data.transform.com) FREE(s->frame_data.transform.com);
    if (s->frame_data.eigen.vector) FREE(s->frame_data.eigen.vector);
    if (s->frame_data.eigen.value) FREE(s->frame_data.eigen.value);
    if (s->frame_data.determinant.abs) FREE(s->frame_data.determinant.abs);
    if (s->frame_data.determinant.rel) FREE(s->frame_data.determinant.rel);
}

static void init_structure_data(Structure* s, ID id, int32 num_frames) {
    ASSERT(s);
    free_structure_data(s);
    s->id = id;
    s->num_frames = num_frames;

    s->frame_data.transform.rotation.absolute = (quat*)MALLOC(sizeof(quat) * num_frames);
    s->frame_data.transform.rotation.relative = (quat*)MALLOC(sizeof(quat) * num_frames);
    s->frame_data.transform.rotation.corrected = (quat*)MALLOC(sizeof(quat) * num_frames);
    s->frame_data.transform.com = (vec3*)MALLOC(sizeof(vec3) * num_frames);
    s->frame_data.eigen.vector = (mat3*)MALLOC(sizeof(mat3) * num_frames);
    s->frame_data.eigen.value = (vec3*)MALLOC(sizeof(vec3) * num_frames);
    s->frame_data.determinant.abs = (float*)MALLOC(sizeof(float) * num_frames);
    s->frame_data.determinant.rel = (float*)MALLOC(sizeof(float) * num_frames);
}

static Structure* find_structure(ID id) {
    for (const auto& e : context->entries) {
        if (e.id == id) return e.ptr;
    }
    return nullptr;
}

void initialize() {
    if (!context) {
        context = NEW(Context);
        context->entries.reserve(8);
    }
}

void shutdown() {
    if (context) {
        clear_structures();
        DELETE(context);
    }
}

ID create_structure() {
    ASSERT(context);
    ID id = context->next_hash;
    context->next_hash = hash::crc32(context->next_hash);

    // This is to ensure that no collision happens
    ASSERT(find_structure(id) == nullptr);

    Structure* s = new (MALLOC(sizeof(Structure))) Structure();
    context->entries.push_back({id, s});
    return id;
}

bool remove_structure(ID id) {
    ASSERT(context);
    for (auto& e : context->entries) {
        if (e.id == id) {
            if (e.ptr) free_structure_data(e.ptr);
            FREE(e.ptr);
            context->entries.swap_back_and_pop(&e);
            return true;
        }
    }
    LOG_ERROR("Could not remove structure: it does not exists!");
    return false;
}

void clear_structures() {
    ASSERT(context);
    for (auto& e : context->entries) {
        if (e.ptr) free_structure_data(e.ptr);
    }
    context->entries.clear();
}
#if 0
inline vec3 compute_shape_space_weights(const vec3& eigen_values) {
    const vec3& l = eigen_values;
    const float l_sum = l[0] + l[1] + l[2];
    if (l_sum < 1.0e-6f) return {0, 0, 0};

    const float one_over_denom = 1.0f / l_sum;
    const float c_lin = (l[0] - l[1]) * one_over_denom;
    const float c_pla = 2.0f * (l[1] - l[2]) * one_over_denom;
    const float c_iso = 3.0f * l[2] * one_over_denom;

    return {c_lin, c_pla, c_iso};
}

inline vec3 compute_weights(const vec3& eigen_values) {
    const vec3 w = compute_shape_space_weights(eigen_values);

    const vec3 w_lin = w.x * vec3(0.0f, 0.5f, 0.5f);  // Linear case, emphasize weight to the mid and min eigen directions
    const vec3 w_pla = w.y * vec3(0.0f, 0.0f, 1.0f);  // Planar case, emphasize weight to min eigen directions
    const vec3 w_iso = w.z * vec3(0.0f, 0.0f, 0.0f);  // Isotropic case, no emphasis on any eigen direction

    return w_lin + w_pla + w_iso;
}

/*
void compute_support_frame(SupportFrame& frame, const mat3& eigen_vectors, const vec3& eigen_values, const vec3& com, float tot_mass, const SupportFrame* ref = nullptr) {
    vec3 dir[6];
    float mass[6];
    const vec3 axis_weights = compute_weights(eigen_values) * tot_mass;

    for (int i = 0; i < 3; i++) {
        dir[i * 2 + 0] = eigen_vectors[i];
        dir[i * 2 + 1] = -eigen_vectors[i];
        mass[i * 2 + 0] = axis_weights[i];
        mass[i * 2 + 1] = axis_weights[i];
    }




if (ref) {
    DynamicArray<int> avail_idx = {0, 1, 2, 3, 4, 5};
    for (int i = 0; i < 6; i++) {
        // Find best matching axis in reference support frame
        float max_d = -FLT_MAX;
        int* max_it = 0;
        for (auto& j : avail_idx) {
            const float d = math::dot(dir[i], ref->axis[j].dir);
            if (d > max_d) {
                max_d = d;
                max_it = &j;
            }
                    }
        ASSERT(max_it != 0);
        const auto max_idx = *max_it;
        frame.axis[i].dir = dir[max_idx];
        frame.axis[i].mass = mass[max_idx];
        avail_idx.swap_back_and_pop(max_it);
    }
} else



    {
        for (int i = 0; i < 6; i++) {
            frame.axis[i].dir = dir[i];
            frame.axis[i].mass = mass[i];
        }
    }
}

inline void set_stabilization_point_pos(float* RESTRICT x, float* RESTRICT y, float* RESTRICT z, const SupportFrame& frame, const vec3& com) {
    for (int i = 0; i < 6; i++) {
        const vec3 p = com + frame.axis[i].dir;
        x[i] = p.x;
        y[i] = p.y;
        z[i] = p.z;
    }
}

inline void set_stabilization_point_mass(float* mass, const SupportFrame& frame) {
    for (int i = 0; i < 6; i++) {
        mass[i] = frame.axis[i].mass;
    }
}
*/

#endif

static inline mat3 orthogonalize(const mat3& M) {
    mat3 U, S, V;
    compute_svd(M, U, S, V);
    return U * math::transpose(V);
}

void smooth_rotation_matrices(mat3* RESTRICT R_out, const mat3* RESTRICT R_in, int N) {
    constexpr int window = 2;
    // constexpr float w[3] = {0.7, 0.1, 0.05};
    constexpr float w[3] = {0.2, 0.2, 0.2};

    ASSERT(N >= 4);

    R_out[0] = orthogonalize((w[0] + w[1]) * R_in[0] + w[1] * R_in[1] + w[0] * R_in[2]);
    R_out[1] = orthogonalize((w[1] + w[2]) * R_in[0] + w[0] * R_in[1] + w[1] * R_in[2] + w[2] * R_in[3]);

    for (int i = window; i < N - window; i++) {
        R_out[i] = orthogonalize(w[2] * R_in[i - 2] + w[1] * R_in[i - 1] + w[0] * R_in[i] + w[1] * R_in[i + 1] + w[2] * R_in[i + 2]);
    }

    R_out[N - 2] = orthogonalize(w[2] * R_in[N - 4] + w[1] * R_in[N - 3] + w[0] * R_in[N - 2] + (w[1] + w[2]) * R_in[N - 1]);
    R_out[N - 1] = orthogonalize(w[2] * R_in[N - 3] + w[1] * R_in[N - 2] + (w[0] + w[1]) * R_in[N - 1]);
}

bool compute_trajectory_transform_data(ID id, Bitfield atom_mask, const MoleculeDynamic& dynamic) {
    ASSERT(context);

    const int32 num_frames = (int32)dynamic.trajectory.num_frames;

    Structure* s = find_structure(id);
    if (s == nullptr) {
        LOG_ERROR("Could not compute tracking data, supplied id is not valid.");
        return false;
    }

    const int num_points = (int)bitfield::number_of_bits_set(atom_mask);
    if (num_points == 0) {
        LOG_ERROR("Supplied atom mask is empty.");
        return false;
    }

    // Allocate memory for all data
    init_structure_data(s, id, num_frames);

    // Scratch data
    const auto mem_size = sizeof(float) * num_points * 13;
    void* mem = TMP_MALLOC(mem_size);
    defer { TMP_FREE(mem); };
    memset(mem, 0, mem_size);

    float* cur_x = (float*)mem;
    float* cur_y = cur_x + num_points;
    float* cur_z = cur_y + num_points;

    float* ref_x = cur_z + num_points;
    float* ref_y = ref_x + num_points;
    float* ref_z = ref_y + num_points;

    float* prv_x = ref_z + num_points;
    float* prv_y = prv_x + num_points;
    float* prv_z = prv_y + num_points;

    float* cor_x = prv_z + num_points;
    float* cor_y = cor_x + num_points;
    float* cor_z = cor_y + num_points;

    float* mass = cor_z + num_points;

    bitfield::extract_data_from_mask(cur_x, get_trajectory_position_x(dynamic.trajectory, 0).data(), atom_mask);
    bitfield::extract_data_from_mask(cur_y, get_trajectory_position_y(dynamic.trajectory, 0).data(), atom_mask);
    bitfield::extract_data_from_mask(cur_z, get_trajectory_position_z(dynamic.trajectory, 0).data(), atom_mask);
    bitfield::extract_data_from_mask(ref_x, get_trajectory_position_x(dynamic.trajectory, 0).data(), atom_mask);
    bitfield::extract_data_from_mask(ref_y, get_trajectory_position_y(dynamic.trajectory, 0).data(), atom_mask);
    bitfield::extract_data_from_mask(ref_z, get_trajectory_position_z(dynamic.trajectory, 0).data(), atom_mask);

    bitfield::extract_data_from_mask(mass, dynamic.molecule.atom.mass, atom_mask);

    /*
    float tot_mass = 0.0f;
    for (int i = 0; i < num_points; i++) {
        tot_mass += mass[i];
    }
    */

    const vec3 ref_com = compute_com(ref_x, ref_y, ref_z, mass, num_points);
    vec3 prv_com = ref_com;
    vec3 cor_com = ref_com;
    vec3 cur_com = {0, 0, 0};

    quat q_relative = {1, 0, 0, 0};
    quat q_corrected = {1, 0, 0, 0};

    const mat3 ref_cov_mat = compute_mass_weighted_covariance_matrix(ref_x, ref_y, ref_z, mass, num_points, ref_com);
    mat3 ref_eigen_vectors;
    vec3 ref_eigen_values;
    compute_eigen(ref_cov_mat, (vec3(&)[3])ref_eigen_vectors, (float(&)[3])ref_eigen_values);

    memcpy(prv_x, ref_x, num_points * sizeof(float) * 3);

    // Set first frame explicitly
    s->frame_data.transform.rotation.absolute[0] = {1, 0, 0, 0};
    s->frame_data.transform.rotation.relative[0] = {1, 0, 0, 0};
    s->frame_data.transform.rotation.corrected[0] = {1, 0, 0, 0};
    s->frame_data.transform.com[0] = ref_com;
    s->frame_data.eigen.vector[0] = ref_eigen_vectors;
    s->frame_data.eigen.value[0] = ref_eigen_values;
    s->frame_data.determinant.abs[0] = 1.0f;
    s->frame_data.determinant.rel[0] = 1.0f;

    for (int32 i = 1; i < num_frames; i++) {
        // Copy previous frame data
        memcpy(prv_x, cur_x, num_points * sizeof(float) * 3);
        prv_com = cur_com;

        // Fetch current
        bitfield::extract_data_from_mask(cur_x, get_trajectory_position_x(dynamic.trajectory, i).data(), atom_mask);
        bitfield::extract_data_from_mask(cur_y, get_trajectory_position_y(dynamic.trajectory, i).data(), atom_mask);
        bitfield::extract_data_from_mask(cur_z, get_trajectory_position_z(dynamic.trajectory, i).data(), atom_mask);

        cur_com = compute_com(cur_x, cur_y, cur_z, mass, num_points);

        const mat3 abs_mat = compute_mass_weighted_cross_covariance_matrix(cur_x, cur_y, cur_z, ref_x, ref_y, ref_z, mass, num_points, cur_com, ref_com);
        const mat3 rel_mat = compute_mass_weighted_cross_covariance_matrix(prv_x, prv_y, prv_z, cur_x, cur_y, cur_z, mass, num_points, prv_com, cur_com);
        const mat3 cov_mat = compute_mass_weighted_covariance_matrix(cur_x, cur_y, cur_z, mass, num_points, cur_com);

        const mat3 abs_rot = extract_rotation(abs_mat);
        const mat3 rel_rot = extract_rotation(rel_mat);

        const quat q_del = math::quat_cast(rel_rot);
        const quat q_abs = math::quat_cast(abs_rot);

        q_relative = math::normalize(q_relative * q_del);

        const float abs_det = math::determinant(abs_mat / cov_mat);
        const float rel_det = math::determinant(rel_mat / cov_mat);
        const quat q_rel = math::conjugate(q_relative);
#if 0
        memcpy(cor_x, ref_x, num_points * sizeof(float) * 3);
        cor_com = ref_com;

        {
            const quat old_q = q_corrected;
            const quat new_q = q_corrected * q_del;
            q_corrected = math::dot(old_q, new_q) < 0.0f ? -new_q : new_q;
        }

        {
            const mat4 R = math::mat4_cast(q_corrected);
            const mat4 T1 = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {-cor_com, 1}};
            const mat4 T2 = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {cor_com, 1}};
            const mat4 M = T2 * R * T1;
            transform_positions_ref(cor_x, cor_y, cor_z, num_points, M);
        }

        const mat3 err_mat = compute_mass_weighted_cross_covariance_matrix(cor_x, cor_y, cor_z, cur_x, cur_y, cur_z, mass, num_points, cor_com, cur_com);
        const mat3 err_rot = extract_rotation(err_mat);
        const quat q_err = math::inverse(math::normalize(math::quat_cast(err_rot)));

        printf("err angle: %3.3f\n", math::rad_to_deg(math::angle(q_err)));

        q_corrected = math::normalize(q_err * q_corrected);
        const quat q_cor = q_corrected;
#else
        q_corrected = math::nlerp(q_corrected * q_del, math::conjugate(q_abs), 0.5f);
        const quat q_cor = math::conjugate(q_corrected);
#endif

        s->frame_data.transform.rotation.absolute[i] = math::dot(s->frame_data.transform.rotation.absolute[i - 1], q_abs) > 0.0f ? q_abs : -q_abs;
        s->frame_data.transform.rotation.relative[i] = math::dot(s->frame_data.transform.rotation.relative[i - 1], q_rel) > 0.0f ? q_rel : -q_rel;
        s->frame_data.transform.rotation.corrected[i] = math::dot(s->frame_data.transform.rotation.corrected[i - 1], q_cor) > 0.0f ? q_cor : -q_cor;
        compute_eigen(cov_mat, (vec3(&)[3])s->frame_data.eigen.vector[i], (float(&)[3])s->frame_data.eigen.value[i]);
        s->frame_data.transform.com[i] = cur_com;
        s->frame_data.determinant.abs[i] = abs_det;
        s->frame_data.determinant.rel[i] = rel_det;
    }

    return true;
}

const Array<const vec3> get_com(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return {s->frame_data.transform.com, s->num_frames};
    }
    LOG_ERROR("Supplied id is not valid.");
    return {};
}

const Array<const quat> get_rot_absolute(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return {s->frame_data.transform.rotation.absolute, s->num_frames};
    }
    LOG_ERROR("Supplied id is not valid.");
    return {};
}

const Array<const quat> get_rot_relative(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return {s->frame_data.transform.rotation.relative, s->num_frames};
    }
    LOG_ERROR("Supplied id is not valid.");
    return {};
}

const Array<const quat> get_rot_corrected(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return {s->frame_data.transform.rotation.corrected, s->num_frames};
    }
    LOG_ERROR("Supplied id is not valid.");
    return {};
}

const Array<const mat3> get_eigen_vectors(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return {s->frame_data.eigen.vector, s->num_frames};
    }
    LOG_ERROR("Supplied id is not valid.");
    return {};
}

const Array<const vec3> get_eigen_values(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return {s->frame_data.eigen.value, s->num_frames};
    }
    LOG_ERROR("Supplied id is not valid.");
    return {};
}

const Array<const float> get_abs_det(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return {s->frame_data.determinant.abs, s->num_frames};
    }
    LOG_ERROR("Supplied id is not valid.");
    return {};
}

const Array<const float> get_rel_det(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return {s->frame_data.determinant.rel, s->num_frames};
    }
    LOG_ERROR("Supplied id is not valid.");
    return {};
}

}  // namespace structure_tracking