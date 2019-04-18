#include "structure_tracking.h"

#include <core/log.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#pragma warning(disable : 4127)  // disable warnings about expressions which could be constexpr in Eigen
#include <Eigen/Eigen>

namespace structure_tracking {

struct SupportFrame {
    struct {
        vec3 pos;
        float mass;
    } point[6];
};

struct Structure {
    ID id = 0;
    int32 ref_frame_idx = 0;
    int32 num_points = 0;
    int32 num_frames = 0;

    struct {
        Transform* transform = nullptr;
        struct {
            mat3* vector = nullptr;
            vec3* value = nullptr;
        } eigen;
        struct {
            float* abs = nullptr;
            float* rel = nullptr;
        } determinant;
        SupportFrame* support_frames = nullptr;
    } frame_data;
};

struct Entry {
    ID id;
    Structure* ptr;
};

struct Context {
    uint32 next_hash = 0xdeadb00b;
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
    D[0][0] = sqrtf(D[0][0]);
    D[1][1] = sqrtf(D[1][1]);
    D[2][2] = sqrtf(D[2][2]);
    // @NOTE: Should one zero every non diagonal element???
    *S = math::inverse(Q) * D * Q;
    *R = M * math::inverse(*S);
}

static mat3 compute_covariance_matrix(const float* x0, const float* y0, const float* z0, const float* x1, const float* y1, const float* z1, const float* mass, int64 count, const vec3& com0,
                                      const vec3& com1) {
    mat3 A{0};
    for (int64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float q_x = x0[i] - com0.x;
        const float q_y = y0[i] - com0.y;
        const float q_z = z0[i] - com0.z;

        const float p_x = x1[i] - com1.x;
        const float p_y = y1[i] - com1.y;
        const float p_z = z1[i] - com1.z;

        A[0][0] += mass[i] * p_x * q_x;
        A[0][1] += mass[i] * p_y * q_x;
        A[0][2] += mass[i] * p_z * q_x;
        A[1][0] += mass[i] * p_x * q_y;
        A[1][1] += mass[i] * p_y * q_y;
        A[1][2] += mass[i] * p_z * q_y;
        A[2][0] += mass[i] * p_x * q_z;
        A[2][1] += mass[i] * p_y * q_z;
        A[2][2] += mass[i] * p_z * q_z;
    }

    return A;
}

#include <svd3/svd3.h>
#define MATRIX_ARGUMENTS(M) M[0][0], M[0][1], M[0][2], M[1][0], M[1][1], M[1][2], M[2][0], M[2][1], M[2][2]

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
    svd(MATRIX_ARGUMENTS(M), MATRIX_ARGUMENTS(U), MATRIX_ARGUMENTS(S), MATRIX_ARGUMENTS(V));
    float max_val = math::max(S[0][0], math::max(S[1][1], S[2][2]));
    S = S / max_val;
    const mat3 Ut = glm::transpose(U);

    value[0] = S[0][0];
    value[1] = S[1][1];
    value[2] = S[2][2];

    vectors[0] = Ut[0];
    vectors[1] = Ut[1];
    vectors[2] = Ut[2];
}

static mat3 compute_rotation(const mat3& M) {
    mat3 R, S;
    decompose(M, &R, &S);
    return R;
}

mat3 compute_rotation(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0, const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
                      const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1) {
    const mat3 Apq = compute_covariance_matrix(x0, y0, z0, x1, y1, z1, mass, count, com0, com1) / (float)(count - 1);
    const mat3 Aqq = compute_covariance_matrix(x0, y0, z0, x0, y0, z0, mass, count, com0, com0) / (float)(count - 1);

    const mat3 A = Apq / Aqq;
    // return A; // Return complete linear transform with skewing and all
    return compute_rotation(A);  // Return rotational part
}

static void compute_residual_error(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, const float* RESTRICT src_x, const float* RESTRICT src_y, const float* RESTRICT src_z,
                                   const float* RESTRICT ref_x, const float* RESTRICT ref_y, const float* RESTRICT ref_z, int64 count, const mat4& matrix) {
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

/*
struct SupportFrame {
    struct {
        vec3 pos;
        float mass;
        } point[6];
};
*/

static void free_structure_data(Structure* s) {
    ASSERT(s);
    if (s->frame_data.transform) FREE(s->frame_data.transform);
    if (s->frame_data.eigen.vector) FREE(s->frame_data.eigen.vector);
    if (s->frame_data.eigen.value) FREE(s->frame_data.eigen.value);
    if (s->frame_data.determinant.abs) FREE(s->frame_data.determinant.abs);
    if (s->frame_data.support_frames) FREE(s->frame_data.support_frames);
}

static void init_structure_data(Structure* s, ID id, int32 ref_frame_idx, int32 num_points, int32 num_frames) {
    ASSERT(s);
    free_structure_data(s);
    s->id = id;
    s->ref_frame_idx = ref_frame_idx;
    s->num_frames = num_frames;
    s->num_points = num_points;

    s->frame_data.transform = (Transform*)MALLOC(sizeof(Transform) * num_frames);
    s->frame_data.eigen.vector = (mat3*)MALLOC(sizeof(mat3) * num_frames);
    s->frame_data.eigen.value = (vec3*)MALLOC(sizeof(vec3) * num_frames);

    {
        float* data = (float*)MALLOC(sizeof(float) * num_frames * 2);
        s->frame_data.determinant.abs = data + 0 * num_frames;
        s->frame_data.determinant.rel = data + 1 * num_frames;
    }

    s->frame_data.support_frames = (SupportFrame*)MALLOC(sizeof(SupportFrame) * num_frames);
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

void compute_support_frame(SupportFrame& frame, const mat3& eigen_vectors, const vec3& eigen_values, const vec3& com, float weight, const SupportFrame* ref = nullptr) {
    // clang-format off
	vec3 pos[6];
	float mass[6];

    for (int i = 0; i < 3; i++) {
        pos[i + 0] = com + eigen_vectors[i];
        pos[i + 1] = com - eigen_vectors[i];
		const float w = weight * eigen_values[i];
		mass[i + 0] = w;
		mass[i + 1] = w;
    }

	if (ref) {
		DynamicArray<int> avail_points(6);
		for (int i = 0; i < 6; i++) avail_points[i] = i;

		for (int i = 0; i < 6; i++) {
			// Find closest matching point in ref
			float min_d2 = FLT_MAX;
			int min_idx = 0;
			for (auto j : avail_points) {
				const float d2 = math::distance2(frame.point[i].pos, ref->point[j].pos);
				if (d2 < min_d2) {
					min_d2 = d2;
					min_idx = j;
				}
			}
		}
	}
	else {

	}
    // clang-format on
}

bool compute_trajectory_transform_data(ID id, Bitfield atom_mask, const MoleculeDynamic& dynamic, int32 target_frame_idx) {
    ASSERT(context);

    const int32 num_frames = (int32)dynamic.trajectory.num_frames;
    if (target_frame_idx >= num_frames) {
        LOG_ERROR("Supplied target frame index is out of range.");
        return false;
    }

    Structure* s = find_structure(id);
    if (s == nullptr) {
        LOG_ERROR("Could not compute tracking data, supplied id is not valid.");
        return false;
    }

    const int num_points = bitfield::number_of_bits_set(atom_mask);
    if (num_points == 0) {
        LOG_ERROR("Supplied atom mask is empty.");
        return false;
    }
    const int num_support_points = 6;
    const int tot_points = num_points + num_support_points;

    // Allocate memory for all data
    init_structure_data(s, id, target_frame_idx, num_points, num_frames);

    // Scratch data
    const auto mem_size = sizeof(float) * tot_points * 10;
    void* mem = TMP_MALLOC(mem_size);
    defer { TMP_FREE(mem); };

    memset(mem, 0, mem_size);

    float* cur_x = (float*)mem;
    float* cur_y = cur_x + tot_points;
    float* cur_z = cur_y + tot_points;
    float* ref_x = cur_z + tot_points;
    float* ref_y = ref_x + tot_points;
    float* ref_z = ref_y + tot_points;
    float* prv_x = ref_z + tot_points;
    float* prv_y = prv_x + tot_points;
    float* prv_z = prv_y + tot_points;
    float* mass = prv_z + tot_points;

    bitfield::extract_data_from_mask(cur_x, get_trajectory_position_x(dynamic.trajectory, target_frame_idx).data(), atom_mask);
    bitfield::extract_data_from_mask(cur_y, get_trajectory_position_y(dynamic.trajectory, target_frame_idx).data(), atom_mask);
    bitfield::extract_data_from_mask(cur_z, get_trajectory_position_z(dynamic.trajectory, target_frame_idx).data(), atom_mask);
    bitfield::extract_data_from_mask(ref_x, get_trajectory_position_x(dynamic.trajectory, target_frame_idx).data(), atom_mask);
    bitfield::extract_data_from_mask(ref_y, get_trajectory_position_y(dynamic.trajectory, target_frame_idx).data(), atom_mask);
    bitfield::extract_data_from_mask(ref_z, get_trajectory_position_z(dynamic.trajectory, target_frame_idx).data(), atom_mask);

    bitfield::extract_data_from_mask(mass, dynamic.molecule.atom.mass, atom_mask);

    const vec3 ref_com = compute_com(ref_x, ref_y, ref_z, mass, num_points);

    float ref_support_point_mass[6];
    float tot_mass = 0.0f;
    for (int i = 0; i < num_points; i++) {
        tot_mass += mass[i];
    }

    vec3 prv_com = {0, 0, 0};
    vec3 cur_com = {0, 0, 0};

    const mat3 ref_cov_mat = compute_covariance_matrix(ref_x, ref_y, ref_z, ref_x, ref_y, ref_z, mass, num_points, ref_com, ref_com);
    mat3 ref_eigen_vectors;
    vec3 ref_eigen_values;
    compute_eigen(ref_cov_mat, (vec3(&)[3])ref_eigen_vectors, (float(&)[3])ref_eigen_values);

    // Set target frame explicitly
    s->frame_data.transform[target_frame_idx] = {mat3(1), ref_com};
    s->frame_data.eigen.vector[target_frame_idx] = ref_eigen_vectors;
    s->frame_data.eigen.value[target_frame_idx] = ref_eigen_values;
    s->frame_data.determinant.abs[target_frame_idx] = 1.0f;
    s->frame_data.determinant.rel[target_frame_idx] = 1.0f;

    const float support_w = tot_mass * 0.1f;
    compute_support_frame(s->frame_data.support_frames[target_frame_idx], ref_eigen_vectors, ref_eigen_values, ref_com, support_w);

    for (int32 i = 0; i < num_frames; i++) {
        if (i == target_frame_idx) continue;

        // Copy previous
        memcpy(prv_x, cur_x, num_points * sizeof(float) * 3);
        prv_com = cur_com;

        // Fetch current
        bitfield::extract_data_from_mask(cur_x, get_trajectory_position_x(dynamic.trajectory, i).data(), atom_mask);
        bitfield::extract_data_from_mask(cur_y, get_trajectory_position_y(dynamic.trajectory, i).data(), atom_mask);
        bitfield::extract_data_from_mask(cur_z, get_trajectory_position_z(dynamic.trajectory, i).data(), atom_mask);

        cur_com = compute_com(cur_x, cur_y, cur_z, mass, num_points);

        const mat3 abs_mat = compute_covariance_matrix(cur_x, cur_y, cur_z, ref_x, ref_y, ref_z, mass, num_points, cur_com, ref_com);
        const mat3 rel_mat = compute_covariance_matrix(cur_x, cur_y, cur_z, prv_x, prv_y, prv_z, mass, num_points, cur_com, prv_com);
        const mat3 cov_mat = compute_covariance_matrix(cur_x, cur_y, cur_z, cur_x, cur_y, cur_z, mass, num_points, cur_com, cur_com);

        const mat3 cur_rot = (abs_mat / cov_mat);

        const float abs_det = math::determinant(abs_mat / cov_mat);
        const float rel_det = math::determinant(rel_mat / cov_mat);

        compute_eigen(cov_mat, (vec3(&)[3])s->frame_data.eigen.vector[i], (float(&)[3])s->frame_data.eigen.value[i]);

        s->frame_data.transform[i].rotation = cur_rot;
        s->frame_data.transform[i].com = cur_com;

        s->frame_data.determinant.abs[i] = abs_det;
        s->frame_data.determinant.rel[i] = rel_det;
    }

    return true;
}

const Transform& get_transform_to_target_frame(ID id, int32 source_frame) {
    ASSERT(context);
    static const Transform invalid_transform{};

    Structure* s = find_structure(id);
    if (s == nullptr) {
        LOG_ERROR("Supplied id is not valid.");
        return invalid_transform;
    }

    if (source_frame < 0 || s->num_frames <= source_frame) {
        LOG_ERROR("Supplied frame is out of range");
        return invalid_transform;
    }

    return s->frame_data.transform[source_frame];
}

const Array<const mat3> get_eigen_vectors(ID id) {
    Structure* s = find_structure(id);
    if (s == nullptr) {
        LOG_ERROR("Supplied id is not valid.");
        return {};
    }

    return {s->frame_data.eigen.vector, s->num_frames};
}

const Array<const vec3> get_eigen_values(ID id) {
    Structure* s = find_structure(id);
    if (s == nullptr) {
        LOG_ERROR("Supplied id is not valid.");
        return {};
    }

    return {s->frame_data.eigen.value, s->num_frames};
}

const Array<const float> get_abs_det(ID id) {
    Structure* s = find_structure(id);
    if (s == nullptr) {
        LOG_ERROR("Supplied id is not valid.");
        return {};
    }

    return {s->frame_data.determinant.abs, s->num_frames};
}

const Array<const float> get_rel_det(ID id) {
    Structure* s = find_structure(id);
    if (s == nullptr) {
        LOG_ERROR("Supplied id is not valid.");
        return {};
    }

    return {s->frame_data.determinant.rel, s->num_frames};
}

}  // namespace structure_tracking