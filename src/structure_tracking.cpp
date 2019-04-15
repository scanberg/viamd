#include "structure_tracking.h"

#include <core/log.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#pragma warning( disable : 4127 ) // disable warnings about expressions which could be constexpr in Eigen
#include <Eigen/Eigen>

namespace structure_tracking {

struct Structure {
	ID id = 0;
	int32 ref_frame_idx = 0;
	int32 num_points = 0;
	int32 num_frames = 0;

	struct {
		Transform* transform = nullptr;
		struct {
			struct {
				float* x = nullptr;
				float* y = nullptr;
				float* z = nullptr;
			} vector[3];
			float* value[3] = { nullptr, nullptr, nullptr };
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
static void Diagonalize(const float(&A)[3][3], float(&Q)[3][3], float(&D)[3][3]) {
	// A must be a symmetric matrix.
	// returns Q and D such that
	// Diagonal matrix D = QT * A * Q;  and  A = Q*D*QT
	const int maxsteps = 24;  // certainly wont need that many.
	int k0, k1, k2;
	float o[3], m[3];
	float q[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
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
	const float det = math::determinant(AtA); // For debugging
	D[0][0] = sqrtf(D[0][0]);
	D[1][1] = sqrtf(D[1][1]);
	D[2][2] = sqrtf(D[2][2]);
	// @NOTE: Should one zero every non diagonal element???
	*S = math::inverse(Q) * D * Q;
	*R = M * math::inverse(*S);
}

static mat3 compute_linear_transform(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
								     const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
									 const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1) {
	mat3 Apq{ 0 };
	mat3 Aqq{ 0 };

	for (int64 i = 0; i < count; i++) {
		// @TODO: Vectorize...
		const float q_x = x0[i] - com0.x;
		const float q_y = y0[i] - com0.y;
		const float q_z = z0[i] - com0.z;

		const float p_x = x1[i] - com1.x;
		const float p_y = y1[i] - com1.y;
		const float p_z = z1[i] - com1.z;

		Apq[0][0] += mass[i] * p_x * q_x;
		Apq[0][1] += mass[i] * p_y * q_x;
		Apq[0][2] += mass[i] * p_z * q_x;
		Apq[1][0] += mass[i] * p_x * q_y;
		Apq[1][1] += mass[i] * p_y * q_y;
		Apq[1][2] += mass[i] * p_z * q_y;
		Apq[2][0] += mass[i] * p_x * q_z;
		Apq[2][1] += mass[i] * p_y * q_z;
		Apq[2][2] += mass[i] * p_z * q_z;

		Aqq[0][0] += mass[i] * q_x * q_x;
		Aqq[0][1] += mass[i] * q_y * q_x;
		Aqq[0][2] += mass[i] * q_z * q_x;
		Aqq[1][0] += mass[i] * q_x * q_y;
		Aqq[1][1] += mass[i] * q_y * q_y;
		Aqq[1][2] += mass[i] * q_z * q_y;
		Aqq[2][0] += mass[i] * q_x * q_z;
		Aqq[2][1] += mass[i] * q_y * q_z;
		Aqq[2][2] += mass[i] * q_z * q_z;
	}

	return Apq / Aqq;
}

static mat3 compute_covariance_matrix(const float* x0, const float* y0, const float* z0,
									  const float* x1, const float* y1, const float* z1,
									  const float* mass, int64 count, const vec3& com0, const vec3& com1)
{
	mat3 A{ 0 };
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

static mat3 compute_pq_matrix(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
							  const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
							  const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1)
{
	mat3 Apq{ 0 };

	for (int64 i = 0; i < count; i++) {
		// @TODO: Vectorize...
		const float q_x = x0[i] - com0.x;
		const float q_y = y0[i] - com0.y;
		const float q_z = z0[i] - com0.z;

		const float p_x = x1[i] - com1.x;
		const float p_y = y1[i] - com1.y;
		const float p_z = z1[i] - com1.z;

		Apq[0][0] += mass[i] * p_x * q_x;
		Apq[0][1] += mass[i] * p_y * q_x;
		Apq[0][2] += mass[i] * p_z * q_x;
		Apq[1][0] += mass[i] * p_x * q_y;
		Apq[1][1] += mass[i] * p_y * q_y;
		Apq[1][2] += mass[i] * p_z * q_y;
		Apq[2][0] += mass[i] * p_x * q_z;
		Apq[2][1] += mass[i] * p_y * q_z;
		Apq[2][2] += mass[i] * p_z * q_z;
	}

	return Apq / (float)(count - 1);
}

// Compute weighted covariance matrix from point data (x,y,z,weight) 
static mat3 compute_qq_matrix(const float* RESTRICT x, const float* RESTRICT y, const float* RESTRICT z, const float* RESTRICT w,
						      int64 count, const vec3& com)
{
	mat3 Aqq{ 0 };
	for (int64 i = 0; i < count; i++) {
		const float q_x = x[i] - com.x;
		const float q_y = y[i] - com.y;
		const float q_z = z[i] - com.z;

		Aqq[0][0] += w[i] * q_x * q_x;
		Aqq[0][1] += w[i] * q_y * q_x;
		Aqq[0][2] += w[i] * q_z * q_x;
		Aqq[1][0] += w[i] * q_x * q_y;
		Aqq[1][1] += w[i] * q_y * q_y;
		Aqq[1][2] += w[i] * q_z * q_y;
		Aqq[2][0] += w[i] * q_x * q_z;
		Aqq[2][1] += w[i] * q_y * q_z;
		Aqq[2][2] += w[i] * q_z * q_z;
	}

	return Aqq / (float)(count - 1);
}

#include <svd3/svd3.h>
#define MATRIX_ARGUMENTS(M) M[0][0], M[0][1], M[0][2], M[1][0], M[1][1], M[1][2], M[2][0], M[2][1], M[2][2]

static void compute_eigen(const mat3& M, vec3(&vectors)[3], float(&value)[3]) {
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

static mat3 compute_rotation_SVD(const mat3& M) {

	/*
	Eigen::Matrix3f B_eigen = Eigen::Matrix3f({
					{M[0][0], M[1][0], M[2][0]},
					{M[0][1], M[1][1], M[2][1]},
					{M[0][2], M[1][2], M[2][2]} });

	Eigen::JacobiSVD<Eigen::Matrix3f> svd(B_eigen, Eigen::ComputeFullU | Eigen::ComputeFullV);

	const auto U = svd.matrixU();
	const auto V = svd.matrixV();
	auto s = svd.singularValues();

	s[0] = 1;
	s[1] = 1;
	s[2] = U.determinant() * V.determinant();
	const Eigen::Matrix3f D = Eigen::DiagonalMatrix<float, 3, 3>(s[0], s[1], s[2]).toDenseMatrix();

	const auto A_eigen = U * D * V.transpose();
	mat3 A = {A_eigen(0,0), A_eigen(0,1), A_eigen(0,2),
			  A_eigen(1,0), A_eigen(1,1), A_eigen(1,2),
			  A_eigen(2,0), A_eigen(2,1), A_eigen(2,2) };

			  */
	mat3 U, S, V;
	svd(MATRIX_ARGUMENTS(M), MATRIX_ARGUMENTS(U), MATRIX_ARGUMENTS(S), MATRIX_ARGUMENTS(V));
	vec3 s = { 1, 1, math::determinant(U) * math::determinant(V) };
	mat3 D = mat3(1);
	D[2][2] = math::determinant(U) * math::determinant(V);

	return U * D * math::transpose(V);
}

static mat3 compute_rotation(const mat3& M) {
	mat3 R, S;
	decompose(M, &R, &S);
	return R;
}

mat3 compute_rotation(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
					  const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
					  const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1) {
	const mat3 M = compute_covariance_matrix(x0, y0, z0, x1, y1, z1, mass, count, com0, com1) / (float)(count - 1);
	return compute_rotation(M);
}

static void compute_residual_error(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
								   const float* RESTRICT src_x, const float* RESTRICT src_y, const float* RESTRICT src_z,
								   const float* RESTRICT ref_x, const float* RESTRICT ref_y, const float* RESTRICT ref_z,
								   int64 count, const mat4& matrix) {
	for (int32 i = 0; i < count; i++) {
		// @TODO: Vectorize this...

		const vec4 r = { ref_x[i], ref_y[i], ref_z[i], 1.0f };
		const vec4 v = { src_x[i], src_y[i], src_z[i], 1.0f };
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
	if (s->frame_data.transform) FREE(s->frame_data.transform);
	if (s->frame_data.eigen.vector[0].x) FREE(s->frame_data.eigen.vector[0].x);
	if (s->frame_data.determinant.abs) FREE(s->frame_data.determinant.abs);
}

static void init_structure_data(Structure* s, ID id, int32 ref_frame_idx, int32 num_points, int32 num_frames) {
	ASSERT(s);
	free_structure_data(s);
	s->id = id;
	s->ref_frame_idx = ref_frame_idx;
	s->num_frames = num_frames;
	s->num_points = num_points;

	s->frame_data.transform = (Transform*)MALLOC(sizeof(Transform) * num_frames);

	float* eigen_data = (float*)MALLOC(sizeof(float) * num_frames * 12);
	s->frame_data.eigen.vector[0].x = eigen_data + num_frames * 0;
	s->frame_data.eigen.vector[0].y = eigen_data + num_frames * 1;
	s->frame_data.eigen.vector[0].z = eigen_data + num_frames * 2;
	s->frame_data.eigen.vector[1].x = eigen_data + num_frames * 3;
	s->frame_data.eigen.vector[1].y = eigen_data + num_frames * 4;
	s->frame_data.eigen.vector[1].z = eigen_data + num_frames * 5;
	s->frame_data.eigen.vector[2].x = eigen_data + num_frames * 6;
	s->frame_data.eigen.vector[2].y = eigen_data + num_frames * 7;
	s->frame_data.eigen.vector[2].z = eigen_data + num_frames * 8;
	s->frame_data.eigen.value[0]	= eigen_data + num_frames * 9;
	s->frame_data.eigen.value[1]	= eigen_data + num_frames * 10;
	s->frame_data.eigen.value[2]	= eigen_data + num_frames * 11;

	float* det_data = (float*)MALLOC(sizeof(float) * num_frames * 2);
	s->frame_data.determinant.abs = det_data;
	s->frame_data.determinant.rel = det_data + num_frames;
/*
	s->data.rbf.radial_cutoff = radial_cutoff;
	s->data.rbf.weight.x = (float*)MALLOC(sizeof(float) * num_points * num_frames * 3);
	s->data.rbf.weight.y = s->data.rbf.weight.x + num_points * num_frames;
	s->data.rbf.weight.z = s->data.rbf.weight.y + num_points * num_frames;
	s->data.rbf.pos.x = (float*)MALLOC(sizeof(float) * num_points * 3);
	s->data.rbf.pos.y = s->data.rbf.pos.x + num_points;
	s->data.rbf.pos.z = s->data.rbf.pos.y + num_points;
	*/
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
		context->~Context();
		FREE(context);
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

template <typename T>
void extract_data_from_indices(T* RESTRICT dst_data, const T* RESTRICT src_data, int32* RESTRICT indices, int32 count) {
	for (int32 i = 0; i < count; i++) {
		dst_data[i] = src_data[indices[i]];
	}
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

	int* indices = (int*)TMP_MALLOC(sizeof(int) * atom_mask.size());
	defer { TMP_FREE(indices); };
	int num_points = 0;
	for (int i = 0; i < (int32)atom_mask.size(); i++) {
		if (bitfield::get_bit(atom_mask, i)) {
			indices[num_points] = i;
			num_points++;
		}
	}

	if (num_points == 0) {
		LOG_ERROR("Supplied atom mask is empty.");
		return false;
	}

	// Allocate memory for all data
	init_structure_data(s, id, target_frame_idx, num_points, num_frames);

	// Scratch data
	const auto mem_size = sizeof(float) * num_points * 10;
	void* mem = TMP_MALLOC(mem_size);
	defer{ TMP_FREE(mem); };

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
	float* mass = prv_z + num_points;

	extract_data_from_indices(ref_x, get_trajectory_position_x(dynamic.trajectory, target_frame_idx).data(), indices, num_points);
	extract_data_from_indices(ref_y, get_trajectory_position_y(dynamic.trajectory, target_frame_idx).data(), indices, num_points);
	extract_data_from_indices(ref_z, get_trajectory_position_z(dynamic.trajectory, target_frame_idx).data(), indices, num_points);
	extract_data_from_indices(mass, dynamic.molecule.atom.mass, indices, num_points);
	const vec3 ref_com = compute_com(ref_x, ref_y, ref_z, mass, num_points);

	vec3 prv_com = { 0,0,0 };
	vec3 cur_com = { 0,0,0 };

	// Set target frame explicitly
	s->frame_data.transform[target_frame_idx] = { mat3(1), ref_com };
	s->frame_data.eigen.vector[0].x[target_frame_idx] = 0;
	s->frame_data.eigen.vector[0].y[target_frame_idx] = 0;
	s->frame_data.eigen.vector[0].z[target_frame_idx] = 0;
	s->frame_data.eigen.vector[1].x[target_frame_idx] = 0;
	s->frame_data.eigen.vector[1].y[target_frame_idx] = 0;
	s->frame_data.eigen.vector[1].z[target_frame_idx] = 0;
	s->frame_data.eigen.vector[2].x[target_frame_idx] = 0;
	s->frame_data.eigen.vector[2].y[target_frame_idx] = 0;
	s->frame_data.eigen.vector[2].z[target_frame_idx] = 0;
	s->frame_data.eigen.value[0][target_frame_idx] = 0;
	s->frame_data.eigen.value[1][target_frame_idx] = 0;
	s->frame_data.eigen.value[2][target_frame_idx] = 0;

	s->frame_data.determinant.abs[target_frame_idx] = 1.0f;
	s->frame_data.determinant.rel[target_frame_idx] = 1.0f;

	for (int32 cur_idx = 0; cur_idx < num_frames; cur_idx++) {
		if (cur_idx == target_frame_idx) continue;
		extract_data_from_indices(cur_x, get_trajectory_position_x(dynamic.trajectory, cur_idx).data(), indices, num_points);
		extract_data_from_indices(cur_y, get_trajectory_position_y(dynamic.trajectory, cur_idx).data(), indices, num_points);
		extract_data_from_indices(cur_z, get_trajectory_position_z(dynamic.trajectory, cur_idx).data(), indices, num_points);

		memcpy(prv_x, cur_x, num_points * sizeof(float) * 3);
		prv_com = cur_com;

		cur_com = compute_com(cur_x, cur_y, cur_z, mass, num_points);
		
		//float* rbf_x = s->data.rbf.weight.x + i * num_points;
		//float* rbf_y = s->data.rbf.weight.y + i * num_points;
		//float* rbf_z = s->data.rbf.weight.z + i * num_points;

		// @NOTE: Compute linear transformation matrix between the two sets of points.

		const mat3 abs_mat = compute_covariance_matrix(cur_x, cur_y, cur_z, ref_x, ref_y, ref_z, mass, num_points, cur_com, ref_com);
		const mat3 rel_mat = compute_covariance_matrix(cur_x, cur_y, cur_z, prv_x, prv_y, prv_z, mass, num_points, cur_com, prv_com);
		const mat3 cov_mat = compute_covariance_matrix(cur_x, cur_y, cur_z, cur_x, cur_y, cur_z, mass, num_points, cur_com, cur_com);

		const mat3 cur_rot = (abs_mat / cov_mat);

		const float abs_det = math::determinant(abs_mat / cov_mat);
		const float rel_det = math::determinant(rel_mat / cov_mat);

		vec3 eigen_vectors[3];
		float eigen_values[3];
		compute_eigen(cov_mat, eigen_vectors, eigen_values);
	
		// @NOTE: Compute residual error between the two sets of points.
		//compute_residual_error(err_x, err_y, err_z, cur_x, cur_y, cur_z, ref_x, ref_y, ref_z, num_points, *M);

		// @NOTE: Encode residual error as rbf at reference points
		//compute_rbf_weights(rbf_x, rbf_y, rbf_z, ref_x, ref_y, ref_z, err_x, err_y, err_z, num_points, rbf_radial_cutoff);

		s->frame_data.transform[cur_idx].rotation = cur_rot;
		s->frame_data.transform[cur_idx].com = cur_com;

		s->frame_data.eigen.vector[0].x[cur_idx] = eigen_vectors[0].x;
		s->frame_data.eigen.vector[0].y[cur_idx] = eigen_vectors[0].y;
		s->frame_data.eigen.vector[0].z[cur_idx] = eigen_vectors[0].z;

		s->frame_data.eigen.vector[1].x[cur_idx] = eigen_vectors[1].x;
		s->frame_data.eigen.vector[1].y[cur_idx] = eigen_vectors[1].y;
		s->frame_data.eigen.vector[1].z[cur_idx] = eigen_vectors[1].z;

		s->frame_data.eigen.vector[2].x[cur_idx] = eigen_vectors[2].x;
		s->frame_data.eigen.vector[2].y[cur_idx] = eigen_vectors[2].y;
		s->frame_data.eigen.vector[2].z[cur_idx] = eigen_vectors[2].z;

		s->frame_data.eigen.value[0][cur_idx] = eigen_values[0];
		s->frame_data.eigen.value[1][cur_idx] = eigen_values[1];
		s->frame_data.eigen.value[2][cur_idx] = eigen_values[2];

		s->frame_data.determinant.abs[cur_idx] = abs_det;
		s->frame_data.determinant.rel[cur_idx] = rel_det;
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

const Array<const float> get_eigen_vector_x(ID id, int64 idx) {
	Structure* s = find_structure(id);
	if (s == nullptr) {
		LOG_ERROR("Supplied id is not valid.");
		return {};
	}

	if (idx < 0 || 2 < idx) {
		LOG_ERROR("Supplied idx[%lli] is out of range [0,2]", idx);
		return {};
	}

	return { s->frame_data.eigen.vector[idx].x, s->num_frames };
}

const Array<const float> get_eigen_vector_y(ID id, int64 idx) {
	Structure* s = find_structure(id);
	if (s == nullptr) {
		LOG_ERROR("Supplied id is not valid.");
		return {};
	}

	if (idx < 0 || 2 < idx) {
		LOG_ERROR("Supplied idx[%lli] is out of range [0,2]", idx);
		return {};
	}

	return { s->frame_data.eigen.vector[idx].y, s->num_frames };
}

const Array<const float> get_eigen_vector_z(ID id, int64 idx) {
	Structure* s = find_structure(id);
	if (s == nullptr) {
		LOG_ERROR("Supplied id is not valid.");
		return {};
	}

	if (idx < 0 || 2 < idx) {
		LOG_ERROR("Supplied idx[%lli] is out of range [0,2]", idx);
		return {};
	}

	return { s->frame_data.eigen.vector[idx].z, s->num_frames };
}

const Array<const float> get_eigen_value(ID id, int64 idx) {
	Structure* s = find_structure(id);
	if (s == nullptr) {
		LOG_ERROR("Supplied id is not valid.");
		return {};
	}

	if (idx < 0 || 2 < idx) {
		LOG_ERROR("Supplied idx[%lli] is out of range [0,2]", idx);
		return {};
	}

	return { s->frame_data.eigen.value[idx], s->num_frames };
}

const Array<const float> get_abs_det(ID id) {
	Structure* s = find_structure(id);
	if (s == nullptr) {
		LOG_ERROR("Supplied id is not valid.");
		return {};
	}

	return { s->frame_data.determinant.abs, s->num_frames };
}

const Array<const float> get_rel_det(ID id) {
	Structure* s = find_structure(id);
	if (s == nullptr) {
		LOG_ERROR("Supplied id is not valid.");
		return {};
	}

	return { s->frame_data.determinant.rel, s->num_frames };
}

} // namespace structure_tracking