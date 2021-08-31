#include "structure_tracking.h"

#include <core/log.h>
#include <core/hash.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#include <svd3/svd3.h>

namespace structure_tracking {

struct Structure {
    ID id = 0;
    TrackingData tracking_data;
};

struct Entry {
    ID id;
    Structure* ptr;
};

struct Context {
    u32 next_hash = 0xdeadf00d;
    DynamicArray<Entry> entries{};
};

Context* context = nullptr;

// clang-format off
static mat3 compute_weighted_cross_covariance_matrix(const float* x0, const float* y0, const float* z0,
                                                     const float* x1, const float* y1, const float* z1,
                                                     const float* weight,
                                                     i64 count,
                                                     const vec3& com0 = {0,0,0}, const vec3& com1 = {0,0,0})
// clang-format on
{
    mat3 A{0};
    for (i64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float px = x0[i] - com0.x;
        const float py = y0[i] - com0.y;
        const float pz = z0[i] - com0.z;

        const float qx = x1[i] - com1.x;
        const float qy = y1[i] - com1.y;
        const float qz = z1[i] - com1.z;

        const float w = weight[i];

        A[0][0] += w * px * qx;
        A[0][1] += w * px * qy;
        A[0][2] += w * px * qz;
        A[1][0] += w * py * qx;
        A[1][1] += w * py * qy;
        A[1][2] += w * py * qz;
        A[2][0] += w * pz * qx;
        A[2][1] += w * pz * qy;
        A[2][2] += w * pz * qz;
    }

    return A;
}

// clang-format off
static mat3 compute_weighted_covariance_matrix(const float* x, const float* y, const float* z,
                                               const float* weight,
                                               i64 count, const vec3& com = {0,0,0})
// clang-format on
{
    mat3 A{0};
    for (i64 i = 0; i < count; i++) {
        // @TODO: Vectorize...
        const float qx = x[i] - com.x;
        const float qy = y[i] - com.y;
        const float qz = z[i] - com.z;
        const float w = weight[i];

        A[0][0] += w * qx * qx;
        A[0][1] += w * qx * qy;
        A[0][2] += w * qx * qz;
        A[1][0] += w * qy * qx;
        A[1][1] += w * qy * qy;
        A[1][2] += w * qy * qz;
        A[2][0] += w * qz * qx;
        A[2][1] += w * qz * qy;
        A[2][2] += w * qz * qz;
    }

    return A;
}

#define ARGS(M) M[0][0], M[1][0], M[2][0], M[0][1], M[1][1], M[2][1], M[0][2], M[1][2], M[2][2]

static void compute_eigen(const mat3& M, vec3 (&vectors)[3], float (&value)[3]) {
    mat3 U, S, V;
    svd(ARGS(M), ARGS(U), ARGS(S), ARGS(V));
    float max_val = math::max(S[0][0], math::max(S[1][1], S[2][2]));

    const float e_val[] = {S[0][0] / max_val, S[1][1] / max_val, S[2][2] / max_val};
    const vec3 e_vec[] = {U[0], U[1], U[2]};
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
//static void compute_svd(const dmat3& A, dmat3& U, dmat3& S, dmat3& V) { svd(ARGS(A), ARGS(U), ARGS(S), ARGS(V)); }

#undef ARGS

vec3 compute_eigen_values(const float* RESTRICT x, const float* RESTRICT y, const float* RESTRICT z, const float* RESTRICT mass, i64 count) {
    const vec3 com = compute_com(x, y, z, mass, count);
    const mat3 A = compute_weighted_covariance_matrix(x, y, z, mass, count, com);
    vec3 vecs[3];
    vec3 vals;
    compute_eigen(A, vecs, (float(&)[3])vals);
    return vals;
}

static mat3 extract_rotation(const mat3& M) {
    mat3 U, S, V;
    compute_svd(M, U, S, V);

    const mat3 Ut = math::transpose(U);
    const float d = math::determinant(V * Ut);
    const mat3 D = {1, 0, 0, 0, 1, 0, 0, 0, d};
    const mat3 R = V * D * Ut;
    return R;
}

// clang-format off
mat3 compute_rotation(const float* x0, const float* y0, const float* z0,
                      const float* x1, const float* y1, const float* z1,
                      const float* mass, i64 count, const vec3& com0, const vec3& com1)
// clang-format on
{
    const mat3 Apq = compute_weighted_cross_covariance_matrix(x0, y0, z0, x1, y1, z1, mass, count, com0, com1);
    return extract_rotation(Apq);  // Return rotational part
}

#if 0
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
#endif

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
    if (s->tracking_data.transform.rotation.absolute) FREE(s->tracking_data.transform.rotation.absolute);
    if (s->tracking_data.transform.rotation.relative) FREE(s->tracking_data.transform.rotation.relative);
    if (s->tracking_data.transform.rotation.hybrid) FREE(s->tracking_data.transform.rotation.hybrid);
    if (s->tracking_data.transform.com) FREE(s->tracking_data.transform.com);
    if (s->tracking_data.eigen.vectors) FREE(s->tracking_data.eigen.vectors);
    if (s->tracking_data.eigen.values) FREE(s->tracking_data.eigen.values);
    if (s->tracking_data.average_structure.x) FREE(s->tracking_data.average_structure.x);
    if (s->tracking_data.average_structure.y) FREE(s->tracking_data.average_structure.y);
    if (s->tracking_data.average_structure.z) FREE(s->tracking_data.average_structure.z);
    *s = {};
}

static void init_structure_data(Structure* s, ID id, i32 num_frames, i32 num_atoms) {
    ASSERT(s);
    free_structure_data(s);
    s->id = id;
    s->tracking_data.count = num_frames;

    s->tracking_data.transform.rotation.absolute = (quat*)MALLOC(sizeof(quat) * num_frames);
    s->tracking_data.transform.rotation.relative = (quat*)MALLOC(sizeof(quat) * num_frames);
    s->tracking_data.transform.rotation.hybrid = (quat*)MALLOC(sizeof(quat) * num_frames);
    s->tracking_data.transform.com = (vec3*)MALLOC(sizeof(vec3) * num_frames);
    s->tracking_data.eigen.vectors = (mat3*)MALLOC(sizeof(mat3) * num_frames);
    s->tracking_data.eigen.values = (vec3*)MALLOC(sizeof(vec3) * num_frames);
    s->tracking_data.average_structure.atom_count = num_atoms;
    s->tracking_data.average_structure.x = (float*)MALLOC(sizeof(float) * num_atoms);
    s->tracking_data.average_structure.y = (float*)MALLOC(sizeof(float) * num_atoms);
    s->tracking_data.average_structure.z = (float*)MALLOC(sizeof(float) * num_atoms);
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

    Structure* s = NEW(Structure);
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

#endif

bool compute_trajectory_transform_data(ID id, const MoleculeDynamic& dynamic, Bitfield atom_mask, i64 mask_offset) {
    ASSERT(context);

    const i32 num_frames = (i32)dynamic.trajectory.num_frames;

    Structure* s = find_structure(id);
    if (s == nullptr) {
        LOG_ERROR("Could not compute tracking data, supplied id is not valid.");
        return false;
    }

    const int num_atoms = (int)bitfield::number_of_bits_set(atom_mask);
    if (num_atoms == 0) {
        LOG_ERROR("Supplied atom mask is empty.");
        return false;
    }

    // Allocate memory for all data
    init_structure_data(s, id, num_frames, num_atoms);

    // Scratch data
    const auto flt_mem_size = sizeof(float) * num_atoms * 10;
    void* flt_mem = TMP_MALLOC(flt_mem_size);
    defer { TMP_FREE(flt_mem); };
    memset(flt_mem, 0, flt_mem_size);

    float* cur_x = (float*)flt_mem + 0 * num_atoms;
    float* cur_y = (float*)flt_mem + 1 * num_atoms;
    float* cur_z = (float*)flt_mem + 2 * num_atoms;
    float* ref_x = (float*)flt_mem + 3 * num_atoms;
    float* ref_y = (float*)flt_mem + 4 * num_atoms;
    float* ref_z = (float*)flt_mem + 5 * num_atoms;
    float* int_x = (float*)flt_mem + 6 * num_atoms;
    float* int_y = (float*)flt_mem + 7 * num_atoms;
    float* int_z = (float*)flt_mem + 8 * num_atoms;
    float* mass  = (float*)flt_mem + 9 * num_atoms;

    mat3 box = {};

    bitfield::gather_masked(ref_x, get_trajectory_position_x(dynamic.trajectory, 0).data(), atom_mask, mask_offset);
    bitfield::gather_masked(ref_y, get_trajectory_position_y(dynamic.trajectory, 0).data(), atom_mask, mask_offset);
    bitfield::gather_masked(ref_z, get_trajectory_position_z(dynamic.trajectory, 0).data(), atom_mask, mask_offset);
    bitfield::gather_masked(mass, dynamic.molecule.atom.mass, atom_mask, mask_offset);
    box = get_trajectory_frame(dynamic.trajectory, 0).box;
    apply_pbc({ref_x, ref_y, ref_z},  num_atoms, box);

    memcpy(int_x, ref_x, num_atoms * sizeof(float) * 3);

    const vec3 ref_com = compute_com(ref_x, ref_y, ref_z, mass, num_atoms);
    vec3 int_com = ref_com;
    vec3 cur_com = {0, 0, 0};

    float* avg_x = s->tracking_data.average_structure.x;
    float* avg_y = s->tracking_data.average_structure.y;
    float* avg_z = s->tracking_data.average_structure.z;

    memcpy(avg_x, ref_x, num_atoms * sizeof(float));
    memcpy(avg_y, ref_y, num_atoms * sizeof(float));
    memcpy(avg_z, ref_z, num_atoms * sizeof(float));
    translate({avg_x, avg_y, avg_z}, num_atoms, -ref_com);

    // Compute average structure
    for (int i = 1; i < num_frames; i++) {
        bitfield::gather_masked(cur_x, dynamic.trajectory.frame_buffer[i].atom_position.x, atom_mask, mask_offset);
        bitfield::gather_masked(cur_y, dynamic.trajectory.frame_buffer[i].atom_position.y, atom_mask, mask_offset);
        bitfield::gather_masked(cur_z, dynamic.trajectory.frame_buffer[i].atom_position.z, atom_mask, mask_offset);
        cur_com = compute_com(cur_x, cur_y, cur_z, num_atoms);
        translate({cur_x, cur_y, cur_z}, num_atoms, -cur_com);

        const mat3 R = extract_rotation(compute_weighted_cross_covariance_matrix(cur_x, cur_y, cur_z, avg_x, avg_y, avg_z, mass, num_atoms));
        transform({cur_x, cur_y, cur_z}, num_atoms, R);

        const float w = 1.0f / (float)(i + 1);
        for (int j = 0; j < num_atoms; j++) {
            // Moving average
            avg_x[j] += (cur_x[j] - avg_x[j]) * w;
            avg_y[j] += (cur_y[j] - avg_y[j]) * w;
            avg_z[j] += (cur_z[j] - avg_z[j]) * w;
        }
    }

    // const mat3 AVG_PCA = compute_eigen_frame(avg_x, avg_y, avg_z, mass, num_atoms).vectors;
    // transform(avg_x, avg_y, avg_z, num_atoms, AVG_PCA);

    quat q_relative = {1, 0, 0, 0};
    quat q_hybrid = {1, 0, 0, 0};

    // memcpy(prv_x, ref_x, num_atoms * sizeof(float) * 3);

    // Set first frame explicitly
    s->tracking_data.transform.rotation.absolute[0] = {1, 0, 0, 0};
    s->tracking_data.transform.rotation.relative[0] = {1, 0, 0, 0};
    s->tracking_data.transform.com[0] = ref_com;

    compute_eigen(compute_weighted_covariance_matrix(ref_x, ref_y, ref_z, mass, num_atoms, ref_com), (vec3(&)[3])s->tracking_data.eigen.vectors[0],
                  (float(&)[3])s->tracking_data.eigen.values[0]);

    if (num_atoms > 1) {
        for (i32 i = 1; i < num_frames; ++i) {
            // Fetch current
            bitfield::gather_masked(cur_x, dynamic.trajectory.frame_buffer[i].atom_position.x, atom_mask, mask_offset);
            bitfield::gather_masked(cur_y, dynamic.trajectory.frame_buffer[i].atom_position.y, atom_mask, mask_offset);
            bitfield::gather_masked(cur_z, dynamic.trajectory.frame_buffer[i].atom_position.z, atom_mask, mask_offset);
            box = get_trajectory_frame(dynamic.trajectory, i).box;
            apply_pbc({cur_x, cur_y, cur_z}, num_atoms, box);

            // cur_com = compute_com(cur_x, cur_y, cur_z, mass, num_atoms);

            vec3 sum_pos = {0, 0, 0};
            float sum_mass = 0;
            for (i64 j = 0; j < num_atoms; ++j) {
                const float m = mass[j];
                sum_pos.x += cur_x[j] * m;
                sum_pos.y += cur_y[j] * m;
                sum_pos.z += cur_z[j] * m;
                sum_mass += m;
            }
            cur_com = sum_pos / sum_mass;

            const mat3 abs_mat =
                compute_weighted_cross_covariance_matrix(ref_x, ref_y, ref_z, cur_x, cur_y, cur_z, mass, num_atoms, ref_com, cur_com);
            // const mat3 rel_mat = compute_weighted_cross_covariance_matrix(prv_x, prv_y, prv_z, cur_x, cur_y, cur_z, mass, num_atoms, prv_com,
            // cur_com);
            const mat3 cov_mat = compute_weighted_covariance_matrix(cur_x, cur_y, cur_z, mass, num_atoms, cur_com);

            const mat3 abs_rot = extract_rotation(abs_mat);
            // const mat3 rel_rot = extract_rotation(rel_mat);

            quat q_abs = math::normalize(math::quat_cast(abs_rot));

            // Concatenate delta to relative transform
            // const quat q_del = math::normalize(math::quat_cast(rel_rot));
            // q_relative = math::normalize(q_relative * q_del);
            // quat q_rel = q_relative;

            // Make sure we take shortest path from previous orientation
            q_abs = math::dot(s->tracking_data.transform.rotation.absolute[i - 1], q_abs) > 0.0f ? q_abs : -q_abs;
            // q_rel = math::dot(s->tracking_data.transform.rotation.relative[i - 1], q_rel) > 0.0f ? q_rel : -q_rel;

            mat3 eigen_vectors;
            vec3 eigen_values;
            compute_eigen(cov_mat, (vec3(&)[3])eigen_vectors, (float(&)[3])eigen_values);
            {
                // @NOTE: ...

#if 0
            const auto find_max_proj = [](const vec3& ref, Array<vec3> vecs) -> int {
                int idx = 0;
                float max = -1.0f;
                for (int i = 0; i < vecs.size(); ++i) {
                    const float proj = math::abs(math::dot(ref, vecs[i]));
                    if (proj > max) {
                        max = proj;
                        idx = i;
                    }
                }
                return idx;
            };

            const mat3 pca = s->tracking_data.eigen.vectors[0];
            vec3 abs_pca[3] = {math::normalize(math::conjugate(q_abs) * pca[0]), math::normalize(math::conjugate(q_abs) * pca[1]), math::normalize(math::conjugate(q_abs) * pca[2])};
            /*

            int maj_idx = find_max_proj(eigen_vectors[0], {abs_pca, 3});
            vec3 maj = abs_pca[maj_idx];
            abs_pca[maj_idx] = abs_pca[2]; // @NOTE: Swap back
            
            int mid_idx = find_max_proj(eigen_vectors[1], {abs_pca, 2});
            vec3 mid = abs_pca[mid_idx];
            abs_pca[mid_idx] = abs_pca[1];  // @NOTE: Swap back

            vec3 min = abs_pca[0];
            */

            eigen_vectors[0] = abs_pca[0];
            eigen_vectors[1] = abs_pca[1];
            eigen_vectors[2] = abs_pca[2];
#endif
            }

#if 0
        // Full correction of relative path
        memcpy(cor_x, ref_x, num_points * sizeof(float) * 3);
        cor_com = ref_com;

        {
            const quat old_q = q_hybrid;
            const quat new_q = q_hybrid * q_del;
            q_hybrid = math::dot(old_q, new_q) < 0.0f ? -new_q : new_q;
        }

        {
            const mat4 R = math::mat4_cast(q_hybrid);
            const mat4 T1 = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {-cor_com, 1}};
            const mat4 T2 = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {cor_com, 1}};
            const mat4 M = T2 * R * T1;
            transform_positions_ref(cor_x, cor_y, cor_z, num_points, M);
        }

        const mat3 err_mat = compute_mass_weighted_cross_covariance_matrix(cor_x, cor_y, cor_z, cur_x, cur_y, cur_z, mass, num_atoms, cor_com, cur_com);
        const mat3 err_rot = extract_rotation(err_mat);
        const quat q_err = math::inverse(math::normalize(math::quat_cast(err_rot)));

        printf("err angle: %3.3f\n", math::rad_to_deg(math::angle(q_err)));

        q_hybrid = math::normalize(q_err * q_hybrid);
        const quat q_cor = q_hybrid;
#elif 0
            // Fixed ratio Slerp 80% relative, 20% absolute
            q_hybrid = math::normalize(math::slerp(q_hybrid * q_del, math::conjugate(q_abs), 0.2f));
            const quat q_cor = math::conjugate(q_hybrid);
#elif 1
            // @NOTE: Relative approach from Chevrot et al.
            // Avoids accumulative error from concatenation of matrices by storing an internal copy of the structure
            // And modifies that in an iterative fashion
            const mat3 rel_mat =
                compute_weighted_cross_covariance_matrix(int_x, int_y, int_z, cur_x, cur_y, cur_z, mass, num_atoms, int_com, cur_com);
            const mat3 rel_rot = extract_rotation(rel_mat);

            const quat q = math::normalize(math::quat_cast(rel_rot));
            const quat q_rel = math::dot(s->tracking_data.transform.rotation.relative[i - 1], q) > 0.0f ? q : -q;

            // @NOTE: Update internal representation
            for (int j = 0; j < num_atoms; j++) {
                const vec3 t = cur_com;
                const vec3 r = {cur_x[j] - t.x, cur_y[j] - t.y, cur_z[j] - t.z};

#if 1
                // Matrix multiply
                int_x[j] = rel_rot[0][0] * r.x + rel_rot[1][0] * r.y + rel_rot[2][0] * r.z + t.x;
                int_y[j] = rel_rot[0][1] * r.x + rel_rot[1][1] * r.y + rel_rot[2][1] * r.z + t.y;
                int_z[j] = rel_rot[0][2] * r.x + rel_rot[1][2] * r.y + rel_rot[2][2] * r.z + t.z;
#else
                // Quaternion multiply
                const vec3 p = q_relative * v + t;
                int_x[j] = p.x;
                int_y[j] = p.y;
                int_z[j] = p.z;
#endif
            }
            int_com = cur_com;

#elif 0
            // Partial correction based on geometric anisotropy
            // Align relative to absolute PCA axes based on shape.
            // If linear    -> Align with PCA[0]
            // If planar    -> Align with PCA[0] and PCA[1]
            // If spherical -> Align with Absolute
            const mat3 ref_pca = s->tracking_data.eigen.vectors[0];
            const mat3 cur_pca = eigen_vectors;

            const vec3 ref_ev = s->tracking_data.eigen.values[0];
            const vec3 cur_ev = eigen_values;

            const vec3 ev = cur_ev;
            const float denom = 1.0f / (ev[0] + ev[1] + ev[2]);
            const float cl = (ev[0] - ev[1]) * denom;
            const float cp = 2.0f * (ev[1] - ev[2]) * denom;
            const float cs = 3.0f * ev[2] * denom;

            quat q = q_hybrid * q_del;

            // const float wl = math::pow(1.0f - cs, 8.0f) * math::pow(1.0f - cp, 2.0f);

            float wl = cl * cp * 4.0f;
            // wl = wl * pow(bc[0], 4.0) * 16.0;

            const float w0 = wl;
            const float w1 = cp;
            const float w2 = 1.0f - w0;

            {
                quat q_pca = q;

                {
                    // Align with Major
                    const vec3 target = ref_pca[0];
                    const vec3 v_src = q_pca * cur_pca[0];
                    const vec3 v_dst = math::dot(target, v_src) > 0.0f ? target : -target;
                    const quat q_cor = math::two_direction_vectors(math::normalize(v_src), math::normalize(v_dst));
                    q_pca = math::slerp(q_pca, q_cor * q_pca, 1.0f);
                }

                {
                    // Align with Mid
                    const vec3 target = ref_pca[1];
                    const vec3 v_src = q_pca * cur_pca[1];
                    const vec3 v_dst = math::dot(target, v_src) > 0.0f ? target : -target;
                    const quat q_cor = math::two_direction_vectors(math::normalize(v_src), math::normalize(v_dst));
                    q_pca = math::slerp(q_pca, q_cor * q_pca, 1.0f);
                }

                math::normalize(q_pca);

                q = math::slerp(q, q_pca, 1.0f);
            }

            // Align with Absolute
            //{ q = math::slerp(q, q_abs, w2); }

            q = math::normalize(q);
            const quat q_hyb = math::dot(s->tracking_data.transform.rotation.hybrid[i - 1], q) > 0.0f ? q : -q;
            q_hybrid = q_hyb;
#else
            // Always go with Shortest path (absolute or prediction)
            // Result -> Bad
            const quat q_pred = q_hybrid * q_del;
            const float d_abs = math::dot(q_hybrid, q_abs);
            const float d_pred = math::dot(q_hybrid, q_pred);

            if (math::abs(d_abs) > math::abs(d_pred)) {
                const float sign = math::sign(d_abs);
                q_hybrid = math::normalize(sign * q_abs);
            } else {
                const float sign = math::sign(d_pred);
                q_hybrid = math::normalize(sign * q_pred);
            }
            const quat q_cor = q_hybrid;
#endif
            // Dynamic ratio Slerp
            // absolute contributes with a factor based on the cosine of the angle between the absolute and relative orientation
            const float d = math::dot(q_rel, q_abs);
            const float t = math::pow(d, 8.0f);
            q_hybrid = math::normalize(math::slerp(q_rel, q_abs, t));
            const quat q_hyb = math::dot(s->tracking_data.transform.rotation.hybrid[i - 1], q_hybrid) > 0.0f ? q_hybrid : -q_hybrid;
            q_hybrid = q_hyb;

            // Store results
            s->tracking_data.transform.rotation.absolute[i] = q_abs;
            s->tracking_data.transform.rotation.relative[i] = q_rel;
            s->tracking_data.transform.rotation.hybrid[i] = q_hyb;
            s->tracking_data.transform.com[i] = cur_com;
            s->tracking_data.eigen.vectors[i] = eigen_vectors;
            s->tracking_data.eigen.values[i] = eigen_values;

            // Copy previous frame data
            // memcpy(prv_x, cur_x, num_atoms * sizeof(float) * 3);
            // prv_com = cur_com;
        }
    } else if (num_atoms == 1) {
        const i64 atom_idx = bitfield::find_first_bit_set(atom_mask) + mask_offset;
        for (int i = 1; i < num_frames; i++) {
            const float x = dynamic.trajectory.frame_buffer[i].atom_position.x[atom_idx];
            const float y = dynamic.trajectory.frame_buffer[i].atom_position.y[atom_idx];
            const float z = dynamic.trajectory.frame_buffer[i].atom_position.z[atom_idx];

            s->tracking_data.transform.rotation.absolute[i] = {1, 0, 0, 0};
            s->tracking_data.transform.rotation.relative[i] = {1, 0, 0, 0};
            s->tracking_data.transform.rotation.hybrid[i] = {1, 0, 0, 0};
            s->tracking_data.transform.com[i] = {x, y, z};
            s->tracking_data.eigen.vectors[i] = mat3(1);
            s->tracking_data.eigen.values[i] = {1, 1, 1};
        }
    }

    const auto ext = get_trajectory_frame(dynamic.trajectory, 0).box * vec3(1, 1, 1);
    int i[3] = {0, 1, 2};

    const auto swap = [](int& x, int& y) {
        int tmp = x;
        x = y;
        y = tmp;
    };

    const mat3 I = mat3(1);
    if (ext[i[0]] < ext[i[1]]) swap(i[0], i[1]);
    if (ext[i[1]] < ext[i[2]]) swap(i[1], i[2]);
    if (ext[i[0]] < ext[i[1]]) swap(i[0], i[1]);

    const mat3 M = mat3(I[i[0]], I[i[1]], I[i[2]]);
    s->tracking_data.simulation_box_aligned_pca = M * math::transpose(s->tracking_data.eigen.vectors[0]);
    s->tracking_data.pca = s->tracking_data.eigen.vectors[0];

    return true;
}  // namespace structure_tracking

void transform_to_internal_frame(MoleculeDynamic& dynamic, Bitfield atom_mask, i64 mask_offset) {
    // Scratch data
    const int64_t num_atoms = bitfield::number_of_bits_set(atom_mask);
    const int64_t num_frames = dynamic.trajectory.num_frames;

    const auto mem_size = sizeof(float) * num_atoms * 7;
    void* mem = TMP_MALLOC(mem_size);
    defer { TMP_FREE(mem); };
    memset(mem, 0, mem_size);

    float* cur_x = (float*)mem + 0 * num_atoms;
    float* cur_y = (float*)mem + 1 * num_atoms;
    float* cur_z = (float*)mem + 2 * num_atoms;
    float* int_x = (float*)mem + 3 * num_atoms;
    float* int_y = (float*)mem + 4 * num_atoms;
    float* int_z = (float*)mem + 5 * num_atoms;
    float* mass = (float*)mem + 6 * num_atoms;

    bitfield::gather_masked(int_x, get_trajectory_position_x(dynamic.trajectory, 0).data(), atom_mask, mask_offset);
    bitfield::gather_masked(int_y, get_trajectory_position_y(dynamic.trajectory, 0).data(), atom_mask, mask_offset);
    bitfield::gather_masked(int_z, get_trajectory_position_z(dynamic.trajectory, 0).data(), atom_mask, mask_offset);
    bitfield::gather_masked(mass, dynamic.molecule.atom.mass, atom_mask, mask_offset);

    vec3 int_com = compute_com(int_x, int_y, int_z, num_atoms);

    for (i32 i = 1; i < num_frames; i++) {
        float* pos_x = dynamic.trajectory.frame_buffer[i].atom_position.x;
        float* pos_y = dynamic.trajectory.frame_buffer[i].atom_position.y;
        float* pos_z = dynamic.trajectory.frame_buffer[i].atom_position.z;

        // Fetch current
        bitfield::gather_masked(cur_x, pos_x, atom_mask, mask_offset);
        bitfield::gather_masked(cur_y, pos_y, atom_mask, mask_offset);
        bitfield::gather_masked(cur_z, pos_z, atom_mask, mask_offset);
        const vec3 cur_com = compute_com(cur_x, cur_y, cur_z, mass, num_atoms);

        const mat3 m1 = compute_weighted_cross_covariance_matrix(int_x, int_y, int_z, cur_x, cur_y, cur_z, mass, num_atoms, int_com, cur_com);
        mat3 r1 = extract_rotation(m1);

        for (i64 j = 0; j < num_atoms; j++) {
            const vec3 t = cur_com;
            const vec3 v = {cur_x[j] - t.x, cur_y[j] - t.y, cur_z[j] - t.z};
            int_x[j] = r1[0][0] * v.x + r1[1][0] * v.y + r1[2][0] * v.z + t.x;
            int_y[j] = r1[0][1] * v.x + r1[1][1] * v.y + r1[2][1] * v.z + t.y;
            int_z[j] = r1[0][2] * v.x + r1[1][2] * v.y + r1[2][2] * v.z + t.z;
        }
        int_com = cur_com;

        bitfield::scatter_masked(pos_x, int_x, atom_mask);
        bitfield::scatter_masked(pos_y, int_y, atom_mask);
        bitfield::scatter_masked(pos_z, int_z, atom_mask);
    }
}

const TrackingData* get_tracking_data(ID id) {
    ASSERT(context);
    if (auto* s = find_structure(id)) {
        return &s->tracking_data;
    }
    LOG_ERROR("Supplied id is not valid.");
    return nullptr;
}

}  // namespace structure_tracking