#pragma once

#include <core/types.h>
#include <core/bitfield.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/hash.h>

struct MoleculeDynamic;

namespace structure_tracking {
typedef uint32 ID;

void initialize();
void shutdown();

ID create_structure();
bool remove_structure(ID structure_id);
void clear_structures();

struct TrackingData {
    int64 count = 0;
    struct {
        vec3* com = nullptr;
        struct {
            quat* absolute = nullptr;
            quat* relative = nullptr;
            quat* hybrid = nullptr;
        } rotation;
    } transform;
    struct {
        mat3* vector = nullptr;
        vec3* value = nullptr;
    } eigen;
};

/*
vec3 compute_eigen_values(const float* RESTRICT x, const float* RESTRICT y, const float* RESTRICT z, const float* RESTRICT mass, int64 count);

mat3 compute_rotation(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
                                          const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
                                          const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1);
*/

bool compute_trajectory_transform_data(ID structure_id, Bitfield atom_mask, const MoleculeDynamic& dynamic);

const TrackingData* get_tracking_data(ID structure_id);

/*
const Array<const vec3> get_com(ID structure_id);
const Array<const quat> get_rot_absolute(ID structure_id);
const Array<const quat> get_rot_relative(ID structure_id);
const Array<const quat> get_rot_hybrid(ID structure_id);
*/
/*
const Array<const mat3> get_eigen_vectors(ID structure_id);
const Array<const vec3> get_eigen_values(ID structure_id);
*/

// const Array<const SupportFrame> get_support_frames(ID structure_id);

}  // namespace structure_tracking
