#pragma once

#include <core/types.h>
#include <core/bitfield.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/hash.h>

struct MoleculeDynamic;

struct Transform {
    mat3 rotation = {};
	vec3 com = {};
};

/*
typedef int SupportAxis;

enum {
	SupportAxis_PosX,
	SupportAxis_NegX,
	SupportAxis_PosY,
	SupportAxis_NegY,
	SupportAxis_PosZ,
	SupportAxis_NegZ
};

struct SupportFrame {
    struct {
        vec3 dir;
        float mass;
    } axis[6];
};
*/

namespace structure_tracking {
typedef uint32 ID;

void initialize();
void shutdown();

ID create_structure();
bool remove_structure(ID structure_id);

void clear_structures();

vec3 compute_eigen_values(const float* RESTRICT x, const float* RESTRICT y, const float* RESTRICT z, const float* RESTRICT mass, int64 count);

mat3 compute_rotation(const float* RESTRICT x0, const float* RESTRICT y0, const float* RESTRICT z0,
					  const float* RESTRICT x1, const float* RESTRICT y1, const float* RESTRICT z1,
					  const float* RESTRICT mass, int64 count, const vec3& com0, const vec3& com1);

bool compute_trajectory_transform_data(ID structure_id, Bitfield atom_mask, const MoleculeDynamic& dynamic, int32 target_frame_idx = 0);

const Transform& get_transform_to_target_frame(ID structure_id, int32 source_frame);

const Array<const mat3> get_eigen_vectors(ID structure_id);
const Array<const vec3> get_eigen_values(ID structure_id);

const Array<const float> get_abs_det(ID structure_id);
const Array<const float> get_rel_det(ID structure_id);

//const Array<const SupportFrame> get_support_frames(ID structure_id);

}  // namespace structure_tracking
