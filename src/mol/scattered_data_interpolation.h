#include <core/types.h>
#include <core/array.h>

struct RadialBasis {
    DynamicArray<vec3> control_points = {};
    DynamicArray<vec3> weights = {};
};

void compute_displacement_basis(RadialBasis* dst, Array<const vec3> orig_points, Array<const vec3> disp_points);
RadialBasis compute_displacement_basis(Array<const vec3> orig_points, Array<const vec3> disp_points);

vec3 evaluate_displacement(const RadialBasis& basis, const vec3& x);
