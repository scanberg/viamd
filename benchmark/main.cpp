#include <core/types.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <core/camera.h>
#include <core/camera_utils.h>
#include <core/string_utils.h>
#include <core/bitfield.h>
#include <core/spatial_hash.h>

#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/filter.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>

#include <glm/gtx/io.hpp>

#include <chrono>
#define TIME() std::chrono::high_resolution_clock::now()
#define NANOSEC(x, y) std::chrono::duration<double, std::nano>(y - x).count()
#define MILLISEC(x, y) std::chrono::duration<double, std::milli>(y - x).count()
#define SECONDS(x, y) std::chrono::duration<double>(y - x).count()

#define SECTION(x) printf("\n########   %s   ########\n", x);

#define DATASET VIAMD_DATA_DIR"/alanine/two4REP-CH3_450K.pdb"

extern void cubic_interpolation_pbc_scalar(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z,
										   const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                                           const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1,
										   const float* RESTRICT in_x2, const float* RESTRICT in_y2, const float* RESTRICT in_z2,
										   const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3,
										   int64 count, float t, const mat3& sim_box);

extern void cubic_interpolation_pbc_128(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                                        const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1, const float* RESTRICT in_x2, const float* RESTRICT in_y2,
                                        const float* RESTRICT in_z2, const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3, int64 count, float t, const mat3& sim_box);

extern void cubic_interpolation_pbc_256(float* RESTRICT out_x, float* RESTRICT out_y, float* RESTRICT out_z, const float* RESTRICT in_x0, const float* RESTRICT in_y0, const float* RESTRICT in_z0,
                                           const float* RESTRICT in_x1, const float* RESTRICT in_y1, const float* RESTRICT in_z1, const float* RESTRICT in_x2, const float* RESTRICT in_y2,
                                           const float* RESTRICT in_z2, const float* RESTRICT in_x3, const float* RESTRICT in_y3, const float* RESTRICT in_z3, int64 count, float t,
                                           const mat3& sim_box);

int main() {
	MoleculeDynamic md;
	
	SECTION("Load Data")
	const auto t0_load = TIME();
	bool res_mol = pdb::load_molecule_from_file(&md.molecule, DATASET);
	const auto t1_load = TIME();
	bool res_traj = pdb::load_trajectory_from_file(&md.trajectory, DATASET);
	const auto t2_load = TIME();

	if (!res_mol || !res_traj) {
		printf("Could not load dataset: %s\n", DATASET);
		return 1;
	}

	const int32 num_atoms = (int32)md.molecule.atom.count;
	const int32 num_frames = (int32)md.trajectory.num_frames;

	printf("Dataset loaded successfully: %s\nnum_atoms: %i\nnum_frames: %i\n", DATASET, num_atoms, num_frames);
	printf("Time to load molecule: %.2fms\n", MILLISEC(t0_load, t1_load));
	printf("Time to load trajectory: %.2fms\n", MILLISEC(t1_load, t2_load));
	printf("Loaded trajectory size in RAM: %.2fMB\n", (double)(num_atoms * num_frames * sizeof(float) * 3) / MEGABYTES(1));

	// AOS layout for comparison
	vec3* xyz = (vec3*)TMP_MALLOC(num_atoms * num_frames * sizeof(vec3));
	defer{ TMP_FREE(xyz); };
	for (int32 i = 0; i < num_frames * num_atoms; i++) {
		xyz[i] = { md.trajectory.position_data.x[i], md.trajectory.position_data.y[i], md.trajectory.position_data.z[i] };
	}

	SECTION("MEMSET REFERENCE")
	{
		const auto num_iter = 100;
		const auto t0 = TIME();

		const auto x = md.trajectory.position_data.x;
		const auto y = md.trajectory.position_data.y;
		const auto z = md.trajectory.position_data.z;
		const auto size = num_atoms * num_frames;
		const auto size_in_bytes = size * sizeof(float) * 3 * num_iter;

		for (int32 i = 0; i < num_iter; i++) {
			memset(x, 0, size * sizeof(float));
			memset(y, 0, size * sizeof(float));
			memset(z, 0, size * sizeof(float));
		}
		const auto t1 = TIME();
		const auto time_naive = MILLISEC(t0, t1) / (double)num_iter;
		const auto throughput_naive = size_in_bytes / SECONDS(t0, t1);
		printf("Time: %.2fms\n", time_naive);
		printf("Throughput: %.2fMB/s\n", throughput_naive / (double)MEGABYTES(1));
	}

	SECTION("CUBIC INTERPOLATION PBC");
	{
		const auto x0 = get_trajectory_position_x(md.trajectory, 100).data();
		const auto y0 = get_trajectory_position_y(md.trajectory, 100).data();
		const auto z0 = get_trajectory_position_z(md.trajectory, 100).data();
		const auto x1 = get_trajectory_position_x(md.trajectory, 101).data();
		const auto y1 = get_trajectory_position_y(md.trajectory, 101).data();
		const auto z1 = get_trajectory_position_z(md.trajectory, 101).data();
		const auto x2 = get_trajectory_position_x(md.trajectory, 102).data();
		const auto y2 = get_trajectory_position_y(md.trajectory, 102).data();
		const auto z2 = get_trajectory_position_z(md.trajectory, 102).data();
		const auto x3 = get_trajectory_position_x(md.trajectory, 103).data();
		const auto y3 = get_trajectory_position_y(md.trajectory, 103).data();
		const auto z3 = get_trajectory_position_z(md.trajectory, 103).data();

		auto x = md.molecule.atom.position.x;
		auto y = md.molecule.atom.position.y;
		auto z = md.molecule.atom.position.z;

		const auto box = get_trajectory_frame(md.trajectory, 101).box;

		const int32 num_iter = 10000;
		const float t = 0.5f;
		const auto size_in_bytes = num_atoms * 12 * sizeof(float) * num_iter;

		const auto t0 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			cubic_interpolation_pbc_scalar(x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, num_atoms, t, box);
		}
		const auto t1 = TIME();
		const auto time_cubic_pbc_scalar = MILLISEC(t0, t1) / (double)num_iter;
		const auto throughput_scalar = size_in_bytes / SECONDS(t0, t1);
		printf("Time (scalar): %.2fms\n", time_cubic_pbc_scalar);
		printf("Throughput (scalar): %.2fMB/s\n", throughput_scalar / MEGABYTES(1));

		const auto t2 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			cubic_interpolation_pbc_128(x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, num_atoms, t, box);
		}
		const auto t3 = TIME();
		const auto time_cubic_pbc_128 = MILLISEC(t2, t3) / (double)num_iter;
		const auto throughput_128 = size_in_bytes / SECONDS(t2, t3);
		printf("Time (128-bit): %.2fms (%.1fx) speedup\n", time_cubic_pbc_128, time_cubic_pbc_scalar / time_cubic_pbc_128);
		printf("Throughput (128-bit): %.2fMB/s\n", throughput_128 / MEGABYTES(1));


#ifdef __AVX__
		const auto t4 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			cubic_interpolation_pbc_256(x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, num_atoms, t, box);
		}
		const auto t5 = TIME();
		const auto time_cubic_pbc_256 = MILLISEC(t4, t5) / (double)num_iter;
		const auto throughput_256 = size_in_bytes / SECONDS(t4, t5);
		printf("Time (256-bit): %.2fms (%.1fx) speedup\n", time_cubic_pbc_256, time_cubic_pbc_scalar / time_cubic_pbc_256);
		printf("Throughput (256-bit): %.2fMB/s\n", throughput_256 / MEGABYTES(1));

#endif
	}

	SECTION("TRANSLATE")
	{
		const auto num_iter = 100;

		const auto x = md.trajectory.position_data.x;
		const auto y = md.trajectory.position_data.y;
		const auto z = md.trajectory.position_data.z;
		const auto size = num_atoms * num_frames;
		const auto size_in_bytes = size * sizeof(float) * 3 * num_iter;

		const vec3 t = { 1,2,3 };

		const auto t0 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			for (int32 j = 0; j < size; j++) {
				x[j] += t.x;
				y[j] += t.y;
				z[j] += t.z;
			}
		}
		const auto t1 = TIME();
		const auto time_naive = MILLISEC(t0, t1) / (double)num_iter;
		const auto throughput_naive = size_in_bytes / SECONDS(t0, t1);
		printf("Time (naive): %.2fms\n", time_naive);
		printf("Throughput (naive): %.2fMB/s\n", throughput_naive / (double)MEGABYTES(1));

		const auto t2 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			translate_positions(x, y, z, size, t);
		}
		const auto t3 = TIME();
		const auto time_vec = MILLISEC(t2, t3) / (double)num_iter;
		const auto throughput_vec = size_in_bytes / SECONDS(t2, t3);

		printf("Time (vectorized): %.2fms (%.1fx) speedup\n", time_vec, time_naive / time_vec);
		printf("Throughput (vectorized): %.2fMB/s\n", throughput_vec / (double)MEGABYTES(1));
	}

	SECTION("MATRIX TRANSFORM")
	{
		const auto num_iter = 100;
		mat4 M = mat4(1);
		
		const auto x = md.trajectory.position_data.x;
		const auto y = md.trajectory.position_data.y;
		const auto z = md.trajectory.position_data.z;
		const auto size = num_atoms * num_frames;
		const auto size_in_bytes = size * sizeof(float) * 3 * num_iter;

		const auto t0 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			for (int32 j = 0; j < size; j++) {
				const float v_x = x[j];
				const float v_y = y[j];
				const float v_z = z[j];
				const float v_w = 1.0f;

				x[j] = v_x * M[0][0] + v_y * M[1][0] + v_z * M[2][0] + v_w * M[3][0];
				y[j] = v_x * M[0][1] + v_y * M[1][1] + v_z * M[2][1] + v_w * M[3][1];
				z[j] = v_x * M[0][2] + v_y * M[1][2] + v_z * M[2][2] + v_w * M[3][2];
			}
		}
		const auto t1 = TIME();
		const auto time_naive = MILLISEC(t0, t1) / (double)num_iter;
		const auto throughput_naive = size_in_bytes / SECONDS(t0, t1);
		printf("Time (naive): %.2fms\n", time_naive);
		printf("Throughput (naive): %.2fMB/s\n", throughput_naive / (double)MEGABYTES(1));

		const auto t2 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			transform_positions(x, y, z, size, M);
		}
		const auto t3 = TIME();
		const auto time_vec = MILLISEC(t2, t3) / (double)num_iter;
		const auto throughput_vec = size_in_bytes / SECONDS(t2, t3);

		printf("Time (vectorized): %.2fms (%.1fx) speedup\n", time_vec, time_naive / time_vec);
		printf("Throughput (vectorized): %.2fMB/s\n", throughput_vec / (double)MEGABYTES(1));

		const auto t4 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			for (int32 j = 0; j < size; j++) {
				vec4 v = { xyz[i], 1.0f };
				v = M * v;
				xyz[i] = v;
			}
		}
		const auto t5 = TIME();
		const auto time_aos = MILLISEC(t4, t5) / (double)num_iter;
		const auto throughput_aos = size_in_bytes / SECONDS(t4, t5);

		printf("Time (aos): %.2fms (%.1fx) speedup\n", time_aos, time_naive / time_aos);
		printf("Throughput (aos): %.2fMB/s\n", throughput_aos / (double)MEGABYTES(1));
	}

	SECTION("COM")
	{
		const auto num_iter = 100;

		const auto x = md.trajectory.position_data.x;
		const auto y = md.trajectory.position_data.y;
		const auto z = md.trajectory.position_data.z;
		const auto size = num_atoms * num_frames;
		const auto size_in_bytes = size * sizeof(float) * 3 * num_iter;

		const auto t0 = TIME();
		float com_x = 0.0f;
		float com_y = 0.0f;
		float com_z = 0.0f;
		for (int32 i = 0; i < num_iter; i++) {
			com_x = 0.0f;
			com_y = 0.0f;
			com_z = 0.0f;
			for (int32 j = 0; j < size; j++) {
				com_x += x[j];
				com_y += y[j];
				com_z += z[j];
			}
			com_x /= (float)size;
			com_y /= (float)size;
			com_z /= (float)size;
		}
		printf("com (%.2f, %.2f, %.2f)\n", com_x, com_y, com_z);

		const auto t1 = TIME();
		const auto time_naive = MILLISEC(t0, t1) / (double)num_iter;
		const auto throughput_naive = size_in_bytes / SECONDS(t0, t1);
		printf("Time (naive): %.2fms\n", time_naive);
		printf("Throughput (naive): %.2fMB/s\n", throughput_naive / (double)MEGABYTES(1));

		const auto t2 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			vec3 com = compute_com(x, y, z, size);
		}
		const auto t3 = TIME();
		const auto time_vec = MILLISEC(t2, t3) / (double)num_iter;
		const auto throughput_vec = size_in_bytes / SECONDS(t2, t3);

		printf("Time (vectorized): %.2fms (%.1fx) speedup\n", time_vec, time_naive / time_vec);
		printf("Throughput (vectorized): %.2fMB/s\n", throughput_vec / (double)MEGABYTES(1));
	}

	SECTION("BOUNDING BOX")
	{
		const auto num_iter = 100;

		const auto x = md.trajectory.position_data.x;
		const auto y = md.trajectory.position_data.y;
		const auto z = md.trajectory.position_data.z;
		const auto size = num_atoms * num_frames;

		const auto size_in_bytes = size * sizeof(float) * 3 * num_iter;

		const auto t0 = TIME();
		vec3 min, max;
		for (int32 i = 0; i < num_iter; i++) {
			min = { x[0], y[0], z[0] };
			max = { x[0], y[0], z[0] };
			for (int32 j = 1; j < size; j++) {
				const float v_x = x[j];
				const float v_y = y[j];
				const float v_z = z[j];

				min.x = v_x < min.x ? v_x : min.x;
				min.y = v_y < min.y ? v_y : min.y;
				min.z = v_z < min.z ? v_z : min.z;
				max.x = v_x > max.x ? v_x : max.x;
				max.y = v_y > max.y ? v_y : max.y;
				max.z = v_z > max.z ? v_z : max.z;
			}
		}
		printf("min (%.2f, %.2f, %.2f)\n", min.x, min.y, min.z);
		printf("max (%.2f, %.2f, %.2f)\n", max.x, max.y, max.z);

		const auto t1 = TIME();
		const auto time_naive = MILLISEC(t0, t1) / (double)num_iter;
		const auto throughput_naive = size_in_bytes / SECONDS(t0, t1);
		printf("Time (naive): %.2fms\n", time_naive);
		printf("Throughput (naive): %.2fMB/s\n", throughput_naive / (double)MEGABYTES(1));

		const auto t2 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			vec3 min_box, max_box;
			compute_bounding_box(min_box, max_box, x, y, z, size);
		}
		const auto t3 = TIME();
		const auto time_vec = MILLISEC(t2, t3) / (double)num_iter;
		const auto throughput_vec = size_in_bytes / SECONDS(t2, t3);

		printf("Time (vectorized): %.2fms (%.1fx) speedup\n", time_vec, time_naive / time_vec);
		printf("Throughput (vectorized): %.2fMB/s\n", throughput_vec / (double)MEGABYTES(1));
	}

	return 0;
}