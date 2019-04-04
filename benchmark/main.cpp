#include <core/gl.h>
#include <core/types.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/math_utils.h>
#include <core/camera.h>
#include <core/camera_utils.h>
#include <core/string_utils.h>
#include <core/volume.h>
#include <core/bitfield.h>

#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/trajectory_utils.h>
#include <mol/molecule_utils.h>
#include <mol/hydrogen_bond.h>
#include <mol/filter.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>
#include <mol/spatial_hash.h>

#include <glm/gtx/io.hpp>

#include <chrono>
#define TIME() std::chrono::high_resolution_clock::now()
#define MILLISEC(x, y) std::chrono::duration<double, std::milli>(y - x).count()

#define SECTION(x) printf("\n########   %s   ########\n", x);

#define DATASET VIAMD_DATA_DIR"/alanine/two4REP-CH3_450K.pdb"

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

	printf("Dataset loaded successfully: %s\n -num_atoms: %i\n -num_frames: %i\n", DATASET, (int32)md.molecule.atom.count, (int32)md.trajectory.num_frames);
	printf("Time to load molecule: %.2fms\n", MILLISEC(t0_load, t1_load));
	printf("Time to load trajectory: %.2fms\n", MILLISEC(t1_load, t2_load));

	SECTION("Interpolation");
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

		const int32 num_atoms = (int32)md.molecule.atom.count;
		const int32 num_iter = 1000;
		const float t = 0.5f;
		const auto t0 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			cubic_interpolation_pbc_scalar(x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, num_atoms, t, box);
		}
		const auto t1 = TIME();
		const auto time_cubic_pbc_scalar = MILLISEC(t0, t1) / (double)num_iter;
		printf("Time to interpolate cubic pbc (scalar): %.2fms\n", time_cubic_pbc_scalar);

		const auto t2 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			cubic_interpolation_pbc(x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, num_atoms, t, box);
		}
		const auto t3 = TIME();
		const auto time_cubic_pbc_128 = MILLISEC(t2, t3) / (double)num_iter;
		printf("Time to interpolate cubic pbc (128 wide): %.2fms (%.1fx) speedup\n", time_cubic_pbc_128, time_cubic_pbc_scalar / time_cubic_pbc_128);

#ifdef __AVX__
		const auto t4 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			cubic_interpolation_pbc_256(x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, num_atoms, t, box);
		}
		const auto t5 = TIME();
		const auto time_cubic_pbc_256 = MILLISEC(t4, t5) / (double)num_iter;
		printf("Time to interpolate cubic pbc (256 wide): %.2fms (%.1fx) speedup\n", time_cubic_pbc_256, time_cubic_pbc_scalar / time_cubic_pbc_256);
#endif
	}

	SECTION("MATRIX TRANSFORM")
	{
		const auto num_iter = 1000;
		mat4 M = mat4(1);
		
		const auto x = md.trajectory.position_data.x;
		const auto y = md.trajectory.position_data.y;
		const auto z = md.trajectory.position_data.z;
		const auto size = md.trajectory.num_atoms * md.trajectory.num_frames;

		const auto t0 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			for (int32 j = 0; j < size; j++) {
				vec4 v = { x[j], y[j], z[j], 1.0f };
				v = M * v;
				x[j] = v.x;
				y[j] = v.y;
				z[j] = v.z;
			}
		}
		const auto t1 = TIME();
		const auto time_glm = MILLISEC(t0, t1) / (double)num_iter;
		printf("Time to transform points (glm): %.2fms\n", time_glm);

		const auto t2 = TIME();
		for (int32 i = 0; i < num_iter; i++) {
			transform_positions(x, y, z, size, M);
		}
		const auto t3 = TIME();
		const auto time_vec = MILLISEC(t2, t3) / (double)num_iter;
		printf("Time to transform points (vectorized): %.2fms (%.1fx) speedup\n", time_vec, time_glm / time_vec);

	}

	return 0;
}