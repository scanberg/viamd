#pragma once

#include <mol/pdb.h>
#include <mol/trajectory.h>
#include <core/common.h>
#include <core/string_utils.h>

struct PdbResult {
	PdbData* mol;
	Trajectory* traj;
	int num_models;
};

enum PdbLoadParams {
	PDB_READ_ATOM = BIT(0),
	PDB_READ_HETATM = BIT(1),
	PDB_TREAT_MODELS_AS_FRAMES = BIT(2),
	PDB_DEFAULT = 0xFFFFFFFF
};

PdbResult load_pdb_from_file(const char* filename, PdbLoadParams params = PDB_DEFAULT);
PdbResult parse_pdb_from_string(CString string, PdbLoadParams params = PDB_DEFAULT);