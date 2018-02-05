#pragma once

#include <mol/pdb.h>
#include <core/common.h>
#include <core/array.h>

struct PdbResult {
	bool success = false;
	PdbStructure pdb = {};

	operator bool() { return success; }
}

enum PdbLoadParams {
	PDB_IGNORE_HETATM = BIT(1);
	PDB_DEFAULT = PDB_IGNORE_HETATM;
}

PdbResult load_pdb_from_file(const char* filename, PdbLoadParams params = PDB_DEFAULT);
PdbResult parse_pdb_from_string(String string, PdbLoadParams params = PDB_DEFAULT);