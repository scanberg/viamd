#pragma once

#include <mol/pdb.h>
#include <core/common.h>

struct PdbResult {
	bool success = false;
	PdbStructure pdb = {};
}

enum PdbLoadParams {
	PDB_IGNORE_HETATM = BIT(1);
	PDB_DEFAULT = PDB_IGNORE_HETATM;
}

void load_pdb_from_file(const char* filename, PdbLoadParams params = PDB_IGNORE_HETATM);