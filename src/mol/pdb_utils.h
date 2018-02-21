#pragma once

#include <mol/pdb.h>
#include <core/common.h>
#include <core/allocator.h>
#include <core/array.h>

struct PdbResult {
	bool success = false;
	PdbStructure pdb = {};

	operator bool() { return success; }
};

enum PdbLoadParams {
	PDB_READ_ATOM = BIT(1),
	PDB_READ_HETATM = BIT(2),
	PDB_DEFAULT = PDB_READ_ATOM
};

PdbResult load_pdb_from_file(const char* filename, PdbLoadParams params = PDB_DEFAULT, Allocator* alloc = &default_alloc);
PdbResult parse_pdb_from_string(CString string, PdbLoadParams params = PDB_DEFAULT, Allocator* alloc = &default_alloc);