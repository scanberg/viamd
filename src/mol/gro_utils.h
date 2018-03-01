#pragma once

#include <mol/molecule.h>
#include <core/string_utils.h>

MoleculeStructure* allocate_and_load_gro_from_file(const char* filename);
MoleculeStructure* allocate_and_parse_gro_from_string(CString string);