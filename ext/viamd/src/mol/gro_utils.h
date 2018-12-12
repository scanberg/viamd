#pragma once

#include <mol/molecule_structure.h>
#include <core/string_utils.h>

bool allocate_and_load_gro_from_file(MoleculeStructure* mol, const char* filename);
bool allocate_and_parse_gro_from_string(MoleculeStructure* mol, CString string);