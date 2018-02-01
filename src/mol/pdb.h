#pragma once

#include <mol/molecule.h>

struct PdbStructure : MoleculeStructure {
	Array<char>    alt_loc;
	Array<float>   occupancy;
	Array<float>   temp_factor;
	Array<char>	   charge;

	Array<MoleculeStructure> models;
};

