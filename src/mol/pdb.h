#pragma once

#include <mol/molecule.h>

// Interface to get Pdb data
struct PdbStructure : MoleculeStructure {
	Array<char>    alt_loc;
	Array<float>   occupancy;
	Array<float>   temp_factor;
	Array<char>	   charge;
};

