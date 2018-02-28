#pragma once

#include <mol/molecule.h>

// Interface to get Pdb data
struct PdbInterface : MoleculeInterface {
    Array<char> alt_loc;
    Array<float> occupancy;
    Array<float> temp_factor;
    Array<char> charge;
};

struct PdbData : MoleculeData {
    DynamicArray<char> alt_loc_data;
    DynamicArray<float> occupancy_data;
    DynamicArray<float> temp_factor_data;
    DynamicArray<char> charge_data;

	/*
    operator PdbInterface() {
        return {atom_position_data, atom_element_data, atom_label_data, atom_residue_index_data, bond_data,  residue_data,
                chain_data,         alt_loc_data,      occupancy_data,  temp_factor_data,        charge_data};
    }
	*/
};