#pragma once

#include <core/types.h>
#include <core/array.h>
#include <string.h>

using Element = char;

struct Label {
	static constexpr int MaxLength = 16;
	Label() = default;

	template<int32 N>
	Label(const char (&cstr)[N]) {
		int32 len = N < MaxLength ? N, MaxLength - 1;
		strncpy_s(data, other.data, len);
	}

	Label(const Label& other) {
		strncpy_s(data, other.data, other.length);
	}

	char data[MaxLength] = {};
	int32 length = 0;
};

struct Residue {
	Label id;
	int32 beg_atom_idx;
	int32 end_atom_idx;
};

struct Chain {
	Label id;
	int32 beg_res_idx;
	int32 end_res_idx;
};

struct Molecule {
	Array<vec3>		atom_position;
	Array<Element>	atom_element;
	Array<Label>	atom_label;
	Array<int32>	atom_residue_index;

	Array<Residue>	residues;
	Array<Chain>	chains;
};