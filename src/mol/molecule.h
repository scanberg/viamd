#pragma once

#include <core/types.h>
#include <core/array.h>
#include <string.h>

using Element = char;

struct Bond {
	int32 atom_idx_a;
	int32 atom_idx_b;
};

struct Label {
	static constexpr int MaxLength = 8;
	Label() = default;

	template<int32 N>
	Label(const char (&cstr)[N]) {
		int32 len = N < MaxLength ? N : MaxLength - 1;
		strncpy(data, cstr.data, len);
	}

	Label(const Label& other) {
		ASSERT(other.length < MaxLength);
		strncpy(data, other.data, other.length);
		length = other.length;
	}

	Label& operator =(const Label& other) {
		if (this != &other) {
			ASSERT(other.length < MaxLength);
			strncpy(data, other.data, other.length);
			length = other.length;
		}
		return *this;
	}

	template<int32 N>
	Label& operator =(const char (&cstr)[N]) {
		if (data != cstr) {
			int32 len = N < MaxLength ? N : MaxLength - 1;
			strncpy(data, cstr.data, len);
		}
		return *this;
	}

	operator const char*() const { return data; }
	char* begin() { return data; }
	char* beg() { return data; }
	char* end() { return data + length; }

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

struct Backbone {
	Array<int32> atom_ca_idx;
	// TODO: SHOULD ONE USE HA FOR ANCHOR DIRECTION VECTOR???
	Array<int32> atom_ha_idx;
};

struct MoleculeStructure {
	Array<vec3>		atom_positions;
	Array<Element>	atom_elements;
	Array<Label>	atom_labels;
	Array<Bond>		atom_bonds;
	Array<int32>	atom_residue_indices;

	Array<Residue>	residues;
	Array<Chain>	chains;
};