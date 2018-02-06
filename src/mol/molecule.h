#pragma once

#include <core/types.h>
#include <core/array.h>
#include <mol/element.h>
#include <string.h>

#ifdef WIN32
#pragma warning(disable:4996) // strncpy instead of strncpy_s (windows)
#endif

struct Bond {
	int32 atom_idx_a;
	int32 atom_idx_b;
};

struct Label {
	static constexpr int MAX_LENGTH = 8;
	Label() = default;

	template<int32 N>
	Label(const char (&cstr)[N]) {
		int32 len = N < MAX_LENGTH ? N : MAX_LENGTH - 1;
		strncpy(data, cstr.data, len);
		length = len;
		data[length] = '\0';
	}

	Label(const char* cstr) {
		// TODO: IS THIS ENOUGH?
		int32 len = (int)strnlen(cstr, MAX_LENGTH - 1);
		strncpy(data, cstr, len);
		length = len;
		data[length] = '\0';
	}

	Label(const Label& other) {
		ASSERT(other.length < MAX_LENGTH);
		memcpy(data, other.data, MAX_LENGTH);
		length = other.length;
		data[length] = '\0';
	}

	Label(const CString& cstr) {
		int32 len = cstr.count < MAX_LENGTH ? (int)cstr.count : MAX_LENGTH - 1;
		strncpy(data, cstr.data, len);
		length = len;
		data[length] = '\0';
	}

	Label& operator =(const Label& other) {
		if (this != &other) {
			ASSERT(other.length < MAX_LENGTH);
			memcpy(data, other.data, MAX_LENGTH);
			length = other.length;
		}
		return *this;
	}

	Label& operator =(const CString& cstr) {
		int32 len = cstr.count < MAX_LENGTH ? (int)cstr.count : MAX_LENGTH - 1;
		strncpy(data, cstr.data, len);
		length = len;
		data[length] = '\0';
		return *this;
	}

	template<int32 N>
	Label& operator =(const char (&cstr)[N]) {
		if (data != cstr) {
			int32 len = N < MAX_LENGTH ? N : MAX_LENGTH - 1;
			strncpy(data, cstr.data, len);
			length = len;
			data[length] = '\0';
		}
		return *this;
	}

	operator String() { return String(data, length); }
	operator CString() const { return CString(data, length); }
	operator const char*() const { return data; }
	char* begin() { return data; }
	char* beg() { return data; }
	char* end() { return data + length; }

	char data[MAX_LENGTH] = {};
	int32 length = 0;
};

struct Residue {
	Label id;
	int32 beg_atom_idx;
	int32 end_atom_idx;
};

struct Chain {
	char id;
	int32 beg_res_idx;
	int32 end_res_idx;
};

struct Backbone {
	Array<int32> atom_ca_idx;
	// Use C -> O vector as orientation vector for ribbons
	Array<int32> atom_o_idx;
	Array<int32> atom_c_idx;
};

// Interface to access molecular data
struct MoleculeStructure {
	Array<vec3>		atom_positions;
	Array<Element>	atom_elements;
	Array<Label>	atom_labels;
	Array<int32>	atom_residue_indices;

	Array<Bond>		bonds;
	Array<Residue>	residues;
	Array<Chain>	chains;
};

// This is the actual owner of the data which is exposed through MoleculeStructure Interface.
struct MoleculeData : MoleculeStructure {
	void* data_block = nullptr;
};

enum MoleculeStructureAllocationFlags {
	MOL_POSITIONS = BIT(0),
	MOL_ELEMENTS = BIT(1),
	MOL_LABELS = BIT(2),
	MOL_RESIDUE_INDICES = BIT(3),
	MOL_ALL = 0xffffffff
};

MoleculeStructure allocate_molecule_structure(int atom_count, int bond_count, int residue_count, int chain_count, MoleculeStructureAllocationFlags alloc_flags = MOL_ALL, Allocator& allocator = default_alloc);

// This is a bit risky
void free_molecule_structure(MoleculeStructure& mol, Allocator& allocator = default_alloc);