#pragma once

#include <core/types.h>
#include <core/string_utils.h>


struct MoleculeDynamic;
struct MoleculeStructure;
struct Residue;

namespace stats {

typedef uint64 ID;
constexpr ID INVALID_ID = 0;

typedef bool (*PropertyComputeFunc)(void* data, const Array<CString> args, const MoleculeDynamic* dynamic, int res_idx);
typedef bool (*ResidueMatchFunc)(const Array<CString> args, const MoleculeStructure* mol, int res_idx);

enum struct PropertyType { FLOAT32, UNKNOWN };

inline int32 get_stride(PropertyType type) {
	switch (type) {
	case PropertyType::FLOAT32:
		return 4;
	default:
		return 4;
	}
}

bool compute_stats(MoleculeDynamic* dynamic);
void clear_stats();
void store_stats(CString filename);
void load_stats(CString filename);

ID      create_group(CString name, CString cmd, CString args);
void	remove_group(ID group_id);

ID      get_group(CString name);
ID      get_group(int32 idx);
int32   get_group_count();

ID      get_property(ID group_id, CString name);
ID      get_property(ID group_id, int32 idx);
int32   get_property_count(ID group_id);

ID		create_property(ID group_id, CString name, CString cmd, CString args);
void	remove_property(ID prop_id);

void*        get_property_data(ID prop_id, int32 residue_idx);
void*        get_property_avg_data(ID prop_id);
int32        get_property_data_count(ID prop_id);
PropertyType get_property_type(ID prop_id);
CString	     get_property_name(ID prop_id);

void	  register_property_command(CString command, PropertyType type, PropertyComputeFunc func);
void	  register_group_command(CString command, ResidueMatchFunc func);

void initialize();
void shutdown();

}  // namespace stats