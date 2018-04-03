#pragma once

#include <core/types.h>
#include <core/string_utils.h>


struct MoleculeDynamic;
struct MoleculeStructure;
struct Residue;

namespace stats {

typedef uint64 ID;
constexpr ID INVALID_ID = 0;

typedef bool (*PropertyComputeFunc)(float* data, const Array<CString> args, const MoleculeDynamic* dynamic, int res_idx);
typedef bool (*ResidueMatchFunc)(const Array<CString> args, const MoleculeStructure* mol, int res_idx);

typedef vec2 Range;

struct Histogram {
    DynamicArray<float> bins = {};
    Range val_range = {};
    Range bin_range = {};
    int32 num_samples = 0;
};

struct PropertyCommandDescriptor {
	CString cmd;
	PropertyComputeFunc func;
	float min_val = 0.f;
	float max_val = 1.f;
	bool periodic = false;
	CString unit = "";
};

// HISTOGRAM
Histogram compute_histogram(int32 num_bins, Array<float> data);
Histogram compute_histogram(int32 num_bins, Array<float> data, float min_val, float max_val);

void compute_histogram(Histogram* hist, int32 num_bins, Array<float> data);
void compute_histogram(Histogram* hist, int32 num_bins, Array<float> data, float min_val, float max_val);

// STATS
bool compute_stats(MoleculeDynamic* dynamic);
void clear_stats();
void store_stats(CString filename);
void load_stats(CString filename);

void register_property_command(CString command, PropertyComputeFunc func);
void register_group_command(CString command, ResidueMatchFunc func);

void initialize();
void shutdown();

// GROUP
ID      create_group(CString name, CString cmd_and_args);
void	remove_group(ID group_id);

ID      get_group(CString name);
ID      get_group(int32 idx);
int32   get_group_count();

int32   get_property_count(ID group_id);
ID      get_property(ID group_id, CString name);
ID      get_property(ID group_id, int32 idx);

Array<ID> get_groups_with_property(CString property_name);

// PROPERTY
ID	 create_property(ID group_id, CString name, CString cmd_and_args);
void remove_property(ID prop_id);

CString			 get_property_name(ID prop_id);
CString			 get_property_unit(ID prop_id);
float			 get_property_min_val(ID prop_id);
float			 get_property_max_val(ID prop_id);
void*			 get_property_data(ID prop_id, int32 residue_idx);
void*			 get_property_avg_data(ID prop_id);
int32			 get_property_data_count(ID prop_id);

Histogram*		 get_property_histogram(ID prop_id, int32 residue_idx);
Histogram*		 get_property_avg_histogram(ID prop_id);

// PROPERTY FILTER
void set_property_filter_min(ID prop_id, float min_val);
void set_property_filter_max(ID prop_id, float max_val);
void set_property_filter(ID prop_id, float min_val, float max_val);


}  // namespace stats