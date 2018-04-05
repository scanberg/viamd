#pragma once

#include <core/types.h>
#include <core/string_utils.h>


struct MoleculeDynamic;
struct MoleculeStructure;

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
	CString command_keyword;
	PropertyComputeFunc compute_function;
	float min_val = FLT_MAX;
	float max_val = FLT_MAX;
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

void register_property_command(PropertyCommandDescriptor cmd_desc);
void register_group_command(CString command, ResidueMatchFunc func);

void initialize();
void shutdown();

// GROUP
ID      create_group(CString name, CString cmd_and_args);
void	remove_group(ID group_id);

ID      get_group(CString name);
ID      get_group(int32 idx);
int32   get_group_count();

// PROPERTY
ID		create_property(CString name, CString cmd_and_args);
void	remove_property(ID prop_id);

ID		get_property(CString name);
ID		get_property(int32 idx);
int32	get_property_count();

CString	get_property_name(ID prop_id);
CString	get_property_unit(ID prop_id);
bool	get_property_periodic(ID prop_id);
float	get_property_min_val(ID prop_id);
float	get_property_max_val(ID prop_id);

void	set_property_filter_min(ID prop_id, float min_val);
void	set_property_filter_max(ID prop_id, float max_val);
void	set_property_filter(ID prop_id, float min_val, float max_val);

// PROPERTY DATA
Array<float> get_property_data(ID prop_id, int32 instance_idx);
Array<float> get_property_avg_data(ID prop_id);

//Histogram*	 get_property_histogram(ID prop_id, int32 residue_idx);
//Histogram*	 get_property_avg_histogram(ID prop_id);

// PROPERTY FILTER



}  // namespace stats