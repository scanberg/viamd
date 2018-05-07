#pragma once

#include <core/types.h>
#include <core/string_utils.h>

struct MoleculeDynamic;
struct MoleculeStructure;

namespace stats {

typedef uint64 ID;
constexpr ID INVALID_ID = 0;

typedef vec2 Range;

struct Histogram {
    Array<float> bins = {};
    Range value_range = {};
    Range bin_range = {};
    int32 num_samples = 0;
};

struct Structure {
	int32 beg_idx;
	int32 end_idx;
};

enum AggregationStrategy {
	COM,
	NONE
};

struct StructureData {
	DynamicArray<Structure> structures{};
	AggregationStrategy strategy = NONE;
};

struct InstanceData {
	InstanceData(int32 size) { data = allocate_array<float>(size); }
	~InstanceData() { free_array(&data); }

	Array<float> data;
};

struct Property {
	StringBuffer<32>  name{};
	StringBuffer<256> args{};
	StringBuffer<32>  unit{};
	StringBuffer<128> error_msg{};

	bool valid = false;
	bool periodic = false;
	bool visualize = false;

	Range filter{ 0,0 };
	Range data_range{ 0,0 };
	DynamicArray<float> data{};
	Histogram histogram;

	Array<InstanceData> instance_data{};
	Array<StructureData> structure_data{};
};

// Helper functions
void set_error_message(CString msg);

void init_instance_data(Array<InstanceData>* instance_data, int32 num_instances, int32 num_frames);
void free_instance_data(Array<InstanceData>* instance_data);

void init_structure_data(Array<StructureData>* structure_data, int32 count);
bool sync_structure_data_length(Array<StructureData> data);
bool extract_args_structures(Array<StructureData> data, Array<CString> arg, const MoleculeStructure& structure);

Array<const vec3> extract_positions(Structure structure, Array<const vec3> atom_positions);
vec3 compute_com(Array<const vec3> positions);

typedef bool (*PropertyComputeFunc)(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic);
typedef bool (*PropertyVisualizeFunc)(const Property& prop, const MoleculeDynamic& dynamic);

// HISTOGRAM
void init_histogram(Histogram* hist, int32 num_bins, Range value_range);
void free_histogram(Histogram* hist);

//Histogram compute_histogram(int32 num_bins, Array<float> data);
//Histogram compute_histogram(int32 num_bins, Array<float> data, float min_val, float max_val);

//void compute_histogram(Histogram* hist, Array<float> data);
void compute_histogram(Histogram* hist, Array<float> data, float min_val, float max_val);
void clear_histogram(Histogram* hist);

// STATS
void initialize();
void shutdown();

bool compute_stats(const MoleculeDynamic& dynamic);
void visualize(const MoleculeDynamic& dynamic);
void update_property(Property* prop);

bool register_property_command(CString cmd_keyword, PropertyComputeFunc compute_func, PropertyVisualizeFunc visualize_func);

Array<CString> get_property_commands();
Array<CString> get_structure_commands();
Array<CString> get_property_names();

// PROPERTY
Property* create_property(CString name = {}, CString args = {});
void	  remove_property(Property* prop);
void	  remove_all_properties();

// Keep property, but remove the generated data
void	  clear_property(Property* prop);
void      clear_all_properties();

int32	  get_property_count();
Property* get_property(int32 idx);
Property* get_property(CString name);


}  // namespace stats
