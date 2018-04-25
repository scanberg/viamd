#pragma once

#include <core/types.h>
#include <core/string_utils.h>

struct MoleculeDynamic;
struct MoleculeStructure;

namespace stats {

typedef uint64 ID;
constexpr ID INVALID_ID = 0;

typedef vec2 Range;

struct Structure {
    int32 beg_atom_idx = 0;
    int32 end_atom_idx = 0;
};

struct Property {
    CString name{};
    CString unit{};
    bool valid = false;
    CString error_message{};
    Range filter {0,0};
};

struct Histogram {
    Array<float> bins = {};
    Range value_range = {};
    Range bin_range = {};
    int32 num_samples = 0;
};

// Helper function for matching structures from arguments
DynamicArray<Structure> match_structures(CString arg);

typedef bool (*PropertyComputeFunc)(Array<float> data, const Array<CString> args, const MoleculeDynamic& dynamic);
//typedef DynamicArray<Structure> (*StructureExtractFunc)(const Array<CString> args, const MoleculeStructure& mol);

// HISTOGRAM
void init_histogram(Histogram* hist, int32 num_bins);
void free_histogram(Histogram* hist);

Histogram compute_histogram(int32 num_bins, Array<float> data);
Histogram compute_histogram(int32 num_bins, Array<float> data, float min_val, float max_val);

void compute_histogram(Histogram* hist, int32 num_bins, Array<float> data);
void compute_histogram(Histogram* hist, int32 num_bins, Array<float> data, float min_val, float max_val);

// STATS
void initialize();
void shutdown();

bool validate_properties();
bool validate_arguments();
bool compute_stats(const MoleculeDynamic& dynamic);
void clear();

void register_property_command(CString cmd_keyword, PropertyComputeFunc func, Range value_range, bool periodic, CString unit);

Array<CString> get_property_commands();
Array<CString> get_structure_commands();

// PROPERTY
ID      create_property(CString name = {}, CString cmd_and_args = {});
void	remove_property(ID prop_id);

void	clear_property(ID prop_id);
void    clear_properties();

ID		get_property(CString name);
ID		get_property(int32 idx);
int32	get_property_count();
 
StringBuffer<32>*  get_property_name_buf(ID prop_id);
StringBuffer<128>* get_property_args_buf(ID prop_id);

// PROPERTY DATA
Range        get_property_data_range(ID prop_id);
Array<float> get_property_data(ID prop_id);
Histogram*   get_property_histogram(ID prop_id);

void		 clear_property_data();


}  // namespace stats
