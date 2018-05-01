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

/*

Examples of what to support and syntax
#name #args
prop1 distance com(resname(ala)) (resname ala, resatom 1)
prop2 distance (residue 2)
a2	  angle    com(resname(PFT)) atom 1
d1    dihedral resatom(resname(PFT), 1 2 3 4)
d2    dihedral (resname PFT, resatom 1) (resname PFT, resatom 2) (resname PFT, resatom 3) (resname PFT, resatom 4)
d3    dihedral (resname PFT, resatom 1) (resname PFT, resatom 2) (resname PFT, resatom 3) (resname PFT, resatom 4)
d4    dihedral (resname PFT, resatom 1) (resname PFT, resatom 2) (resname PFT, resatom 3) (resname PFT, resatom 4)
rdf   rdf	   com(resname(PFT))
den   density  
rmsd  rmsd
plan  planarity d1 d2 d3 d4

prop5 coulombic resname(ala) 
prop7 lj12		resname(pft) protein

*/

struct Structure {
	int32 beg_idx;
	int32 end_idx;
};

enum AggregationStrategy {
	COM,
	NONE
};

struct StructureData {
	DynamicArray<Structure> structures;
	AggregationStrategy strategy;
};

struct StructData {
	struct ArgData {
		Array<Structure> structures;
		AggregationStrategy strategy;
	};

	DynamicArray<Structure> structures;
	DynamicArray<ArgData> arg_structures;
};

// Helper functions
bool extract_args_structures(Array<StructureData> data, Array<CString> arg, const MoleculeStructure& structure);
void compute_frame_positions(Array<vec3> dst, const StructureData& data, const MoleculeDynamic& dynamic, int frame_idx);

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

bool compute_stats(const MoleculeDynamic& dynamic);
//void clear();

void register_property_command(CString cmd_keyword, PropertyComputeFunc func, Range value_range, bool periodic, CString unit);

Array<CString> get_property_commands();
Array<CString> get_structure_commands();

// PROPERTY
ID      create_property(CString name = {}, CString args = {});
void	remove_property(ID prop_id);

void	clear_property(ID prop_id);
void    clear_properties();

ID		get_property(CString name);
ID		get_property(int32 idx);

int32	get_property_count();

bool			   get_property_valid(ID prop_id); 
CString			   get_property_error_message(ID prop_id);
StringBuffer<32>*  get_property_name_buf(ID prop_id);
StringBuffer<256>* get_property_args_buf(ID prop_id);
CString			   get_property_name(ID prop_id);
Range			   get_property_data_range(ID prop_id);
Array<float>	   get_property_data(ID prop_id);
bool			   get_property_periodic(ID prop_id);
Histogram*		   get_property_histogram(ID prop_id);


//void		 clear_property_data();


}  // namespace stats
