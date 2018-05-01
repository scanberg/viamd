#include "stats.h"
#include <core/math_utils.h>
#include <core/hash.h>
#include <mol/molecule_dynamic.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#include <ctype.h>

#define COMPUTE_ID(x) (hash::crc64(x))

namespace stats {

typedef DynamicArray<Structure>(*StructureMatchFunc)(const Array<CString> args, const MoleculeStructure& molecule);

struct PropertyCommand {
    ID id = INVALID_ID;
	CString keyword = "";
    PropertyComputeFunc func = nullptr;
    Range val_range {};
    bool periodic = false;
    CString unit = "";
};

struct StructureCommand {
	ID id = INVALID_ID;
	CString keyword = "";
};

struct Property {
	ID id = INVALID_ID;
	StringBuffer<32>  name{};
	StringBuffer<256> args{};
	StringBuffer<32>  unit{};
	StringBuffer<128> error_msg{};
	bool valid = false;
	bool periodic = false;
	Range filter{ 0,0 };
	Range data_range{ 0,0 };
	Array<float> data{};
	Histogram hist;
};

struct PropertyEntry {
	ID id = INVALID_ID;
	Property* ptr = nullptr;
};

struct StatisticsContext {
    ID next_id = 1;
    DynamicArray<ID> free_ids {};
    DynamicArray<String> string_buffer {};

    DynamicArray<PropertyEntry> property_entries{};

    DynamicArray<PropertyCommand> property_commands {};
    DynamicArray<StructureCommand> structure_commands {};
};

static StatisticsContext ctx;

static Property* find_property(ID id) {
	for (auto& e : ctx.property_entries) {
		if (id == e.id) {
			return e.ptr;
		}
	}
	return nullptr;
}

template <typename T>
static T* find_id(Array<T> data, ID id) {
	for (auto& item : data) {
		if (item.id == id) return &item;
	}
	return nullptr;
}

static ID create_id() {
	if (ctx.free_ids.count > 0) return ctx.free_ids.pop_back();
	return ctx.next_id++;
}

static CString alloc_string(CString str) {
	char* data = (char*)MALLOC(str.count);
	ctx.string_buffer.push_back({ data, str.count });
	copy(ctx.string_buffer.back(), str);
	return ctx.string_buffer.back();
}

static void free_string(CString str) {
	for (auto int_str : ctx.string_buffer) {
		if (int_str.count == str.count && compare(int_str, str)) {
			ctx.string_buffer.remove(&int_str);
			return;
		}
	}
}

// HISTOGRAMS
void init_histogram(Histogram* hist, int32 num_bins) {
	ASSERT(hist);
	if (hist->bins.data) {
		FREE(hist->bins.data);
	}
	hist->bins.data = (float*)MALLOC(num_bins * sizeof(float));
	hist->bins.count = num_bins;

	ASSERT(hist->bins.data);
}

void free_histogram(Histogram* hist) {
	ASSERT(hist);
	if (hist->bins.data) {
		FREE(hist->bins.data);
		hist->bins.data = nullptr;
		hist->bins.count = 0;
	}
}

Histogram compute_histogram(int32 num_bins, Array<float> data) {
	Histogram hist;
	compute_histogram(&hist, num_bins, data);
	return hist;
}

Histogram compute_histogram(int32 num_bins, Array<float> data, float min_val, float max_val) {
	Histogram hist;
	compute_histogram(&hist, num_bins, data, min_val, max_val);
	return hist;
}

void compute_histogram(Histogram* hist, int32 num_bins, Array<float> data) {
	ASSERT(hist);
	if (data.count == 0) return;
	float min_val = FLT_MAX;
	float max_val = -FLT_MAX;
	for (const auto& d : data) {
		min_val = math::min(min_val, d);
		max_val = math::max(max_val, d);
	}
	compute_histogram(hist, num_bins, data, min_val, max_val);
}

void compute_histogram(Histogram* hist, int32 num_bins, Array<float> data, float min_val, float max_val) {
	ASSERT(hist);
	ASSERT(num_bins > 0);

	init_histogram(hist, num_bins);
	memset(hist->bins.data, 0, hist->bins.count * sizeof(float));

	const float scl = num_bins / (max_val - min_val);
	hist->bin_range = { 0,0 };
	for (auto v : data) {
		int32 bin_idx = math::clamp((int32)((v - min_val) * scl), 0, num_bins - 1);
		hist->bins[bin_idx]++;
		hist->bin_range.y = math::max(hist->bin_range.y, hist->bins[bin_idx]);
	}
	hist->value_range = { min_val, max_val };
}

static void set_error_message(CString msg) {
    printf("%s\n", msg.beg());
}

typedef bool(*StructureFunc)(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule);


bool structure_match_resname(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule) {
	ASSERT(data);

	// Expect args.count to be > 0
	if (args.count == 0) {
		set_error_message("Expects one or more arguments for resname");
		return false;
	}

	for (const auto& res : molecule.residues) {
		for (const auto& arg : args) {
			if (compare(res.name, arg)) {
				data->structures.push_back({ res.beg_atom_idx, res.end_atom_idx });
				break;
			}
		}
	}

	return true;
}

bool structure_match_resid(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule) {
	ASSERT(data);

	// Expect args to be  > 0
	if (args.count == 0) {
		set_error_message("Expects one or more arguments for resid");
		return false;
	}

	for (const auto& arg : args) {
		auto id = to_int(arg);
		if (!id.success) {
			set_error_message("Failed to parse argument for resid");
			return false;
		}
		for (const auto& res : molecule.residues) {
			if (res.id == id) {
				data->structures.push_back({ res.beg_atom_idx, res.end_atom_idx });
				break;
			}
		}
	}

	return true;
}

bool structure_match_residue(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule) {
	ASSERT(data);

	// Expect args to be  > 0
	if (args.count == 0) {
		set_error_message("Expects one or more arguments for residue");
		return false;
	}

	for (const auto& arg : args) {
		auto id = to_int(arg);
		if (!id.success) {
			set_error_message("Failed to parse argument for residue");
			return false;
		}
		auto idx = id - 1;
		if (idx < 0 || molecule.residues.count <= idx) {
			set_error_message("Index for residue is out of bounds");
			return false;
		}
		const auto& res = molecule.residues[idx];
		data->structures.push_back({ res.beg_atom_idx, res.end_atom_idx });
	}

	return true;
}

bool structure_match_atom(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule) {
	ASSERT(data);

	// Expect args to be  > 0
	if (args.count == 0) {
		set_error_message("Expects one or more arguments for atom");
		return false;
	}

	for (const auto& arg : args) {
		auto id = to_int(arg);
		if (!id.success) {
			set_error_message("Failed to parse argument for atom");
			return false;
		}
		auto idx = id - 1;
		if (idx < 0 || molecule.atom_positions.count <= idx) {
			set_error_message("Index for atom is out of bounds");
			return false;
		}
		data->structures.push_back({ idx, idx + 1 });
	}

	return true;
}

bool structure_extract_resatom(StructureData* data, const Array<CString> args, const MoleculeStructure&) {
	ASSERT(data);
	if (args.count == 0) {
		set_error_message("resatom requires 1 argument");
		return false;
	}

	auto res = to_int(args[0]);
	if (!res.success) {
		set_error_message("resatom: failed to parse argument");
		return false;
	}
	int offset = res - 1;

	for (auto& s : data->structures) {
		int count = s.end_idx - s.beg_idx;
		if (offset >= count) {
			set_error_message("restom: Index is out of range for structure");
			return false;
		}
		if (count > 1) {
			int new_idx = s.beg_idx + offset;
			s = { new_idx, new_idx + 1 };
		}
	}
    return true;
}

bool structure_apply_aggregation_strategy_com(StructureData* data, const Array<CString>, const MoleculeStructure&) {
	ASSERT(data);
	data->strategy = AggregationStrategy::COM;
	return true;
}

// Helper funcs
inline static vec3 compute_com(Array<const vec3> positions) {
	if (positions.count == 0) return { 0,0,0 };
	if (positions.count == 1) return positions[0];
	
	vec3 com{ 0 };
	for (const auto& p : positions) {
		com += p;
	}

	return com / (float)positions.count;
}

void compute_frame_positions(Array<vec3> dst, const StructureData& data, const MoleculeDynamic& dynamic, int frame_index) {
	ASSERT(dst.data);
	ASSERT(dst.count >= data.structures.count);
	
	Array<vec3> positions = get_trajectory_positions(dynamic.trajectory, frame_index);
	switch (data.strategy) {
	case AggregationStrategy::COM:
		for (int32 i = 0; i < data.structures.count; i++) {
			dst[i] = compute_com(positions.sub_array(data.structures[i].beg_idx, data.structures[i].end_idx - data.structures[i].beg_idx));
		}
		break;
	default:
		for (int32 i = 0; i < data.structures.count; i++) {
			dst[i] = positions[data.structures[i].beg_idx];
		}
		break;
	}
}

DynamicArray<vec3> compute_frame_positions(const StructureData& data, const MoleculeDynamic& dynamic, int frame_index) {
	DynamicArray<vec3> positions(data.structures.count);
	compute_frame_positions(positions, data, dynamic, frame_index);
	return positions;
}

static bool compute_distance(Array<float32> data, const Array<CString> args, const MoleculeDynamic& dynamic) {
	ASSERT(data);
	if (args.count != 2) {
		set_error_message("distance expects 2 arguments");
		return false;
	}

	StructureData arg_structure_data[2];
	extract_args_structures(arg_structure_data, args, dynamic.molecule);

	int32 count = (int32)arg_structure_data[0].structures.count;
	DynamicArray<vec3> pos[2] { {count}, {count} };

	for (int32 i = 0; i < dynamic.trajectory.num_frames; i++) {
		compute_frame_positions(pos[0], arg_structure_data[0], dynamic, i);
		compute_frame_positions(pos[1], arg_structure_data[1], dynamic, i);

        float32 dist = 0.f;
        for (int32 j = 0; j < count; j++) {
		    dist += math::distance(pos[0][j], pos[1][j]);
        }
		data[i] = dist / (float32)count;
	}

	return true;
}

static bool compute_angle(Array<float32> data, const Array<CString> args, const MoleculeDynamic& dynamic) {
	ASSERT(data);
	if (args.count != 3) {
		set_error_message("angle expects 2 arguments");
		return false;
	}

	StructureData arg_structure_data[3];
	if (!extract_args_structures(arg_structure_data, args, dynamic.molecule)) {
		return false;
	}

	int32 count = (int32)arg_structure_data[0].structures.count;
	DynamicArray<vec3> pos[3]{ { count }, { count }, { count } };

	for (int32 i = 0; i < dynamic.trajectory.num_frames; i++) {
		compute_frame_positions(pos[0], arg_structure_data[0], dynamic, i);
		compute_frame_positions(pos[1], arg_structure_data[1], dynamic, i);
		compute_frame_positions(pos[2], arg_structure_data[2], dynamic, i);

		float32 angle = 0.f;
		for (int32 j = 0; j < count; j++) {
			vec3 a = pos[0][j] - pos[1][j];
			vec3 b = pos[2][j] - pos[1][j];
			angle += math::angle(a, b);
		}
		data[i] = angle / (float32)count;
	}

	return true;
}

void initialize() {

}

void shutdown() {

}

static CString extract_parentheses(CString str) {
	const char* beg = str.beg();

	while (beg < str.end() && *beg != '(') beg++;
	if (beg < str.end()) beg++;
	if (beg == str.end()) return { beg, str.end() };

	const char* end = beg;
	int count = 1;
	while (end != str.end()) {
		if (*end == '(') count++;
		else if (*end == ')' && --count == 0) break;
		end++;
	}

	return { beg, end };
}

static const char* find_character(CString str, char c) {
	const char* ptr = str.beg();
	while (ptr < str.end() && *ptr != c) ptr++;
	return ptr;
}

static DynamicArray<CString> extract_arguments(CString str) {
	DynamicArray<CString> args;

	const char* beg = str.beg();
	const char* end = str.beg();
	int32 count = 0;

	while (end < str.end()) {
		if (*end == '(')
			count++;
		else if (*end == ')')
			count--;
		else if (*end == ',') {
			if (count == 0) {
				args.push_back(trim({ beg, end }));
				beg = end + 1;
			}
		}
		end++;
	}
	if (beg != end) args.push_back(trim({ beg, end }));

	return args;
}

static CString extract_command(CString str) {
	str = trim(str);
	const char* ptr = str.beg();
	while (ptr != str.end() && *ptr != '(' && !isspace(*ptr)) ptr++;
	return { str.beg(), ptr };
}

static StructureFunc get_structure_func(CString cmd) {
	if (compare(cmd, "resname")) {
		return structure_match_resname;
	}
	else if (compare(cmd, "residue")) {
		return structure_match_residue;
	}
	else if (compare(cmd, "resid")) {
		return structure_match_resid;
	}
	else if (compare(cmd, "atom")) {
		return structure_match_atom;
	}
	else if (compare(cmd, "resatom")) {
		return structure_extract_resatom;
	}
	else if (compare(cmd, "com")) {
		return structure_apply_aggregation_strategy_com;
	}
	return nullptr;
}

bool extract_structures(StructureData* data, CString arg, const MoleculeStructure& molecule) {
	CString cmd = extract_command(arg);
	auto func = get_structure_func(cmd);

	if (!func) {
		char buf[64];
		snprintf(buf, 64, "Could not identify command '%s'", cmd.beg());
		set_error_message(buf);
		return false;
	}

	CString outer = extract_parentheses(arg);
	DynamicArray<CString> cmd_args = extract_arguments(outer);

	// @NOTE: ONLY ALLOW RECURSION FOR FIRST ARGUMENT?
	if (cmd_args.count > 0 && find_character(cmd_args[0], '(') != cmd_args[0].end()) {
		if (!extract_structures(data, cmd_args[0], molecule)) return false;
		cmd_args = cmd_args.sub_array(1);
	}
	if (!func(data, cmd_args, molecule)) return false;

	return true;
}

bool extract_args_structures(Array<StructureData> data, Array<CString> args, const MoleculeStructure& molecule) {
	ASSERT(data.count == args.count);

	for (int i = 0; i < data.count; i++) {
		if (!extract_structures(&data[i], args[i], molecule)) return false;
	}

	int64 multi_idx = -1;
	int64 multi_count = 0;
	for (int64 i = 0; i < data.count; i++) {
		if (data[i].structures.count == 0) {
			set_error_message("One argument did not match any structures");
			return false;
		}

		if (data[i].structures.count > 1) {
			if (multi_idx != -1 && multi_count != data[i].structures.count) {
				set_error_message("Multiple structures found for more than one argument, but their count did not match");
				return false;
			}
			multi_idx = i;
			multi_count = data[i].structures.count;
		}
	}

	// Extend and copy data from first element to match multi_count
	for (int64 i = 0; i < data.count; i++) {
		while (data[i].structures.count < multi_count) {
			data[i].structures.push_back(data[i].structures.front());
		}
	}

	return true;
}

bool balanced_parentheses(CString str) {
	int count = 0;
	const char* ptr = str.beg();
	while (ptr != str.end()) {
		if (*ptr == '(') count++;
		else if (*ptr == ')') count--;
		ptr++;
	}
	return count == 0;
}

static PropertyComputeFunc get_property_func(CString cmd) {
	if (compare(cmd, "distance")) {
		return compute_distance;
	}
	else if (compare(cmd, "angle")) {
		return compute_angle;
	}
	return nullptr;
}

static Range compute_range(Array<float> data) {
	if (data.count == 0) {
		return { 0,0 };
	}
	Range range{ FLT_MAX, -FLT_MAX };
	for (float v : data) {
		range.x = math::min(range.x, v);
		range.y = math::max(range.y, v);
	}
	return range;
}

bool compute_stats(const MoleculeDynamic& dynamic) {
	int num_frames = dynamic.trajectory.num_frames;
	
	for (auto& pe : ctx.property_entries) {
		auto& prop = *pe.ptr;
		if (prop.data.count != num_frames) {
			if (prop.data.data) FREE(prop.data.data);
			prop.data = { (float*)MALLOC(num_frames * sizeof(float)), num_frames };
		}
		memset(prop.data.data, 0, prop.data.count * sizeof(float));

		if (!balanced_parentheses(prop.args)) {
			snprintf(prop.error_msg.beg(), prop.error_msg.MAX_LENGTH, "Unbalanced parantheses!");
			prop.valid = false;
			continue;
		}

		DynamicArray<CString> args;

		// Extract big argument chunks
		const char* beg = prop.args.beg();
		const char* end = prop.args.beg();
		int count = 0;

		// Use space separation unless we are inside a parenthesis
		while (end != prop.args.end()) {
			if (*end == '(') count++;
			else if (*end == ')') count--;
			else if (isspace(*end) && count == 0) {
				args.push_back(trim({ beg, end }));
				beg = end + 1;
			}
			end++;
		}
		if (beg != end) args.push_back(trim({ beg,end }));

		if (args.count == 0) {
			prop.valid = false;
			continue;
		}

		CString cmd = args[0];
		args = args.sub_array(1);
		
		auto func = get_property_func(cmd);
		if (!func) {
			snprintf(prop.error_msg.beg(), prop.error_msg.MAX_LENGTH, "Could not recognize command '%s'", cmd);
			prop.valid = false;
			continue;
		}

		if (func(prop.data, args, dynamic)) {
			prop.data_range = compute_range(prop.data);
			prop.valid = true;
		}
	}

	return true;
}

ID create_property(CString name, CString args) {
	Property* prop = (Property*)MALLOC(sizeof(Property));
	new(prop) Property();
	prop->id = create_id();
	prop->name = name;
	prop->args = args;
	prop->valid = false;

	ctx.property_entries.push_back({ prop->id, prop });
	return prop->id;
}

void remove_property(ID prop_id) {
	for (auto& pe : ctx.property_entries) {
		if (pe.id == prop_id) {
			auto prop = pe.ptr;
			if (prop->data.data) FREE(prop->data.data);
			ctx.property_entries.swap_back_and_pop(&pe);
		}
	}
}

ID get_property(CString name) {
	for (const auto &pe : ctx.property_entries) {
		if (compare(pe.ptr->name, name)) return pe.id;
	}
	return INVALID_ID;
}

ID get_property(int32 idx) {
	if (-1 < idx && idx < ctx.property_entries.count) return ctx.property_entries[idx].id;
	return INVALID_ID;
}

int32 get_property_count() {
	return (int32)ctx.property_entries.count;
}

Range get_property_data_range(ID prop_id) {
	auto prop = find_property(prop_id);
	if (prop) {
		return prop->data_range;
	}
	return {};
}

Array<float> get_property_data(ID prop_id) {
    auto prop = find_property(prop_id);
    if (prop) {
		return prop->data;
    }
    return {};
}

Histogram* get_property_histogram(ID prop_id) {
	auto prop = find_property(prop_id);
	if (prop) {
		return &prop->hist;
	}
	return nullptr;
}

StringBuffer<32>* get_property_name_buf(ID prop_id) {
    auto prop = find_property(prop_id);
    if (prop) {
        return &prop->name;
    }
    return nullptr;
}
    
StringBuffer<256>* get_property_args_buf(ID prop_id) {
    auto prop = find_property(prop_id);
    if (prop) {
        return &prop->args;
    }
    return nullptr;
}

CString get_property_name(ID prop_id) {
	auto prop = find_property(prop_id);
	if (prop) {
		return prop->name;
	}
	return {};
}
    
bool get_property_valid(ID prop_id) {
    auto prop = find_property(prop_id);
    if (prop) {
        return prop->valid;
    }
    return false;
}

CString get_property_error_message(ID prop_id) {
	auto prop = find_property(prop_id);
	if (prop) {
		return prop->error_msg;
	}
	return "";
}

CString get_property_unit(ID prop_id) {
	auto prop = find_property(prop_id);
	if (prop) {
		return prop->unit;
	}
	return "";
}

bool get_property_periodic(ID prop_id) {
	auto prop = find_property(prop_id);
	if (prop) {
		return false;
	}
	return false;
}

}  // namespace stats
