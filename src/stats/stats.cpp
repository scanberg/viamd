#include "stats.h"
#include <core/math_utils.h>
#include <core/hash.h>
#include <mol/molecule_dynamic.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <ctype.h>

#define HASH(x) (hash::crc64(x))

namespace stats {

typedef bool(*StructureFunc)(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule);

struct PropertyFuncEntry {
	ID hash = INVALID_ID;
	PropertyComputeFunc compute_func = nullptr;
	PropertyVisualizeFunc visualize_func = nullptr;
};

struct StructureFuncEntry {
	ID hash = INVALID_ID;
	StructureFunc func = nullptr;
};

struct StatisticsContext {
	DynamicArray<PropertyFuncEntry> property_func_entries{};
	DynamicArray<StructureFuncEntry> structure_func_entries{};
	DynamicArray<Property*> properties {};
    DynamicArray<String> string_buffer {};
};

static StatisticsContext ctx;

static PropertyFuncEntry* find_property_func_entry(ID hash) {
	for (auto& e : ctx.property_func_entries) {
		if (hash == e.hash) {
			return &e;
		}
	}
	return nullptr;
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

void clear_histogram(Histogram* hist) {
	ASSERT(hist);
	if (hist->bins.data) {
		memset(hist->bins.data, 0, hist->bins.count * sizeof(float));
	}
}

static void set_error_message(CString msg) {
    printf("%s\n", msg.beg());
}

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

void compute_positions(Array<vec3> dst, const StructureData& data, Array<const vec3> atom_positions) {
	ASSERT(dst.data);
	ASSERT(dst.count >= data.structures.count);
	
	switch (data.strategy) {
	case AggregationStrategy::COM:
		for (int32 i = 0; i < data.structures.count; i++) {
			dst[i] = compute_com(atom_positions.sub_array(data.structures[i].beg_idx, data.structures[i].end_idx - data.structures[i].beg_idx));
		}
		break;
	default:
		for (int32 i = 0; i < data.structures.count; i++) {
			dst[i] = atom_positions[data.structures[i].beg_idx];
		}
		break;
	}
}

DynamicArray<vec3> compute_positions(const StructureData& data, Array<const vec3> atom_positions) {
	DynamicArray<vec3> positions(data.structures.count);
	compute_positions(positions, data, atom_positions);
	return positions;
}

static bool compute_distance(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic) {
	ASSERT(prop);
	if (args.count != 2) {
		set_error_message("distance expects 2 arguments");
		return false;
	}

	init_structure_data(&prop->structure_data, 2);
	extract_args_structures(prop->structure_data, args, dynamic.molecule);

	const int32 count = (int32)prop->structure_data[0].structures.count;
	void* tmp_data = TMP_MALLOC(count * 2 * sizeof(vec3));
	Array<vec3> pos[2] = {
		{ (vec3*)tmp_data + count * 0, count },
		{ (vec3*)tmp_data + count * 1, count } };

	const float32 scl = 1.f / (float32)count;
	for (int32 i = 0; i < dynamic.trajectory.num_frames; i++) {
		Array<const vec3> atom_positions = get_trajectory_positions(dynamic.trajectory, i);
		compute_positions(pos[0], prop->structure_data[0], atom_positions);
		compute_positions(pos[1], prop->structure_data[1], atom_positions);

        float32 dist = 0.f;
        for (int32 j = 0; j < count; j++) {
		    dist += math::distance(pos[0][j], pos[1][j]);
        }
		prop->data[i] = dist * scl;
	}

	TMP_FREE(tmp_data);

	return true;
}

static bool compute_angle(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic) {
	ASSERT(prop);
	if (args.count != 3) {
		set_error_message("angle expects 3 arguments");
		return false;
	}

	init_structure_data(&prop->structure_data, 3);
	if (!extract_args_structures(prop->structure_data, args, dynamic.molecule)) {
		return false;
	}

	const int32 count = (int32)prop->structure_data[0].structures.count;
	void* tmp_data = TMP_MALLOC(count * 3 * sizeof(vec3));
	Array<vec3> pos[3] = {
		{ (vec3*)tmp_data + count * 0, count },
		{ (vec3*)tmp_data + count * 1, count },
		{ (vec3*)tmp_data + count * 2, count } };

	const float32 scl = 1.f / (float32)count;
	for (int32 i = 0; i < dynamic.trajectory.num_frames; i++) {
		Array<const vec3> atom_positions = get_trajectory_positions(dynamic.trajectory, i);
		compute_positions(pos[0], prop->structure_data[0], atom_positions);
		compute_positions(pos[1], prop->structure_data[1], atom_positions);
		compute_positions(pos[2], prop->structure_data[2], atom_positions);

		float32 angle = 0.f;
		for (int32 j = 0; j < count; j++) {
			vec3 a = pos[0][j] - pos[1][j];
			vec3 b = pos[2][j] - pos[1][j];
			angle += math::angle(a, b);
		}
		prop->data[i] = angle * scl;
	}

	TMP_FREE(tmp_data);

	return true;
}

static bool compute_dihedral(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic) {
	ASSERT(prop);
	if (args.count != 4) {
		set_error_message("dihedral expects 4 arguments");
		return false;
	}

	init_structure_data(&prop->structure_data, 4);
	if (!extract_args_structures(prop->structure_data, args, dynamic.molecule)) {
		return false;
	}

	const int32 count = (int32)prop->structure_data[0].structures.count;
	void* tmp_data = TMP_MALLOC(count * 4 * sizeof(vec3));
	Array<vec3> pos[4] = {
		{ (vec3*)tmp_data + count * 0, count },
		{ (vec3*)tmp_data + count * 1, count }, 
		{ (vec3*)tmp_data + count * 2, count }, 
		{ (vec3*)tmp_data + count * 3, count } };

	const float32 scl = 1.f / (float32)count;
	for (int32 i = 0; i < dynamic.trajectory.num_frames; i++) {
		Array<const vec3> atom_positions = get_trajectory_positions(dynamic.trajectory, i);
		compute_positions(pos[0], prop->structure_data[0], atom_positions);
		compute_positions(pos[1], prop->structure_data[1], atom_positions);
		compute_positions(pos[2], prop->structure_data[2], atom_positions);
		compute_positions(pos[3], prop->structure_data[3], atom_positions);

		float32 angle = 0.f;
		for (int32 j = 0; j < count; j++) {
			angle += math::dihedral_angle(pos[0][j], pos[1][j], pos[2][j], pos[3][j]);
		}
		prop->data[i] = angle * scl;
	}

	TMP_FREE(tmp_data);

	return true;
}

static bool visualize_structures(const Property& prop, const MoleculeDynamic& dynamic) {
	if (prop.structure_data.count == 0) return false;

	if (prop.structure_data.count == 1) {
		auto pos = compute_positions(prop.structure_data[0], dynamic.molecule.atom_positions);
		for (const auto& p : pos) {
			immediate::draw_point(p);
		}
	}
	else {
		int32 count = (int32)prop.structure_data[0].structures.count;
		DynamicArray<vec3> pos_prev(count);
		DynamicArray<vec3> pos_next(count);

		pos_prev = compute_positions(prop.structure_data[0], dynamic.molecule.atom_positions);
		for (const auto& p : pos_prev) {
			immediate::draw_point(p);
		}
		for (int32 i = 1; i < prop.structure_data.count; i++) {
			pos_next = compute_positions(prop.structure_data[i], dynamic.molecule.atom_positions);
			for (const auto& p : pos_next) {
				immediate::draw_point(p);
			}
			for (int32 j = 0; j < pos_next.count; j++) {
				immediate::draw_line(pos_prev[j], pos_next[j]);
			}
			pos_prev = pos_next;
		}
	}

	return true;
}

void initialize() {
	ctx.property_func_entries.push_back({ HASH("distance"), compute_distance, visualize_structures });
	ctx.property_func_entries.push_back({ HASH("angle"),	compute_angle,	  visualize_structures });
	ctx.property_func_entries.push_back({ HASH("dihedral"), compute_dihedral, visualize_structures });

	ctx.structure_func_entries.push_back({ HASH("resname"), structure_match_resname });
	ctx.structure_func_entries.push_back({ HASH("resid"),   structure_match_resid });
	ctx.structure_func_entries.push_back({ HASH("residue"), structure_match_residue });
	ctx.structure_func_entries.push_back({ HASH("atom"),    structure_match_atom });
	ctx.structure_func_entries.push_back({ HASH("resatom"), structure_extract_resatom });
	ctx.structure_func_entries.push_back({ HASH("com"),		structure_apply_aggregation_strategy_com });
}

void shutdown() {

}

bool register_property_command(CString cmd_keyword, PropertyComputeFunc compute_func, PropertyVisualizeFunc visualize_func) {
	if (!cmd_keyword.data || cmd_keyword.count == 0) {
		printf("Error! Property command cannot be an empty string!\n");
		return false;
	}
	
	if (contains_whitespace(cmd_keyword)) {
		printf("Error! Property command cannot contain whitespace!\n");
		return false;
	}

	if (!compute_func) {
		printf("Error! Property command must have compute function!\n");
		return false;
	}
	
	ID hash = HASH(cmd_keyword);
	if (find_property_func_entry(hash) != nullptr) {
		printf("Error! Property command already registered!\n");
		return false;
	}

	ctx.property_func_entries.push_back({ hash, compute_func, visualize_func });
	return true;
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
	int32 num_frames = dynamic.trajectory.num_frames;
	
	for (auto p : ctx.properties) {
		auto& prop = *p;
		if (prop.data.count != num_frames) {
			prop.data.resize(num_frames);
		}
		prop.data.set_mem_to_zero();
		clear_histogram(&prop.hist);

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
			snprintf(prop.error_msg.beg(), prop.error_msg.MAX_LENGTH, "Could not recognize command '%s'", cmd.cstr());
			prop.valid = false;
			continue;
		}

		if (func(&prop, args, dynamic)) {
			constexpr int32 NUM_BINS = 64;
			prop.data_range = compute_range(prop.data);
			compute_histogram(&prop.hist, NUM_BINS, prop.data, prop.data_range.x, prop.data_range.y);
			prop.valid = true;
		}
	}

	return true;
}

void visualize(const MoleculeDynamic& dynamic) {
	for (auto p : ctx.properties) {
		if (!p->visualize) continue;
		CString cmd = extract_command(p->args);
		auto entry = find_property_func_entry(HASH(cmd));
		if (entry && entry->visualize_func) {
			entry->visualize_func(*p, dynamic);
		}
	}
}

Property* create_property(CString name, CString args) {
	Property* prop = (Property*)MALLOC(sizeof(Property));
	new(prop) Property();
	prop->name = name;
	prop->args = args;
	prop->valid = false;

	ctx.properties.push_back(prop);
	return prop;
}

void remove_property(Property* prop) {
	for (auto p : ctx.properties) {
		if (p == prop) {
			free_histogram(&p->hist);
			ctx.properties.swap_back_and_pop(&p);
		}
	}
}

Property* get_property(CString name) {
	for (auto p : ctx.properties) {
		if (compare(p->name, name)) return p;
	}
	return nullptr;
}

Property* get_property(int32 idx) {
	if (-1 < idx && idx < ctx.properties.count) return ctx.properties[idx];
	return nullptr;
}

int32 get_property_count() {
	return (int32)ctx.properties.count;
}

}  // namespace stats
