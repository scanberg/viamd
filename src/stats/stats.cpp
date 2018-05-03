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

typedef bool (*StructureFunc)(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule);

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
    DynamicArray<Property*> properties{};
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

static PropertyComputeFunc find_property_compute_func(CString cmd) {
	auto e = find_property_func_entry(HASH(cmd));
	if (e) return e->compute_func;

	return nullptr;
}

static PropertyVisualizeFunc find_property_visualize_func(CString cmd) {
	auto e = find_property_func_entry(HASH(cmd));
	if (e) return e->visualize_func;

	return nullptr;
}


static StructureFuncEntry* find_structure_func_entry(ID hash) {
	for (auto& e : ctx.structure_func_entries) {
		if (hash == e.hash) {
			return &e;
		}
	}
	return nullptr;
}

static StructureFunc find_structure_func(CString cmd) {
	auto e = find_structure_func_entry(HASH(cmd));
	if (e) return e->func;

	return nullptr;
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
    hist->bin_range = {0, 0};
    for (auto v : data) {
        int32 bin_idx = math::clamp((int32)((v - min_val) * scl), 0, num_bins - 1);
        hist->bins[bin_idx]++;
        hist->bin_range.y = math::max(hist->bin_range.y, hist->bins[bin_idx]);
    }
    hist->value_range = {min_val, max_val};
}

void clear_histogram(Histogram* hist) {
    ASSERT(hist);
    if (hist->bins.data) {
        memset(hist->bins.data, 0, hist->bins.count * sizeof(float));
    }
}

static void set_error_message(CString msg) { printf("%s\n", msg.beg()); }

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
                data->structures.push_back({res.beg_atom_idx, res.end_atom_idx});
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
                data->structures.push_back({res.beg_atom_idx, res.end_atom_idx});
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
        data->structures.push_back({res.beg_atom_idx, res.end_atom_idx});
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
        data->structures.push_back({idx, idx + 1});
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
            s = {new_idx, new_idx + 1};
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
    if (positions.count == 0) return {0, 0, 0};
    if (positions.count == 1) return positions[0];

    vec3 com{0};
    for (const auto& p : positions) {
        com += p;
    }

    return com / (float)positions.count;
}

static inline int32 structure_index_count(Structure s) { return s.end_idx - s.beg_idx; }

static inline int32 structures_index_count(Array<const Structure> structures) {
    int32 count = 0;
    for (const auto& s : structures) {
        count += structure_index_count(s);
    }
    return count;
}

static inline int32 compute_positions(Array<vec3> positions, Structure structure, AggregationStrategy strategy, Array<const vec3> atom_positions) {
    switch (strategy) {
        case NONE: {
            int32 count = structure_index_count(structure);
			ASSERT(positions.count >= count);
            for (int32 i = 0; i < count; i++) {
                positions[i] = atom_positions[structure.beg_idx + i];
            }
			return count;
        }
        case COM: {
			ASSERT(positions.count >= 1);
            vec3 com(0);
            for (int32 i = structure.beg_idx; i < structure.end_idx; i++) {
                com += atom_positions[i];
            }
            positions[0] = com / (float)structure_index_count(structure);
			return 1;
        }
    }
	return 0;
}

/*
static inline void compute_positions(Array<vec3> dst, const StructureData& data, Array<const vec3> atom_positions) {
	ASSERT(dst.count >= (data.strategy == COM) ? data.structures.count : structures_index_count(data.structures));

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
*/

static Array<const vec3> extract_positions(Structure structure, Array<const vec3> atom_positions) {
	return atom_positions.sub_array(structure.beg_idx, structure.end_idx - structure.beg_idx);
}

static inline float multi_distance(Array<const vec3> pos_a, Array<const vec3> pos_b) {
	if (pos_a.count == 0 || pos_b.count == 0) return 0.f;
	else if (pos_a.count == 1 && pos_b.count == 1) return math::distance(pos_a[0], pos_b[0]);
	else {
		float dist = 0.f;
		float count = 0.f;

		if (pos_a.count > 1) {
			for (const auto& p : pos_a) dist += math::distance(p, pos_b[0]);
			count = (float)pos_a.count;
		}
		else if (pos_b.count > 1) {
			for (const auto& p : pos_b) dist += math::distance(pos_a[0], p);
			count = (float)pos_b.count;
		}
		else {
			// ERROR
			return 0.f;
		}
		return dist / count;
	}
	return 0.f;
}

static inline float multi_angle(Array<const vec3> pos_a, Array<const vec3> pos_b, Array<const vec3> pos_c) {
	if (pos_a.count == 0 || pos_b.count == 0 || pos_c.count == 0) return 0.f;
	else if (pos_a.count == 1 && pos_b.count == 1 && pos_c.count == 1) return math::angle(pos_a[0], pos_b[0], pos_c[0]);
	else {
		float angle = 0.f;
		float count;
		if (pos_a.count > 1.f) {
			for (const auto& p : pos_a) angle += math::angle(p, pos_b[0], pos_c[0]);
			count = (float)pos_a.count;
		}
		else if (pos_b.count > 1.f) {
			for (const auto& p : pos_b) angle += math::angle(pos_a[0], p, pos_c[0]);
			count = (float)pos_b.count;
		}
		else if (pos_c.count > 1.f) {
			for (const auto& p : pos_c) angle += math::angle(pos_a[0], pos_b[0], p);
			count = (float)pos_c.count;
		}
		else {
			// ERROR
			return 0.f;
		}
		return angle / count;
	}
}

static inline float multi_dihedral(Array<const vec3> pos_a, Array<const vec3> pos_b, Array<const vec3> pos_c, Array<const vec3> pos_d) {
	if (pos_a.count == 0 || pos_b.count == 0 || pos_c.count == 0 || pos_d.count == 0) return 0.f;
	else if (pos_a.count == 1 && pos_b.count == 1 && pos_c.count == 1 && pos_d.count == 1) return math::dihedral_angle(pos_a[0], pos_b[0], pos_c[0], pos_d[0]);
	else {
		float angle = 0.f;
		float count;
		if (pos_a.count > 1.f) {
			for (const auto& p : pos_a) angle += math::dihedral_angle(p, pos_b[0], pos_c[0], pos_d[0]);
			count = (float)pos_a.count;
		}
		else if (pos_b.count > 1.f) {
			for (const auto& p : pos_b) angle += math::dihedral_angle(pos_a[0], p, pos_c[0], pos_d[0]);
			count = (float)pos_b.count;
		}
		else if (pos_c.count > 1.f) {
			for (const auto& p : pos_c) angle += math::dihedral_angle(pos_a[0], pos_b[0], p, pos_d[0]);
			count = (float)pos_c.count;
		}
		else if (pos_c.count > 1.f) {
			for (const auto& p : pos_d) angle += math::dihedral_angle(pos_a[0], pos_b[0], pos_c[0], p);
			count = (float)pos_d.count;
		}
		else {
			// ERROR
			return 0.f;
		}
		return angle / count;
	}
}

static bool compute_distance(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic) {
    ASSERT(prop);
    if (args.count != 2) {
        set_error_message("distance expects 2 arguments");
        return false;
    }

	// Allocate structure data to hold two data for two arguments
    init_structure_data(&prop->structure_data, 2);

	// Extract structures for the arguments
	if (!extract_args_structures(prop->structure_data, args, dynamic.molecule)) {
		return false;
	}

	// Sync the count of structures for between arguments
	if (!sync_structure_count(prop->structure_data)) {
		return false;
	}
    const int32 structure_count = (int32)prop->structure_data[0].structures.count;

    const float32 scl = 1.f / (float32)structure_count;
	Array<const vec3> pos[2];
	vec3 com[2];
    for (int32 i = 0; i < dynamic.trajectory.num_frames; i++) {
		float val = 0.f;
		for (int32 j = 0; j < structure_count; j++) {
			pos[0] = extract_positions(prop->structure_data[0].structures[j], get_trajectory_positions(dynamic.trajectory, i));
			pos[1] = extract_positions(prop->structure_data[1].structures[j], get_trajectory_positions(dynamic.trajectory, i));
			if (prop->structure_data[0].strategy == COM) {
				com[0] = compute_com(pos[0]);
				pos[0] = { &com[0], 1 };
			}
			if (prop->structure_data[1].strategy == COM) {
				com[1] = compute_com(pos[1]);
				pos[1] = { &com[1], 1 };
			}

			val += multi_distance(pos[0], pos[1]);
		}

        prop->data[i] = val * scl;
    }

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

	// Sync the count of structures for between arguments
	if (!sync_structure_count(prop->structure_data)) {
		return false;
	}
	const int32 structure_count = (int32)prop->structure_data[0].structures.count;

	const float32 scl = 1.f / (float32)structure_count;
	Array<const vec3> pos[3];
	vec3 com[3];
	for (int32 i = 0; i < dynamic.trajectory.num_frames; i++) {
		float val = 0.f;
		for (int32 j = 0; j < structure_count; j++) {
			pos[0] = extract_positions(prop->structure_data[0].structures[j], get_trajectory_positions(dynamic.trajectory, i));
			pos[1] = extract_positions(prop->structure_data[1].structures[j], get_trajectory_positions(dynamic.trajectory, i));
			pos[2] = extract_positions(prop->structure_data[2].structures[j], get_trajectory_positions(dynamic.trajectory, i));

			if (prop->structure_data[0].strategy == COM) {
				com[0] = compute_com(pos[0]);
				pos[0] = { &com[0], 1 };
			}
			if (prop->structure_data[1].strategy == COM) {
				com[1] = compute_com(pos[1]);
				pos[1] = { &com[1], 1 };
			}
			if (prop->structure_data[2].strategy == COM) {
				com[2] = compute_com(pos[2]);
				pos[2] = { &com[2], 1 };
			}

			val += multi_angle(pos[0], pos[1], pos[2]);
		}

		prop->data[i] = val * scl;
	}

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

	// Sync the count of structures for between arguments
	if (!sync_structure_count(prop->structure_data)) {
		return false;
	}
	const int32 structure_count = (int32)prop->structure_data[0].structures.count;

	const float32 scl = 1.f / (float32)structure_count;
	Array<const vec3> pos[4];
	vec3 com[4];
	for (int32 i = 0; i < dynamic.trajectory.num_frames; i++) {
		float val = 0.f;
		for (int32 j = 0; j < structure_count; j++) {
			pos[0] = extract_positions(prop->structure_data[0].structures[j], get_trajectory_positions(dynamic.trajectory, i));
			pos[1] = extract_positions(prop->structure_data[1].structures[j], get_trajectory_positions(dynamic.trajectory, i));
			pos[2] = extract_positions(prop->structure_data[2].structures[j], get_trajectory_positions(dynamic.trajectory, i));
			pos[3] = extract_positions(prop->structure_data[3].structures[j], get_trajectory_positions(dynamic.trajectory, i));

			if (prop->structure_data[0].strategy == COM) {
				com[0] = compute_com(pos[0]);
				pos[0] = { &com[0], 1 };
			}
			if (prop->structure_data[1].strategy == COM) {
				com[1] = compute_com(pos[1]);
				pos[1] = { &com[1], 1 };
			}
			if (prop->structure_data[2].strategy == COM) {
				com[2] = compute_com(pos[2]);
				pos[2] = { &com[2], 1 };
			}
			if (prop->structure_data[3].strategy == COM) {
				com[3] = compute_com(pos[3]);
				pos[3] = { &com[3], 1 };
			}

			val += multi_dihedral(pos[0], pos[1], pos[2], pos[3]);
		}

		prop->data[i] = val * scl;
	}
    return true;
}

static bool visualize_structures(const Property& prop, const MoleculeDynamic& dynamic) {
    if (prop.structure_data.count == 0) return false;

    if (prop.structure_data.count == 1) {
		for (const auto& s : prop.structure_data[0].structures) {
			Array<const vec3> pos = extract_positions(s, dynamic.molecule.atom_positions);
			if (prop.structure_data[0].strategy == COM) {
				immediate::draw_point(compute_com(pos));
			}
			else {
				for (const auto& p : pos) {
					immediate::draw_point(p);
				}
			}
		}
    } else {
        int32 count = (int32)prop.structure_data[0].structures.count;
		Array<const vec3> pos_prev;
		Array<const vec3> pos_next;
		vec3 com_prev(0);
		vec3 com_next(0);

		for (int32 i = 0; i < count; i++) {
			pos_prev = extract_positions(prop.structure_data[0].structures[i], dynamic.molecule.atom_positions);
			if (prop.structure_data[0].strategy == COM) {
				com_prev = compute_com(pos_prev);
				pos_prev = { &com_prev, 1 };
			}
			for (const auto& p : pos_prev) {
				immediate::draw_point(p);
			}
			for (int32 j = 1; j < prop.structure_data.count; j++) {
				pos_next = extract_positions(prop.structure_data[j].structures[i], dynamic.molecule.atom_positions);
				if (prop.structure_data[j].strategy == COM) {
					com_next = compute_com(pos_next);
					pos_next = { &com_next, 1 };
				}
				for (const auto& p : pos_next) {
					immediate::draw_point(p);
				}
				if (pos_prev.count == 1 && pos_next.count == 1) {
					immediate::draw_line(pos_prev[0], pos_next[0]);
				}
				if (pos_prev.count > 1) {
					for (int32 k = 0; k < pos_prev.count; k++) {
						immediate::draw_line(pos_prev[k], pos_next[0]);
					}
				}
				else if (pos_next.count > 1) {
					for (int32 k = 0; k < pos_next.count; k++) {
						immediate::draw_line(pos_prev[0], pos_next[k]);
					}
				}

				pos_prev = pos_next;
				com_prev = com_next;
			}
		}
    }

    return true;
}

void initialize() {
    ctx.property_func_entries.push_back({HASH("distance"), compute_distance, visualize_structures});
    ctx.property_func_entries.push_back({HASH("angle"), compute_angle, visualize_structures});
    ctx.property_func_entries.push_back({HASH("dihedral"), compute_dihedral, visualize_structures});

    ctx.structure_func_entries.push_back({HASH("resname"), structure_match_resname});
    ctx.structure_func_entries.push_back({HASH("resid"), structure_match_resid});
    ctx.structure_func_entries.push_back({HASH("residue"), structure_match_residue});
    ctx.structure_func_entries.push_back({HASH("atom"), structure_match_atom});
    ctx.structure_func_entries.push_back({HASH("resatom"), structure_extract_resatom});
    ctx.structure_func_entries.push_back({HASH("com"), structure_apply_aggregation_strategy_com});
}

void shutdown() {}

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

    ctx.property_func_entries.push_back({hash, compute_func, visualize_func});
    return true;
}

static CString extract_parentheses(CString str) {
    const char* beg = str.beg();

    while (beg < str.end() && *beg != '(') beg++;
    if (beg < str.end()) beg++;
    if (beg == str.end()) return {beg, str.end()};

    const char* end = beg;
    int count = 1;
    while (end != str.end()) {
        if (*end == '(')
            count++;
        else if (*end == ')' && --count == 0)
            break;
        end++;
    }

    return {beg, end};
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
                args.push_back(trim({beg, end}));
                beg = end + 1;
            }
        }
        end++;
    }
    if (beg != end) args.push_back(trim({beg, end}));

    return args;
}

static CString extract_command(CString str) {
    str = trim(str);
    const char* ptr = str.beg();
    while (ptr != str.end() && *ptr != '(' && !isspace(*ptr)) ptr++;
    return {str.beg(), ptr};
}

bool extract_structures(StructureData* data, CString arg, const MoleculeStructure& molecule) {
    CString cmd = extract_command(arg);
    auto func = find_structure_func(cmd);

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

    int32 max_count = 0;
    for (const auto& s : data) {
		int32 count = (int32)s.structures.count;
        if (count == 0) {
            set_error_message("One argument did not match any structures");
            return false;
        }

		max_count = math::max(max_count, count);
        if (count > 1 && max_count > 1 && count != max_count) {
			set_error_message("Multiple structures found for more than one argument, but the structure count did not match");
            return false;
        }
    }

    return true;
}

bool sync_structure_count(Array<StructureData> data) {
	int32 max_count = 0;
	for (const auto& s : data) {
		max_count = math::max(max_count, (int32)s.structures.count);
	}

	// Extend and copy data from first element to match multi_count
	for (auto& s : data) {
		while (s.structures.count < max_count) {
			s.structures.push_back(s.structures.front());
		}
	}

	if (data.count > 0) {
		for (int32 i = 0; i < data[0].structures.count; i++) {
			int32 max_structure_count = 0;
			for (const auto& s : data) {
				int32 c = structure_index_count(s.structures[i]);
				max_structure_count = math::max(max_structure_count, c);
				if (c > 1 && c != max_structure_count) {
					set_error_message("Structures matched has different sizes in different arguments, this is not supported");
					return false;
				}
			}
		}
	}

	return true;
}


bool balanced_parentheses(CString str) {
    int count = 0;
    const char* ptr = str.beg();
    while (ptr != str.end()) {
        if (*ptr == '(')
            count++;
        else if (*ptr == ')')
            count--;
        ptr++;
    }
    return count == 0;
}

static Range compute_range(Array<float> data) {
    if (data.count == 0) {
        return {0, 0};
    }
    Range range{FLT_MAX, -FLT_MAX};
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
            if (*end == '(')
                count++;
            else if (*end == ')')
                count--;
            else if (isspace(*end) && count == 0) {
                args.push_back(trim({beg, end}));
                beg = end + 1;
            }
            end++;
        }
        if (beg != end) args.push_back(trim({beg, end}));

        if (args.count == 0) {
            prop.valid = false;
            continue;
        }

        CString cmd = args[0];
        args = args.sub_array(1);

		auto func = find_property_compute_func(cmd);
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
    new (prop) Property();
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

int32 get_property_count() { return (int32)ctx.properties.count; }

}  // namespace stats
