#include "stats.h"
#include <core/math_utils.h>
#include <core/hash.h>
#include <mol/molecule_dynamic.h>
#include <mol/molecule_utils.h>

#define COMPUTE_ID(x) (hash::crc64(x))

namespace stats {

struct PropertyCommand {
    ID id = INVALID_ID;
	CString keyword = "";
    PropertyComputeFunc func = nullptr;
    Range val_range {};
    PropertyType type;
    bool periodic = false;
    CString unit = "";
};

struct GroupCommand {
    ID id = INVALID_ID;
	CString keyword = "";
    StructureExtractFunc func = nullptr;
};

struct Property {
    ID id = INVALID_ID;
    ID data_avg_id = INVALID_ID;
    ID data_beg_id = INVALID_ID;
    int32 data_count = 0;

	float filter_min = 0.f;
	float filter_max = 1.f;

    StringBuffer<32> name_buf{};
    StringBuffer<64> args_buf{};
    CString name{};
    CString cmd{};
    CString args{};
    ID cmd_id = INVALID_ID;
    bool valid = true;
};
    
struct PropertyData {
    ID id = INVALID_ID;
    ID property_id = INVALID_ID;
	ID instance_id = INVALID_ID;
    Array<float> data {};
};

struct Group {
    ID id = INVALID_ID;
	ID instance_beg_id = INVALID_ID;
	int32 instance_count = 0;

    StringBuffer<32> name_buf{};
    StringBuffer<64> args_buf{};
    CString name{};
    CString cmd{};
    CString args{};
    ID cmd_id = INVALID_ID;
    bool valid;
};
    
struct GroupInstance {
    ID id = INVALID_ID;
    ID group_id = INVALID_ID;
    Structure structure {};
};

struct StatisticsContext {
    ID next_id = 1;
    DynamicArray<ID> free_ids {};
    DynamicArray<String> string_buffer {};

    DynamicArray<Property> properties{};
    DynamicArray<PropertyData> property_data{};
    DynamicArray<Group> groups{};
	DynamicArray<GroupInstance> group_instances{};

    DynamicArray<PropertyCommand> property_commands {};
    DynamicArray<GroupCommand> group_commands {};
};

static StatisticsContext ctx;

static bool compute_atomic_distance(float* data, const Array<CString> args, const MoleculeDynamic& dynamic, Structure group_struct) {
	if (args.count != 2) return false;
    if (group_struct.beg_atom_idx == group_struct.end_atom_idx) return false;

	auto int_a = to_int32(args[0]);
	auto int_b = to_int32(args[1]);

	if (!int_a.success || !int_b.success) return false;
	int32 atom_a = group_struct.beg_atom_idx + int_a;
	int32 atom_b = group_struct.beg_atom_idx + int_b;

    int32 count = dynamic.trajectory->num_frames;
    for (int32 i = 0; i < count; i++) {
        vec3 pos_a = dynamic.trajectory->frame_buffer[i].atom_positions[atom_a];
        vec3 pos_b = dynamic.trajectory->frame_buffer[i].atom_positions[atom_b];
        data[i] = math::distance(pos_a, pos_b);
    }

    return true;
}

static bool compute_atomic_angle(float* data, const Array<CString> args, const MoleculeDynamic& dynamic, Structure group_struct) {
	if (args.count != 3) return false;
    if (group_struct.beg_atom_idx == group_struct.end_atom_idx) return false;

	auto int_a = to_int32(args[0]);
	auto int_b = to_int32(args[1]);
	auto int_c = to_int32(args[2]);

	if (!int_a.success || !int_b.success || !int_c.success) return false;
	int32 atom_a = group_struct.beg_atom_idx + int_a;
	int32 atom_b = group_struct.beg_atom_idx + int_b;
	int32 atom_c = group_struct.beg_atom_idx + int_c;

	int32 count = dynamic.trajectory->num_frames;
	for (int32 i = 0; i < count; i++) {
		vec3 pos_a = dynamic.trajectory->frame_buffer[i].atom_positions[atom_a];
		vec3 pos_b = dynamic.trajectory->frame_buffer[i].atom_positions[atom_b];
		vec3 pos_c = dynamic.trajectory->frame_buffer[i].atom_positions[atom_c];

		data[i] = math::angle(pos_a - pos_b, pos_c - pos_b);
	}

	return true;
}

static bool compute_atomic_dihedral(float* data, const Array<CString> args, const MoleculeDynamic& dynamic, Structure group_struct) {
	if (args.count != 4) return false;
    if (group_struct.beg_atom_idx == group_struct.end_atom_idx) return false;

	auto int_a = to_int32(args[0]);
	auto int_b = to_int32(args[1]);
	auto int_c = to_int32(args[2]);
	auto int_d = to_int32(args[2]);

	if (!int_a.success || !int_b.success || !int_c.success || !int_d.success) return false;
    int32 atom_a = group_struct.beg_atom_idx + int_a;
    int32 atom_b = group_struct.beg_atom_idx + int_b;
    int32 atom_c = group_struct.beg_atom_idx + int_c;
    int32 atom_d = group_struct.beg_atom_idx + int_d;

	int32 count = dynamic.trajectory->num_frames;
	for (int32 i = 0; i < count; i++) {
		vec3 pos_a = dynamic.trajectory->frame_buffer[i].atom_positions[atom_a];
		vec3 pos_b = dynamic.trajectory->frame_buffer[i].atom_positions[atom_b];
		vec3 pos_c = dynamic.trajectory->frame_buffer[i].atom_positions[atom_c];
		vec3 pos_d = dynamic.trajectory->frame_buffer[i].atom_positions[atom_d];

		data[i] = dihedral_angle(pos_a, pos_b, pos_c, pos_d);
	}

	return true;
}

static DynamicArray<Structure> match_by_resname(const Array<CString> args, const MoleculeStructure& mol) {
    DynamicArray<Structure> result;
    for (const auto& res : mol.residues) {
        for (const auto& arg : args) {
            if (compare(res.name, arg)) {
                result.push_back({res.beg_atom_idx, res.end_atom_idx});
            }
        }
    }
    return result;
}

void initialize() {
    ctx.property_commands.push_back({ COMPUTE_ID("dist"),     "dist",	  compute_atomic_distance, {0, FLT_MAX}, INTRA, false, "å" });
	ctx.property_commands.push_back({ COMPUTE_ID("bond"),     "bond",     compute_atomic_distance, {0, FLT_MAX}, INTRA, false, "å" });
	ctx.property_commands.push_back({ COMPUTE_ID("angle"),    "angle",    compute_atomic_angle,    {0, math::PI}, INTRA, true, "°" });
	ctx.property_commands.push_back({ COMPUTE_ID("dihedral"), "dihedral", compute_atomic_dihedral, {-math::PI, math::PI}, INTRA, true, "°" });

    //ctx.group_commands.push_back({ COMPUTE_ID("resid"), match_by_resname});
    ctx.group_commands.push_back({ COMPUTE_ID("resname"), "resname", match_by_resname});
}

void shutdown() {

}

void clear() {
    for (auto str : ctx.string_buffer) {
        FREE(&str.data);
    }
    ctx.string_buffer.clear();
    ctx.properties.clear();
    ctx.property_data.clear();
    ctx.groups.clear();
	ctx.group_instances.clear();

    //ctx.property_commands.clear();
    //ctx.group_commands.clear();
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
    ctx.string_buffer.push_back({data, str.count});
	copy(ctx.string_buffer.back(), str);
    return ctx.string_buffer.back();
}

static void free_string(CString str) {
    for (auto int_str : ctx.string_buffer)  {
        if (int_str.count == str.count && compare(int_str, str)) {
            ctx.string_buffer.remove(&int_str);
            return;
        }
    }
}

// HISTOGRAMS
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

	hist->bins.resize(num_bins);
	memset(hist->bins.data, 0, hist->bins.count * sizeof(float));

	const float scl = num_bins / (max_val - min_val);
	for (auto v : data) {
		int32 bin_idx = math::clamp((int32)((v - min_val) * scl), 0, num_bins - 1);
		hist->bins[bin_idx]++;
	}
}

bool compute_stats(const MoleculeDynamic& dynamic) {
    if (!dynamic.molecule) {
        printf("ERROR! Computing statistics: molecule is not set");
        return false;
    }
    if (!dynamic.trajectory) {
        printf("ERROR! Computing statistics: trajectory is not set");
        return false;
    }

    // Find uninitialized groups and compute their instances
    for (auto& group : ctx.groups) {
		if (group.valid && group.instance_count == 0) {
            GroupCommand* group_cmd = find_id(ctx.group_commands, group.cmd_id);
			ASSERT(group_cmd);
            DynamicArray<CString> args = ctokenize(group.args);
            auto matching_structures = group_cmd->func(args, *dynamic.molecule);

			if (matching_structures.count == 0) {
				printf("WARNING! group '%s' did not match any structures.\n", group.name.beg());
				continue;
			}

			for (int32 i = 0; i < matching_structures.count; i++) {
				char buf[64];
				snprintf(buf, 64, "%s.%i", group.name.beg(), i);
				GroupInstance instance;
				instance.id = COMPUTE_ID(buf);
				instance.group_id = group.id;
				instance.structure = matching_structures[i];
				ctx.group_instances.push_back(instance);

				if (i == 0) {
					group.instance_beg_id = instance.id;
					group.instance_count = 0;
				}
				group.instance_count++;
			}
        }
    }

    int32 count = dynamic.trajectory->num_frames;
    for (auto& prop : ctx.properties) {
        if (prop.valid && prop.data_beg_id == INVALID_ID) {
            // NEED TO COMPUTE PROPERTY DATA
            PropertyCommand* prop_cmd = find_id(ctx.property_commands, prop.cmd_id);
			ASSERT(prop_cmd);
            DynamicArray<CString> args = ctokenize(prop.args);

			if (prop_cmd->type == INTRA) {
				if (args == 0) {
					printf("WARNING! Property '%s': Missing arguments!", prop.name.beg());
					continue;
				}
                CString group_name = args[0];
				ID group_id = get_group(group_name);
				Group* group = find_id(ctx.groups, group_id);

                args = args.sub_array(1);

				if (group_id == INVALID_ID) {
                    StringBuffer<32> buf = group_name;
					printf("WARNING! Property '%s': could not find group with name '%s'\n", prop.name.beg(), buf.buffer);
					continue;
				}

				// DATA AVERAGE
				StringBuffer<64> prop_data_name;
				snprintf(prop_data_name.beg(), 64, "%s.%s.avg", group->name.beg(), prop.name.beg());

				PropertyData prop_avg_data;
				prop_avg_data.id = COMPUTE_ID(prop_data_name.operator CString());
				prop_avg_data.instance_id = INVALID_ID;
				prop_avg_data.property_id = prop.id;
				prop_avg_data.data = { (float*)CALLOC(count, sizeof(float)), count };
				prop.data_avg_id = prop_avg_data.id;
				ctx.property_data.push_back(prop_avg_data);

				// DATA
				for (int32 i = 0; i < group->instance_count; i++) {
					ID instance_id = get_group_instance(group_id, i);
					auto inst = find_id(ctx.group_instances, instance_id);
					float* data = (float*)CALLOC(count, sizeof(float));
					prop_cmd->func(data, args, dynamic, inst->structure);

					for (int32 j = 0; j < count; j++) {
						prop_avg_data.data[j] += data[j] / (float)group->instance_count;
					}

					snprintf(prop_data_name.beg(), 64, "%s.%s.%i", group->name.beg(), prop.name.beg(), i);

					PropertyData prop_data;
					prop_data.id = COMPUTE_ID(prop_data_name);
					prop_data.instance_id = instance_id;
					prop_data.property_id = prop.id;
					prop_data.data = { data, count };

					ctx.property_data.push_back(prop_data);
					if (i == 0) {
						prop.data_beg_id = prop_data.id;
						prop.data_count = 0;
					}
					prop.data_count++;
				}
			}
			else if (prop_cmd->type == INTER) {
                // @TODO: Implement
			}

        }
    }

    return true;
}

void register_property_command(CString command, PropertyCommandDescriptor desc) {
    ID id = COMPUTE_ID(command);
    auto cmd = find_id(ctx.property_commands, id);
    if (cmd != nullptr) {
        printf("ERROR: PROPERTY COMMAND %s ALREADY REGISTERED!", command.beg());
        return;
    }

	auto keyword = alloc_string(command);
    auto unit = alloc_string(desc.unit);
    ctx.property_commands.push_back({id, keyword, desc.compute_function, desc.val_range, desc.type, desc.periodic, unit});
}

void register_group_command(CString command, StructureExtractFunc func) {
    ID id = COMPUTE_ID(command);
    auto cmd = find_id(ctx.group_commands, id);
    if (cmd != nullptr) {
        printf("ERROR: GROUP COMMAND %s ALREADY REGISTERED!", command.beg());
        return;
    }

	auto keyword = alloc_string(command);
    ctx.group_commands.push_back({id, keyword, func});
}

int32 get_property_command_count() {
	return (int32)ctx.property_commands.count;
}
CString get_property_command_keyword(int32 idx) {
	return ctx.property_commands[idx].keyword;
}

int32 get_group_command_count() {
	return (int32)ctx.group_commands.count;
}

CString get_group_command_keyword(int32 idx) {
	return ctx.group_commands[idx].keyword;
}
    
bool validate_group(Group* group) {
    ASSERT(group);
    group->valid = false;
    
    group->name = CString(group->name_buf.beg());
    auto tokens = ctokenize(group->args_buf);
    
    group->cmd_id = INVALID_ID;
    group->cmd = {};
    group->args = {};
    
    if (tokens.count > 0) {
        group->cmd = tokens[0];
        if (tokens.count > 1) {
            group->args = {tokens[1].beg(), tokens.back().end()};
        }
    }
    
    if (group->name) {
        for (const auto& g : ctx.groups) {
            if (g.id != group->id && compare(g.name, group->name)) {
                return false;
            }
        }
    } else {
        return false;
    }
    
    if (group->cmd) {
        for (const auto& c : ctx.group_commands) {
            if (compare(c.keyword, group->cmd)) {
                group->cmd_id = c.id;
                break;
            }
        }
        if (group->cmd_id == INVALID_ID) {
            return false;
        }
    } else {
        return false;
    }
    
    group->valid = true;
    return true;
}
    
bool validate_group(ID group_id) {
    auto group = find_id(ctx.groups, group_id);
    if (!group) return false;
    return validate_group(group);
}

ID create_group(CString name, CString cmd_and_args) {
    Group group;
    group.id = create_id();
    group.instance_beg_id = INVALID_ID;
    group.instance_count = 0;
    group.name_buf = name;
    group.args_buf = cmd_and_args;
    group.valid = false;
    
    ctx.groups.push_back(group);
    validate_group(&group);
    
    return group.id;
}

void remove_group(ID group_id) {
	Group* group = find_id(ctx.groups, group_id);
	if (!group) {
		printf("ERROR: COULD NOT FIND GROUP!\n");
		return;
	}

	free_string(group->name);
	free_string(group->args);

    ctx.groups.remove(group);
}

ID get_group(CString name) {
    for (const auto& g : ctx.groups) {
        if (compare(name, g.name)) return g.id;
    }
    return INVALID_ID;
}

ID get_group(int32 idx) {
    if (idx < ctx.groups.count) {
        return ctx.groups[idx].id;
    }
    return INVALID_ID;
}

int32 get_group_count() {
    return (int32)ctx.groups.count;
}
    
StringBuffer<32>* get_group_name_buf(ID group_id) {
    auto group = find_id(ctx.groups, group_id);
    if (group) {
        return &group->name_buf;
    }
    return nullptr;
}

StringBuffer<64>* get_group_args_buf(ID group_id) {
    auto group = find_id(ctx.groups, group_id);
    if (group) {
        return &group->args_buf;
    }
    return nullptr;
}
    
bool get_group_valid(ID group_id) {
    auto group = find_id(ctx.groups, group_id);
    if (group) {
        return group->valid;
    } else {
        printf("ERROR! Could not find group id!\n");
    }
    return false;
}
    
CString get_group_name(ID group_id) {
    auto group = find_id(ctx.groups, group_id);
    if (group) {
        return group->name;
    } else {
        printf("ERROR! Could not find group id!\n");
    }
    return {};
}
    
bool set_group_name(ID group_id, CString name) {
    if (name.count == 0) return false;
    
    Group* group = find_id(ctx.groups, group_id);
    if (!group) {
        printf("ERROR! Could not find group id!\n");
        return false;
    }
    
    if (group->name) free_string(group->name);
    group->name = alloc_string(name);
    
    return true;
}

ID get_group_instance(ID group_id, int32 idx) {
	int32 count = 0;
	for (const auto &inst : ctx.group_instances) {
		if (inst.group_id == group_id) {
			if (count == idx) return inst.id;
			count++;
		}
	}
	return INVALID_ID;
}

int32 get_group_instance_count(ID group_id) {
	auto group = find_id(ctx.groups, group_id);
	if (group) {
		return group->instance_count;
	}
	return 0;
}

ID get_property(CString name) {
    for (const auto &p : ctx.properties) {
        if (compare(p.name, name)) return p.id;
    }
    return INVALID_ID;
}

ID get_property(int32 idx) {
	if (-1 < idx && idx < ctx.properties.count) return ctx.properties[idx].id;
    return INVALID_ID;
}

int32 get_property_count() {
    return (int32)ctx.properties.count;
}
    
bool validate_property(Property* prop) {
    ASSERT(prop);
    prop->valid = false;
    
    prop->name = CString(prop->name_buf.beg());
    auto tokens = ctokenize(prop->args_buf);
    
    prop->cmd_id = INVALID_ID;
    prop->cmd = {};
    prop->args = {};

    if (tokens.count > 0) {
        prop->cmd = tokens[0];
        if (tokens.count > 1) {
            prop->args = {tokens[1].beg(), tokens.back().end()};
        }
    }
    
    if (prop->name) {
        for (const auto& p : ctx.properties) {
            if (p.id != prop->id && compare(p.name, prop->name)) {
                return false;
            }
        }
    } else {
        return false;
    }
    
    if (prop->cmd) {
        for (const auto& c : ctx.property_commands) {
            if (compare(c.keyword, prop->cmd)) {
                prop->cmd_id = c.id;
                break;
            }
        }
        if (prop->cmd_id == INVALID_ID) {
            return false;
        }
    } else {
        return false;
    }

    prop->valid = true;
    return true;
}
    
bool validate_property(ID prop_id) {
    auto prop = find_id(ctx.properties, prop_id);
    if (!prop) return false;
    return validate_property(prop);
}

ID create_property(CString name, CString cmd_and_args) {
    Property prop;
    prop.id = create_id();
    prop.data_beg_id = INVALID_ID;
    prop.data_count = 0;
    prop.name_buf = name;
    prop.args_buf = cmd_and_args;
    prop.valid = false;
    
    ctx.properties.push_back(prop);
    validate_property(&prop);
    
	return prop.id;
}

void remove_property(ID prop_id) {
    Property* prop = find_id(ctx.properties, prop_id);
    if (prop == nullptr) {
        printf("ERROR: PROPERTY NOT FOUND\n");
        return;
    }

    for (PropertyData* pd = ctx.property_data.beg(); pd != ctx.property_data.end(); pd++) {
        if (pd->property_id == prop_id) ctx.property_data.remove(pd);
    }

    ctx.properties.remove(prop);
}

Array<float> get_property_data(ID prop_id, int32 idx) {
    auto prop = find_id(ctx.properties, prop_id);
    if (prop && prop->data_beg_id != INVALID_ID && idx < prop->data_count) {
        auto prop_data = find_id(ctx.property_data, prop->data_beg_id);
        return prop_data[idx].data;
    }
    return {};
}

Array<float> get_property_avg_data(ID prop_id) {
    auto prop = find_id(ctx.properties, prop_id);
    if (prop && prop->data_avg_id != INVALID_ID) {
        auto prop_data = find_id(ctx.property_data, prop->data_avg_id);
        return prop_data->data;
    }
    return {};
}

int32 get_property_data_count(ID prop_id) {
    auto prop = find_id(ctx.properties, prop_id);
    if (prop) {
        return prop->data_count;
    }
	return 0;
}
    
StringBuffer<32>* get_property_name_buf(ID prop_id) {
    auto prop = find_id(ctx.properties, prop_id);
    if (prop) {
        return &prop->name_buf;
    }
    return nullptr;
}
    
StringBuffer<64>* get_property_args_buf(ID prop_id) {
    auto prop = find_id(ctx.properties, prop_id);
    if (prop) {
        return &prop->args_buf;
    }
    return nullptr;
}

CString get_property_name(ID prop_id) {
	auto prop = find_id(ctx.properties, prop_id);
	if (prop) {
		return prop->name;
	}
	return {};
}
    
bool get_property_valid(ID prop_id) {
    auto prop = find_id(ctx.properties, prop_id);
    if (prop) {
        return prop->valid;
    }
    return false;
}

}  // namespace stats
