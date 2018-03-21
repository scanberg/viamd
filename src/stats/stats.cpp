#include "stats.h"
#include <core/math_utils.h>
#include <core/hash.h>
#include <mol/molecule_dynamic.h>
#include <mol/molecule_utils.h>

#define COMPUTE_ID(x) (hash::crc64(x))

namespace stats {

struct PropertyCommand {
    ID id;
    PropertyType type;
    PropertyComputeFunc func;
};

struct GroupCommand {
    ID id;
    ResidueMatchFunc func;
};

struct Property {
    ID id = INVALID_ID;
    ID group_id = INVALID_ID;
    ID data_avg_id = INVALID_ID;
    ID data_beg_id = INVALID_ID;
    int32 data_count = 0;

    ID cmd_id = INVALID_ID;
    CString name;
    CString args;

    PropertyType type = PropertyType::FLOAT32;
};

struct PropertyData {
    ID id = INVALID_ID;
    ID group_id = INVALID_ID;
    ID property_id = INVALID_ID;

    int32 residue_idx;
    void* data;
    int32 count;
};

struct Group {
    ID id = INVALID_ID;
    ID property_beg_id = INVALID_ID;
    int32 property_count = 0;

    ID cmd_id = INVALID_ID;
    CString name;
    CString args;

    Array<int32> residues {};
};

struct StatisticsContext {
    DynamicArray<String> string_buffer {};

    DynamicArray<Property> properties{};
    DynamicArray<PropertyData> property_data{};
    DynamicArray<Group> groups{};

    DynamicArray<PropertyCommand> property_commands;
    DynamicArray<GroupCommand> group_commands;
};

static StatisticsContext ctx;

static bool compute_atomic_distance(void* data, const Array<CString> args, const MoleculeDynamic* dynamic, int res_idx) {
	if (args != 2) return false;

    auto res = dynamic->molecule->residues[res_idx];
	auto int_a = to_int32(args[0]);
	auto int_b = to_int32(args[1]);

	if (!int_a.success || !int_b.success) return false;
	int32 atom_a = res.beg_atom_idx + int_a;
	int32 atom_b = res.beg_atom_idx + int_b;

    int32 count = dynamic->trajectory->num_frames;
    float* f_data = (float*)data;
    for (int32 i = 0; i < count; i++) {
        vec3 pos_a = dynamic->trajectory->frame_buffer[i].atom_positions[atom_a];
        vec3 pos_b = dynamic->trajectory->frame_buffer[i].atom_positions[atom_b];
        f_data[i] = math::distance(pos_a, pos_b);
    }

    return true;
}

static bool compute_atomic_angle(void* data, const Array<CString> args, const MoleculeDynamic* dynamic, int res_idx) {
	if (args != 3) return false;

	auto res = dynamic->molecule->residues[res_idx];
	auto int_a = to_int32(args[0]);
	auto int_b = to_int32(args[1]);
	auto int_c = to_int32(args[2]);

	if (!int_a.success || !int_b.success || !int_c.success) return false;
	int32 atom_a = res.beg_atom_idx + int_a;
	int32 atom_b = res.beg_atom_idx + int_b;
	int32 atom_c = res.beg_atom_idx + int_c;

	int32 count = dynamic->trajectory->num_frames;
	float* f_data = (float*)data;
	for (int32 i = 0; i < count; i++) {
		vec3 pos_a = dynamic->trajectory->frame_buffer[i].atom_positions[atom_a];
		vec3 pos_b = dynamic->trajectory->frame_buffer[i].atom_positions[atom_b];
		vec3 pos_c = dynamic->trajectory->frame_buffer[i].atom_positions[atom_c];

		f_data[i] = math::angle(pos_a - pos_b, pos_c - pos_b);
	}

	return true;
}

static bool compute_atomic_dihedral(void* data, const Array<CString> args, const MoleculeDynamic* dynamic, int res_idx) {
	if (args != 4) return false;

	auto res = dynamic->molecule->residues[res_idx];
	auto int_a = to_int32(args[0]);
	auto int_b = to_int32(args[1]);
	auto int_c = to_int32(args[2]);
	auto int_d = to_int32(args[2]);

	if (!int_a.success || !int_b.success || !int_c.success, !int_d.success) return false;
	int32 atom_a = res.beg_atom_idx + int_a;
	int32 atom_b = res.beg_atom_idx + int_b;
	int32 atom_c = res.beg_atom_idx + int_c;
	int32 atom_d = res.beg_atom_idx + int_d;

	int32 count = dynamic->trajectory->num_frames;
	float* f_data = (float*)data;
	for (int32 i = 0; i < count; i++) {
		vec3 pos_a = dynamic->trajectory->frame_buffer[i].atom_positions[atom_a];
		vec3 pos_b = dynamic->trajectory->frame_buffer[i].atom_positions[atom_b];
		vec3 pos_c = dynamic->trajectory->frame_buffer[i].atom_positions[atom_c];
		vec3 pos_d = dynamic->trajectory->frame_buffer[i].atom_positions[atom_d];

		f_data[i] = dihedral_angle(pos_a, pos_b, pos_c, pos_d);
	}

	return true;
}

static bool match_by_resid(const Array<CString> args, const MoleculeStructure* mol, int32 res_idx) {
    const auto& res = mol->residues[res_idx];
    for (const auto& arg : args) {
        if (compare(res.name, arg)) return true;
    }
    return false;
}

void initialize() {
    ctx.property_commands.push_back({ COMPUTE_ID("dist"), PropertyType::FLOAT32, compute_atomic_distance});
	ctx.property_commands.push_back({ COMPUTE_ID("angle"), PropertyType::FLOAT32, compute_atomic_angle });
	ctx.property_commands.push_back({ COMPUTE_ID("dihedral"), PropertyType::FLOAT32, compute_atomic_dihedral });

    ctx.group_commands.push_back({ COMPUTE_ID("resid"), match_by_resid});
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

bool compute_stats(MoleculeDynamic* dynamic) {
    ASSERT(dynamic);
    ASSERT(dynamic->molecule);
    ASSERT(dynamic->trajectory);

    for (auto& group : ctx.groups) {
        if (!group.residues) {
            // MATCH RESIDUE INDICES
            GroupCommand* group_cmd = find_id(ctx.group_commands, group.cmd_id);
            DynamicArray<CString> args = ctokenize(group.args);
            DynamicArray<int32> residue_indices;
            for (int32 res_idx = 0; res_idx < dynamic->molecule->residues.count; res_idx++) {
                if (group_cmd->func(args, dynamic->molecule, res_idx)) {
                    residue_indices.push_back(res_idx);
                }
            }
            group.residues.data = (int32*)MALLOC(residue_indices.count * sizeof(int32));
            group.residues.count = residue_indices.count;
            memcpy(group.residues.data, residue_indices.data, residue_indices.count * sizeof(int32));
        }
    }

    int32 count = dynamic->trajectory->num_frames;
    for (auto& prop : ctx.properties) {
        if (prop.data_beg_id == INVALID_ID) {
            // NEED TO COMPUTE PROPERTY DATA FOR RESIDUES
            PropertyCommand* prop_cmd = find_id(ctx.property_commands, prop.cmd_id);
            DynamicArray<CString> args = ctokenize(prop.args);
            Group* group = find_id(ctx.groups, prop.group_id);
            auto byte_size = count * get_stride(prop.type);

            // DATA AVERAGE
            StringBuffer<64> prop_data_name;
            snprintf(prop_data_name.beg(), 64, "%s.%s.avg", group->name.beg(), prop.name.beg());

            PropertyData prop_avg_data;
            prop_avg_data.id = COMPUTE_ID(prop_data_name.operator CString());
            prop_avg_data.group_id = group->id;
            prop_avg_data.property_id = prop.id;
            prop_avg_data.residue_idx = -1;
            prop_avg_data.data = MALLOC(byte_size);
            prop_avg_data.count = count;
            memset(prop_avg_data.data, 0, byte_size);

            prop.data_avg_id = prop_avg_data.id;
			ctx.property_data.push_back(prop_avg_data);

            // DATA
            for (int32 i = 0; i < (int32)group->residues.count; i++) {
                int32 res_idx = group->residues[i];
                void* data = MALLOC(byte_size);
                prop_cmd->func(data, args, dynamic, res_idx);

                for (int32 j = 0; j < count; j++) {
                    switch(prop.type) {
                        case PropertyType::FLOAT32:
                        ((float*)prop_avg_data.data)[j] += ((float*)data)[j] / (float)group->residues.count;
                    }
                }

                StringBuffer<64> prop_data_name;
                snprintf(prop_data_name.beg(), 64, "%s.%s.%i", group->name.beg(), prop.name.beg(), i);

                PropertyData prop_data;
                prop_data.id = COMPUTE_ID(prop_data_name.operator CString());
                prop_data.group_id = group->id;
                prop_data.property_id = prop.id;
                prop_data.residue_idx = res_idx;
                prop_data.data = data;
                prop_data.count = count;

                ctx.property_data.push_back(prop_data);
                if (i == 0) {
                    prop.data_beg_id = prop_data.id;
                    prop.data_count = 1;
                } else {
                    prop.data_count++;
                }
            }
        }
    }

    return true;
}

void register_property_command(CString command, PropertyType type, PropertyComputeFunc func) {
    ID id = COMPUTE_ID(command);
    auto cmd = find_id(ctx.property_commands, id);
    if (cmd != nullptr) {
        printf("ERROR: PROPERTY COMMAND %s ALREADY REGISTERED!", command.beg());
        return;
    }

    ctx.property_commands.push_back({id, type, func});
}

void register_group_command(CString command, ResidueMatchFunc func) {
    ID id = COMPUTE_ID(command);
    auto cmd = find_id(ctx.group_commands, id);
    if (cmd != nullptr) {
        printf("ERROR: GROUP COMMAND %s ALREADY REGISTERED!", command.beg());
        return;
    }

    ctx.group_commands.push_back({id, func});
}

ID create_group(CString name, CString cmd, CString args) {
	ID grp_id = COMPUTE_ID(name);
	Group* grp = find_id(ctx.groups, grp_id);
	if (grp != nullptr) {
		StringBuffer<32> buf = name;
		printf("ERROR: GROUP '%s' ALREADY REGISTERED!", buf.beg());
		return INVALID_ID;
	}

	ID grp_cmd_id = COMPUTE_ID(cmd);
	GroupCommand* grp_cmd = find_id(ctx.group_commands, grp_cmd_id);
	if (grp_cmd == nullptr) {
		StringBuffer<32> buf = cmd;
		printf("ERROR: UNIDENTIFIED GROUP COMMAND '%s'!", buf.beg());
		return INVALID_ID;
	}

	Group group;
	group.id = grp_id;
	group.property_beg_id = INVALID_ID;
	group.property_count = 0;
	group.name = alloc_string(name);
    group.args = alloc_string(args);
	group.cmd_id = grp_cmd_id;
    group.residues = {};
    //group.residue_indices = nullptr;
    //group.residue_count = 0;

	ctx.groups.push_back(group);

	return group.id;
}

void remove_group(ID group_id) {
	Group* group = find_id(ctx.groups, group_id);
	if (!group) {
		printf("ERROR: COULD NOT FIND GROUP!\n");
		return;
	}

    for (PropertyData* pd = ctx.property_data.beg(); pd != ctx.property_data.end(); pd++) {
        if (pd->group_id == group_id) ctx.property_data.remove(pd);
    }

	for (Property* p = ctx.properties.beg(); p != ctx.properties.end(); p++) {
		if (p->group_id == group_id) ctx.properties.remove(p);
	}

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
    return ctx.groups.count;
}

ID get_property(ID group_id, CString name) {
    for (const auto &p : ctx.properties) {
        if (p.group_id == group_id && compare(p.name, name)) return p.id;
    }
    return INVALID_ID;
}

ID get_property(ID group_id, int32 idx) {
    // @TODO: Perhaps do this in another way...
    int32 counter = 0;
    for (const auto& p : ctx.properties) {
        if (p.group_id == group_id && counter++ == idx) return p.id;
    }
    return INVALID_ID;
}

int32 get_property_count(ID group_id) {
    auto group = find_id(ctx.groups, group_id);
    if (group) {
        return group->property_count;
    }
    return 0;
}

ID create_property(ID group_id, CString name, CString cmd, CString args) {
    Group* group = find_id(ctx.groups, group_id);
    if (group == nullptr) {
        StringBuffer<32> buf = name;
        printf("ERROR: GROUP '%s' NOT FOUND!", buf.beg());
        return INVALID_ID;
    }

    StringBuffer<32> name_buf = name;
    StringBuffer<32> buf;
    snprintf(buf.beg(), 32, "%s.%s", group->name.beg(), name_buf.beg());
    ID prop_id = COMPUTE_ID(buf.operator CString());
    Property* old_prop = find_id(ctx.properties, prop_id);
    if (old_prop != nullptr) {
        printf("ERROR: PROPERTY '%s' ALREADY EXISTS!", buf.beg());
        return INVALID_ID;
    }

    ID prop_cmd_id = COMPUTE_ID(cmd);
    PropertyCommand* prop_cmd = find_id(ctx.property_commands, prop_cmd_id);
    if (prop_cmd == nullptr) {
        StringBuffer<32> cmd_buf = cmd;
        printf("ERROR: UNIDENTIFIED PROPERTY COMMAND '%s'!", cmd_buf.beg());
        return INVALID_ID;
    }

    Property prop;
    prop.id = prop_id;
    prop.group_id = group_id;
    prop.data_beg_id = INVALID_ID;
    prop.data_count = 0;

    prop.cmd_id = prop_cmd_id;
    prop.type = prop_cmd->type;
    prop.name = alloc_string(name);
    prop.args = alloc_string(args);

    Property* group_beg_prop = find_id(ctx.properties, group->property_beg_id);
    if (group_beg_prop != nullptr) {
        ctx.properties.insert(group_beg_prop + group->property_count, prop);
    }
    else {
        ctx.properties.push_back(prop);
    }
    group->property_count++;

	return prop.id;
}

void remove_property(ID prop_id) {
    // @TODO remove from instances
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

void* get_property_data(ID prop_id, int32 residue_idx) {
	if (prop_id != INVALID_ID) {
		for (const auto& prop_data : ctx.property_data) {
			if (prop_data.property_id == prop_id && prop_data.residue_idx == residue_idx)
				return prop_data.data;
		}
	}
	return nullptr;
}

void* get_property_avg_data(ID prop_id) {
	if (prop_id != INVALID_ID) {
		for (const auto& prop_data : ctx.property_data) {
			if (prop_data.property_id == prop_id && prop_data.residue_idx == -1)
				return prop_data.data;
		}
	}
	return nullptr;
}

int32 get_property_data_count(ID prop_id) {
	if (prop_id != INVALID_ID) {
		for (const auto& prop_data : ctx.property_data) {
			if (prop_data.property_id == prop_id)
				return prop_data.count;
		}
	}
	return 0;
}

PropertyType get_property_type(ID prop_id) {
	auto prop = find_id(ctx.properties, prop_id);
	if (prop) {
		return prop->type;
	}
	return PropertyType::UNKNOWN;
}

CString	get_property_name(ID prop_id) {
	auto prop = find_id(ctx.properties, prop_id);
	if (prop) {
		return prop->name;
	}
	return {};
}

}  // namespace stats