#include "stats.h"
#include <core/math_utils.h>
#include <core/hash.h>
#include <mol/molecule_dynamic.h>

#define COMPUTE_ID(x) (hash::crc64(x))

namespace stats {

struct PropertyCommand {
    ID id;
    PropertyType data_type;
    PropertyComputeFunc func;
};

struct GroupCommand {
    ID id;
    ResidueMatchFunc func;
};

struct PropertyRecipe {
	StringBuffer<32> name {};
	StringBuffer<32> cmd {};
	StringBuffer<64> args {};
};

struct Recipe {
	StringBuffer<32> name{};
	StringBuffer<64> args{};
	ID cmd_id;
};

struct Property {
    ID id = INVALID_ID;
    ID instance_id = INVALID_ID;
    ID group_id = INVALID_ID;

    void* data = nullptr;
    int32 count = 0;
    PropertyType type = PropertyType::FLOAT32;

	// Recipe for computing property values
	Recipe recipe;
};

struct Instance {
    ID id = INVALID_ID;
    ID group_id = INVALID_ID;
    ID property_beg_id = INVALID_ID;
    int32 property_count = 0;

    int32 mol_res_idx = -1;
};

struct Group {
    ID id = INVALID_ID;
    ID instance_avg_id = INVALID_ID;
    ID instance_beg_id = INVALID_ID;
    int32 instance_count = 0;

	// Recipe for matching residues
	Recipe recipe;
};

struct StatisticsContext {
	//DynamicArray<GroupRecipe> group_recipes;
	//DynamicArray<PropertyRecipe> property_recipes;

    DynamicArray<Property> properties{};
    DynamicArray<Instance> instances{};
    DynamicArray<Group> groups{};
    DynamicArray<PropertyCommand> property_commands;
    DynamicArray<GroupCommand> group_commands;
};

static StatisticsContext ctx;

static bool compute_atomic_distance(void* data, const Array<CString> args, MoleculeDynamic* dynamic, int res_idx) {
    ASSERT(args.count == 2);

    auto res = dynamic->molecule->residues[res_idx];
    int32 atom_a = res.beg_atom_idx + to_int32(args[0]);
    int32 atom_b = res.beg_atom_idx + to_int32(args[1]);

    int32 count = dynamic->trajectory->num_frames;
    float* f_data = (float*)data;
    for (int32 i = 0; i < count; i++) {
        vec3 pos_a = dynamic->trajectory->frame_buffer[i].atom_positions[atom_a];
        vec3 pos_b = dynamic->trajectory->frame_buffer[i].atom_positions[atom_b];
        f_data[i] = math::distance(pos_a, pos_b);
    }

    return true;
}

static bool match_by_resid(const Array<CString> args, MoleculeStructure* mol, Residue* res) {
    for (const auto& arg : args) {
        if (compare(res->id, arg)) return true;
    }
    return false;
}

void initialize() {
    ctx.property_commands.push_back({ COMPUTE_ID("dist"), PropertyType::FLOAT32, compute_atomic_distance});

    ctx.group_commands.push_back({ COMPUTE_ID("resid"), match_by_resid});
}

template <typename T>
static T* find_id(Array<T> data, ID id) {
    for (auto& item : data) {
        if (item.id == id) return &item;
    }
    return nullptr;
}

/*
bool compute_stats(StatisticsContext* dst, MoleculeDynamic* dynamic, Array<GroupRecipe> group_recipes) {
    ASSERT(dst);
    ASSERT(dynamic);
    ASSERT(dynamic->molecule);
    ASSERT(dynamic->trajectory);

    dst->properties.clear();
    dst->instances.clear();
    dst->groups.clear();

    for (const auto& group_recipe : group_recipes) {
        ID group_id = COMPUTE_ID(group_recipe.name);
        ID group_cmd_id = COMPUTE_ID(group_recipe.cmd);

        GroupCommand* group_cmd = find_id(group_commands, group_cmd_id);
        if (!group_cmd) {
            printf("Error! could not find group command '%s'\n", group_recipe.cmd.beg());
            return false;
        }

        // Assert that all property commands exist
        for (const auto& prop_recipe : group_recipe.properties) {
            ID prop_cmd_id = COMPUTE_ID(prop_recipe.cmd);
            PropertyCommand* prop_cmd = find_id(property_commands, prop_cmd_id);
            if (!prop_cmd) {
                printf("Error! could not find property command '%s'\n", prop_recipe.cmd.beg());
                return false;
            }
        }

        int32 group_idx = dst->groups.size();
        Group group;
        group.id = group_id;
        group.instance_beg_idx = dst->instances.size();
        group.instance_end_idx = group.instance_beg_idx;
        copy(group.name, group_recipe.name);
        auto res_indices = group_cmd->func(group_recipe.args, dynamic->molecule);

        for (int32 i = 0; i < res_indices.size(); i++) {
            int32 instance_idx = group.instance_end_idx++;
            int32 res_idx = res_indices[i];

            // Create instance
            Instance instance;
            instance.group_idx = group_idx;
            instance.property_beg_idx = dst->properties.size();
            instance.property_end_idx = instance.property_beg_idx;
            instance.mol_res_idx = res_indices[i];
            StringBuffer<32> buf;
            int buf_len = snprintf(buf.beg(), 32, "%s.%i", group.name.beg(), i);
            instance.id = COMPUTE_ID(buf.operator CString());

            // Create and compute properties for instance
            for (const auto& prop_recipe : group_recipe.properties) {
                int32 prop_idx = instance.property_end_idx++;
                ID prop_cmd_id = COMPUTE_ID(prop_recipe.cmd);
                PropertyCommand* prop_cmd = find_id(property_commands, prop_cmd_id);
                auto prop_result = prop_cmd->func(prop_recipe.args, dynamic, res_idx);
                Property prop;
                prop.data = prop_result.data;
                prop.count = prop_result.count;
                StringBuffer<32> prop_name = prop_recipe.name;
                snprintf(buf.beg() + buf_len, 32 - buf_len, ".%s", prop_name);
                prop.id = COMPUTE_ID(buf.operator CString());
                prop.instance_idx = instance_idx;
                prop.group_idx = group_idx;
                copy(prop.name, prop_recipe.name);
                prop.type = prop_cmd->data_type;

                dst->properties.push_back(prop);
            }
            dst->instances.push_back(instance);
        }

        StringBuffer<32> buf;
        snprintf(buf.beg(), 32, "%s.avg", group.name.beg());
        Instance avg_instance;
        avg_instance.group_idx = group_idx;
        avg_instance.id = COMPUTE_ID(buf.operator CString());
        avg_instance.mol_res_idx = -1;
        avg_instance.property_beg_idx = dst->properties.size();
        avg_instance.property_end_idx = avg_instance.property_beg_idx + 1;

        group.instance_avg_idx = dst->instances.size();

        // Create and Compute average properties for instance
        int property_count = group_recipe.properties.count;
        int instance_count = group.instance_end_idx - group.instance_beg_idx;

        for (int p_idx = 0; p_idx < property_count; p_idx++) {
            Property avg_prop;
            avg_prop.group_idx = group_idx;
            avg_prop.instance_idx = group.instance_avg_idx;
            StringBuffer<32> prop_name = group_recipe.properties[p_idx].name;
            snprintf(avg_prop.name.beg(), 32, "%s.%s", buf, prop_name);
            avg_prop.id = COMPUTE_ID(avg_prop.name.operator CString());

            for (int i_idx = group.instance_beg_idx; i_idx < group.instance_end_idx; i_idx++) {
                const Instance& i = dst->instances[i_idx];
                const Property& p = dst->properties[i.property_beg_idx + p_idx];
                if (p_idx == 0) {
                    int32 byte_size = p.count * get_stride(p.type);
                    avg_prop.data = MALLOC(byte_size);
                    avg_prop.count = p.count;
                    avg_prop.type = p.type;
                } else {
                    switch (avg_prop.type) {
                        case (PropertyType::FLOAT32):
                            float* dst = (float*)avg_prop.data;
                            float* src = (float*)p.data;
                            for (int i = 0; i < avg_prop.count; i++) {
                                dst[i] += src[i] / (float)instance_count;
                            }
                            break;
                    }
                }
            }
            dst->properties.push_back(avg_prop);
        }
        dst->instances.push_back(avg_instance);

        dst->groups.push_back(group);
    }

    return true;
}
*/

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
	group.instance_avg_id = INVALID_ID;
	group.instance_beg_id = INVALID_ID;
	group.instance_count = 0;
	copy(group.recipe.name, name);
	copy(group.recipe.args, args);
	group.recipe.cmd_id = grp_cmd_id;

	ctx.groups.push_back(group);

	return group.id;
}

void remove_group(ID group_id) {
	Group* group = find_id(ctx.groups, group_id);
	if (!group) {
		printf("ERROR: COULD NOT FIND GROUP!\n");
		return;
	}
	
	/*
	if (group->instance_avg_id != INVALID_ID) {
		Instance* inst = find_id(ctx.instances, group->instance_avg_id);
		if (inst) {
			Property* prop = find_id(ctx.properties, inst->property_beg_id);
			if (prop) {
				ctx.properties.remove(prop, inst->property_count);
			}
			ctx.instances.remove(inst);
		}
	}

	if (group->instance_beg_id != INVALID_ID) {
		Instance* inst = find_id(ctx.instances, group->instance_beg_id);
		if (inst) {
			ctx.instances.remove(inst, group->instance_count);
		}
	}
	*/

	for (Instance* i = ctx.instances.beg(); i != ctx.instances.end(); i++) {
		if (i->group_id == group_id) ctx.instances.remove(i);
	}

	for (Property* p = ctx.properties.beg(); p != ctx.properties.end(); p++) {
		if (p->group_id == group_id) ctx.properties.remove(p);
	}
}

ID create_property(ID group_id, CString name, CString args) {
	return INVALID_ID;
}

void remove_property(ID prop_id) {

}

}  // namespace stats