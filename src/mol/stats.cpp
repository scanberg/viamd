#include "stats.h"
#include <core/math_utils.h>
#include <core/hash.h>
#include <core/log.h>
#include <mol/molecule.h>
#include <mol/molecule_dynamic.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <ctype.h>
#include <stdarg.h>
#include <new>
#include <tinyexpr.h>

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

    Property* current_property = nullptr;
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

void free_structure_data(Array<StructureData>* structure_data) {
    ASSERT(structure_data);
    if (structure_data) {
        for (auto& s : *structure_data) {
            s.StructureData::~StructureData();
        }
        free_array(structure_data);
    }
}

void init_instance_data(Array<InstanceData>* instance_data, int32 num_instances, int32 num_frames) {
    ASSERT(instance_data);
    free_instance_data(instance_data);

    *instance_data = allocate_array<InstanceData>(num_instances);
    for (auto& inst : *instance_data) {
        new (&inst) InstanceData(num_frames);
    }
}

void free_instance_data(Array<InstanceData>* instance_data) {
    ASSERT(instance_data);
    if (instance_data->data) {
        for (auto inst : *instance_data) {
            inst.~InstanceData();
        }
    }
    free_array(instance_data);
}

void init_structure_data(Array<StructureData>* structure_data, int32 count) {
    ASSERT(count > 0);
    ASSERT(structure_data);
    free_structure_data(structure_data);
    *structure_data = allocate_array<StructureData>(count);
    for (auto& s : *structure_data) {
        new (&s) StructureData();
    }
}

// HISTOGRAMS
void init_histogram(Histogram* hist, int32 num_bins, Range value_range) {
    ASSERT(hist);
    if (hist->bins.data) {
        FREE(hist->bins.data);
    }
    hist->bins.data = (float*)MALLOC(num_bins * sizeof(float));
    hist->bins.count = num_bins;
    hist->value_range = value_range;
    hist->bin_range = {0, 0};

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

void compute_histogram(Histogram* hist, Array<const float> data) {
    ASSERT(hist);
    const int32 num_bins = (int32)hist->bins.count;
    const float scl = num_bins / (hist->value_range.y - hist->value_range.x);
    hist->bin_range = {0, 0};
    for (auto v : data) {
        int32 bin_idx = math::clamp((int32)((v - hist->value_range.x) * scl), 0, num_bins - 1);
        hist->bins[bin_idx]++;
        hist->bin_range.y = math::max(hist->bin_range.y, hist->bins[bin_idx]);
    }
}

void clear_histogram(Histogram* hist) {
    ASSERT(hist);
    if (hist->bins.data) {
        memset(hist->bins.data, 0, hist->bins.count * sizeof(float));
    }
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

static Range compute_range(const Property& prop) {
    Range range{0, 0};
    if (prop.instance_data) {
        for (int32 i = 0; i < prop.instance_data.count; i++) {
            Range r = compute_range(prop.instance_data[i].data);
            if (i == 0)
                range = r;
            else {
                range.x = math::min(range.x, r.x);
                range.y = math::max(range.y, r.y);
            }
        }
    } else {
        range = compute_range(prop.data);
    }
    return range;
}

void set_error_message(const char* fmt, ...) {
    if (ctx.current_property) {
        char* buf = ctx.current_property->error_msg.buffer;
        int32 len = ctx.current_property->error_msg.MAX_LENGTH;
        va_list ap;
        va_start(ap, fmt);
        vsnprintf(buf, len, fmt, ap);
        va_end(ap);
        LOG_ERROR("Error when evaluating property '%s': %s", ctx.current_property->name.cstr(), buf);
    }
}

// Parse arithmetic expression

/*
// Operations
namespace op {
        void negate(Array<float> data) {
                for (auto& v : data) v = -v;
        }

        void add(Array<float> dst, Array<const float> val) {
                ASSERT(dst.count == val.count);
                for (int64 i = 0; i < dst.count; i++) {
                        dst[i] += val[i];
                }
        }

        void sub(Array<float> dst, Array<const float> val) {
                ASSERT(dst.count == val.count);
                for (int64 i = 0; i < dst.count; i++) {
                        dst[i] -= val[i];
                }
        }

        void mul(Array<float> dst, Array<const float> val) {
                ASSERT(dst.count == val.count);
                for (int64 i = 0; i < dst.count; i++) {
                        dst[i] *= val[i];
                }
        }

        void div(Array<float> dst, Array<const float> val) {
                ASSERT(dst.count == val.count);
                for (int64 i = 0; i < dst.count; i++) {
                        dst[i] /= val[i];
                }
        }

        typedef void(*ExprFunc)(Array<float> dst, Array<const float> val);

        void abs(Array<float> dst, Array<const float>) {
                for (auto& v : dst) {
                        v = math::abs(v);
                }
        }

        void pi(Array<float> dst, Array<const float>) {

        }
}


// parsing
char peek(CString expr) { return expr[0]; }

char get(CString* expr) { return *(expr->data)++; }

DynamicArray<float> number(CString* expr, int32 count) {
    const char* beg = expr->beg();
    const char* end = beg + 1;
    while (end != expr->end() && (isdigit(*end) || *end == '.')) end++;
    CString float_expr(beg, end);
        *expr = { end, expr->end() };

    auto res = to_float(float_expr);
    if (!res) {
        set_error_message("Could not parse number given in expression: '%s'", make_tmp_str(float_expr).cstr());
    }
    return DynamicArray<float>(count, res.value);
}

DynamicArray<float> variable(CString* expr, Array<Property*> properties, int32 count) {
    const char* beg = expr->beg();
    const char* end = beg + 1;
    while (end != expr->end() && (isalpha(*end) || isdigit(*end))) end++;
    CString prop_name(beg, end);
        *expr = { end, expr->end() };

    for (auto p : properties) {
                if (compare(prop_name, p->name)) {
                        return p->data;
                }
    }
        set_error_message("Could not find property: '%s'", make_tmp_str(prop_name).cstr());
    return DynamicArray<float>(count, 0);
}

DynamicArray<float> expression(CString* expr, Array<Property*> properties, int32 count);



DynamicArray<float> factor(CString* expr, Array<Property*> properties, int32 count) {
        struct Command {
                const char* cmd = "";
                ExprFunc func = nullptr;
        };

        const Command cmd_list[] = { {"abs", abs}, {"PI", pi} };

    if (isdigit(peek(*expr)))
        return number(expr, count);
    else if (isalpha(peek(*expr)))
        return variable(expr, properties, count);
    else if (peek(*expr) == '(') {
        get(expr);  // '('
        auto result = expression(expr, properties, count);
        get(expr);  // ')'
        return result;
    } else if (peek(*expr) == '-') {
        get(expr);
        auto result = factor(expr, properties, count);
        negate(result);
        return result;
    }
        set_error_message("Some thing went wrong when parsing the expression");
        return DynamicArray<float>(count, 0.f);
}

DynamicArray<float> term(CString* expr, Array<Property*> properties, int32 count) {
    auto result = factor(expr, properties, count);
    while (peek(*expr) == '*' || peek(*expr) == '/')
        if (get(expr) == '*')
            mul(result, factor(expr, properties, count));
        else
            div(result, factor(expr, properties, count));
    return result;
}

DynamicArray<float> expression(CString* expr, Array<Property*> properties, int32 count) {
    auto result = term(expr, properties, count);
    while (peek(*expr) == '+' || peek(*expr) == '-')
        if (get(expr) == '+')
            add(result, term(expr, properties, count));
        else
            sub(result, term(expr, properties, count));
    return result;
}
*/

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
        int32 first = 0;
        int32 last = 0;

        if (is_range(arg)) {
            if (!extract_range(&first, &last, arg)) {
                set_error_message("Failed to parse range in argument for residue");
                return false;
            }
            if (first == -1) first = 1;
            if (last == -1) last = molecule.residues.count;
        } else {
            auto id = to_int(arg);
            if (!id.success) {
                set_error_message("Failed to parse argument for residue");
                return false;
            }
            first = last = id;
        }

        if (first < 1 || molecule.residues.count < last) {
            set_error_message("Index for residue is out of bounds");
            return false;
        }
        for (int32 i = first - 1; i < last; i++) {
            const auto& res = molecule.residues[i];
            data->structures.push_back({res.beg_atom_idx, res.end_atom_idx});
        }
    }

    return true;
}

bool structure_match_chainid(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule) {
    ASSERT(data);

    // Expect args.count to be > 0
    if (args.count == 0) {
        set_error_message("Expects one or more arguments for chainid");
        return false;
    }

    for (const auto& chain : molecule.chains) {
        for (const auto& arg : args) {
            if (compare(chain.id, arg)) {
                data->structures.push_back({get_atom_beg_idx(molecule, chain), get_atom_end_idx(molecule, chain)});
                break;
            }
        }
    }

    return true;
}

bool structure_match_chain(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule) {
    ASSERT(data);

    // Expect args to be  > 0
    if (args.count == 0) {
        set_error_message("Expects one or more arguments for chain");
        return false;
    }

    for (const auto& arg : args) {
        int32 first = 0;
        int32 last = 0;

        if (is_range(arg)) {
            if (!extract_range(&first, &last, arg)) {
                set_error_message("Failed to parse range in argument for chain");
                return false;
            }
            if (first == -1) first = 1;
            if (last == -1) last = molecule.chains.count;
        } else {
            auto id = to_int(arg);
            if (!id.success) {
                set_error_message("Failed to parse argument for chain");
                return false;
            }
            first = last = id;
        }

        if (first < 1 || molecule.chains.count < last) {
            set_error_message("Index for chain is out of bounds");
            return false;
        }
        for (int32 i = first - 1; i <= last - 1; i++) {
            const auto& chain = molecule.chains[i];
            data->structures.push_back({get_atom_beg_idx(molecule, chain), get_atom_end_idx(molecule, chain)});
        }
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
        int32 first = 0;
        int32 last = 0;

        if (is_range(arg)) {
            if (!extract_range(&first, &last, arg)) {
                set_error_message("Failed to parse range in argument for atom");
                return false;
            }
            if (first == -1) first = 1;
            if (last == -1) last = molecule.chains.count;
        } else {
            auto id = to_int(arg);
            if (!id.success) {
                set_error_message("Failed to parse argument for atom");
                return false;
            }
            first = last = id;
        }

        if (first < 1 || molecule.atom_positions.count < last) {
            set_error_message("Index for atom is out of bounds");
            return false;
        }

        data->structures.push_back({first - 1, last});
    }

    return true;
}

bool structure_extract_resatom(StructureData* data, const Array<CString> args, const MoleculeStructure&) {
    ASSERT(data);
    if (args.count != 1) {
        set_error_message("resatom requires exactly 1 argument");
        return false;
    }

    int32 first = 0;
    int32 last = 0;

    if (is_range(args[0])) {
        if (!extract_range(&first, &last, args[0])) {
            set_error_message("Failed to parse range in argument for resatom");
            return false;
        }
    } else {
        auto id = to_int(args[0]);
        if (!id.success) {
            set_error_message("Failed to parse argument for resatom");
            return false;
        }
        first = last = id;
    }

    for (auto& s : data->structures) {
        int32 count = s.end_idx - s.beg_idx;
        if (first == -1) first = 1;
        if (last == -1) last = count;

        if (count < 0 || first < 1 || count < last) {
            set_error_message("restom: Index is out of range for structure");
            return false;
        }

        int new_beg = s.beg_idx + first - 1;
        int new_end = s.beg_idx + last;
        s = {new_beg, new_end};
    }

    return true;
}

bool structure_apply_aggregation_strategy_com(StructureData* data, const Array<CString>, const MoleculeStructure&) {
    ASSERT(data);
    data->strategy = AggregationStrategy::COM;
    return true;
}

// Helper funcs
vec3 compute_com(Array<const vec3> positions) {
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

Array<const vec3> extract_positions(Structure structure, Array<const vec3> atom_positions) {
    return atom_positions.sub_array(structure.beg_idx, structure.end_idx - structure.beg_idx);
}

static inline float multi_distance(Array<const vec3> pos_a, Array<const vec3> pos_b) {
    if (pos_a.count == 0 || pos_b.count == 0)
        return 0.f;
    else if (pos_a.count == 1 && pos_b.count == 1)
        return math::distance(pos_a[0], pos_b[0]);
    else {
        float dist = 0.f;
        float count = 0.f;

        if (pos_a.count > 1) {
            for (const auto& p : pos_a) dist += math::distance(p, pos_b[0]);
            count = (float)pos_a.count;
        } else if (pos_b.count > 1) {
            for (const auto& p : pos_b) dist += math::distance(pos_a[0], p);
            count = (float)pos_b.count;
        } else {
            // ERROR
            return 0.f;
        }
        return dist / count;
    }
}

static inline float multi_angle(Array<const vec3> pos_a, Array<const vec3> pos_b, Array<const vec3> pos_c) {
    if (pos_a.count == 0 || pos_b.count == 0 || pos_c.count == 0)
        return 0.f;
    else if (pos_a.count == 1 && pos_b.count == 1 && pos_c.count == 1)
        return math::angle(pos_a[0], pos_b[0], pos_c[0]);
    else {
        float angle = 0.f;
        float count;
        if (pos_a.count > 1.f) {
            for (const auto& p : pos_a) angle += math::angle(p, pos_b[0], pos_c[0]);
            count = (float)pos_a.count;
        } else if (pos_b.count > 1.f) {
            for (const auto& p : pos_b) angle += math::angle(pos_a[0], p, pos_c[0]);
            count = (float)pos_b.count;
        } else if (pos_c.count > 1.f) {
            for (const auto& p : pos_c) angle += math::angle(pos_a[0], pos_b[0], p);
            count = (float)pos_c.count;
        } else {
            // ERROR
            return 0.f;
        }
        return angle / count;
    }
}

static inline float multi_dihedral(Array<const vec3> pos_a, Array<const vec3> pos_b, Array<const vec3> pos_c, Array<const vec3> pos_d) {
    if (pos_a.count == 0 || pos_b.count == 0 || pos_c.count == 0 || pos_d.count == 0)
        return 0.f;
    else if (pos_a.count == 1 && pos_b.count == 1 && pos_c.count == 1 && pos_d.count == 1)
        return math::dihedral_angle(pos_a[0], pos_b[0], pos_c[0], pos_d[0]);
    else {
        float angle = 0.f;
        float count;
        if (pos_a.count > 1.f) {
            for (const auto& p : pos_a) angle += math::dihedral_angle(p, pos_b[0], pos_c[0], pos_d[0]);
            count = (float)pos_a.count;
        } else if (pos_b.count > 1.f) {
            for (const auto& p : pos_b) angle += math::dihedral_angle(pos_a[0], p, pos_c[0], pos_d[0]);
            count = (float)pos_b.count;
        } else if (pos_c.count > 1.f) {
            for (const auto& p : pos_c) angle += math::dihedral_angle(pos_a[0], pos_b[0], p, pos_d[0]);
            count = (float)pos_c.count;
        } else if (pos_c.count > 1.f) {
            for (const auto& p : pos_d) angle += math::dihedral_angle(pos_a[0], pos_b[0], pos_c[0], p);
            count = (float)pos_d.count;
        } else {
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

    // Sync the number of structures between arguments
    if (!sync_structure_data_length(prop->structure_data)) {
        return false;
    }
    const int32 structure_count = (int32)prop->structure_data[0].structures.count;

    init_instance_data(&prop->instance_data, structure_count, dynamic.trajectory.num_frames);

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
                pos[0] = {&com[0], 1};
            }
            if (prop->structure_data[1].strategy == COM) {
                com[1] = compute_com(pos[1]);
                pos[1] = {&com[1], 1};
            }

            prop->instance_data[j].data[i] = multi_distance(pos[0], pos[1]);
            val += prop->instance_data[j].data[i];
        }

        prop->data[i] = val * scl;
    }

    prop->data_range = compute_range(*prop);
    prop->periodic = false;
    prop->unit = "?";

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

    // Sync the number of structures between arguments
    if (!sync_structure_data_length(prop->structure_data)) {
        return false;
    }
    const int32 structure_count = (int32)prop->structure_data[0].structures.count;

    init_instance_data(&prop->instance_data, structure_count, dynamic.trajectory.num_frames);

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
                pos[0] = {&com[0], 1};
            }
            if (prop->structure_data[1].strategy == COM) {
                com[1] = compute_com(pos[1]);
                pos[1] = {&com[1], 1};
            }
            if (prop->structure_data[2].strategy == COM) {
                com[2] = compute_com(pos[2]);
                pos[2] = {&com[2], 1};
            }

            prop->instance_data[j].data[i] = multi_angle(pos[0], pos[1], pos[2]);
            val += prop->instance_data[j].data[i];
        }

        prop->data[i] = val * scl;
    }

    prop->data_range = {0, math::PI};
    prop->periodic = true;
    prop->unit = "?";

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

    // Sync the number of structures between arguments
    if (!sync_structure_data_length(prop->structure_data)) {
        return false;
    }
    const int32 structure_count = (int32)prop->structure_data[0].structures.count;

    init_instance_data(&prop->instance_data, structure_count, dynamic.trajectory.num_frames);

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
                pos[0] = {&com[0], 1};
            }
            if (prop->structure_data[1].strategy == COM) {
                com[1] = compute_com(pos[1]);
                pos[1] = {&com[1], 1};
            }
            if (prop->structure_data[2].strategy == COM) {
                com[2] = compute_com(pos[2]);
                pos[2] = {&com[2], 1};
            }
            if (prop->structure_data[3].strategy == COM) {
                com[3] = compute_com(pos[3]);
                pos[3] = {&com[3], 1};
            }

            prop->instance_data[j].data[i] = multi_dihedral(pos[0], pos[1], pos[2], pos[3]);
            val += prop->instance_data[j].data[i];
        }

        prop->data[i] = val * scl;
    }

    prop->data_range = {-math::PI, math::PI};
    prop->periodic = true;
    prop->unit = "?";

    return true;
}

static bool compute_expression(Property* prop, const Array<CString> args, const MoleculeDynamic&) {
    ASSERT(prop);
    if (args.count == 0) {
        set_error_message("expression expects 1 or more arguments");
        return false;
    }

    // Concatenate all arguments
    auto expr_str = make_tmp_str(CString(args.front().beg(), args.back().end()));

    if (!balanced_parentheses(expr_str)) {
        set_error_message("Expression contains unbalanced parentheses!");
        return false;
    }

    // Extract which properties preceedes this property
    Array<Property*> properties = ctx.properties;
    properties.count = 0;
    for (int i = 0; i < ctx.properties.count; i++) {
        if (prop == ctx.properties[i]) {
            properties.count = i;
            break;
        }
    }

    DynamicArray<double> values(properties.count, 0);
    DynamicArray<te_variable> vars;
    for (int32 i = 0; i < properties.count; i++) {
        vars.push_back({properties[i]->name.cstr(), &values[i]});
    }

    int err;
    te_expr* expr = te_compile(expr_str.cstr(), vars.data, vars.count, &err);

    if (expr) {
        int32 max_instance_count = 0;
        for (auto* p : properties) {
            max_instance_count = math::max(max_instance_count, (int32)p->instance_data.count);
        }

        init_instance_data(&prop->instance_data, max_instance_count, (int32)prop->data.count);

        float scl = 1.f / (float)max_instance_count;
        for (int32 frame = 0; frame < prop->data.count; frame++) {
            float val = 0.f;
            for (int32 i = 0; i < max_instance_count; i++) {
                for (int32 j = 0; j < values.count; j++) {
                    if (properties[j]->instance_data.count == max_instance_count) {
                        values[j] = properties[j]->instance_data[i].data[frame];
                    } else {
                        values[j] = properties[j]->instance_data[0].data[frame];
                    }
                }
                prop->instance_data[i].data[frame] = (float)te_eval(expr);
                val += prop->instance_data[i].data[frame];
            }
            prop->data[frame] = val * scl;
        }

        te_free(expr);
    } else {
        set_error_message("Malformed expression:\n%s\n%*s^\nError near here", expr_str.cstr(), err - 1, "");
        return false;
    }

    prop->data_range = compute_range(*prop);
    prop->periodic = false;
    prop->unit = "";

    return true;
}

static bool visualize_structures(const Property& prop, const MoleculeDynamic& dynamic) {
    if (prop.structure_data.count == 0) return false;

    if (prop.structure_data.count == 1) {
        for (const auto& s : prop.structure_data[0].structures) {
            Array<const vec3> pos = extract_positions(s, dynamic.molecule.atom_positions);
            if (prop.structure_data[0].strategy == COM) {
                immediate::draw_point(compute_com(pos));
            } else {
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

        const int32 NUM_COLORS = 4;
        const uint32 COLORS[NUM_COLORS]{0xffe3cea6, 0xffb4781f, 0xff8adfb2, 0xff2ca033};
        const uint32 LINE_COLOR = 0xffcccccc;

        for (int32 i = 0; i < count; i++) {
            pos_prev = extract_positions(prop.structure_data[0].structures[i], dynamic.molecule.atom_positions);
            if (prop.structure_data[0].strategy == COM) {
                com_prev = compute_com(pos_prev);
                pos_prev = {&com_prev, 1};
            }
            for (const auto& p : pos_prev) {
                immediate::draw_point(p, COLORS[0]);
            }
            for (int32 j = 1; j < prop.structure_data.count; j++) {
                pos_next = extract_positions(prop.structure_data[j].structures[i], dynamic.molecule.atom_positions);
                if (prop.structure_data[j].strategy == COM) {
                    com_next = compute_com(pos_next);
                    pos_next = {&com_next, 1};
                }
                for (const auto& p : pos_next) {
                    immediate::draw_point(p, COLORS[j % NUM_COLORS]);
                }
                if (pos_prev.count == 1 && pos_next.count == 1) {
                    immediate::draw_line(pos_prev[0], pos_next[0], LINE_COLOR);
                }
                if (pos_prev.count > 1) {
                    for (int32 k = 0; k < pos_prev.count; k++) {
                        immediate::draw_line(pos_prev[k], pos_next[0], LINE_COLOR);
                    }
                } else if (pos_next.count > 1) {
                    for (int32 k = 0; k < pos_next.count; k++) {
                        immediate::draw_line(pos_prev[0], pos_next[k], LINE_COLOR);
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
    ctx.property_func_entries.push_back({HASH("expression"), compute_expression, nullptr});

    ctx.structure_func_entries.push_back({HASH("resname"), structure_match_resname});
    ctx.structure_func_entries.push_back({HASH("resid"), structure_match_resid});
    ctx.structure_func_entries.push_back({HASH("residue"), structure_match_residue});
    ctx.structure_func_entries.push_back({HASH("chainid"), structure_match_chainid});
    ctx.structure_func_entries.push_back({HASH("chain"), structure_match_chain});
    ctx.structure_func_entries.push_back({HASH("atom"), structure_match_atom});
    ctx.structure_func_entries.push_back({HASH("resatom"), structure_extract_resatom});
    ctx.structure_func_entries.push_back({HASH("com"), structure_apply_aggregation_strategy_com});
}

void shutdown() {}

bool register_property_command(CString cmd_keyword, PropertyComputeFunc compute_func, PropertyVisualizeFunc visualize_func) {
    if (!cmd_keyword.data || cmd_keyword.count == 0) {
        LOG_ERROR("Property command cannot be an empty string!");
        return false;
    }

    if (contains_whitespace(cmd_keyword)) {
        LOG_ERROR("Property command cannot contain whitespace!");
        return false;
    }

    if (!compute_func) {
        LOG_ERROR("Property command must have compute function!");
        return false;
    }

    ID hash = HASH(cmd_keyword);
    if (find_property_func_entry(hash) != nullptr) {
        LOG_ERROR("Property command already registered!");
        return false;
    }

    ctx.property_func_entries.push_back({hash, compute_func, visualize_func});
    return true;
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
        set_error_message("Could not identify command: '%s'", make_tmp_str(cmd).cstr());
        return false;
    }

    DynamicArray<CString> cmd_args;
    CString outer = extract_parentheses(arg);
    if (outer) {
        outer = {outer.beg() + 1, outer.end() - 1};
        cmd_args = extract_arguments(outer);
    }

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

bool sync_structure_data_length(Array<StructureData> data) {
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

bool compute_stats(const MoleculeDynamic& dynamic) {
    int32 num_frames = dynamic.trajectory.num_frames;

    for (auto p : ctx.properties) {
        auto& prop = *p;
        ctx.current_property = p;
        if (prop.data.count != num_frames) {
            prop.data.resize(num_frames);
        }
        prop.data.set_mem_to_zero();
        clear_histogram(&prop.histogram);

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
            set_error_message("Could not recognize command '%s'", make_tmp_str(cmd).cstr());
            prop.valid = false;
            continue;
        }

        if (!func(&prop, args, dynamic)) {
            prop.valid = false;
            continue;
        }

		update_property(&prop);
        prop.valid = true;
    }
    ctx.current_property = nullptr;

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

void update_property(Property* prop, Range frame_range) {
	ASSERT(prop);
	for (auto& p : ctx.properties) {
		if (p == prop) {
			constexpr int32 NUM_BINS = 64;
			init_histogram(&p->histogram, NUM_BINS, p->data_range);
			clear_histogram(&p->histogram);

			if (frame_range.x == frame_range.y && frame_range.x == 0.f) {
				frame_range.x = 0;
				frame_range.y = prop->data.count;
			}

			int32 beg_idx = math::clamp((int32)frame_range.x, 0, (int32)prop->data.count);
			int32 end_idx = math::clamp((int32)frame_range.y, 0, (int32)prop->data.count);
			int32 offset = beg_idx;
			int32 count = end_idx - beg_idx;

			if (p->instance_data) {
				for (const auto& inst : p->instance_data) {
					compute_histogram(&p->histogram, inst.data.sub_array(offset, count));
				}
			}
			else {
				compute_histogram(&p->histogram, p->data.sub_array(offset, count));
			}
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

static void free_property_data(Property* prop) {
    ASSERT(prop);
    free_histogram(&prop->histogram);
    free_structure_data(&prop->structure_data);
    free_instance_data(&prop->instance_data);
}

void remove_property(Property* prop) {
    ASSERT(prop);
    for (auto& p : ctx.properties) {
        if (p == prop) {
            free_property_data(prop);
            ctx.properties.remove(&p);
        }
    }
}

void remove_all_properties() {
    for (auto prop : ctx.properties) {
        free_property_data(prop);
    }
    ctx.properties.clear();
}

void move_property_up(Property* prop) {
	if (ctx.properties.count <= 1) return;
	for (int32 i = 1; i < (int32)ctx.properties.count; i++) {
		if (ctx.properties[i] == prop) {
			// swap
			Property* tmp = ctx.properties[i - 1];
			ctx.properties[i - 1] = ctx.properties[i];
			ctx.properties[i] = tmp;
			break;
		}
	}
}

void move_property_down(Property* prop) {
	if (ctx.properties.count <= 1) return;
	for (int32 i = 0; i < (int32)ctx.properties.count - 1; i++) {
		if (ctx.properties[i] == prop) {
			// swap
			Property* tmp = ctx.properties[i + 1];
			ctx.properties[i + 1] = ctx.properties[i];
			ctx.properties[i] = tmp;
			break;
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

void clear_property(Property* prop) {
    ASSERT(prop);
    free_structure_data(&prop->structure_data);
    free_instance_data(&prop->instance_data);
    clear_histogram(&prop->histogram);
    prop->data.clear();
}

void clear_all_properties() {
    for (auto p : ctx.properties) {
        clear_property(p);
    }
}

}  // namespace stats
