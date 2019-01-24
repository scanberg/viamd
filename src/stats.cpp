#include "stats.h"
#include <core/math_utils.h>
#include <core/hash.h>
#include <core/log.h>
#include <core/volume.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/molecule_dynamic.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <ctype.h>
#include <stdarg.h>
#include <new>
#include <thread>
#include <mutex>
#include <tinyexpr.h>

#define COMPUTE_ID(x) (hash::crc64(x))
constexpr int32 NUM_BINS = 128;

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
    Volume volume{};
    bool volume_dirty = false;

    VisualizationStyle style;

    //   std::mutex thread_mutex;
    volatile bool thread_running = false;
    volatile bool stop_signal = false;
    volatile float fraction_done = 0.f;
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
    auto e = find_property_func_entry(COMPUTE_ID(cmd));
    if (e) return e->compute_func;

    return nullptr;
}

static PropertyVisualizeFunc find_property_visualize_func(CString cmd) {
    auto e = find_property_func_entry(COMPUTE_ID(cmd));
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
    auto e = find_structure_func_entry(COMPUTE_ID(cmd));
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
    static int32 count = 0;
    count++;
    if (instance_data->ptr) {
        for (auto& inst : *instance_data) {
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
void init_histogram(Histogram* hist, int32 num_bins) {
    ASSERT(hist);
    free_array(&hist->bins);
    hist->bins = allocate_array<float>(num_bins);
    hist->value_range = {0, 0};
    hist->bin_range = {0, 0};
    ASSERT(hist->bins);
}

void free_histogram(Histogram* hist) {
    ASSERT(hist);
    free_array(&hist->bins);
}

void compute_histogram(Histogram* hist, Array<const float> data) {
    ASSERT(hist);
    const int32 num_bins = (int32)hist->bins.count;
    const float scl = num_bins / (hist->value_range.y - hist->value_range.x);
    for (auto v : data) {
        int32 bin_idx = math::clamp((int32)((v - hist->value_range.x) * scl), 0, num_bins - 1);
        hist->bins[bin_idx]++;
    }
    hist->num_samples += (int32)data.count;
}

void compute_histogram(Histogram* hist, Array<const float> data, Range filter) {
    ASSERT(hist);
    const int32 num_bins = (int32)hist->bins.count;
    const float scl = num_bins / (hist->value_range.y - hist->value_range.x);
    for (auto v : data) {
        if (filter.x <= v && v <= filter.y) {
            int32 bin_idx = math::clamp((int32)((v - hist->value_range.x) * scl), 0, num_bins - 1);
            hist->bins[bin_idx]++;
            // hist->num_samples++;
        }
    }
    hist->num_samples += (int32)data.count;
}

void clear_histogram(Histogram* hist) {
    ASSERT(hist);
    if (hist->bins.ptr) {
        memset(hist->bins.ptr, 0, hist->bins.count * sizeof(float));
    }
    hist->bin_range = {0, 0};
    hist->num_samples = 0;
}

void normalize_histogram(Histogram* hist, int32 num_samples) {
    hist->bin_range = {0, 0};
    const float bin_scl = 1.f / (float)num_samples;
    for (auto& b : hist->bins) {
        b *= bin_scl;
        hist->bin_range.y = math::max(hist->bin_range.y, b);
    }
}

inline bool point_in_aabb(const vec3& p, const vec3& min_box, const vec3& max_box) {
    if (p.x < min_box.x || max_box.x < p.x) return false;
    if (p.y < min_box.y || max_box.y < p.y) return false;
    if (p.z < min_box.z || max_box.z < p.z) return false;
    return true;
}

// from here https://stackoverflow.com/questions/6127503/shuffle-array-in-c
/* Arrange the N elements of ARRAY in random order.
Only effective if N is much smaller than RAND_MAX;
if this may not be the case, use a better random
number generator. */
inline void shuffle(int32* arr, size_t n) {
    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            int32 t = arr[j];
            arr[j] = arr[i];
            arr[i] = t;
        }
    }
}

void compute_density_volume(Volume* vol, const mat4& world_to_volume_matrix, const MoleculeTrajectory& traj, Range frame_range) {
    ASSERT(vol);
    if (vol->dim.x == 0 || vol->dim.y == 0 || vol->dim.z == 0) {
        LOG_ERROR("One or more volume dimension are zero...");
        return;
    }

    if (ctx.properties.count == 0) return;
    float num_frames = (float)ctx.properties.front()->avg_data.count;

    if (frame_range.x == 0 && frame_range.y == 0) {
        frame_range.y = num_frames;
    }

    frame_range.x = math::clamp(frame_range.x, 0.f, num_frames);
    frame_range.y = math::clamp(frame_range.y, 0.f, num_frames);

    clear_volume(vol);

    for (auto prop : ctx.properties) {
        if (!prop->enable_volume) continue;
        for (int32 frame_idx = (int32)frame_range.x; frame_idx < (int32)frame_range.y; frame_idx++) {
            const Array<const vec3> atom_positions = get_trajectory_positions(traj, frame_idx);
            for_each_filtered_property_structure_in_frame(prop, frame_idx, [vol, &atom_positions, &world_to_volume_matrix](const Structure& s) {
                for (int32 i = s.beg_idx; i < s.end_idx; i++) {
                    const vec4 tc = world_to_volume_matrix * vec4(atom_positions[i], 1);
                    if (tc.x < 0.f || 1.f < tc.x) continue;
                    if (tc.y < 0.f || 1.f < tc.y) continue;
                    if (tc.z < 0.f || 1.f < tc.z) continue;
                    const ivec3 c = vec3(tc) * (vec3)vol->dim;
                    const int32 voxel_idx = c.z * vol->dim.x * vol->dim.y + c.y * vol->dim.x + c.x;
                    vol->voxel_data[voxel_idx]++;
                    vol->voxel_range.y = math::max(vol->voxel_range.y, vol->voxel_data[voxel_idx]);
                }
            });
        }
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
        range = compute_range(prop.avg_data);
    }
    return range;
}

void set_error_message(Property* prop, const char* fmt, ...) {
    ASSERT(prop);
    auto buf = prop->error_msg_buf.cstr();
    auto len = prop->error_msg_buf.capacity();
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, len, fmt, ap);
    va_end(ap);
    LOG_ERROR("Error when evaluating property '%s': %s", ctx.current_property->name_buf.cstr(), buf);
}

static DynamicArray<CString> extract_arguments(CString str) {
    DynamicArray<CString> args;

    const uint8* beg = str.beg();
    const uint8* end = str.beg();
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
    const uint8* ptr = str.beg();
    while (ptr != str.end() && *ptr != '(' && !isspace(*ptr)) ptr++;
    return {str.beg(), ptr};
}

bool extract_structures(StructureData* data, CString arg, const MoleculeStructure& molecule) {
    CString cmd = extract_command(arg);
    auto func = find_structure_func(cmd);

    if (!func) {
        set_error_message(ctx.current_property, "Could not identify command: '%s'", make_tmp_str(cmd).cstr());
        return false;
    }

    DynamicArray<CString> cmd_args;
    CString outer = extract_parentheses(arg);
    if (outer) {
        CString inner = {outer.beg() + 1, outer.end() - 1};
        cmd_args = extract_arguments(inner);
    }

    // @NOTE: ONLY ALLOW RECURSION FOR FIRST ARGUMENT?
    if (cmd_args.count > 0 && find_character(cmd_args[0], '(') != cmd_args[0].end()) {
        if (!extract_structures(data, cmd_args[0], molecule)) return false;
        cmd_args = cmd_args.subarray(1);
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
            set_error_message(ctx.current_property, "One argument did not match any structures");
            return false;
        }

        max_count = math::max(max_count, count);
        if (count > 1 && max_count > 1 && count != max_count) {
            set_error_message(ctx.current_property, "Multiple structures found for more than one argument, but the structure count did not match");
            return false;
        }
    }

    return true;
}

bool structure_match_resname(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule) {
    ASSERT(data);

    // Expect args.count to be > 0
    if (args.count == 0) {
        set_error_message(ctx.current_property, "Expects one or more arguments for resname");
        return false;
    }

    for (const auto& res : molecule.residues) {
        for (const auto& arg : args) {
            if (compare(res.name, arg)) {
                data->structures.push_back({res.atom_idx.beg, res.atom_idx.end});
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
        set_error_message(ctx.current_property, "Expects one or more arguments for resid");
        return false;
    }

    for (const auto& arg : args) {
        auto id = to_int(arg);
        if (!id.success) {
            set_error_message(ctx.current_property, "Failed to parse argument for resid");
            return false;
        }
        for (const auto& res : molecule.residues) {
            if (res.id == id) {
                data->structures.push_back({res.atom_idx.beg, res.atom_idx.end});
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
        set_error_message(ctx.current_property, "Expects one or more arguments for residue");
        return false;
    }

    for (const auto& arg : args) {
        IntRange range = {-1, -1};

        if (is_range(arg)) {
            if (!extract_range(&range, arg)) {
                set_error_message(ctx.current_property, "Failed to parse range in argument for residue");
                return false;
            }
            if (range.x == -1) range.x = 1;
            if (range.y == -1) range.y = (int32)molecule.residues.count;
        } else {
            auto id = to_int(arg);
            if (!id.success) {
                set_error_message(ctx.current_property, "Failed to parse argument for residue");
                return false;
            }
            range.x = range.y = id;
        }

        if (range.x < 1 || (int32)molecule.residues.count < range.y) {
            set_error_message(ctx.current_property, "Index for residue is out of bounds");
            return false;
        }
        for (int32 i = range.x - 1; i < range.y; i++) {
            const auto& res = molecule.residues[i];
            data->structures.push_back({res.atom_idx.beg, res.atom_idx.end});
        }
    }

    return true;
}

bool structure_match_chainid(StructureData* data, const Array<CString> args, const MoleculeStructure& molecule) {
    ASSERT(data);

    // Expect args.count to be > 0
    if (args.count == 0) {
        set_error_message(ctx.current_property, "Expects one or more arguments for chainid");
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
        set_error_message(ctx.current_property, "Expects one or more arguments for chain");
        return false;
    }

    for (const auto& arg : args) {
        IntRange range{0, 0};

        if (is_range(arg)) {
            if (!extract_range(&range, arg)) {
                set_error_message(ctx.current_property, "Failed to parse range in argument for chain");
                return false;
            }
            if (range.x == -1) range.x = 1;
            if (range.y == -1) range.y = (int32)molecule.chains.count;
        } else {
            auto id = to_int(arg);
            if (!id.success) {
                set_error_message(ctx.current_property, "Failed to parse argument for chain");
                return false;
            }
            range.x = range.y = id;
        }

        if (range.x < 1 || (int32)molecule.chains.count < range.y) {
            set_error_message(ctx.current_property, "Index for chain is out of bounds");
            return false;
        }
        for (int32 i = range.x - 1; i <= range.y - 1; i++) {
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
        set_error_message(ctx.current_property, "Expects one or more arguments for atom");
        return false;
    }

    for (const auto& arg : args) {
        IntRange range{0, 0};

        if (is_range(arg)) {
            if (!extract_range(&range, arg)) {
                set_error_message(ctx.current_property, "Failed to parse range in argument for atom");
                return false;
            }
            if (range.x == -1) range.x = 1;
            if (range.y == -1) range.y = (int32)molecule.atom.count;
        } else {
            auto id = to_int(arg);
            if (!id.success) {
                set_error_message(ctx.current_property, "Failed to parse argument for atom");
                return false;
            }
            range.x = range.y = id;
        }

        if (range.x < 1 || (int32)molecule.atom.count < range.y) {
            set_error_message(ctx.current_property, "Index for atom is out of bounds");
            return false;
        }

        data->structures.push_back({range.x - 1, range.y});
    }

    return true;
}

bool structure_extract_resatom(StructureData* data, const Array<CString> args, const MoleculeStructure&) {
    ASSERT(data);
    if (args.count != 1) {
        set_error_message(ctx.current_property, "resatom requires exactly 1 argument");
        return false;
    }

    IntRange range{0, 0};

    if (is_range(args[0])) {
        if (!extract_range(&range, args[0])) {
            set_error_message(ctx.current_property, "Failed to parse range in argument for resatom");
            return false;
        }
    } else {
        auto id = to_int(args[0]);
        if (!id.success) {
            set_error_message(ctx.current_property, "Failed to parse argument for resatom");
            return false;
        }
        range.x = range.y = id;
    }

    for (auto& s : data->structures) {
        int32 count = s.end_idx - s.beg_idx;

        int32 s_first = (range.x == -1) ? 1 : range.x;
        int32 s_last = (range.y == -1) ? count : range.y;

        if (count < 0 || s_first < 1 || count < s_last) {
            set_error_message(ctx.current_property, "restom: Index is out of range for structure");
            return false;
        }

        int new_beg = s.beg_idx + s_first - 1;
        int new_end = s.beg_idx + s_last;
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
static inline int32 structure_index_count(Structure s) { return s.end_idx - s.beg_idx; }

static inline int32 structures_index_count(Array<const Structure> structures) {
    int32 count = 0;
    for (const auto& s : structures) {
        count += structure_index_count(s);
    }
    return count;
}

inline Array<const vec3> extract_positions(Structure structure, Array<const vec3> atom_positions) { return atom_positions.subarray(structure.beg_idx, structure.end_idx - structure.beg_idx); }

static inline float multi_distance(Array<const vec3> pos_a, Array<const vec3> pos_b, float* variance = nullptr) {
    if (pos_a.count == 0 || pos_b.count == 0)
        return 0.f;
    else if (pos_a.count == 1 && pos_b.count == 1)
        return math::distance(pos_a[0], pos_b[0]);
    else {
        float total = 0.f;
        for (const auto& a : pos_a) {
            for (const auto& b : pos_b) {
                total += math::distance(a, b);
            }
        }
        float count = (float)(pos_a.count * pos_b.count);
        float mean = total / count;

        if (variance) {
            *variance = 0.f;
            for (const auto& a : pos_a) {
                for (const auto& b : pos_b) {
                    float x = math::distance(a, b) - mean;
                    *variance += x * x;
                }
            }
            *variance = *variance / count;
        }

        return mean;
    }
}

static inline float multi_angle(Array<const vec3> pos_a, Array<const vec3> pos_b, Array<const vec3> pos_c, float* variance = nullptr) {
    if (pos_a.count == 0 || pos_b.count == 0 || pos_c.count == 0)
        return 0.f;
    else if (pos_a.count == 1 && pos_b.count == 1 && pos_c.count == 1)
        return math::angle(pos_a[0], pos_b[0], pos_c[0]);
    else {
        float total = 0.f;
        for (const auto& a : pos_a) {
            for (const auto& b : pos_b) {
                for (const auto& c : pos_c) {
                    total += math::angle(a, b, c);
                }
            }
        }

        float count = (float)(pos_a.count * pos_b.count * pos_c.count);
        float mean = total / count;

        if (variance) {
            *variance = 0.f;
            for (const auto& a : pos_a) {
                for (const auto& b : pos_b) {
                    for (const auto& c : pos_c) {
                        float x = math::angle(a, b, c) - mean;
                        *variance += x * x;
                    }
                }
            }
            *variance = *variance / count;
        }

        return mean;

        /*
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
        */
    }
}

static inline float multi_dihedral(Array<const vec3> pos_a, Array<const vec3> pos_b, Array<const vec3> pos_c, Array<const vec3> pos_d, float* variance = nullptr) {
    if (pos_a.count == 0 || pos_b.count == 0 || pos_c.count == 0 || pos_d.count == 0)
        return 0.f;
    else if (pos_a.count == 1 && pos_b.count == 1 && pos_c.count == 1 && pos_d.count == 1)
        return math::dihedral_angle(pos_a[0], pos_b[0], pos_c[0], pos_d[0]);
    else {
        float total = 0.f;
        for (const auto& a : pos_a) {
            for (const auto& b : pos_b) {
                for (const auto& c : pos_c) {
                    for (const auto& d : pos_d) {
                        float x = math::dihedral_angle(a, b, c, d);
                        total += x * x;
                    }
                }
            }
        }

        float count = (float)(pos_a.count * pos_b.count * pos_c.count);
        float mean = total / count;

        // @TODO: This is madness, just use a temporary array and reuse the values
        if (variance) {
            *variance = 0.f;
            for (const auto& a : pos_a) {
                for (const auto& b : pos_b) {
                    for (const auto& c : pos_c) {
                        for (const auto& d : pos_d) {
                            *variance += math::dihedral_angle(a, b, c, d) - mean;
                        }
                    }
                }
            }
            *variance = *variance / count;
        }

        return mean;

        /*
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
        */
    }
}

static bool compute_distance(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic) {
    ASSERT(prop);
    if (args.count != 2) {
        set_error_message(ctx.current_property, "distance expects 2 arguments");
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

    // @IMPORTANT! Use this instead of dynamic.trajectory.num_frames as that can be changing in a different thread!
    const int32 num_frames = (int32)prop->avg_data.count;
    const int32 structure_count = (int32)prop->structure_data[0].structures.count;

    init_instance_data(&prop->instance_data, structure_count, num_frames);

    const float32 scl = 1.f / (float32)structure_count;
    Array<const vec3> pos[2];
    vec3 com[2];
    for (int32 i = 0; i < num_frames; i++) {
        float sum = 0.f;
        float var = 0.f;
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

            float variance = 0.f;
            prop->instance_data[j].data[i] = multi_distance(pos[0], pos[1], &variance);
            sum += prop->instance_data[j].data[i];
            var += variance;
        }

        prop->avg_data[i] = sum * scl;
        if (var > 0.f) prop->std_dev_data.ptr[i] = math::sqrt(var * scl);
    }

    prop->total_data_range = compute_range(*prop);
    prop->avg_data_range = compute_range(prop->avg_data);
    prop->periodic = false;
    prop->unit_buf = "Å";

    return true;
}

static bool compute_angle(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic) {
    ASSERT(prop);
    if (args.count != 3) {
        set_error_message(ctx.current_property, "angle expects 3 arguments");
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

    // @IMPORTANT! Use this instead of dynamic.trajectory.num_frames as that can be changing in a different thread!
    const int32 num_frames = (int32)prop->avg_data.count;
    const int32 structure_count = (int32)prop->structure_data[0].structures.count;

    init_instance_data(&prop->instance_data, structure_count, num_frames);

    const float32 scl = 1.f / (float32)structure_count;
    Array<const vec3> pos[3];
    vec3 com[3];
    for (int32 i = 0; i < num_frames; i++) {
        float sum = 0.f;
        float var = 0.f;
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

            float variance = 0.f;
            prop->instance_data[j].data[i] = multi_angle(pos[0], pos[1], pos[2], &variance);
            sum += prop->instance_data[j].data[i];
            var += variance;
        }

        prop->avg_data[i] = sum * scl;
        if (var > 0.f) prop->std_dev_data[i] = math::sqrt(var * scl);
    }

    prop->total_data_range = {0, math::PI};
    prop->avg_data_range = compute_range(prop->avg_data);
    prop->periodic = true;
    prop->unit_buf = u8"°";

    return true;
}

static bool compute_dihedral(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic) {
    ASSERT(prop);
    if (args.count != 4) {
        set_error_message(ctx.current_property, "dihedral expects 4 arguments");
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

    // @IMPORTANT! Use this instead of dynamic.trajectory.num_frames as that can be changing in a different thread!
    const int32 num_frames = (int32)prop->avg_data.count;
    const int32 structure_count = (int32)prop->structure_data[0].structures.count;

    init_instance_data(&prop->instance_data, structure_count, num_frames);

    const float32 scl = 1.f / (float32)structure_count;
    Array<const vec3> pos[4];
    vec3 com[4];
    for (int32 i = 0; i < num_frames; i++) {
        float sum = 0.f;
        float var = 0.f;
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

            float variance = 0.f;
            prop->instance_data[j].data[i] = multi_dihedral(pos[0], pos[1], pos[2], pos[3], &variance);
            sum += prop->instance_data[j].data[i];
            var += variance;
        }

        prop->avg_data[i] = sum * scl;
        if (var > 0.f) prop->std_dev_data[i] = math::sqrt(var * scl);
    }

    prop->total_data_range = {-math::PI, math::PI};
    prop->avg_data_range = compute_range(prop->avg_data);
    prop->periodic = true;
    prop->unit_buf = u8"°";

    return true;
}

static DynamicArray<Property*> extract_property_dependencies(Array<Property*> properties, CString expression) {
    DynamicArray<Property*> dependencies;
    for (Property* prop : properties) {
        if (!prop->valid) continue;
        CString match = find_first_match(expression, prop->name_buf);
        if (match) {
            if (match.beg() != expression.beg()) {
                char pre_char = *(match.beg() - 1);
                if (isalpha(pre_char)) continue;
                if (isdigit(pre_char)) continue;
            }
            if (match.end() != expression.end()) {
                char post_char = *match.end();
                if (isalpha(post_char)) continue;
                if (isdigit(post_char)) continue;
            }
            dependencies.push_back(prop);
        }
    }
    return dependencies;
}

static bool compute_expression(Property* prop, const Array<CString> args, const MoleculeDynamic&) {
    ASSERT(prop);
    if (args.count == 0) {
        set_error_message(ctx.current_property, "expression expects 1 or more arguments");
        return false;
    }

    // Concatenate all arguments
    auto expr_str = make_tmp_str(CString(args.front().beg(), args.back().end()));

    if (!expr_str) {
        return false;
    }

    if (!balanced_parentheses(expr_str)) {
        set_error_message(ctx.current_property, "Expression contains unbalanced parentheses!");
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

    prop->dependencies = extract_property_dependencies(properties, expr_str);

    DynamicArray<double> values(prop->dependencies.count, 0);
    DynamicArray<te_variable> vars;
    for (int32 i = 0; i < prop->dependencies.count; i++) {
        vars.push_back({prop->dependencies[i]->name_buf.cstr(), &values[i], 0, 0});
    }

    int err;
    te_expr* expr = te_compile(expr_str.cstr(), vars.ptr, (int32)vars.count, &err);

    if (expr) {
        int32 max_instance_count = 0;
        for (auto* p : prop->dependencies) {
            max_instance_count = math::max(max_instance_count, (int32)p->instance_data.count);
        }

        const int32 frame_count = prop->avg_data.count;
        init_instance_data(&prop->instance_data, max_instance_count, frame_count);

        float scl = 1.f / (float)max_instance_count;
        for (int32 frame = 0; frame < frame_count; frame++) {
            float val = 0.f;
            for (int32 i = 0; i < max_instance_count; i++) {
                for (int32 j = 0; j < values.count; j++) {
                    values[j] = 0;
                    if (prop->dependencies[j]->instance_data.count == max_instance_count) {
                        if (frame < prop->dependencies[j]->instance_data[i].data.count) {
                            values[j] = prop->dependencies[j]->instance_data[i].data[frame];
                        }
                    } else {
                        if (frame < prop->dependencies[j]->instance_data[0].data.count) {
                            values[j] = prop->dependencies[j]->instance_data[0].data[frame];
                        }
                    }
                }
                prop->instance_data[i].data[frame] = (float)te_eval(expr);
                val += prop->instance_data[i].data[frame];
            }
            prop->avg_data[frame] = val * scl;
        }

        te_free(expr);
    } else {
        set_error_message(ctx.current_property, "Malformed expression:\n%s\n%*s^\nError near here", expr_str.cstr(), err - 1, "");
        return false;
    }

    prop->total_data_range = compute_range(*prop);
    prop->avg_data_range = compute_range(prop->avg_data);
    prop->periodic = false;
    prop->unit_buf = "";

    return true;
}

static bool visualize_structures(const Property& prop, const MoleculeDynamic& dynamic) {
    if (prop.structure_data.count == 0) return false;

    if (prop.structure_data.count == 1) {
        for (const auto& s : prop.structure_data[0].structures) {
            Array<const vec3> pos = extract_positions(s, get_positions(dynamic.molecule));
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

        for (int32 i = 0; i < count; i++) {
            pos_prev = extract_positions(prop.structure_data[0].structures[i], get_positions(dynamic.molecule));
            if (prop.structure_data[0].strategy == COM) {
                com_prev = compute_com(pos_prev);
                pos_prev = {&com_prev, 1};
            }
            for (const auto& p : pos_prev) {
                immediate::draw_point(p, ctx.style.point_colors[0]);
            }
            for (int32 j = 1; j < prop.structure_data.count; j++) {
                const int32 col_idx = j % VisualizationStyle::NUM_COLORS;
                pos_next = extract_positions(prop.structure_data[j].structures[i], get_positions(dynamic.molecule));
                if (prop.structure_data[j].strategy == COM) {
                    com_next = compute_com(pos_next);
                    pos_next = {&com_next, 1};
                }
                for (const auto& p : pos_next) {
                    immediate::draw_point(p, ctx.style.point_colors[col_idx]);
                }
                if (pos_prev.count == 1 && pos_next.count == 1) {
                    immediate::draw_line(pos_prev[0], pos_next[0], ctx.style.line_color);
                }
                if (pos_prev.count > 1 && pos_next.count == 1) {
                    for (const auto& pp : pos_prev) {
                        immediate::draw_line(pp, pos_next[0]), ctx.style.line_color;
                    }
                } else if (pos_next.count > 1 && pos_prev.count == 1) {
                    for (const auto& pn : pos_next) {
                        immediate::draw_line(pos_prev[0], pn, ctx.style.line_color);
                    }
                } else {
                    // N^2 :'(
                    for (const auto& pp : pos_prev) {
                        for (const auto& pn : pos_next) {
                            immediate::draw_line(pp, pn, ctx.style.line_color);
                        }
                    }
                }

                pos_prev = pos_next;
                com_prev = com_next;
            }
        }
    }

    return true;
}

static bool visualize_dependencies(const Property& prop, const MoleculeDynamic& dynamic) {
    for (Property* p : prop.dependencies) {
        CString cmd = extract_command(p->args_buf);
        PropertyFuncEntry* entry = find_property_func_entry(COMPUTE_ID(cmd));
        if (entry && entry->visualize_func) {
            entry->visualize_func(*p, dynamic);
        }
    }

    return true;
}

void initialize() {
    ctx.property_func_entries.push_back({COMPUTE_ID("distance"), compute_distance, visualize_structures});
    ctx.property_func_entries.push_back({COMPUTE_ID("angle"), compute_angle, visualize_structures});
    ctx.property_func_entries.push_back({COMPUTE_ID("dihedral"), compute_dihedral, visualize_structures});
    ctx.property_func_entries.push_back({COMPUTE_ID("expression"), compute_expression, visualize_dependencies});

    ctx.structure_func_entries.push_back({COMPUTE_ID("resname"), structure_match_resname});
    ctx.structure_func_entries.push_back({COMPUTE_ID("resid"), structure_match_resid});
    ctx.structure_func_entries.push_back({COMPUTE_ID("residue"), structure_match_residue});
    ctx.structure_func_entries.push_back({COMPUTE_ID("chainid"), structure_match_chainid});
    ctx.structure_func_entries.push_back({COMPUTE_ID("chain"), structure_match_chain});
    ctx.structure_func_entries.push_back({COMPUTE_ID("atom"), structure_match_atom});
    ctx.structure_func_entries.push_back({COMPUTE_ID("resatom"), structure_extract_resatom});
    ctx.structure_func_entries.push_back({COMPUTE_ID("com"), structure_apply_aggregation_strategy_com});
}

void shutdown() {}

bool register_property_command(CString cmd_keyword, PropertyComputeFunc compute_func, PropertyVisualizeFunc visualize_func) {
    if (!cmd_keyword.ptr || cmd_keyword.count == 0) {
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

    ID hash = COMPUTE_ID(cmd_keyword);
    if (find_property_func_entry(hash) != nullptr) {
        LOG_ERROR("Property command already registered!");
        return false;
    }

    ctx.property_func_entries.push_back({hash, compute_func, visualize_func});
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
                    set_error_message(ctx.current_property, "Structures matched has different sizes in different arguments, this is not supported");
                    return false;
                }
            }
        }
    }

    return true;
}

static bool compute_property_data(Property* prop, const MoleculeDynamic& dynamic, int32 num_frames) {
    ASSERT(prop);

    ctx.current_property = prop;
    prop->error_msg_buf = "";
    if (prop->avg_data.count != num_frames) {
        prop->avg_data.resize(num_frames);
    }
    if (prop->std_dev_data.count != num_frames) {
        prop->std_dev_data.resize(num_frames);
    }
    if (prop->filter_fraction.count != num_frames) {
        prop->filter_fraction.resize(num_frames);
    }
    prop->avg_data.set_mem_to_zero();
    prop->std_dev_data.set_mem_to_zero();
    prop->filter_fraction.set_mem_to_zero();
    clear_histogram(&prop->full_histogram);
    clear_histogram(&prop->filt_histogram);

    for (const auto p : ctx.properties) {
        if (p != prop && compare(p->name_buf, prop->name_buf)) {
            set_error_message(prop, "A property is already defined with that name!");
            prop->valid = false;
            return false;
        }
        if (p == prop) break;
    }

    if (!balanced_parentheses(prop->args_buf)) {
        set_error_message(prop, "Unbalanced parantheses!");
        prop->valid = false;
        return false;
    }

    DynamicArray<CString> args;

    // Extract big argument chunks
    const uint8* beg = (uint8*)prop->args_buf.beg();
    const uint8* end = (uint8*)prop->args_buf.beg();
    int count = 0;

    // Use space separation unless we are inside a parenthesis
    while (end != prop->args_buf.end()) {
        if (*end == '(')
            count++;
        else if (*end == ')')
            count--;
        else if (count == 0 && isspace(*end)) {
            CString arg = trim({beg, end});
            if (arg.size() > 0) {
                args.push_back(trim({beg, end}));
            }
            beg = end + 1;
        }
        end++;
    }
    if (beg != end) {
        CString arg = trim({beg, end});
        if (arg.size() > 0) {
            args.push_back(trim({beg, end}));
        }
    }

    if (args.count == 0) {
        prop->valid = false;
        return false;
    }

    CString cmd = args[0];
    args = args.subarray(1);

    auto func = find_property_compute_func(cmd);
    if (!func) {
        set_error_message(prop, "Could not recognize command '%.*s'", cmd.length(), cmd.cstr());
        prop->valid = false;
        return false;
    }

    Range pre_range = prop->total_data_range;
    if (!func(prop, args, dynamic)) {
        prop->valid = false;
        return false;
    }

    if (prop->instance_data.count > 1) {
        if (prop->std_dev_data[0] == 0.f) {
            // Compute std dev of instances
            const float scl = 1.f / (float)prop->instance_data.count;
            for (int32 i = 0; i < num_frames; i++) {
                float sum = 0.f;
                for (const auto& inst : prop->instance_data) {
                    float x = inst.data[i] - prop->avg_data[i];
                    sum += x * x;
                }
                prop->std_dev_data[i] = math::sqrt(sum * scl);
            }
        }
    }

    if (pre_range != prop->total_data_range) {
        prop->filter = prop->total_data_range;
    }

    prop->full_histogram.value_range = prop->total_data_range;
    prop->filt_histogram.value_range = prop->total_data_range;
    prop->filter_dirty = true;
    prop->valid = true;

    ctx.current_property = nullptr;

    return true;
}

bool properties_dirty() {
    for (const auto p : ctx.properties) {
        if (p->data_dirty) return true;
        if (p->filter_dirty) return true;
    }

    return false;
}

void async_update(const MoleculeDynamic& dynamic, Range frame_filter, void (*on_finished)(void*), void* usr_data) {
    if (frame_filter.x == 0.f && frame_filter.y == 0.f) {
        frame_filter.y = (float)dynamic.trajectory.num_frames;
    }

    if (!dynamic) return;
    if (dynamic.trajectory.num_frames == 0) return;

    static Range prev_frame_filter{0, 0};

    if ((prev_frame_filter != frame_filter || properties_dirty()) && !ctx.thread_running && !ctx.stop_signal) {
        prev_frame_filter = frame_filter;
        ctx.thread_running = true;
        std::thread([&dynamic, frame_filter, on_finished, usr_data]() {
            Histogram tmp_hist;
            init_histogram(&tmp_hist, NUM_BINS);
            ctx.fraction_done = 0.f;

            // @NOTE IMPORTANT: This is the one 'true' frame count which should be used for properties.
            // It is important that this is used so that all properties data lengths are in sync
            // When dealing with dependencies.
            const int32 num_frames = dynamic.trajectory.num_frames;

            for (int32 i = 0; i < ctx.properties.count; i++) {
                auto p = ctx.properties[i];
                ctx.fraction_done = (i / (float)ctx.properties.count);
                auto filter = p->filter;

                if (p->data_dirty) {
                    compute_property_data(p, dynamic, num_frames);

                    // recompute full histogram
                    clear_histogram(&p->full_histogram);
                    if (p->instance_data) {
                        for (const auto& inst : p->instance_data) {
                            compute_histogram(&p->full_histogram, inst.data);
                        }
                    } else {
                        compute_histogram(&p->full_histogram, p->avg_data);
                    }
                    normalize_histogram(&p->full_histogram, p->full_histogram.num_samples);

                    p->data_dirty = false;
                }

                if (ctx.stop_signal) break;

                if (p->filter_dirty) {
                    clear_histogram(&tmp_hist);
                    tmp_hist.value_range = p->filt_histogram.value_range;

                    int32 beg_idx = math::clamp((int32)frame_filter.x, 0, (int32)p->avg_data.count);
                    int32 end_idx = math::clamp((int32)frame_filter.y, 0, (int32)p->avg_data.count);

                    if (beg_idx != end_idx) {
                        // Since the data is probably showing, perform the operations on tmp data then copy the results
                        if (p->instance_data) {
                            for (const auto& inst : p->instance_data) {
                                compute_histogram(&tmp_hist, inst.data.subarray(beg_idx, end_idx - beg_idx));
                            }
                        } else {
                            compute_histogram(&tmp_hist, p->avg_data.subarray(beg_idx, end_idx - beg_idx));
                        }
                        normalize_histogram(&tmp_hist, tmp_hist.num_samples);
                        p->filt_histogram.bin_range = tmp_hist.bin_range;
                        p->filt_histogram.num_samples = tmp_hist.num_samples;
                        memcpy(p->filt_histogram.bins.ptr, tmp_hist.bins.ptr, p->filt_histogram.bins.size_in_bytes());
                    }

                    // Compute filter fractions for frames
                    if (p->instance_data) {
                        for (int32 j = 0; j < p->filter_fraction.count; j++) {
                            float val = 0.f;
                            for (const auto& inst : p->instance_data) {
                                if (filter.x <= inst.data[j] && inst.data[j] <= filter.y) {
                                    val += 1.f;
                                }
                            }
                            p->filter_fraction[j] = val / (float)p->instance_data.count;
                        }
                    }

                    if (p->filter == filter) p->filter_dirty = false;
                }

                if (ctx.stop_signal) break;
            }

            if (on_finished) {
                on_finished(usr_data);
            }

            free_histogram(&tmp_hist);
            ctx.fraction_done = 1.f;
            ctx.thread_running = false;
            ctx.stop_signal = false;
        }).detach();
    }
}

// void lock_thread_mutex() { ctx.thread_mutex.lock(); }

// void unlock_thread_mutex() { ctx.thread_mutex.unlock(); }

bool thread_running() { return ctx.thread_running; }

void signal_stop() { ctx.stop_signal = true; }

void signal_stop_and_wait() {
    ctx.stop_signal = true;
    while (ctx.thread_running) {
        // possibly sleep 1 ms or so
    }
    ctx.stop_signal = false;
}

void signal_start() { ctx.stop_signal = false; }

float fraction_done() { return ctx.fraction_done; }

void visualize(const MoleculeDynamic& dynamic) {
    if (thread_running()) return;
    for (auto p : ctx.properties) {
        if (!p->enable_visualization) continue;
        CString cmd = extract_command(p->args_buf);
        auto entry = find_property_func_entry(COMPUTE_ID(cmd));
        if (entry && entry->visualize_func) {
            entry->visualize_func(*p, dynamic);
        }
    }
}

const Volume& get_density_volume() { return ctx.volume; }

VisualizationStyle* get_style() { return &ctx.style; }

Property* create_property(CString name, CString args) {
    Property* prop = (Property*)MALLOC(sizeof(Property));
    new (prop) Property();
    prop->name_buf = name;
    prop->args_buf = args;
    prop->valid = false;
    prop->data_dirty = true;
    prop->filter_dirty = true;

    init_histogram(&prop->full_histogram, NUM_BINS);
    clear_histogram(&prop->full_histogram);

    init_histogram(&prop->filt_histogram, NUM_BINS);
    clear_histogram(&prop->filt_histogram);

    ctx.properties.push_back(prop);
    return prop;
}

static void free_property_data(Property* prop) {
    signal_stop_and_wait();
    ASSERT(prop);
    free_histogram(&prop->full_histogram);
    free_histogram(&prop->filt_histogram);
    free_structure_data(&prop->structure_data);
    free_instance_data(&prop->instance_data);
    ctx.stop_signal = false;
}

void remove_property(Property* prop) {
    ASSERT(prop);
    signal_stop_and_wait();
    for (auto& p : ctx.properties) {
        if (p == prop) {
            free_property_data(prop);
            ctx.properties.remove(&p);
        }
    }
    ctx.stop_signal = false;
}

void remove_all_properties() {
    signal_stop_and_wait();
    for (auto prop : ctx.properties) {
        free_property_data(prop);
    }
    ctx.properties.clear();
    ctx.stop_signal = false;
}

void move_property_up(Property* prop) {
    if (ctx.properties.count <= 1) return;
    signal_stop_and_wait();
    for (int32 i = 1; i < (int32)ctx.properties.count; i++) {
        if (ctx.properties[i] == prop) {
            // swap
            Property* tmp = ctx.properties[i - 1];
            ctx.properties[i - 1] = ctx.properties[i];
            ctx.properties[i] = tmp;
            break;
        }
    }
    ctx.stop_signal = false;
}

void move_property_down(Property* prop) {
    if (ctx.properties.count <= 1) return;
    signal_stop_and_wait();
    for (int32 i = 0; i < (int32)ctx.properties.count - 1; i++) {
        if (ctx.properties[i] == prop) {
            // swap
            Property* tmp = ctx.properties[i + 1];
            ctx.properties[i + 1] = ctx.properties[i];
            ctx.properties[i] = tmp;
            break;
        }
    }
    ctx.stop_signal = false;
}

Array<Property*> get_properties() { return ctx.properties; }

Property* find_property(CString name) {
    for (auto p : ctx.properties) {
        if (compare(p->name_buf, name)) return p;
    }
    return nullptr;
}

void clear_property(Property* prop) {
    signal_stop_and_wait();
    ASSERT(prop);
    for (auto p : ctx.properties) {
        if (p == prop) {
            free_structure_data(&prop->structure_data);
            free_instance_data(&prop->instance_data);
            clear_histogram(&prop->full_histogram);
            prop->avg_data.clear();
        }
    }
    ctx.stop_signal = false;
}

void clear_all_properties() {
    signal_stop_and_wait();
    for (auto p : ctx.properties) {
        free_structure_data(&p->structure_data);
        free_instance_data(&p->instance_data);
        clear_histogram(&p->full_histogram);
        p->avg_data.clear();
    }
    ctx.stop_signal = false;
}

void set_all_property_flags(bool data_dirty, bool filter_dirty) {
    for (auto p : ctx.properties) {
        p->data_dirty |= data_dirty;
        p->filter_dirty |= filter_dirty;
    }
}

}  // namespace stats
