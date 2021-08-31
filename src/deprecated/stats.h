#pragma once

#include <core/common.h>
#include <core/types.h>
#include <core/log.h>
#include <core/array_types.h>
#include <core/vector_types.h>
#include <core/math_utils.h>
#include <core/string_utils.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#include "volume.h"

struct MoleculeDynamic;
struct MoleculeStructure;
struct MoleculeTrajectory;

namespace stats {

typedef u64 ID;
constexpr ID INVALID_ID = 0;

struct Histogram {
    Array<float> bins = {};
    Range<float> value_range = {0, 0};
    Range<float> bin_range = {0, 0};
    i32 num_samples = 0;
};

struct Structure {
    i32 beg_idx;
    i32 end_idx;
};

enum AggregationStrategy { COM, NONE };

struct StructureData {
    DynamicArray<Structure> structures{};
    AggregationStrategy strategy = NONE;
};

struct InstanceData {
    InstanceData(i32 size) { data = allocate_array<float>(size); }
    ~InstanceData() { free_array(&data); }

    Array<float> data;
};

/*
typedef int PropertyFlags;
enum PropertyFlags_ {
    PropertyFlags_Valid             = BIT(1);
    PropertyFlags_Periodic          = BIT(2);
    PropertyFlags_ShowOverlay       = BIT(3);
    PropertyFlags_ShowTimeline      = BIT(4);
    PropertyFlags_ShowDistribution  = BIT(5);
    PropertyFlags_ShowVolume        = BIT(5);
    PropertyFlags_DataDirty         = BIT(6);
    PropertyFlags_FilterDirty       = BIT(7);
}
*/

struct Property {
    StringBuffer<32> name_buf{};
    StringBuffer<256> args_buf{};
    StringBuffer<32> unit_buf{};
    StringBuffer<256> error_msg_buf{};

    bool valid = false;
    bool periodic = false;

    bool enable_visualization = false;
    bool enable_timeline = true;
    bool enable_distribution = true;
    bool enable_volume = true;

    bool data_dirty = false;
    bool filter_dirty = false;

    Range<float> filter{0, 0};
    Range<float> total_data_range{0, 0};
    Range<float> avg_data_range{0, 0};
    DynamicArray<float> filter_fraction{};
    DynamicArray<float> avg_data{};
    DynamicArray<float> std_dev_data{};
    Histogram full_histogram{};
    Histogram filt_histogram{};

    Array<InstanceData> instance_data{};
    Array<StructureData> structure_data{};
    DynamicArray<Property*> dependencies{};
};

// This is to use as base colors for overlaying graphics
struct VisualizationStyle {
    static const i32 NUM_COLORS = 4;
    u32 point_colors[NUM_COLORS] = {0xffe3cea6, 0xffb4781f, 0xff8adfb2, 0xff2ca033};
    u32 line_color = 0x55cccccc;
};

// Helper functions
void set_error_message(Property* prop, const char* fmt, ...);

void init_instance_data(Array<InstanceData>* instance_data, i32 num_instances, i32 num_frames);
void free_instance_data(Array<InstanceData>* instance_data);

void init_structure_data(Array<StructureData>* structure_data, i32 count);
bool sync_structure_data_length(Array<StructureData> data);
bool extract_args_structures(Array<StructureData> data, Array<CStringView> arg, const MoleculeStructure& structure);

Array<const vec3> extract_positions(Structure structure, Array<const vec3> atom_positions);

template <typename Callback>
void for_each_filtered_property_structure_in_frame(Property* prop, i32 frame_idx, Callback cb) {
    ASSERT(prop);
    for (i32 inst_idx = 0; inst_idx < prop->instance_data.count; inst_idx++) {
        const float v = prop->instance_data[inst_idx].data[frame_idx];
        if (prop->filter.beg <= v && v <= prop->filter.end) {
            for (const auto& s_data : prop->structure_data) {
                for (const auto& s : s_data.structures) {
                    cb(s);
                }
            }
        }
    }
}

typedef bool (*PropertyComputeFunc)(Property* prop, const Array<CStringView> args, const MoleculeDynamic& dynamic);
typedef bool (*PropertyVisualizeFunc)(const Property& prop, const MoleculeDynamic& dynamic);

// HISTOGRAM
void init_histogram(Histogram* hist, i32 num_bins);
void free_histogram(Histogram* hist);

void compute_histogram(Histogram* hist, Array<const float> data);
void compute_histogram(Histogram* hist, Array<const float> data, Range<float> filter);

void clear_histogram(Histogram* hist);
void normalize_histogram(Histogram* hist);

// STATS
void initialize();
void shutdown();

void async_update(const MoleculeDynamic& dynamic, Range<i32> frame_filter = {0, 0}, void (*on_finished)(void*) = nullptr, void* usr_data = nullptr);

// Async functionality
// void lock_thread_mutex();
// void unlock_thread_mutex();
bool thread_running();
void signal_stop();
void signal_stop_and_wait();
void signal_start();
float fraction_done();

// bool compute_stats(const MoleculeDynamic& dynamic);
void visualize(const MoleculeDynamic& dynamic);

const Volume& get_density_volume();
VisualizationStyle* get_style();

// void compute_property(Property* prop, const MoleculeDynamic& dynamic);
// void compute_property_histograms(Property* prop);
// void compute_property_histograms(Property* prop, Range frame_filter);

bool register_property_command(CStringView cmd_keyword, PropertyComputeFunc compute_func, PropertyVisualizeFunc visualize_func);

Array<CStringView> get_property_commands();
Array<CStringView> get_structure_commands();
Array<CStringView> get_property_names();

// PROPERTY
Property* create_property(CStringView name = {}, CStringView args = {});
void remove_property(Property* prop);
void remove_all_properties();
void move_property_up(Property* prop);
void move_property_down(Property* prop);

// Keep property, but remove the generated data
void clear_property(Property* prop);
void clear_all_properties();

void set_all_property_flags(bool data_dirty, bool filter_dirty);

Array<Property*> get_properties();
Property* find_property(CStringView name);

// DENSITY VOLUME
void compute_density_volume(Volume* vol, const MoleculeTrajectory& traj, Range<i32> frame_range, const mat4& world_to_volume);


// @NOTE: Implement this using proper lightweight std::function wrapper which supports lambda,
// template is kind of obfuscating here
// the expected definition of the function is:
//
// vec4 TransformWorldToVolumeFunc(const vec4& world_coord, int32 frame_idx)

template <typename TransformWorldToVolumeFunc>
void compute_density_volume_with_basis(Volume* vol, const MoleculeTrajectory& traj, Range<i32> frame_range, TransformWorldToVolumeFunc func) {
    ASSERT(vol);
    if (vol->dim.x == 0 || vol->dim.y == 0 || vol->dim.z == 0) {
        LOG_ERROR("One or more volume dimension are zero...");
        return;
    }

    auto properties = get_properties();
    if (properties.count == 0) return;
    const i32 num_frames = (i32)properties.front()->avg_data.count;

    if (frame_range.beg == 0 && frame_range.end == 0) {
        frame_range.beg = num_frames;
    }

    frame_range.beg = math::clamp(frame_range.beg, 0, num_frames);
    frame_range.end = math::clamp(frame_range.end, 0, num_frames);

    clear_volume(vol);

    for (auto prop : get_properties()) {
        if (!prop->enable_volume) continue;
        for (i32 frame_idx = frame_range.beg; frame_idx < frame_range.end; frame_idx++) {
            const Array<const float> pos_x = get_trajectory_position_x(traj, frame_idx);
            const Array<const float> pos_y = get_trajectory_position_y(traj, frame_idx);
            const Array<const float> pos_z = get_trajectory_position_z(traj, frame_idx);
            //const mat4 world_to_volume_matrix = func(frame_idx);
            for_each_filtered_property_structure_in_frame(prop, frame_idx, [vol, frame_idx, &pos_x, &pos_y, &pos_z, &func](const Structure& s) {
                for (i32 i = s.beg_idx; i < s.end_idx; i++) {
                    const vec4 wc = {pos_x[i], pos_y[i], pos_z[i], 1.0f};
                    const vec4 vc = func(wc, frame_idx);
                    const vec3 tc = apply_pbc((vec3)vc);  // PBC
                    const ivec3 c = (ivec3)(tc * (vec3)vol->dim + 0.5f);
                    const i32 voxel_idx = c.z * vol->dim.x * vol->dim.y + c.y * vol->dim.x + c.x;
                    vol->voxel_data[voxel_idx]++;
                    vol->voxel_range.end = math::max(vol->voxel_range.end, vol->voxel_data[voxel_idx]);
                }
            });
        }
    }
}

}  // namespace stats
