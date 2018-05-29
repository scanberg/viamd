#pragma once

#include <core/types.h>
#include <core/string_utils.h>
#include <mutex>

struct MoleculeDynamic;
struct MoleculeStructure;
struct MoleculeTrajectory;
struct Volume;

namespace stats {

typedef uint64 ID;
constexpr ID INVALID_ID = 0;

struct Histogram {
    Array<float> bins = {};
    Range value_range = {};
    Range bin_range = {};
    int32 num_samples = 0;
};

struct Structure {
    int32 beg_idx;
    int32 end_idx;
};

enum AggregationStrategy { COM, NONE };

struct StructureData {
    DynamicArray<Structure> structures{};
    AggregationStrategy strategy = NONE;
};

struct InstanceData {
    InstanceData(int32 size) { data = allocate_array<float>(size); }
    ~InstanceData() { free_array(&data); }

    Array<float> data;
};

struct Property {
    StringBuffer<32> name_buf{};
    StringBuffer<256> args_buf{};
    StringBuffer<32> unit_buf{};
    StringBuffer<256> error_msg_buf{};

    bool valid = false;
    bool periodic = false;

    bool visualize = false;
    bool show_as_timeline = true;
    bool show_as_distribution = true;

    bool data_dirty = false;
    bool full_hist_dirty = false;
    bool filt_hist_dirty = false;

    Range filter{0, 0};
    Range data_range{0, 0};
    DynamicArray<float> filter_fraction{};
    DynamicArray<float> avg_data{};
    DynamicArray<float> std_dev_data{};
    Histogram full_histogram{};
    Histogram filt_histogram{};

    Array<InstanceData> instance_data{};
    Array<StructureData> structure_data{};
    DynamicArray<Property*> dependencies{};
};

// Helper functions
void set_error_message(const char* fmt, ...);

void init_instance_data(Array<InstanceData>* instance_data, int32 num_instances, int32 num_frames);
void free_instance_data(Array<InstanceData>* instance_data);

void init_structure_data(Array<StructureData>* structure_data, int32 count);
bool sync_structure_data_length(Array<StructureData> data);
bool extract_args_structures(Array<StructureData> data, Array<CString> arg, const MoleculeStructure& structure);

Array<const vec3> extract_positions(Structure structure, Array<const vec3> atom_positions);
vec3 compute_com(Array<const vec3> positions);

typedef bool (*PropertyComputeFunc)(Property* prop, const Array<CString> args, const MoleculeDynamic& dynamic);
typedef bool (*PropertyVisualizeFunc)(const Property& prop, const MoleculeDynamic& dynamic);

// HISTOGRAM
void init_histogram(Histogram* hist, int32 num_bins);
void free_histogram(Histogram* hist);

void compute_histogram(Histogram* hist, Array<const float> data);
void compute_histogram(Histogram* hist, Array<const float> data, Range filter);

void clear_histogram(Histogram* hist);
void normalize_histogram(Histogram* hist);

// DENSITY VOLUME
void compute_density_volume(Volume* vol, const mat4& world_to_volume, const MoleculeTrajectory& traj, Range frame_range);

// STATS
void initialize();
void shutdown();

void async_update(const MoleculeDynamic& dynamic, Range frame_filter = {0, 0}, void (*on_finished)(void*) = nullptr, void* usr_data = nullptr);

// ASync functionality
void lock_thread_mutex();
void unlock_thread_mutex();
bool thread_running();
void signal_stop();
void signal_stop_and_wait();
void signal_start();
float fraction_done();

// bool compute_stats(const MoleculeDynamic& dynamic);
void visualize(const MoleculeDynamic& dynamic);

const Volume& get_density_volume();

// void compute_property(Property* prop, const MoleculeDynamic& dynamic);
// void compute_property_histograms(Property* prop);
// void compute_property_histograms(Property* prop, Range frame_filter);

bool register_property_command(CString cmd_keyword, PropertyComputeFunc compute_func, PropertyVisualizeFunc visualize_func);

Array<CString> get_property_commands();
Array<CString> get_structure_commands();
Array<CString> get_property_names();

// PROPERTY
Property* create_property(CString name = {}, CString args = {});
void remove_property(Property* prop);
void remove_all_properties();
void move_property_up(Property* prop);
void move_property_down(Property* prop);

// Keep property, but remove the generated data
void clear_property(Property* prop);
void clear_all_properties();

void set_all_property_flags(bool data_dirty, bool full_hist_dirty, bool filt_hist_dirty);

int32 get_property_count();
Property* get_property(int32 idx);
Property* get_property(CString name);

}  // namespace stats
