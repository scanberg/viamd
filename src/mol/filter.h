#pragma once

#include <core/types.h>
#include <core/array_types.h>
#include <core/string_types.h>
#include <mol/molecule_dynamic.h>

typedef bool (*FilterCommandFunc)(Array<bool> mask, const MoleculeDynamic& dynamic, Array<const CString> args);

struct FilterCommand {
    StringBuffer<16> keyword{};
    FilterCommandFunc func = nullptr;
};

namespace filter {
void initialize();
void shutdown();
bool compute_filter_mask(Array<bool> mask, const MoleculeDynamic& dynamic, CString filter);
void filter_colors(Array<uint32> colors, Array<bool> mask);

template <typename T>
void extract_filtered_data(DynamicArray<T>* dst_data, Array<const T> src_data, Array<const bool> mask) {
    ASSERT(dst_data);
    dst_data->clear();
    for (int32 i = 0; i < mask.count; i++) {
        if (mask[i]) dst_data->push_back(src_data[i]);
    }
}

template <typename T>
DynamicArray<T> extract_filtered_data(Array<const T> data, Array<const bool> mask) {
    DynamicArray<T> result{};
    extract_filtered_data(&result, data, mask);
    return result;
}

}  // namespace filter