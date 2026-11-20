#pragma once

#include <stddef.h>
#include <stdint.h>
#include <vector>

#include <core/md_str.h>
#include <core/md_unit.h>
#include <core/md_vec_math.h>
#include <system/md_attributes.h>

struct md_system_t;

// ATTRIBUTE ACCESS
//
// The few ways viamd reads the system's attribute table (md_attributes.h), the same for every component. mdlib's
// own calls do everything; these name what a reader asks for, by path, with the shape it expects checked:
//
//   attribute_exists     whether the table has the path at all ("does this file carry X")
//   attribute_shaped     an attribute of the expected shape (rank, components, at least min_rows), else NULL
//   attribute_count      the rows of one (shape[0]), 0 when absent
//   attribute_values     its values as doubles IN A UNIT, copied: converted from any stored type and unit, a virtual
//                        attribute computed. What a reader that depends on the unit takes (a length in Angstrom).
//   attribute_view<T>    its stored values IN PLACE, when resident and stored as T: what a panel takes every frame
//                        from what a reader published in the unit it knows, without a copy. Valid until the table
//                        changes (md_attributes_get's rule): never kept past the frame.
//   attribute_scalar     one value in a unit, or a fallback
//   attribute_str        one string, or empty
//
// md_unit_none() as a unit is "as stored" (md_attribute_extract_f64): what a caller passes when the unit is the
// stored one by construction, and says so by passing it.

bool attribute_exists(const md_system_t& sys, md_str_t path);

// Whether attr has the shape: rank, components (0: any) and at least min_rows along its first axis; a number
bool attribute_has_shape(const md_attribute_t* attr, uint32_t rank, uint32_t components, size_t min_rows = 0);
const md_attribute_t* attribute_shaped(const md_system_t& sys, md_str_t path, uint32_t rank, uint32_t components, size_t min_rows = 0);

size_t attribute_count(const md_system_t& sys, md_str_t path);

// All of attr (any shape), as doubles in unit; empty when it cannot be (absent, STR, dimensions that do not match)
std::vector<double> attribute_values(const md_system_t& sys, const md_attribute_t* attr, md_unit_t unit);
// A series (rank 1) of `components` values each, of exactly `rows` of them (0: any), in unit; empty otherwise
std::vector<double> attribute_values(const md_system_t& sys, md_str_t path, md_unit_t unit, uint32_t components = 1, size_t rows = 0);

double   attribute_scalar(const md_system_t& sys, md_str_t path, md_unit_t unit, double fallback = 0.0);
md_str_t attribute_str(const md_system_t& sys, md_str_t path, size_t index = 0);

// The stored values in place (see above). rows = shape[0], cols = shape[1] of a rank 2 attribute (1 of rank 1).
template <typename T>
struct AttributeView {
    const T* data = nullptr;
    size_t   rows = 0;
    size_t   cols = 0;
    explicit operator bool() const { return data != nullptr; }
    // The data, its rows (and columns) written out: for a caller keeping counts of its own
    const T* get(size_t* out_rows, size_t* out_cols = nullptr) const {
        if (out_rows) *out_rows = rows;
        if (out_cols) *out_cols = cols;
        return data;
    }
    // Row r of a rank 2 one; NULL past the rows, or for a row shorter than min_cols (a caller walking it over
    // something of its own, the QM atoms)
    const T* row(size_t r, size_t min_cols = 0) const { return data && r < rows && cols >= min_cols ? data + r * cols : nullptr; }
};

// What T is stored as: its type and components
struct AttributeElement { md_attribute_type_t type; uint32_t components; };
constexpr AttributeElement attribute_element(const double*)       { return {MD_ATTRIBUTE_TYPE_F64, 1}; }
constexpr AttributeElement attribute_element(const float*)        { return {MD_ATTRIBUTE_TYPE_F32, 1}; }
constexpr AttributeElement attribute_element(const uint8_t*)      { return {MD_ATTRIBUTE_TYPE_U8,  1}; }
constexpr AttributeElement attribute_element(const int32_t*)      { return {MD_ATTRIBUTE_TYPE_I32, 1}; }
constexpr AttributeElement attribute_element(const uint32_t*)     { return {MD_ATTRIBUTE_TYPE_U32, 1}; }
constexpr AttributeElement attribute_element(const md_double3_t*) { return {MD_ATTRIBUTE_TYPE_F64, 3}; }

const void* attribute_view_raw(const md_system_t& sys, md_str_t path, AttributeElement e, uint32_t rank, size_t min_rows, size_t* rows, size_t* cols);

// At least min_rows of them (fewer: none), of the given rank
template <typename T>
AttributeView<T> attribute_view(const md_system_t& sys, md_str_t path, uint32_t rank = 1, size_t min_rows = 0) {
    AttributeView<T> v;
    v.data = (const T*)attribute_view_raw(sys, path, attribute_element((const T*)nullptr), rank, min_rows, &v.rows, &v.cols);
    return v;
}
