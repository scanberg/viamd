#include <attribute_access.h>

#include <md_system.h>

bool attribute_exists(const md_system_t& sys, md_str_t path) {
    return md_attributes_find(&sys.attributes, path) != nullptr;
}

bool attribute_has_shape(const md_attribute_t* a, uint32_t rank, uint32_t components, size_t min_rows) {
    if (!a || a->format.rank != rank || !md_attribute_type_is_numeric(a->format.type)) return false;
    if (components && a->format.components != components) return false;
    return rank == 0 || a->format.shape[0] >= min_rows;
}

const md_attribute_t* attribute_shaped(const md_system_t& sys, md_str_t path, uint32_t rank, uint32_t components, size_t min_rows) {
    const md_attribute_t* a = md_attributes_find(&sys.attributes, path);
    return attribute_has_shape(a, rank, components, min_rows) ? a : nullptr;
}

size_t attribute_count(const md_system_t& sys, md_str_t path) {
    const md_attribute_t* a = md_attributes_find(&sys.attributes, path);
    return a && a->format.rank >= 1 ? a->format.shape[0] : 0;
}

std::vector<double> attribute_values(const md_system_t& sys, const md_attribute_t* attr, md_unit_t unit) {
    std::vector<double> v;
    if (!attr || !md_attribute_type_is_numeric(attr->format.type)) return v;
    const size_t n = md_attribute_element_count(&attr->format);
    v.resize(n);
    if (n == 0 || md_attribute_extract_f64(v.data(), n, &sys.attributes, attr, md_attribute_slice_all(), unit) != n) v.clear();
    return v;
}

std::vector<double> attribute_values(const md_system_t& sys, md_str_t path, md_unit_t unit, uint32_t components, size_t rows) {
    const md_attribute_t* a = attribute_shaped(sys, path, 1, components);
    if (!a || (rows && a->format.shape[0] != rows)) return {};
    return attribute_values(sys, a, unit);
}

double attribute_scalar(const md_system_t& sys, md_str_t path, md_unit_t unit, double fallback) {
    // anything not single valued is the fallback: the extract refuses a cap smaller than the attribute
    const md_attribute_t* a = md_attributes_find(&sys.attributes, path);
    double value = 0.0;
    if (!a || md_attribute_extract_f64(&value, 1, &sys.attributes, a, md_attribute_slice_all(), unit) != 1) return fallback;
    return value;
}

md_str_t attribute_str(const md_system_t& sys, md_str_t path, size_t index) {
    const md_attribute_t* a = md_attributes_find(&sys.attributes, path);
    if (!a || a->format.type != MD_ATTRIBUTE_TYPE_STR || index >= md_attribute_value_count(&a->format)) return {};
    return md_attribute_str(&sys.attributes, a, index);
}

const void* attribute_view_raw(const md_system_t& sys, md_str_t path, AttributeElement e, uint32_t rank, size_t min_rows, size_t* rows, size_t* cols) {
    *rows = 0;
    *cols = 0;
    const md_attribute_t* a = md_attributes_find(&sys.attributes, path);
    const void* data = md_attribute_view(a, e.type, e.components, rank);
    if (!data || (rank >= 1 && a->format.shape[0] < min_rows)) return nullptr;
    *rows = rank >= 1 ? a->format.shape[0] : 1;
    *cols = rank >= 2 ? a->format.shape[1] : 1;
    return data;
}
