#pragma once

#include <stddef.h>
#include <stdint.h>

#include <core/md_str.h>
#include <core/md_array.h>
#include <md_attributes.h>

#include <surface_field.h>

struct md_system_t;
struct md_allocator_i;
enum class ElectronicStructureSource;

// FIELDS
//
// The fields over R^3 a Field representation shows (and colours a surface by, and takes a surface from), named by a
// path as attributes are: its segments a hierarchy, its first the source the field comes from, the rest what it is.
//   es/orbital                         (spin, index)
//   es/nto                             (state, component, index)
//   es/transition_density              (state, component)
//   es/density                         (spin)
//   es/density_property/<attribute>    a leaf of the density properties of the file (vlx/density_property/...)
//   es/elf, es/lol
//   es/esp, es/embedding_potential, es/alie, es/leae, es/density_scale_length
//   es/gostshyp/wall_potential, es/gostshyp/effective_potential, es/gostshyp/energy_density
//   cube/<name>                        a grid of a file, the path of its attribute (grid_field.h)   (member)
// A path names a field, or a family of them told apart by parameters (an orbital by its spin and index), as an array
// attribute is one path with an index. The path is what a workspace keeps of what a representation shows; the label
// what the picker shows; a search matches either, and some words more (FieldEntry::keywords).
//
// Each binds to what evaluates it: of the electronic structure, the representation's source (ElectronicStructureSource,
// and a density property's attribute) and through it the kind of surface field (surface_field.h) where it is one.

#define FIELD_PATH_MAX 160

struct FieldEntry {
    char path[FIELD_PATH_MAX] = "";
    md_str_t label = {};                    // the pretty name, searched as the path is; borrowed: static, or the attribute's
    const char* keywords = "";              // more words a search finds it by
    ElectronicStructureSource source = {};
    md_attribute_id_t attribute_key = MD_ATTRIBUTE_INVALID;
    SurfaceFieldKind field_kind = SurfaceFieldKind::Count;  // the surface field it is, if it is one: what can colour
};

// The fields the system has, in the order of the hierarchy, into *out (pushed). Their labels are borrowed from the
// system's attributes: valid while it is.
size_t field_catalog(md_array(FieldEntry)* out, const md_system_t& sys, md_allocator_i* alloc);

// The label of a group of paths, the segments before a leaf ("es" Electronic Structure, "es/gostshyp" GOSTSHYP); the
// segment itself for one without
md_str_t field_group_label(md_str_t group);

// The path of the field a representation's source is (and its density property), and of a surface field
size_t field_path_of_source(char* buf, size_t cap, ElectronicStructureSource source, md_attribute_id_t attribute_key, const md_system_t& sys);
size_t field_path_of_kind(char* buf, size_t cap, SurfaceFieldKind kind);

// What a path is: false for one no field has (a density property's is resolved without the system: its key is its
// attribute's path)
bool field_resolve(ElectronicStructureSource* source, md_attribute_id_t* attribute_key, md_str_t path);
bool field_resolve_kind(SurfaceFieldKind* kind, md_str_t path);

// A combo of the fields: its popup a filter (every word of it in the path, the label or the keywords) and the fields
// under the labels of their groups. `pickable` leaves out what cannot be picked here (NULL: none); `selected` the
// path shown as picked. Returns the index of the one picked, -1 for none.
int field_picker(const char* id, const char* preview, const FieldEntry* entries, size_t count, md_str_t selected, bool (*pickable)(const FieldEntry&));
