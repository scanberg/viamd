#pragma once

#include <stdint.h>
#include <core/md_str.h>
#include <core/md_writer.h>
#include <core/md_vec.h>
#include <core/md_log.h>

struct md_bitset_t;

// Reading and writing the workspace text format:
//
//     [Section]
//     Ident=value
//     Text="""a value spanning
//     several lines"""
//
// A section is a header line in brackets, followed by entries until the next header. An entry is
// one line split at its first '='; both sides are trimmed. A value that would not survive that - one
// spanning lines, or with whitespace at either end - is written between triple quotes, and read back
// verbatim. Lines that are neither (blank lines, the banner, anything starting with '#' or ';') are
// skipped.
//
// A key can be a PATH into what its section is of (see Paths below): Field.Range, Slice.0.Opacity.
//
// Readers are forgiving and writers are exact: numbers are written with the digits it takes to read
// back the same value, an extract_* that cannot parse its argument returns false and leaves the
// destination untouched, and an entry nobody reads is skipped.
//
// The decimal point is '.', in every locale. A vector is read only when it has exactly the components
// asked for: one written with decimal commas has more, and is logged and left unread rather than
// read as other values.

namespace viamd {

struct deserialization_state_t {
	md_str_t filename;
	md_str_t text;
	md_str_t cur_section;
	void*    reader = nullptr;	// What reads it, for what its sections name beyond themselves (load_workspace: the systems)
};

struct serialization_state_t {
	md_str_t filename;
	md_str_writer_t sw;
};

// Advances to the next section header, skipping whatever the current section's reader left unread.
// Skipping goes entry by entry, so a multiline value holding a line that starts with '[' is not taken
// for a header.
bool next_section_header(md_str_t& section, deserialization_state_t& state);

// Use these two when checking the current section and parsing the entries within it
inline md_str_t section_header(deserialization_state_t& state) { return state.cur_section; }
bool  next_entry(md_str_t& ident, md_str_t& arg, deserialization_state_t& state);

void write_section_header(serialization_state_t& state, md_str_t section);
void write_int(serialization_state_t& state, md_str_t ident, int64_t val);
void write_int_vec(serialization_state_t& state, md_str_t ident, const int* elem, size_t len);
void write_flt(serialization_state_t& state, md_str_t ident, float val);
void write_dbl(serialization_state_t& state, md_str_t ident, double val);
void write_flt_vec(serialization_state_t& state, md_str_t ident, const float* elem, size_t len);
void write_dbl_vec(serialization_state_t& state, md_str_t ident, const double* elem, size_t len);
void write_str(serialization_state_t& state, md_str_t ident, md_str_t str);
// Base64 between ### markers: any bit pattern survives the line based format
void write_bitset(serialization_state_t& state, md_str_t ident, const md_bitset_t* bs);

static inline void write_vec3(serialization_state_t& state, md_str_t ident, md_float3_t v)  { write_flt_vec(state, ident, v.elem, 3); }
static inline void write_vec4(serialization_state_t& state, md_str_t ident, md_float4_t v)  { write_flt_vec(state, ident, v.elem, 4); }
static inline void write_quat(serialization_state_t& state, md_str_t ident, md_quat_t q)  { write_flt_vec(state, ident, q.elem, 4); }
static inline void write_bool(serialization_state_t& state, md_str_t ident, bool val) { write_int(state, ident, (int)val); }

// 0/1, and true/false in any case
bool extract_bool(bool& val,  md_str_t arg);
bool extract_int (int&  val,  md_str_t arg);
bool extract_int64(int64_t& val, md_str_t arg);
bool extract_int_vec(int* elem, size_t len, md_str_t arg);
bool extract_dbl (double& val, md_str_t arg);
bool extract_flt (float& val, md_str_t arg);
bool extract_flt_vec (float* elem, size_t len, md_str_t arg);
bool extract_dbl_vec (double* elem, size_t len, md_str_t arg);
bool extract_str (md_str_t& str, md_str_t arg);
bool extract_to_char_buf(char* buf, size_t cap, md_str_t arg);
// The set is allocated with alloc. Reads the byte format of the former md_bitfield masks, and the one for sets
// reaching beyond 2^24 (see md_bitset_serialize)
bool extract_bitset(md_bitset_t* bs, md_str_t arg, struct md_allocator_i* alloc);

// An integer into an enum, accepted only within [0, count)
template <typename T>
static inline bool extract_enum(T& val, md_str_t arg, int count) {
	int i;
	if (extract_int(i, arg) && 0 <= i && i < count) {
		val = (T)i;
		return true;
	}
	return false;
}

static inline bool extract_vec3(md_float3_t& v, md_str_t arg) { return extract_flt_vec(v.elem, 3, arg); }
static inline bool extract_vec4(md_float4_t& v, md_str_t arg) { return extract_flt_vec(v.elem, 4, arg); }
static inline bool extract_quat(md_quat_t& q, md_str_t arg) { return extract_flt_vec(q.elem, 4, arg); }

// A list of any length up to cap: count is how many there were, and an empty value is the empty list. False, with
// nothing written, when one of them is not a number or there are more than cap. Lists are of the paths, which no
// VIAMD wrote with decimal commas, so a ',' always separates.
bool extract_flt_list(float* elem, size_t* count, size_t cap, md_str_t arg);

// ## Paths
//
// A key can be a PATH into what its section is of: steps separated by '.', each a name or the index of an item of an
// array - Field.Range, Slice.0.Opacity. A struct is written by a function of its own, its keys under a prefix
// (PathPrefix), and read by another from the entries of its section (Entries), by path. The keys of a workspace from
// before the paths are single steps: those that only moved are put at their paths as they are read (PathAlias), and
// those whose meaning changed are read by their old names, by whoever reads the struct.

// Where a key of a workspace from before the paths is now. '#' in both stands for the number of a numbered key:
// {"ElectronicStructureSliceOpacity#", "ElectronicStructure.Slice.#.Opacity"}. A '#' is digits and nothing else, so a
// family is never mistaken for a longer one which begins like it (ElectronicStructureSlice# for SliceOpacity3).
struct PathAlias {
	const char* key;
	const char* path;
};

// The path of key, into buf: empty when it is none of the aliases
md_str_t path_alias(char* buf, size_t cap, md_str_t key, const PathAlias* aliases, size_t count);

// The keys of an object: p("Iso") is "<prefix>Iso"; p.sub("Slice", 3) the prefix of the item Slice.3, p.sub("Field")
// of the block Field. A key is valid until the next one of the same prefix.
struct PathPrefix {
	char     buf[192] = {};
	uint32_t len = 0;
	char     key_buf[256] = {};

	PathPrefix() = default;
	explicit PathPrefix(const char* prefix);
	PathPrefix sub(const char* key, int index = -1) const;
	md_str_t operator()(const char* key);
};

// The entries of a section, read whole and in their order: each a path and a value, the keys that only moved put at
// their paths. Read by path (get_*), each one marked as read, so that what nobody read can be told (entries_unread).
struct Entry {
	md_str_t path;
	md_str_t value;
	bool     read;
};

struct Entries {
	Entry* items = nullptr;   // md_array
	size_t count = 0;
	struct md_allocator_i* alloc = nullptr;
};

// The entries up to the next section header, in alloc (the paths an alias made; the values are views of the text)
Entries read_entries(deserialization_state_t& state, const PathAlias* aliases, size_t num_aliases, struct md_allocator_i* alloc);
template <size_t N>
static inline Entries read_entries(deserialization_state_t& state, const PathAlias (&aliases)[N], struct md_allocator_i* alloc) {
	return read_entries(state, aliases, N, alloc);
}
static inline Entries read_entries(deserialization_state_t& state, struct md_allocator_i* alloc) {
	return read_entries(state, nullptr, 0, alloc);
}

// Logs (debug) the entries nobody read, of the section
void entries_unread(const Entries& e, md_str_t section);

// The value at path, marked read: the last there is (a key and an older one put at the same path: the later stands).
// NULL for none.
const md_str_t* get_value(Entries& e, md_str_t path);

// The value at path into val: true when there and read. A value that does not read is left out (and logged, debug).
bool get(Entries& e, md_str_t path, bool& val);
bool get(Entries& e, md_str_t path, int& val);
bool get(Entries& e, md_str_t path, float& val);
bool get(Entries& e, md_str_t path, double& val);
bool get(Entries& e, md_str_t path, md_float2_t& val);
bool get(Entries& e, md_str_t path, md_float3_t& val);
bool get(Entries& e, md_str_t path, md_float4_t& val);
bool get(Entries& e, md_str_t path, float* elem, size_t len);   // exactly len numbers
bool get(Entries& e, md_str_t path, char* buf, size_t cap);
// A view of the text being read: copy what is kept
bool get(Entries& e, md_str_t path, md_str_t& val);

// The names of the values of an enum, as a workspace keeps them: the value's identifier, which does not change when
// what it is shown as does. The enum may grow at its end; its values are not reordered, as an older workspace keeps
// the index.
struct EnumNames {
	const char* const* names;
	int count;
};

template <size_t N>
static inline EnumNames enum_names(const char* const (&names)[N]) { return { names, (int)N }; }

// An enum by its name, or by its index as a workspace from before the names kept it
bool get_enum_int(Entries& e, md_str_t path, int& val, const EnumNames& names);
template <typename E>
static inline bool get_enum(Entries& e, md_str_t path, E& val, const EnumNames& names) {
	int i = (int)val;
	if (!get_enum_int(e, path, i, names)) return false;
	val = (E)i;
	return true;
}
void write_enum(serialization_state_t& state, md_str_t path, int val, const EnumNames& names);

// A value that has a name elsewhere (a colour map, by ImPlot's name): name_of gives the name of a value (NULL for
// none, which is written as the number), index_of the value of a name (-1 for none)
bool get_named(Entries& e, md_str_t path, int& val, int (*index_of)(md_str_t));
void write_named(serialization_state_t& state, md_str_t path, int val, const char* (*name_of)(int));

// A list of count items of 'arity' numbers each, on one line: Isolines=-0.01,0.01. Or, as an older workspace has it,
// an item a key: path.<i>. True when either was there.
bool get_list(Entries& e, md_str_t path, float* elem, int& count, int cap, int arity = 1);
void write_list(serialization_state_t& state, md_str_t path, const float* elem, int count, int arity = 1);

// How many items an array has: its count (path=<n>), at least one past the last item there is (path.<i>....), at
// most cap; unchanged when there is neither. Its items are then read as path.<i>.
bool get_count(Entries& e, md_str_t path, int& count, int cap);

}
