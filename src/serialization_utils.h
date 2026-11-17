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
// A key can be a PATH into what its section is of: steps separated by '.', each a name, or the index of an item of
// an array - Field.Range, Slice.0.Opacity. A section header can be one too: [Representation.ElectronicStructure]
// right after a [Representation] continues it, and its keys are paths within ElectronicStructure. To a reader the
// two are the same path; the writer puts a large block in a section of its own so that its lines are short
// (Visitor::section). The keys of a workspace from before the paths are single steps, and those that moved into a
// block are mapped onto their paths (ElectronicStructureSlice3 onto ElectronicStructure.Slice.3.Plane, PathAlias).
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

// Takes the first step off a path: false when there is none
bool path_next(md_str_t& step, md_str_t& path);

// a.b into buf (b alone when a is empty), empty when it does not fit
md_str_t path_join(char* buf, size_t cap, md_str_t a, md_str_t b);

// Advances to the next section when it continues 'section' ([section.sub]), sub its path within it. Otherwise false,
// and the next section is left for next_section_header: a continuation that does not follow what it continues is
// read by nobody, as any section nobody knows.
bool next_subsection(md_str_t& sub, md_str_t section, deserialization_state_t& state);

// The entries of the current section, and of the sections continuing it, each as its path and value: read(path, arg)
// says whether it was one it knows, and one that was not is logged (debug).
template <typename F>
void read_section(deserialization_state_t& state, md_str_t section, F&& read) {
	md_str_t ident, arg;
	while (next_entry(ident, arg, state)) {
		if (!read(ident, arg)) {
			MD_LOG_DEBUG("Workspace: [" MD_STR_FMT "] has no '" MD_STR_FMT "', it is skipped", MD_STR_ARG(section), MD_STR_ARG(ident));
		}
	}
	md_str_t sub;
	while (next_subsection(sub, section, state)) {
		while (next_entry(ident, arg, state)) {
			char buf[256];
			const md_str_t path = path_join(buf, sizeof(buf), sub, ident);
			if (md_str_empty(path) || !read(path, arg)) {
				MD_LOG_DEBUG("Workspace: [" MD_STR_FMT "] has no '" MD_STR_FMT "', it is skipped", MD_STR_ARG(state.cur_section), MD_STR_ARG(ident));
			}
		}
	}
}

// Where a key of a workspace from before the paths is now. '#' in both stands for the number of a numbered key:
// {"ElectronicStructureSliceOpacity#", "ElectronicStructure.Slice.#.Opacity"}. A '#' is digits and nothing else, so a
// family is never mistaken for a longer one which begins like it (ElectronicStructureSlice# for SliceOpacity3).
struct PathAlias {
	const char* key;
	const char* path;
};

// The path of key, into buf: empty when it is none of the aliases
md_str_t path_alias(char* buf, size_t cap, md_str_t key, const PathAlias* aliases, size_t count);

template <size_t N>
static inline md_str_t path_alias(char* buf, size_t cap, md_str_t key, const PathAlias (&aliases)[N]) {
	return path_alias(buf, cap, key, aliases, N);
}

// ## Visitor
//
// One function per struct lists its fields, and it is called both to write the struct and to read each entry back
// into it: what is written and what is read are the same list, and a struct used in several places (a colour scale)
// has its list once.
//
//     void visit(viamd::Visitor& v, Thing& t) {
//         if (v.field("Source", t.source, source_names)) {
//             thing_set_source_defaults(&t);           // only when it was read
//         }
//         v.field("Iso", t.iso);
//         if (v.when(t.has_field) && v.block("Field")) {
//             visit(v, t.field_scale);                 // Field.Colormap, Field.Range, ...
//             v.end();
//         }
//         for (int i; v.each("Slice", i, t.num_slices, MAX_SLICES); ) {
//             v.field("Opacity", t.slices[i].opacity); // Slice.0.Opacity, ...
//         }
//     }
//
// Writing goes through the list in its order. Reading takes one entry at a time, so the fields are read in the order
// of the file, which is the order they were written in: a field which sets others up (a source, and the defaults of
// its spin) is read before them, as long as it comes before them in the list.
//
// A field returns true only when reading, only for the field the entry was for, and only when its value read: what
// follows reading a value - a clamp, what it sets up - goes in that branch. Writing, it returns false. when(),
// block(), section() and each() say what is visited, writing and reading alike. Writing does not change the struct,
// so a writer may hand in one it was given as const.

// The names of the values of an enum, as a workspace keeps them: the value's identifier, which does not change when
// the enum is reordered or what it is shown as is
struct EnumNames {
	const char* const* names;
	int count;
};

template <size_t N>
static inline EnumNames enum_names(const char* const (&names)[N]) { return { names, (int)N }; }

struct Visitor {
	// Writes the fields of the struct (an object) into the current section, which is 'section': a block written in a
	// section of its own is [section.Name]
	static Visitor writer(serialization_state_t& state, md_str_t section);
	// Reads one entry, by its path within the object
	static Visitor reader(md_str_t path, md_str_t arg);

	bool writing() const { return ser != nullptr; }
	bool reading() const { return ser == nullptr; }
	// What is written only when cond holds is read whenever the file has it
	bool when(bool cond) const { return reading() || cond; }
	// Reading: whether the entry was one of the fields (also when its value did not read)
	bool consumed() const { return taken; }

	bool field(const char* key, bool& val);
	bool field(const char* key, int& val);
	bool field(const char* key, float& val);
	bool field(const char* key, double& val);
	bool field(const char* key, md_float2_t& val);
	bool field(const char* key, md_float3_t& val);
	bool field(const char* key, md_float4_t& val);
	// Exactly len numbers
	bool field(const char* key, float* elem, size_t len);
	// Read, a view of the text being read
	bool field(const char* key, md_str_t& val);
	bool field(const char* key, char* buf, size_t cap);

	// An enum, by its name. Its index is read too, which is how a value was kept before it had a name.
	template <typename E>
	bool field(const char* key, E& val, const EnumNames& names) {
		int i = (int)val;
		if (!field_enum(key, i, names)) return false;
		val = (E)i;
		return true;
	}
	// An enum by its index alone, in [0, count): for the keys an older VIAMD reads too, which knows no names
	template <typename E>
	bool field_index(const char* key, E& val, int count) {
		int i = (int)val;
		if (!field_index_int(key, i, count)) return false;
		val = (E)i;
		return true;
	}
	// A value that has a name elsewhere (a colour map, by ImPlot's name): name_of gives the name of a value (NULL for
	// none, which is written as the number), index_of the value of a name (-1 for none). A number is read as the
	// value it is when it has a name.
	bool field_named(const char* key, int& val, const char* (*name_of)(int), int (*index_of)(md_str_t));

	// A list of count items of 'arity' numbers each, on one line: Isolines=-0.01,0.01. Reading one replaces the list;
	// an item of it, key.<i>, is read too (the numbered keys of older workspaces), the count then at least i + 1.
	bool list(const char* key, float* elem, int& count, int cap, int arity = 1);

	// A block of fields within the object (key.field): true when its fields are to be visited, then end()
	bool block(const char* key);
	// A block in a section of its own ([section.key]), for one with many fields. It comes after every other field of
	// the object, which would otherwise be written into it. Then end().
	bool section(const char* key);
	// Each item of an array of blocks (key.<i>.field), i its index: for (int i; v.each("Slice", i, count, cap); ) {}.
	// The count is written first (key=<count>), so an array is what it says even where it starts with items of its
	// own; reading an item makes the count at least its index + 1.
	bool each(const char* key, int& i, int& count, int cap);
	void end();

private:
	enum FrameKind : uint8_t { Frame_Block, Frame_Section, Frame_Each };
	struct Frame {
		const char* key;
		int         index;
		uint32_t    prefix_len;
		FrameKind   kind;
	};

	bool field_enum(const char* key, int& val, const EnumNames& names);
	bool field_index_int(const char* key, int& val, int count);
	// Reading: whether the entry is key, within the blocks entered (and takes it)
	bool take(const char* key);
	// Reading: whether the entry is within key, the rest after it
	bool enter(const char* key, md_str_t* rest) const;
	// Writing: the path of key, in the blocks entered
	md_str_t key_path(const char* key);
	void push(FrameKind kind, const char* key, int index);
	bool done() const { return taken || finished; }

	// Writing
	serialization_state_t* ser = nullptr;
	md_str_t section_name = {};
	char     prefix[192] = {};
	uint32_t prefix_len = 0;
	bool     sectioned = false;   // a section was written: the object's own fields are behind it
	char     key_buf[256] = {};

	// Reading
	md_str_t rest = {};           // what is left of the entry's path, past the blocks entered
	md_str_t arg = {};
	bool     taken = false;       // the entry was a field
	bool     finished = false;    // the entry is within a block that was visited: nothing after it can be it

	Frame frames[8] = {};
	int   depth = 0;
};

}
