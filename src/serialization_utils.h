#pragma once

#include <stdint.h>
#include <core/md_str.h>
#include <core/md_str_builder.h>
#include <core/md_vec_math.h>

struct md_bitfield_t;

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
// Readers are forgiving and writers are exact: numbers are written with the digits it takes to read
// back the same value, an extract_* that cannot parse its argument returns false and leaves the
// destination untouched, and an entry nobody reads is skipped.

namespace viamd {

struct deserialization_state_t {
	str_t filename;
	str_t text;
	str_t cur_section;
};

struct serialization_state_t {
	str_t filename;
	md_strb_t sb;
};

// Advances to the next section header, skipping whatever the current section's reader left unread.
// Skipping goes entry by entry, so a multiline value holding a line that starts with '[' is not taken
// for a header.
bool next_section_header(str_t& section, deserialization_state_t& state);

// Use these two when checking the current section and parsing the entries within it
inline str_t section_header(deserialization_state_t& state) { return state.cur_section; }
bool  next_entry(str_t& ident, str_t& arg, deserialization_state_t& state);

void write_section_header(serialization_state_t& state, str_t section);
void write_int(serialization_state_t& state, str_t ident, int64_t val);
void write_int_vec(serialization_state_t& state, str_t ident, const int* elem, size_t len);
void write_flt(serialization_state_t& state, str_t ident, float val);
void write_dbl(serialization_state_t& state, str_t ident, double val);
void write_flt_vec(serialization_state_t& state, str_t ident, const float* elem, size_t len);
void write_str(serialization_state_t& state, str_t ident, str_t str);
// Base64 between ### markers: any bit pattern survives the line based format
void write_bitfield(serialization_state_t& state, str_t ident, const md_bitfield_t* bf);

static inline void write_vec3(serialization_state_t& state, str_t ident, vec3_t v)  { write_flt_vec(state, ident, v.elem, 3); }
static inline void write_vec4(serialization_state_t& state, str_t ident, vec4_t v)  { write_flt_vec(state, ident, v.elem, 4); }
static inline void write_quat(serialization_state_t& state, str_t ident, quat_t q)  { write_flt_vec(state, ident, q.elem, 4); }
static inline void write_bool(serialization_state_t& state, str_t ident, bool val) { write_int(state, ident, (int)val); }

// 0/1, and true/false in any case
bool extract_bool(bool& val,  str_t arg);
bool extract_int (int&  val,  str_t arg);
bool extract_int64(int64_t& val, str_t arg);
bool extract_int_vec(int* elem, size_t len, str_t arg);
bool extract_dbl (double& val, str_t arg);
bool extract_flt (float& val, str_t arg);
bool extract_flt_vec (float* elem, size_t len, str_t arg);
bool extract_str (str_t& str, str_t arg);
bool extract_to_char_buf(char* buf, size_t cap, str_t arg);
bool extract_bitfield(md_bitfield_t* bf, str_t arg);

// An integer into an enum, accepted only within [0, count)
template <typename T>
static inline bool extract_enum(T& val, str_t arg, int count) {
	int i;
	if (extract_int(i, arg) && 0 <= i && i < count) {
		val = (T)i;
		return true;
	}
	return false;
}

static inline bool extract_vec3(vec3_t& v, str_t arg) { return extract_flt_vec(v.elem, 3, arg); }
static inline bool extract_vec4(vec4_t& v, str_t arg) { return extract_flt_vec(v.elem, 4, arg); }
static inline bool extract_quat(quat_t& q, str_t arg) { return extract_flt_vec(q.elem, 4, arg); }

}
