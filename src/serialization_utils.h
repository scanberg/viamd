#pragma once

#include <stdint.h>
#include <core/md_str.h>
#include <core/md_str_builder.h>
#include <core/md_vec_math.h>

struct md_bitfield_t;

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

bool next_section_header(str_t& section, deserialization_state_t& state);

// Use these two when checking the current section and parsing the entries within it
inline str_t section_header(deserialization_state_t& state) { return state.cur_section; }
bool  next_entry(str_t& ident, str_t& arg, deserialization_state_t& state);

void write_section_header(serialization_state_t& state, str_t section);
void write_int(serialization_state_t& state, str_t ident, int64_t val);
void write_dbl(serialization_state_t& state, str_t ident, double val);
void write_flt_vec(serialization_state_t& state, str_t ident, const float* elem, size_t len);
void write_str(serialization_state_t& state, str_t ident, str_t str);
void write_bitfield(serialization_state_t& state, str_t ident, const md_bitfield_t* bf);

static inline void write_vec3(serialization_state_t& state, str_t ident, vec3_t v)  { write_flt_vec(state, ident, v.elem, 3); }
static inline void write_vec4(serialization_state_t& state, str_t ident, vec4_t v)  { write_flt_vec(state, ident, v.elem, 4); }
static inline void write_quat(serialization_state_t& state, str_t ident, quat_t q)  { write_flt_vec(state, ident, q.elem, 4); }
static inline void write_bool(serialization_state_t& state, str_t ident, bool val) { write_int(state, ident, (int)val); }
static inline void write_flt(serialization_state_t& state, str_t ident, float val) { write_dbl(state, ident, val); }

bool extract_bool(bool& val,  str_t arg);
bool extract_int (int&  val,  str_t arg);
bool extract_dbl (double& val, str_t arg);
bool extract_flt (float& val, str_t arg);
bool extract_flt_vec (float* elem, size_t len, str_t arg);
bool extract_str (str_t& str, str_t arg);
bool extract_bitfield(md_bitfield_t* bf, str_t arg);
}
