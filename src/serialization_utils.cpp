#include <serialization_utils.h>
#include <core/md_allocator.h>
#include <core/md_parse.h>
#include <core/md_bitset.h>
#include <core/md_base64.h>
#include <core/md_log.h>

#include <string.h>
#include <stdio.h>

namespace viamd {

static const md_str_t esc = MD_STR_INIT("\"\"\"");

// Numbers are written with a '.' whatever the locale. A ',' in one is from a workspace written by a VIAMD that
// followed a decimal comma locale, whose numbers cannot be told from the separators between them: the value is not
// read rather than read wrong, and the log says why.
static void log_decimal_comma(md_str_t arg) {
	MD_LOG_ERROR("Workspace: '" MD_STR_FMT "' is not read, its numbers were written with a decimal comma", MD_STR_ARG(arg));
}

// "0,5" or "1,5e-07": a single ',' where a '.' makes it a number
static bool is_decimal_comma_number(md_str_t arg) {
	char buf[64];
	if (arg.len == 0 || arg.len >= sizeof(buf)) return false;
	size_t commas = 0;
	for (size_t i = 0; i < arg.len; ++i) {
		buf[i] = arg.ptr[i];
		if (buf[i] == ',') {
			buf[i] = '.';
			commas += 1;
		}
	}
	return commas == 1 && md_parse_is_float(md_str_t{buf, arg.len});
}

static bool is_comment(md_str_t line) {
	return md_str_begins_with(line, MD_STR_LIT("#")) || md_str_begins_with(line, MD_STR_LIT(";"));
}

bool next_section_header(md_str_t& section, deserialization_state_t& state) {
	// Whatever the previous section's reader did not consume
	md_str_t ident, arg;
	while (next_entry(ident, arg, state)) {}

	md_str_t line;
	while (md_parse_next_line(&line, &state.text)) {
		line = md_str_trim(line);
		if (md_str_begins_with(line, MD_STR_LIT("[")) &&
			md_str_ends_with(line, MD_STR_LIT("]")))
		{
			section = md_str_trim(md_str_substr(line, 1, md_str_len(line) - 2));
			state.cur_section = section;
			return true;
		}
	}
	state.cur_section = {};
	return false;
}

bool next_entry(md_str_t& ident, md_str_t& arg, deserialization_state_t& state) {
	md_str_t line;
	while (md_parse_peek_line(&line, &state.text)) {
		line = md_str_trim(line);
		if (md_str_begins_with(line, MD_STR_LIT("["))) {
			return false;
		}
		md_parse_skip_line(&state.text);
		if (is_comment(line)) {
			continue;
		}

		if (md_str_cut(&ident, &arg, line, md_str_find_char(line, '='))) {
			ident = md_str_trim(ident);
			arg   = md_str_trim(arg);
			if (md_str_begins_with(arg, esc)) {
				// Multiline string: from the opening quotes in the text itself (the line was only a
				// view of its first line) to the matching closing ones
				const char* beg = md_str_beg(arg) + md_str_len(esc);
				const char* end = md_str_end(state.text);
				md_str_t haystack = {beg, (size_t)(end-beg)};
				const size_t loc = md_str_find_str(haystack, esc);
				if (loc != MD_STR_NPOS) {
					arg = md_str_substr(haystack, 0, loc);
					state.text = md_str_substr(haystack, loc + md_str_len(esc));
					// The rest of the closing line
					md_parse_skip_line(&state.text);
				} else {
					MD_LOG_ERROR("Workspace: unbalanced \"\"\" in the value of '" MD_STR_FMT "', the rest of the file is skipped", MD_STR_ARG(ident));
					state.text = {};
					return false;
				}
			}
			return true;
		}
	}
	return false;
}

void write_section_header(serialization_state_t& state, md_str_t section) {
	md_str_writer_fmt(&state.sw, "\n[" MD_STR_FMT "]\n", MD_STR_ARG(section));
}

void write_int(serialization_state_t& state, md_str_t ident, int64_t val) {
	md_str_writer_fmt(&state.sw, MD_STR_FMT "=%lld\n", MD_STR_ARG(ident), (long long)val);
}

void write_int_vec(serialization_state_t& state, md_str_t ident, const int* elem, size_t len) {
	md_str_writer_fmt(&state.sw, MD_STR_FMT "=", MD_STR_ARG(ident));
	for (size_t i = 0; i < len; ++i) {
		md_str_writer_fmt(&state.sw, "%i", elem[i]);
		if (i < len - 1) {
			md_str_writer_push_char(&state.sw, ',');
		}
	}
	md_str_writer_push_char(&state.sw, '\n');
}

// Shortest round trip is not on offer from printf; these digit counts are the ones that always
// read back the same value (FLT_DECIMAL_DIG and DBL_DECIMAL_DIG). %f would print 1e-7 as 0.000000.
void write_flt(serialization_state_t& state, md_str_t ident, float val) {
	md_str_writer_fmt(&state.sw, MD_STR_FMT "=%.9g\n", MD_STR_ARG(ident), (double)val);
}

void write_dbl(serialization_state_t& state, md_str_t ident, double val) {
	md_str_writer_fmt(&state.sw, MD_STR_FMT "=%.17g\n", MD_STR_ARG(ident), val);
}

void write_flt_vec(serialization_state_t& state, md_str_t ident, const float* elem, size_t len) {
	md_str_writer_fmt(&state.sw, MD_STR_FMT "=", MD_STR_ARG(ident));
	for (size_t i = 0; i < len; ++i) {
		md_str_writer_fmt(&state.sw, "%.9g", (double)elem[i]);
		if (i < len - 1) {
			md_str_writer_push_char(&state.sw, ',');
		}
	}
	md_str_writer_push_char(&state.sw, '\n');
}

void write_str(serialization_state_t& state, md_str_t ident, md_str_t str) {
	// Quoted whenever a plain value would not come back as it went: a line break ends it, and the
	// reader trims whitespace at either end
	const bool multiline  = md_str_contains_char(str, '\n') || md_str_contains_char(str, '\r');
	const bool padded     = str.len > 0 && (md_ascii_is_whitespace(str.ptr[0]) || md_ascii_is_whitespace(str.ptr[str.len - 1]));
	const bool looks_esc  = md_str_begins_with(str, esc);
	if (multiline || padded || looks_esc) {
		if (md_str_contains_str(str, esc)) {
			MD_LOG_ERROR("Workspace: the value of '" MD_STR_FMT "' contains \"\"\" and cannot be stored as it is", MD_STR_ARG(ident));
		}
		md_str_writer_fmt(&state.sw, MD_STR_FMT "=" MD_STR_FMT MD_STR_FMT MD_STR_FMT "\n", MD_STR_ARG(ident), MD_STR_ARG(esc), MD_STR_ARG(str), MD_STR_ARG(esc));
	} else {
		md_str_writer_fmt(&state.sw, MD_STR_FMT "=" MD_STR_FMT "\n", MD_STR_ARG(ident), MD_STR_ARG(str));
	}
}

void write_bitset(serialization_state_t& state, md_str_t ident, const md_bitset_t* bs) {
    md_temp_scope_t temp = md_temp_begin();
    MD_DEFER { md_temp_end(temp); };

	if (md_bitset_empty(*bs)) return;
	void*  serialized_data = md_temp_alloc(temp, md_bitset_serialize_size_in_bytes(*bs));
	size_t serialized_size = md_bitset_serialize(serialized_data, *bs);
	if (!serialized_size) {
		// Not silently: the workspace would load without it
		MD_LOG_ERROR("Failed to write the set '" MD_STR_FMT "' to the workspace", MD_STR_ARG(ident));
		return;
	}
	{
		char*  base64_data = (char*)md_temp_alloc(temp, md_base64_encode_size_in_bytes(serialized_size));
		size_t base64_size = md_base64_encode(base64_data, serialized_data, serialized_size);
		if (base64_size) {
			md_str_t base64 = {base64_data, base64_size};
			// The trailing newline is not cosmetic. Every other writer here terminates its line, and
			// next_entry() splits a line at its first '='. Without it the following entry is glued
			// onto the end of this one and silently swallowed - which is how every group label
			// after the first went missing from a saved workspace.
			md_str_writer_fmt(&state.sw, MD_STR_FMT "=###" MD_STR_FMT "###\n", MD_STR_ARG(ident), MD_STR_ARG(base64));
		}
	}
}

bool extract_bool(bool& val, md_str_t arg) {
	if (md_parse_is_int(arg)) {
		int64_t integer = md_parse_int(arg);
		if (integer == 0 || integer == 1) {
			val = (bool)integer;
			return true;
		}
	} else if (md_str_eq_cstr_ignore_case(arg, "true")) {
		val = true;
		return true;
	} else if (md_str_eq_cstr_ignore_case(arg, "false")) {
		val = false;
		return true;
	}
	return false;
}

bool extract_int(int& val, md_str_t arg) {
	if (md_parse_is_int(arg)) {
		const int64_t i = md_parse_int(arg);
		if (INT32_MIN <= i && i <= INT32_MAX) {
			val = (int)i;
			return true;
		}
	}
	return false;
}

bool extract_int64(int64_t& val, md_str_t arg) {
	if (md_parse_is_int(arg)) {
		val = md_parse_int(arg);
		return true;
	}
	return false;
}

bool extract_int_vec(int* elem, size_t len, md_str_t arg) {
	// Exactly len integers, or nothing is written
	int tmp[16];
	if (len > MD_COUNTOF(tmp)) return false;
	md_str_t tok;
	size_t count = 0;
	while (md_parse_next_token_delim(&tok, &arg, ',')) {
		tok = md_str_trim(tok);
		if (count == len || !md_parse_is_int(tok)) return false;
		tmp[count++] = (int)md_parse_int(tok);
	}
	if (count != len) return false;
	md_memcpy(elem, tmp, len * sizeof(int));
	return true;
}

bool extract_dbl (double& val, md_str_t arg) {
	if (md_parse_is_float(arg)) {
		val = md_parse_float(arg);
		return true;
	}
	if (is_decimal_comma_number(arg)) log_decimal_comma(arg);
	return false;
}

bool extract_flt (float& val, md_str_t arg) {
	if (md_parse_is_float(arg)) {
		val = (float)md_parse_float(arg);
		return true;
	}
	if (is_decimal_comma_number(arg)) log_decimal_comma(arg);
	return false;
}

bool extract_flt_vec (float* elem, size_t len, md_str_t arg) {
	// Exactly len numbers, parsed into a copy and written back only then. More than len is not cut down to the
	// first len: that is what a vector written with decimal commas looks like ("0,5,0,25,1,1" for 0.5 0.25 1 1),
	// and its first len would be values nobody wrote.
	float tmp[16];
	if (len > MD_COUNTOF(tmp)) return false;
	const md_str_t all = arg;
	md_str_t tok;
	size_t count = 0;
	size_t extra = 0;
	while (md_parse_next_token_delim(&tok, &arg, ',')) {
		tok = md_str_trim(tok);
		if (!(md_parse_is_float(tok) || md_parse_is_int(tok))) return false;
		if (count < len) {
			tmp[count++] = (float)md_parse_float(tok);
		} else {
			extra += 1;
		}
	}
	if (extra) log_decimal_comma(all);
	if (count != len || extra) return false;
	md_memcpy(elem, tmp, len * sizeof(float));
	return true;
}

bool extract_str(md_str_t& str, md_str_t arg) {
	// next_entry already strips the quotes of a multiline value; this covers a value handed in
	// with them still on
	if (arg.len >= 2 * esc.len && md_str_begins_with(arg, esc) && md_str_ends_with(arg, esc)) {
		str = md_str_substr(arg, esc.len, arg.len - 2 * esc.len);
	} else {
		str = arg;
	}
	return true;
}

bool extract_to_char_buf(char* buf, size_t cap, md_str_t arg) {
	md_str_t str;
	if (extract_str(str, arg)) {
		md_str_copy_to_char_buf(buf, cap, str);
		return true;
	}
	return false;
}

bool extract_bitset(md_bitset_t* bs, md_str_t arg, md_allocator_i* alloc) {
	// Not alloc's arena: the result must outlive this scope
	md_temp_scope_t temp = md_temp_begin_avoid(alloc);
	MD_DEFER { md_temp_end(temp); };

	// Bitfield starts with ###
	// and ends with ###
	md_str_t token = MD_STR_INIT("###");
	if (!md_str_begins_with(arg, token)) {
		MD_LOG_ERROR("Malformed start token for bitfield");
		return false;
	}
	arg = md_str_substr(arg, md_str_len(token));

	size_t loc;
	if ((loc = md_str_find_str(arg, token)) == MD_STR_NPOS) {
		MD_LOG_ERROR("Malformed end token for bitfield");
		return false;
	}

	arg = md_str_substr(arg, 0, loc);
	const size_t raw_cap = md_base64_decode_size_in_bytes(md_str_len(arg));
	void* raw_ptr = md_temp_alloc(temp, raw_cap);

	const size_t raw_len = md_base64_decode(raw_ptr, md_str_ptr(arg), md_str_len(arg));
	if (!raw_len || !md_bitset_deserialize(bs, raw_ptr, raw_len, alloc)) {
		MD_LOG_ERROR("Failed to deserialize bitset");
		*bs = md_bitset_empty_set();
		return false;
	}

	return true;
}

bool extract_flt_list(float* elem, size_t* count, size_t cap, md_str_t arg) {
	// Checked whole before anything is written
	md_str_t rest = arg;
	md_str_t tok;
	size_t n = 0;
	while (md_parse_next_token_delim(&tok, &rest, ',')) {
		tok = md_str_trim(tok);
		if (!(md_parse_is_float(tok) || md_parse_is_int(tok))) return false;
		n += 1;
	}
	if (n > cap) return false;
	rest = arg;
	for (size_t i = 0; i < n; ++i) {
		md_parse_next_token_delim(&tok, &rest, ',');
		elem[i] = (float)md_parse_float(md_str_trim(tok));
	}
	*count = n;
	return true;
}

// ## Paths

bool path_next(md_str_t& step, md_str_t& path) {
	if (md_str_empty(path)) return false;
	const size_t dot = md_str_find_char(path, '.');
	step = md_str_substr(path, 0, dot);
	path = dot == MD_STR_NPOS ? md_str_t{} : md_str_substr(path, dot + 1, SIZE_MAX);
	return true;
}

md_str_t path_join(char* buf, size_t cap, md_str_t a, md_str_t b) {
	if (md_str_empty(a)) return b;
	const size_t len = a.len + 1 + b.len;
	if (len >= cap) return {};
	memcpy(buf, a.ptr, a.len);
	buf[a.len] = '.';
	memcpy(buf + a.len + 1, b.ptr, b.len);
	buf[len] = '\0';
	return {buf, len};
}

bool next_subsection(md_str_t& sub, md_str_t section, deserialization_state_t& state) {
	// What the reader of the current section left unread
	md_str_t ident, arg;
	while (next_entry(ident, arg, state)) {}

	// next_entry stops at the line of the next header, without taking it
	md_str_t line;
	if (!md_parse_peek_line(&line, &state.text)) return false;
	line = md_str_trim(line);
	if (!md_str_begins_with(line, MD_STR_LIT("[")) || !md_str_ends_with(line, MD_STR_LIT("]"))) return false;

	const md_str_t name = md_str_trim(md_str_substr(line, 1, md_str_len(line) - 2));
	if (name.len <= section.len + 1 || name.ptr[section.len] != '.' || !md_str_begins_with(name, section)) return false;

	md_parse_skip_line(&state.text);
	state.cur_section = name;
	sub = md_str_substr(name, section.len + 1, SIZE_MAX);
	return true;
}

// The number a '#' stands for: 1 to 3 digits and nothing else
static bool is_key_number(md_str_t s) {
	if (s.len == 0 || s.len > 3) return false;
	for (size_t i = 0; i < s.len; ++i) {
		if (s.ptr[i] < '0' || s.ptr[i] > '9') return false;
	}
	return true;
}

md_str_t path_alias(char* buf, size_t cap, md_str_t key, const PathAlias* aliases, size_t count) {
	for (size_t i = 0; i < count; ++i) {
		const md_str_t old_key = md_str_from_cstr(aliases[i].key);
		const md_str_t path    = md_str_from_cstr(aliases[i].path);
		const size_t hash = md_str_find_char(old_key, '#');
		if (hash == MD_STR_NPOS) {
			if (md_str_eq(key, old_key)) {
				if (path.len >= cap) return {};
				md_str_copy_to_char_buf(buf, cap, path);
				return {buf, path.len};
			}
			continue;
		}
		// The families end in their number
		MD_ASSERT(hash == old_key.len - 1);
		const md_str_t family = md_str_substr(old_key, 0, hash);
		if (!md_str_begins_with(key, family)) continue;
		const md_str_t num = md_str_substr(key, family.len, SIZE_MAX);
		if (!is_key_number(num)) continue;

		const size_t at = md_str_find_char(path, '#');
		MD_ASSERT(at != MD_STR_NPOS);
		const size_t len = path.len - 1 + num.len;
		if (len >= cap) return {};
		memcpy(buf, path.ptr, at);
		memcpy(buf + at, num.ptr, num.len);
		memcpy(buf + at + num.len, path.ptr + at + 1, path.len - at - 1);
		buf[len] = '\0';
		return {buf, len};
	}
	return {};
}

// ## Visitor

Visitor Visitor::writer(serialization_state_t& state, md_str_t section) {
	Visitor v;
	v.ser = &state;
	v.section_name = section;
	return v;
}

Visitor Visitor::reader(md_str_t path, md_str_t arg) {
	Visitor v;
	v.rest = path;
	v.arg = arg;
	return v;
}

// The index of an item: digits, below cap. -1 otherwise.
static int item_index(md_str_t step, int cap) {
	if (!is_key_number(step)) return -1;
	const int i = (int)md_parse_int(step);
	return i < cap ? i : -1;
}

bool Visitor::take(const char* key) {
	MD_ASSERT(reading());
	if (done() || !md_str_eq_cstr(rest, key)) return false;
	taken = true;
	return true;
}

bool Visitor::enter(const char* key, md_str_t* after) const {
	MD_ASSERT(reading());
	if (done()) return false;
	md_str_t path = rest;
	md_str_t step;
	if (!path_next(step, path) || !md_str_eq_cstr(step, key)) return false;
	*after = path;
	return true;
}

md_str_t Visitor::key_path(const char* key) {
	MD_ASSERT(writing());
	// Behind a section, a field of the object itself would land in it
	MD_ASSERT(!(sectioned && depth == 0));
	const size_t len = strlen(key);
	if (prefix_len + len >= sizeof(key_buf)) {
		MD_LOG_ERROR("Workspace: the path of '%s' is too long to be written", key);
		return {};
	}
	memcpy(key_buf, prefix, prefix_len);
	memcpy(key_buf + prefix_len, key, len);
	key_buf[prefix_len + len] = '\0';
	return {key_buf, prefix_len + len};
}

void Visitor::push(FrameKind kind, const char* key, int index) {
	MD_ASSERT(depth < (int)MD_COUNTOF(frames));
	frames[depth++] = { key, index, prefix_len, kind };
	if (writing() && kind != Frame_Section) {
		const int n = index < 0 ? snprintf(prefix + prefix_len, sizeof(prefix) - prefix_len, "%s.", key)
		                        : snprintf(prefix + prefix_len, sizeof(prefix) - prefix_len, "%s.%d.", key, index);
		MD_ASSERT(n > 0 && prefix_len + (uint32_t)n < sizeof(prefix));
		prefix_len += (uint32_t)n;
	}
}

#define VISITOR_FIELD(write_call, read_call)       \
	if (writing()) {                              \
		const md_str_t k = key_path(key);         \
		if (!md_str_empty(k)) write_call;         \
		return false;                             \
	}                                             \
	return take(key) && (read_call)

bool Visitor::field(const char* key, bool& val)        { VISITOR_FIELD(write_bool(*ser, k, val), extract_bool(val, arg)); }
bool Visitor::field(const char* key, int& val)         { VISITOR_FIELD(write_int (*ser, k, val), extract_int (val, arg)); }
bool Visitor::field(const char* key, float& val)       { VISITOR_FIELD(write_flt (*ser, k, val), extract_flt (val, arg)); }
bool Visitor::field(const char* key, double& val)      { VISITOR_FIELD(write_dbl (*ser, k, val), extract_dbl (val, arg)); }
bool Visitor::field(const char* key, md_float2_t& val) { VISITOR_FIELD(write_flt_vec(*ser, k, val.elem, 2), extract_flt_vec(val.elem, 2, arg)); }
bool Visitor::field(const char* key, md_float3_t& val) { VISITOR_FIELD(write_vec3(*ser, k, val), extract_vec3(val, arg)); }
bool Visitor::field(const char* key, md_float4_t& val) { VISITOR_FIELD(write_vec4(*ser, k, val), extract_vec4(val, arg)); }
bool Visitor::field(const char* key, float* elem, size_t len) { VISITOR_FIELD(write_flt_vec(*ser, k, elem, len), extract_flt_vec(elem, len, arg)); }
bool Visitor::field(const char* key, md_str_t& val)    { VISITOR_FIELD(write_str(*ser, k, val), extract_str(val, arg)); }
bool Visitor::field(const char* key, char* buf, size_t cap) { VISITOR_FIELD(write_str(*ser, k, md_str_from_cstr(buf)), extract_to_char_buf(buf, cap, arg)); }

#undef VISITOR_FIELD

bool Visitor::field_enum(const char* key, int& val, const EnumNames& names) {
	if (writing()) {
		const md_str_t k = key_path(key);
		if (md_str_empty(k)) return false;
		if (0 <= val && val < names.count) {
			write_str(*ser, k, md_str_from_cstr(names.names[val]));
		} else {
			write_int(*ser, k, val);
		}
		return false;
	}
	if (!take(key)) return false;
	for (int i = 0; i < names.count; ++i) {
		if (md_str_eq_cstr(arg, names.names[i])) {
			val = i;
			return true;
		}
	}
	int i;
	if (extract_int(i, arg) && 0 <= i && i < names.count) {
		val = i;
		return true;
	}
	MD_LOG_DEBUG("Workspace: '%s' has no value '" MD_STR_FMT "', it is skipped", key, MD_STR_ARG(arg));
	return false;
}

bool Visitor::field_index_int(const char* key, int& val, int count) {
	if (writing()) {
		const md_str_t k = key_path(key);
		if (!md_str_empty(k)) write_int(*ser, k, val);
		return false;
	}
	if (!take(key)) return false;
	int i;
	if (extract_int(i, arg) && 0 <= i && i < count) {
		val = i;
		return true;
	}
	return false;
}

bool Visitor::field_named(const char* key, int& val, const char* (*name_of)(int), int (*index_of)(md_str_t)) {
	if (writing()) {
		const md_str_t k = key_path(key);
		if (md_str_empty(k)) return false;
		const char* name = name_of(val);
		if (name) {
			write_str(*ser, k, md_str_from_cstr(name));
		} else {
			write_int(*ser, k, val);
		}
		return false;
	}
	if (!take(key)) return false;
	const int idx = index_of(arg);
	if (idx >= 0) {
		val = idx;
		return true;
	}
	int i;
	if (extract_int(i, arg) && name_of(i)) {
		val = i;
		return true;
	}
	MD_LOG_DEBUG("Workspace: '%s' has no value '" MD_STR_FMT "', it is skipped", key, MD_STR_ARG(arg));
	return false;
}

bool Visitor::list(const char* key, float* elem, int& count, int cap, int arity) {
	MD_ASSERT(arity > 0 && count >= 0 && count <= cap);
	if (writing()) {
		const md_str_t k = key_path(key);
		if (!md_str_empty(k)) write_flt_vec(*ser, k, elem, (size_t)(count * arity));
		return false;
	}
	md_str_t after;
	if (!enter(key, &after)) return false;
	taken = true;
	if (md_str_empty(after)) {
		// The whole list
		float tmp[64];
		MD_ASSERT((size_t)(cap * arity) <= MD_COUNTOF(tmp));
		size_t n = 0;
		if (!extract_flt_list(tmp, &n, (size_t)(cap * arity), arg) || n % (size_t)arity != 0) return false;
		md_memcpy(elem, tmp, n * sizeof(float));
		count = (int)(n / (size_t)arity);
		return true;
	}
	// One item of it
	md_str_t step;
	path_next(step, after);
	const int i = md_str_empty(after) ? item_index(step, cap) : -1;
	if (i < 0 || !extract_flt_vec(elem + (size_t)i * (size_t)arity, (size_t)arity, arg)) return false;
	count = MD_MAX(count, i + 1);
	return true;
}

bool Visitor::block(const char* key) {
	if (writing()) {
		push(Frame_Block, key, -1);
		return true;
	}
	md_str_t after;
	if (!enter(key, &after) || md_str_empty(after)) return false;
	rest = after;
	push(Frame_Block, key, -1);
	return true;
}

bool Visitor::section(const char* key) {
	if (writing()) {
		// Of the object itself, not of a block within it
		MD_ASSERT(depth == 0);
		char buf[256];
		const md_str_t name = path_join(buf, sizeof(buf), section_name, md_str_from_cstr(key));
		if (md_str_empty(name)) return false;
		write_section_header(*ser, name);
		push(Frame_Section, key, -1);
		sectioned = true;
		return true;
	}
	return block(key);
}

bool Visitor::each(const char* key, int& i, int& count, int cap) {
	MD_ASSERT(cap > 0);
	// The call after an item: the item ends
	const bool after_item = depth > 0 && frames[depth - 1].kind == Frame_Each && strcmp(frames[depth - 1].key, key) == 0;

	if (writing()) {
		int next = 0;
		if (after_item) {
			depth -= 1;
			prefix_len = frames[depth].prefix_len;
			next = frames[depth].index + 1;
		} else {
			// How many: an array that is empty says so
			const md_str_t k = key_path(key);
			if (md_str_empty(k)) return false;
			write_int(*ser, k, MD_CLAMP(count, 0, cap));
		}
		if (next >= MD_MIN(count, cap)) return false;
		push(Frame_Each, key, next);
		i = next;
		return true;
	}

	if (after_item) {
		// The entry was in that item: nothing after it can be it
		depth -= 1;
		finished = true;
		return false;
	}
	md_str_t after;
	if (!enter(key, &after)) return false;
	if (md_str_empty(after)) {
		// key=<count>
		taken = true;
		int n;
		if (extract_int(n, arg)) count = MD_CLAMP(n, 0, cap);
		return false;
	}
	md_str_t step;
	path_next(step, after);
	const int idx = item_index(step, cap);
	if (idx < 0 || md_str_empty(after)) {
		// Within the array, but no field of an item of it
		finished = true;
		return false;
	}
	rest = after;
	count = MD_MAX(count, idx + 1);
	push(Frame_Each, key, idx);
	i = idx;
	return true;
}

void Visitor::end() {
	MD_ASSERT(depth > 0);
	depth -= 1;
	if (writing()) {
		prefix_len = frames[depth].prefix_len;
	} else {
		// The entry was within the block: nothing after it can be it
		finished = true;
	}
}

}
