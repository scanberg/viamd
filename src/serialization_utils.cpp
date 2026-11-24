#include <serialization_utils.h>
#include <core/md_allocator.h>
#include <core/md_parse.h>
#include <core/md_bitset.h>
#include <core/md_base64.h>
#include <core/md_log.h>
#include <core/md_array.h>

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

void write_dbl_vec(serialization_state_t& state, md_str_t ident, const double* elem, size_t len) {
	md_str_writer_fmt(&state.sw, MD_STR_FMT "=", MD_STR_ARG(ident));
	for (size_t i = 0; i < len; ++i) {
		md_str_writer_fmt(&state.sw, "%.17g", elem[i]);
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

bool extract_dbl_vec (double* elem, size_t len, md_str_t arg) {
	// Exactly len numbers, parsed into a copy and written back only then. More than len is not cut down to the
	// first len: that is what a vector written with decimal commas looks like ("0,5,0,25,1,1" for 0.5 0.25 1 1),
	// and its first len would be values nobody wrote.
	double tmp[16];
	if (len > MD_COUNTOF(tmp)) return false;
	const md_str_t all = arg;
	md_str_t tok;
	size_t count = 0;
	size_t extra = 0;
	while (md_parse_next_token_delim(&tok, &arg, ',')) {
		tok = md_str_trim(tok);
		if (!(md_parse_is_float(tok) || md_parse_is_int(tok))) return false;
		if (count < len) {
			tmp[count++] = md_parse_float(tok);
		} else {
			extra += 1;
		}
	}
	if (extra) log_decimal_comma(all);
	if (count != len || extra) return false;
	md_memcpy(elem, tmp, len * sizeof(double));
	return true;
}

bool extract_flt_vec (float* elem, size_t len, md_str_t arg) {
	double tmp[16];
	if (len > MD_COUNTOF(tmp) || !extract_dbl_vec(tmp, len, arg)) return false;
	for (size_t i = 0; i < len; ++i) {
		elem[i] = (float)tmp[i];
	}
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

PathPrefix::PathPrefix(const char* prefix) {
	const int n = snprintf(buf, sizeof(buf), "%s", prefix);
	MD_ASSERT(n >= 0 && (size_t)n < sizeof(buf));
	len = (uint32_t)n;
}

PathPrefix PathPrefix::sub(const char* key, int index) const {
	PathPrefix p = *this;
	const int n = index < 0 ? snprintf(p.buf + len, sizeof(p.buf) - len, "%s.", key)
	                        : snprintf(p.buf + len, sizeof(p.buf) - len, "%s.%d.", key, index);
	MD_ASSERT(n > 0 && len + (uint32_t)n < sizeof(p.buf));
	p.len = len + (uint32_t)n;
	return p;
}

md_str_t PathPrefix::operator()(const char* key) {
	const int n = snprintf(key_buf, sizeof(key_buf), "%.*s%s", (int)len, buf, key);
	MD_ASSERT(n > 0 && (size_t)n < sizeof(key_buf));
	return {key_buf, (size_t)n};
}

// ## Entries

Entries read_entries(deserialization_state_t& state, const PathAlias* aliases, size_t num_aliases, md_allocator_i* alloc) {
	Entries e = {};
	e.alloc = alloc;
	md_str_t ident, arg;
	while (next_entry(ident, arg, state)) {
		Entry entry = { ident, arg, false };
		// A key from before the paths is one step: where it is now, when it only moved
		if (num_aliases && md_str_find_char(ident, '.') == MD_STR_NPOS) {
			char buf[256];
			const md_str_t moved = path_alias(buf, sizeof(buf), ident, aliases, num_aliases);
			if (!md_str_empty(moved)) entry.path = md_alloc_str(alloc, moved);
		}
		md_array_push(e.items, entry, alloc);
	}
	e.count = md_array_size(e.items);
	return e;
}

void entries_unread(const Entries& e, md_str_t section) {
	for (size_t i = 0; i < e.count; ++i) {
		if (!e.items[i].read) {
			MD_LOG_DEBUG("Workspace: [" MD_STR_FMT "] '" MD_STR_FMT "' is read by nothing, it is skipped", MD_STR_ARG(section), MD_STR_ARG(e.items[i].path));
		}
	}
}

const md_str_t* get_value(Entries& e, md_str_t path) {
	const md_str_t* val = nullptr;
	for (size_t i = 0; i < e.count; ++i) {
		if (md_str_eq(e.items[i].path, path)) {
			e.items[i].read = true;
			val = &e.items[i].value;
		}
	}
	return val;
}

static bool not_read(md_str_t path, md_str_t value) {
	MD_LOG_DEBUG("Workspace: '" MD_STR_FMT "' has a value that does not read, '" MD_STR_FMT "': it is skipped", MD_STR_ARG(path), MD_STR_ARG(value));
	return false;
}

#define ENTRIES_GET(extract_call)                             \
	const md_str_t* v = get_value(e, path);                  \
	if (!v) return false;                                     \
	const md_str_t arg = *v;                                  \
	return (extract_call) || not_read(path, arg)

bool get(Entries& e, md_str_t path, bool& val)        { ENTRIES_GET(extract_bool(val, arg)); }
bool get(Entries& e, md_str_t path, int& val)         { ENTRIES_GET(extract_int(val, arg)); }
bool get(Entries& e, md_str_t path, float& val)       { ENTRIES_GET(extract_flt(val, arg)); }
bool get(Entries& e, md_str_t path, double& val)      { ENTRIES_GET(extract_dbl(val, arg)); }
bool get(Entries& e, md_str_t path, md_float2_t& val) { ENTRIES_GET(extract_flt_vec(val.elem, 2, arg)); }
bool get(Entries& e, md_str_t path, md_float3_t& val) { ENTRIES_GET(extract_vec3(val, arg)); }
bool get(Entries& e, md_str_t path, md_float4_t& val) { ENTRIES_GET(extract_vec4(val, arg)); }
bool get(Entries& e, md_str_t path, float* elem, size_t len) { ENTRIES_GET(extract_flt_vec(elem, len, arg)); }
bool get(Entries& e, md_str_t path, char* buf, size_t cap)   { ENTRIES_GET(extract_to_char_buf(buf, cap, arg)); }
bool get(Entries& e, md_str_t path, md_str_t& val)    { ENTRIES_GET(extract_str(val, arg)); }

#undef ENTRIES_GET

bool get_enum_int(Entries& e, md_str_t path, int& val, const EnumNames& names) {
	const md_str_t* v = get_value(e, path);
	if (!v) return false;
	for (int i = 0; i < names.count; ++i) {
		if (md_str_eq_cstr(*v, names.names[i])) {
			val = i;
			return true;
		}
	}
	int i;
	if (extract_int(i, *v) && 0 <= i && i < names.count) {
		val = i;
		return true;
	}
	return not_read(path, *v);
}

void write_enum(serialization_state_t& state, md_str_t path, int val, const EnumNames& names) {
	if (0 <= val && val < names.count) {
		write_str(state, path, md_str_from_cstr(names.names[val]));
	} else {
		write_int(state, path, val);
	}
}

bool get_named(Entries& e, md_str_t path, int& val, int (*index_of)(md_str_t)) {
	const md_str_t* v = get_value(e, path);
	if (!v) return false;
	const int i = index_of(*v);
	if (i >= 0) {
		val = i;
		return true;
	}
	return extract_int(val, *v) || not_read(path, *v);
}

void write_named(serialization_state_t& state, md_str_t path, int val, const char* (*name_of)(int)) {
	if (const char* name = name_of(val)) {
		write_str(state, path, md_str_from_cstr(name));
	} else {
		write_int(state, path, val);
	}
}

// The index an entry is an item of, path.<i> or path.<i>.<...>: -1 for none
static int item_index(md_str_t entry, md_str_t path) {
	if (entry.len <= path.len + 1 || !md_str_begins_with(entry, path) || entry.ptr[path.len] != '.') return -1;
	md_str_t rest = md_str_substr(entry, path.len + 1, SIZE_MAX);
	const size_t dot = md_str_find_char(rest, '.');
	const md_str_t num = md_str_substr(rest, 0, dot);
	if (num.len == 0 || num.len > 4) return -1;
	int i = 0;
	for (size_t k = 0; k < num.len; ++k) {
		if (num.ptr[k] < '0' || num.ptr[k] > '9') return -1;
		i = i * 10 + (num.ptr[k] - '0');
	}
	return i;
}

bool get_list(Entries& e, md_str_t path, float* elem, int& count, int cap, int arity) {
	MD_ASSERT(arity > 0 && cap >= 0);
	bool found = false;
	if (const md_str_t* v = get_value(e, path)) {
		float tmp[256];
		MD_ASSERT((size_t)(cap * arity) <= MD_COUNTOF(tmp));
		size_t n = 0;
		if (extract_flt_list(tmp, &n, (size_t)(cap * arity), *v) && n % (size_t)arity == 0) {
			md_memcpy(elem, tmp, sizeof(float) * n);
			count = (int)(n / (size_t)arity);
			found = true;
		} else {
			not_read(path, *v);
		}
	}
	// An older workspace's items, a key each
	for (size_t i = 0; i < e.count; ++i) {
		const int idx = item_index(e.items[i].path, path);
		if (idx < 0) continue;
		e.items[i].read = true;
		if (idx >= cap || !extract_flt_vec(elem + idx * arity, (size_t)arity, e.items[i].value)) continue;
		count = MD_MAX(count, idx + 1);
		found = true;
	}
	return found;
}

void write_list(serialization_state_t& state, md_str_t path, const float* elem, int count, int arity) {
	// An empty list is a value too: no items
	if (count == 0) {
		write_str(state, path, MD_STR_LIT(""));
		return;
	}
	write_flt_vec(state, path, elem, (size_t)(count * arity));
}

bool get_count(Entries& e, md_str_t path, int& count, int cap) {
	bool found = false;
	int n = 0;
	if (get(e, path, n)) {
		count = MD_CLAMP(n, 0, cap);
		found = true;
	}
	for (size_t i = 0; i < e.count; ++i) {
		const int idx = item_index(e.items[i].path, path);
		if (idx < 0 || idx >= cap) continue;
		count = MD_MAX(count, idx + 1);
		found = true;
	}
	return found;
}

}
