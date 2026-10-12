#include <serialization_utils.h>
#include <core/md_allocator.h>
#include <core/md_parse.h>
#include <core/md_bitfield.h>
#include <core/md_base64.h>
#include <core/md_log.h>

#include <string.h>

namespace viamd {

static const str_t esc = STR_INIT("\"\"\"");

static bool is_comment(str_t line) {
	return str_begins_with(line, STR_LIT("#")) || str_begins_with(line, STR_LIT(";"));
}

bool next_section_header(str_t& section, deserialization_state_t& state) {
	// Whatever the previous section's reader did not consume
	str_t ident, arg;
	while (next_entry(ident, arg, state)) {}

	str_t line;
	while (str_extract_line(&line, &state.text)) {
		line = str_trim(line);
		if (str_begins_with(line, STR_LIT("[")) &&
			str_ends_with(line, STR_LIT("]")))
		{
			section = str_trim(str_substr(line, 1, str_len(line) - 2));
			state.cur_section = section;
			return true;
		}
	}
	state.cur_section = {};
	return false;
}

bool next_entry(str_t& ident, str_t& arg, deserialization_state_t& state) {
	str_t line;
	while (str_peek_line(&line, &state.text)) {
		line = str_trim(line);
		if (str_begins_with(line, STR_LIT("["))) {
			return false;
		}
		str_skip_line(&state.text);
		if (is_comment(line)) {
			continue;
		}

		size_t loc;
		if (str_find_char(&loc, line, '=')) {
			ident = str_trim(str_substr(line, 0, loc));
			arg   = str_trim(str_substr(line, loc + 1, SIZE_MAX));
			if (str_begins_with(arg, esc)) {
				// Multiline string: from the opening quotes in the text itself (the line was only a
				// view of its first line) to the matching closing ones
				const char* beg = str_beg(arg) + str_len(esc);
				const char* end = str_end(state.text);
				str_t haystack = {beg, (size_t)(end-beg)};
				if (str_find_str(&loc, haystack, esc)) {
					arg = {beg, loc};
					state.text = str_substr(haystack, loc + str_len(esc));
					// The rest of the closing line
					str_skip_line(&state.text);
				} else {
					MD_LOG_ERROR("Workspace: unbalanced \"\"\" in the value of '" STR_FMT "', the rest of the file is skipped", STR_ARG(ident));
					state.text = {};
					return false;
				}
			}
			return true;
		}
	}
	return false;
}

void write_section_header(serialization_state_t& state, str_t section) {
	md_strb_fmt(&state.sb, "\n[" STR_FMT "]\n", STR_ARG(section));
}

void write_int(serialization_state_t& state, str_t ident, int64_t val) {
	md_strb_fmt(&state.sb, STR_FMT "=%lld\n", STR_ARG(ident), (long long)val);
}

void write_int_vec(serialization_state_t& state, str_t ident, const int* elem, size_t len) {
	md_strb_fmt(&state.sb, STR_FMT "=", STR_ARG(ident));
	for (size_t i = 0; i < len; ++i) {
		md_strb_fmt(&state.sb, "%i", elem[i]);
		if (i < len - 1) {
			md_strb_push_char(&state.sb, ',');
		}
	}
	md_strb_push_char(&state.sb, '\n');
}

// Shortest round trip is not on offer from printf; these digit counts are the ones that always
// read back the same value (FLT_DECIMAL_DIG and DBL_DECIMAL_DIG). %f would print 1e-7 as 0.000000.
void write_flt(serialization_state_t& state, str_t ident, float val) {
	md_strb_fmt(&state.sb, STR_FMT "=%.9g\n", STR_ARG(ident), (double)val);
}

void write_dbl(serialization_state_t& state, str_t ident, double val) {
	md_strb_fmt(&state.sb, STR_FMT "=%.17g\n", STR_ARG(ident), val);
}

void write_flt_vec(serialization_state_t& state, str_t ident, const float* elem, size_t len) {
	md_strb_fmt(&state.sb, STR_FMT "=", STR_ARG(ident));
	for (size_t i = 0; i < len; ++i) {
		md_strb_fmt(&state.sb, "%.9g", (double)elem[i]);
		if (i < len - 1) {
			md_strb_push_char(&state.sb, ',');
		}
	}
	md_strb_push_char(&state.sb, '\n');
}

void write_str(serialization_state_t& state, str_t ident, str_t str) {
	// Quoted whenever a plain value would not come back as it went: a line break ends it, and the
	// reader trims whitespace at either end
	const bool multiline  = str_find_char(NULL, str, '\n') || str_find_char(NULL, str, '\r');
	const bool padded     = str.len > 0 && (is_whitespace(str.ptr[0]) || is_whitespace(str.ptr[str.len - 1]));
	const bool looks_esc  = str_begins_with(str, esc);
	if (multiline || padded || looks_esc) {
		if (str_find_str(NULL, str, esc)) {
			MD_LOG_ERROR("Workspace: the value of '" STR_FMT "' contains \"\"\" and cannot be stored as it is", STR_ARG(ident));
		}
		md_strb_fmt(&state.sb, STR_FMT "=" STR_FMT STR_FMT STR_FMT "\n", STR_ARG(ident), STR_ARG(esc), STR_ARG(str), STR_ARG(esc));
	} else {
		md_strb_fmt(&state.sb, STR_FMT "=" STR_FMT "\n", STR_ARG(ident), STR_ARG(str));
	}
}

void write_bitfield(serialization_state_t& state, str_t ident, const md_bitfield_t* bf) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

	void*  serialized_data = md_temp_alloc(temp, md_bitfield_serialize_size_in_bytes(bf));
	size_t serialized_size = md_bitfield_serialize(serialized_data, bf);
	if (serialized_size) {
		char*  base64_data = (char*)md_temp_alloc(temp, md_base64_encode_size_in_bytes(serialized_size));
		size_t base64_size = md_base64_encode(base64_data, serialized_data, serialized_size);
		if (base64_size) {
			str_t base64 = {base64_data, base64_size};
			// The trailing newline is not cosmetic. Every other writer here terminates its line, and
			// next_entry() splits a line at its first '='. Without it the following entry is glued
			// onto the end of this one and silently swallowed - which is how every group label
			// after the first went missing from a saved workspace.
			md_strb_fmt(&state.sb, STR_FMT "=###" STR_FMT "###\n", STR_ARG(ident), STR_ARG(base64));
		}
	}
}

bool extract_bool(bool& val, str_t arg) {
	if (is_int(arg)) {
		int64_t integer = parse_int(arg);
		if (integer == 0 || integer == 1) {
			val = (bool)integer;
			return true;
		}
	} else if (str_eq_cstr_ignore_case(arg, "true")) {
		val = true;
		return true;
	} else if (str_eq_cstr_ignore_case(arg, "false")) {
		val = false;
		return true;
	}
	return false;
}

bool extract_int(int& val, str_t arg) {
	if (is_int(arg)) {
		const int64_t i = parse_int(arg);
		if (INT32_MIN <= i && i <= INT32_MAX) {
			val = (int)i;
			return true;
		}
	}
	return false;
}

bool extract_int64(int64_t& val, str_t arg) {
	if (is_int(arg)) {
		val = parse_int(arg);
		return true;
	}
	return false;
}

bool extract_int_vec(int* elem, size_t len, str_t arg) {
	str_t tok;
	size_t count = 0;
	while (count < len && extract_token_delim(&tok, &arg, ',')) {
		tok = str_trim(tok);
		if (is_int(tok)) {
			elem[count++] = (int)parse_int(tok);
		}
	}
	return count == len;
}

bool extract_dbl (double& val, str_t arg) {
	if (is_float(arg)) {
		val = parse_float(arg);
		return true;
	}
	return false;
}

bool extract_flt (float& val, str_t arg) {
	if (is_float(arg)) {
		val = (float)parse_float(arg);
		return true;
	}
	return false;
}

bool extract_flt_vec (float* elem, size_t len, str_t arg) {
	// Parsed into a copy and written back only when every component was there
	float tmp[32];
	if (len > ARRAY_SIZE(tmp)) return false;
	str_t tok;
	size_t count = 0;
	while (count < len && extract_token_delim(&tok, &arg, ',')) {
		tok = str_trim(tok);
		if (is_float(tok) || is_int(tok)) {
			tmp[count++] = (float)parse_float(tok);
		}
	}
	if (count != len) return false;
	// Anything left over means this is not the vector asked for (e.g. a longer one)
	while (extract_token_delim(&tok, &arg, ',')) {
		if (str_trim(tok).len > 0) return false;
	}
	MEMCPY(elem, tmp, len * sizeof(float));
	return true;
}

bool extract_str(str_t& str, str_t arg) {
	// next_entry already strips the quotes of a multiline value; this covers a value handed in
	// with them still on
	if (arg.len >= 2 * esc.len && str_begins_with(arg, esc) && str_ends_with(arg, esc)) {
		str = str_substr(arg, esc.len, arg.len - 2 * esc.len);
	} else {
		str = arg;
	}
	return true;
}

bool extract_to_char_buf(char* buf, size_t cap, str_t arg) {
	str_t str;
	if (extract_str(str, arg)) {
		str_copy_to_char_buf(buf, cap, str);
		return true;
	}
	return false;
}

bool extract_bitfield(md_bitfield_t* bf, str_t arg) {
	md_temp_scope_t temp = md_temp_begin();
	defer { md_temp_end(temp); };

	// Bitfield starts with ###
	// and ends with ###
	str_t token = STR_INIT("###");
	if (!str_eq_n(arg, token, str_len(token))) {
		MD_LOG_ERROR("Malformed start token for bitfield");
		return false;
	}
	arg = str_substr(arg, str_len(token));

	size_t loc;
	if (!str_find_str(&loc, arg, token)) {
		MD_LOG_ERROR("Malformed end token for bitfield");
		return false;
	}

	arg = str_substr(arg, 0, loc);
	const size_t raw_cap = md_base64_decode_size_in_bytes(str_len(arg));
	void* raw_ptr = md_temp_alloc(temp, raw_cap);

	const size_t raw_len = md_base64_decode(raw_ptr, str_ptr(arg), str_len(arg));
	if (!raw_len || !md_bitfield_deserialize(bf, raw_ptr, raw_len)) {
		MD_LOG_ERROR("Failed to deserialize bitfield");
		md_bitfield_clear(bf);
		return false;
	}

	return true;
}

}
