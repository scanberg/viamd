#include <serialization_utils.h>
#include <core/md_parse.h>
#include <core/md_bitfield.h>
#include <core/md_base64.h>
#include <core/md_log.h>

namespace viamd {

static const str_t esc = STR_LIT("\"\"\"");

bool next_section_header(str_t& section, deserialization_state_t& state) {
	str_t line;
	while (str_extract_line(&line, &state.text)) {
		line = str_trim(line);
		if (str_begins_with(line, STR_LIT("[")) &&
			str_ends_with(line, STR_LIT("]")))
		{
			section = str_substr(line, 1, str_len(line) - 2);
			state.cur_section = section;
			return true;
		}
	}
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

		size_t loc;
		if (str_find_char(&loc, line, '=')) {
			ident = str_trim(str_substr(line, 0, loc));
			arg   = str_trim(str_substr(line, loc + 1, SIZE_MAX));
			if (str_begins_with(arg, esc)) {
				// Multiline string, find matching espace sequence
				const char* beg = str_beg(arg) + str_len(esc);
				const char* end = str_end(state.text);
				str_t haystack = {beg, (size_t)(end-beg)};
				if (str_find_str(&loc, haystack, esc)) {
					arg = {beg, loc};
					state.text = str_substr(state.text, loc + str_len(esc));
				} else {
					// Error
					MD_LOG_ERROR("Unbalanced escape sequence in multiline string");
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
	md_strb_fmt(&state.sb, STR_FMT "=%i\n", STR_ARG(ident), (int)val);
}

void write_dbl(serialization_state_t& state, str_t ident, double val) {
	md_strb_fmt(&state.sb, STR_FMT "=%f\n", STR_ARG(ident), val);
}

void write_flt_vec(serialization_state_t& state, str_t ident, const float* elem, size_t len) {
	md_strb_fmt(&state.sb, STR_FMT "=", STR_ARG(ident));
	for (size_t i = 0; i < len; ++i) {
		md_strb_fmt(&state.sb, "%f", elem[i]);
		if (i < len - 1) {
			md_strb_push_char(&state.sb, ',');
		}
	}
	md_strb_push_char(&state.sb, '\n');
}

void write_str(serialization_state_t& state, str_t ident, str_t str) {
	if (str_find_char(NULL, str, '\n')) {
		md_strb_fmt(&state.sb, STR_FMT "=" STR_FMT STR_FMT STR_FMT "\n", STR_ARG(ident), STR_ARG(esc), STR_ARG(str), STR_ARG(esc));
	} else {
		md_strb_fmt(&state.sb, STR_FMT "=" STR_FMT "\n", STR_ARG(ident), STR_ARG(str));
	}
}

void write_bitfield(serialization_state_t& state, str_t ident, const md_bitfield_t* bf) {
	size_t tmp_pos = md_temp_get_pos();
	defer { md_temp_set_pos_back(tmp_pos); };

	void*  serialized_data = md_temp_push(md_bitfield_serialize_size_in_bytes(bf));
	size_t serialized_size = md_bitfield_serialize(serialized_data, bf);
	if (serialized_size) {
		char*  base64_data = (char*)md_temp_push(md_base64_encode_size_in_bytes(serialized_size));
		size_t base64_size = md_base64_encode(base64_data, serialized_data, serialized_size);
		if (base64_size) {
			str_t base64 = {base64_data, base64_size};
			md_strb_fmt(&state.sb, STR_FMT "=###" STR_FMT "###", STR_ARG(ident), STR_ARG(base64));
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
	}
	return false;
}

bool extract_int(int& val, str_t arg) {
	if (is_int(arg)) {
		val = (int)parse_int(arg);
		return true;
	}
	return false;
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
	str_t tok;
	size_t count = 0;
	while (count < len && extract_token_delim(&tok, &arg, ',')) {
		if (is_float(arg) || is_int(arg)) {
			elem[count++] = (float)parse_float(arg);
		}
	}
	if (count < len && extract_token(&tok, &arg)) {
		if (is_float(arg) || is_int(arg)) {
			elem[count++] = (float)parse_float(arg);
		}
	}
	return count == len;
}

bool extract_str(str_t& str, str_t arg) {
	if (str_begins_with(arg, esc) && str_ends_with(arg, esc)) {
		str = str_substr(arg, 3, arg.len - 6);
	}
	str = arg;
	return true;
}

bool extract_bitfield(md_bitfield_t* bf, str_t arg) {
	size_t tmp_pos = md_temp_get_pos();
	defer { md_temp_set_pos_back(tmp_pos); };

	// Bitfield starts with ###
	// and ends with ###
	str_t token = STR_LIT("###");
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
	void* raw_ptr = md_temp_push(raw_cap);

	const size_t raw_len = md_base64_decode(raw_ptr, str_ptr(arg), str_len(arg));
	if (!raw_len || !md_bitfield_deserialize(bf, raw_ptr, raw_len)) {
		MD_LOG_ERROR("Failed to deserialize bitfield");
		md_bitfield_clear(bf);
		return false;
	}

	return true;
}

}
