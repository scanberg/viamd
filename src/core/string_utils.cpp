#include "string_utils.h"
#include <core/common.h>
#include <fstream>
#include <ctype.h>

#ifdef WIN32
#pragma warning(disable:4996) // Unsafe strncpy
#endif

#define MIN(x,y) ((x < y) ? (x) : (y))
#define MAX(x,y) ((x > y) ? (x) : (y))

static inline bool internal_compare(const char* str_a, const char* str_b, int64 len, bool ignore_case) {
	if (ignore_case) {
		for (int64 i = 0; i < len; i++) {
			if (tolower(str_a[i]) != tolower(str_b[i])) return false;
		}
	}
	else {
		for (int64 i = 0; i < len; i++) {
			if (str_a[i] != str_b[i]) return false;
		}
	}
	return true;
}

bool compare(CString str_a, CString str_b, bool ignore_case) {
	//int64 len = MIN(str_a.count, str_b.count);
	if (str_a.count != str_b.count) return false;
	if (str_a.count == 0) return false;
	return internal_compare(str_a, str_b, str_a.count, ignore_case);
}

bool compare_n(CString str_a, CString str_b, int64 num_chars, bool ignore_case) {
	int64 len = MIN(str_a.count, MIN(str_b.count, num_chars));
	if (len == 0) return false;
	return internal_compare(str_a, str_b, len, ignore_case);
}

void copy(String dst, CString src) {
	ASSERT(dst.data != 0);
	ASSERT(src.data != 0);
	auto len = MIN(dst.count, src.count);
	memcpy(dst.data, src.data, len);
	dst.data[src.count] = '\0';
}

void copy_n(String dst, CString src, int64 num_chars) {
	ASSERT(dst.data != 0);
	ASSERT(src.data != 0);
	auto len = MIN(dst.count, src.count);
	len = MIN(len, num_chars);
	memcpy(dst.data, src.data, len);
	dst.data[num_chars] = '\0';
}

String allocate_string(CString str) {
	if (str.count == 0) return {};
	char* data = (char*)MALLOC(str.count + 1);
	strncpy(data, str.data, str.count + 1);
	return { data, str.count };
}

String allocate_string(int32 length) {
	if (length == 0) return {};
	char* data = (char*)MALLOC(length);
	return { data, length };
}

void free_string(String* str) {
	if (str->data) {
		FREE(str->data);
		str->data = nullptr;
		str->count = 0;
	}
}

bool extract_line(CString& line, CString& str) {
	const char* str_beg = str.data;
	const char* str_end = str.data + str.count;

	if (str_beg == str_end) {
		line = {};
		return false;
	}

	const char* line_beg = str_beg;
	const char* line_end = line_beg;

	// Find return or new line character
	while (line_end < str_end && *line_end != '\r' && *line_end != '\n') ++line_end;

	// Step over return or new line characters
	str_beg = MIN(line_end + 1, str_end);
	while (str_beg < str_end && (*str_beg == '\r' || *str_beg == '\n')) ++str_beg;

	line.count = line_end - line_beg;
	line.data = line_beg;

	str.data = str_beg;
	str.count = str_end - str_beg;

	return true;
}

bool copy_line(String& line, CString& str) {
	const char* str_beg = str.data;
	const char* str_end = str.data + str.count;

	if (str_beg == str_end) {
		line = {};
		return false;
	}

	const char* line_beg = str_beg;
	const char* line_end = line_beg;

	// Find return or new line character
	while (line_end < str_end && *line_end != '\r' && *line_end != '\n') ++line_end;

	// Step over return or new line characters
	str_beg = MIN(line_end + 1, str_end);
	while (str_beg < str_end && (*str_beg == '\r' || *str_beg == '\n')) ++str_beg;

	// @NOTE: Do not modify line.count, its value contains the length of the buffer its pointing to 
	auto count = MIN(line_end - line_beg, line.count - 1);
	line.data = (char*)memcpy(line.data, line_beg, count);
	line.data[count] = '\0';

	str.data = str_beg;
	str.count = str_end - str_beg;

	return true;
}

ConversionResult<float32> to_float32(CString str) {
	// Make sure that the string passed into atof is zero-terminated
	StringBuffer<32> buf = str;
	char* end = nullptr;
	float32 val = strtof(buf, &end);
	return { val, end != buf.beg() };
}

ConversionResult<float64> to_float64(CString str) {
	// Make sure that the string passed into atof is zero-terminated
	StringBuffer<32> buf = str;
	char* end = nullptr;
	float64 val = strtod(buf, &end);
	return { val, end != buf.beg() };
}

ConversionResult<int32> to_int32(CString str) {
	// Make sure that the string passed into atof is zero-terminated
	StringBuffer<32> buf = str;
	char* end = nullptr;
	int32 val = strtol(buf, &end, 10);
	return { val, end != buf.beg() };
}

ConversionResult<int64> to_int64(CString str) {
	// Make sure that the string passed into atof is zero-terminated
	StringBuffer<32> buf = str;
	char* end = nullptr;
	int64 val = strtoll(buf, &end, 10);
	return { val, end != buf.beg() };
}

CString trim(CString str) {
	const char* beg = str.data;
	const char* end = str.data + str.count;

	while (beg < end && isspace(*beg)) ++beg;
	while (end > beg && (isspace(*(end-1)) || *(end-1) == '\0')) --end;

	return CString(beg, end - beg);
}

String trim(String str) {
	char* beg = str.data;
	char* end = str.data + str.count;

	while (beg < end && isspace(*beg)) ++beg;
	while (end > beg && isspace(*(end-1))) --end;

	return String(beg, end - beg);
}

String allocate_and_read_textfile(CString filename) {
	std::ifstream file(filename);
	if (!file) return {};

	file.seekg(0, std::ios::end);
	int64 file_size = file.tellg();
	file.seekg(0, std::ios::beg);

	char* data = (char*)MALLOC(file_size + 1);
	file.read(data, file_size);
	data[file_size] = '\0';

	return {data, file_size + 1};
}

CString get_directory(CString url) {
    if (url.count == 0) {
        return url;
    }
    
    url = trim(url);
    
    const char* beg = url.begin();
    const char* end = url.end() - 1;
    
    while (end != beg && *end != '\\' && *end != '/') {
        end--;
    }
    
    return CString(beg, end - beg);
}

CString get_file(CString url) {
    if (url.count == 0) {
        return url;
    }

    url = trim(url);
    
    const char* beg = url.end() - 1;
    const char* end = url.end();
    
    while (beg != url.begin() && *beg != '\\' && *beg != '/') {
        beg--;
    }
	if (*beg == '\\' || *beg == '/') beg++;
    
    return CString(beg, end - beg);
}

CString get_file_without_extension(CString url) {
    if (url.count == 0) {
        return url;
    }
    
    url = trim(url);
    
    const char* beg = url.end() - 1;
    const char* end = url.end();
    
    while (beg != url.begin() && *beg != '\\' && *beg != '/') beg--;
    while (end != beg && *beg != '.') end--;
    
    return CString(beg, end - beg);
}

CString get_file_extension(CString url) {
    if (url.count == 0) {
        return url;
    }
    
    url = trim(url);
    
    const char* beg = url.end() - 1;
    const char* end = url.end();
    
    while (beg != url.begin() && *beg != '.' && *beg != '\\' && *beg != '/') beg--;
    
    if (beg == url.begin() || *beg == '\\' || *beg == '/') {
        return CString();
    }
    
    return CString(beg + 1, end - beg);
}

inline static bool char_in_string(char c, CString str) {
	for (int64 i = 0; i < str.count; i++) {
		if (c == str[i]) return true;
	}
	return false;
}

StringBuffer<256> get_relative_path(CString from, CString to) {
	const char* c_from = from.beg();
	const char* c_to = to.beg();
	while (c_from != from.end() && c_to != to.end() && *c_from == *c_to) {
		c_from++;
		c_to++;
	}

	// If they have nothing in common. Return absolute path of to
	if (c_to == to.beg()) {
		return to;
	}

	int dir_count = 0;
	for (const char* c = c_from; c != from.end(); c++ /* <- LOL! */)  {
		if (*c == '\\' || *c == '/') dir_count++;
	}

	StringBuffer<256> res;
	int offset = 0;
	for (int i = 0; i < dir_count; i++) {
		offset += snprintf(res.buffer + offset, res.MAX_LENGTH, "../");
	}
	
	StringBuffer<256> to_buf = CString(c_to, to.end());
	snprintf(res.buffer + offset, res.MAX_LENGTH, "%s", to_buf.beg());

	return res;
}

StringBuffer<256> get_absolute_path(CString absolute_reference, CString relative_file) {
	CString abs_dir = get_directory(absolute_reference);
	if (relative_file.count < 3) return {};

	StringBuffer<256> res;
	// If relative path is really an absolute path, just return that
	if (relative_file[0] == '/' || relative_file[1] == ':') {
		res = relative_file;
		return res;
	}

	int dir_count = 0;
	for (const char* c = relative_file.beg(); c < relative_file.end(); c += 3) {
		if (c[0] == '.' && (c + 1) != relative_file.end() && c[1] == '.' && (c + 2) != relative_file.end() && (c[2] == '/' || c[2] == '\\')) dir_count++;
	}

	const char* c = abs_dir.end() - 1;
	while (c > abs_dir.beg() && dir_count > 0) {
		if (*c == '/' || *c == '\\') {
			if (--dir_count == 0) break;
		}
		c--;
	}
	if (dir_count > 0 || c == abs_dir.beg()) return res;

	CString base_dir(abs_dir.beg(), c + 1);
	res = base_dir;
	StringBuffer<128> file = get_file(relative_file);
	snprintf(res.buffer + base_dir.count, res.MAX_LENGTH - base_dir.count, "/%s", file.beg());

	return res;
}

void convert_backslashes(String str) {
	for (char* c = str.beg(); c != str.end(); c++) {
		if (*c == '\\') *c = '/';
	}
}

bool contains_whitespace(CString str) {
	for (const char* c = str.beg(); c != str.end(); c++) {
		if (isspace(*c)) return true;
	}
	return false;
}

DynamicArray<String> tokenize(String str, char delimiter) {
	DynamicArray<String> tokens;

	char* beg = str.beg();
	char* end = str.beg();

	while (end != str.end() && *end != '\0') {
		while (end != str.end() && *end != '\0' && *end != delimiter) end++;
		tokens.push_back(String(beg, end));
		beg = end;
		while (beg != str.end() && *end != '\0' && *beg == delimiter) beg++;
		end = beg;
	}

	return tokens;
}

DynamicArray<String> tokenize(String str, CString delimiter) {
	DynamicArray<String> tokens;

	char* beg = str.beg();
	char* end = str.beg();

	while (end != str.end() && *end != '\0') {
		while (end != str.end() && *end != '\0' && !char_in_string(*end, delimiter)) end++;
		tokens.push_back(String(beg, end));
		beg = end;
		while (beg != str.end() && *end != '\0' && char_in_string(*beg, delimiter)) beg++;
		end = beg;
	}

	return tokens;
}

DynamicArray<CString> ctokenize(CString str, char delimiter) {
	DynamicArray<CString> tokens;

	const char* beg = str.beg();
	const char* end = str.beg();

	while (end != str.end() && *end != '\0') {
		while (end != str.end() && *end != '\0' && *end != delimiter) end++;
		tokens.push_back(CString(beg, end));
		beg = end;
		while (beg != str.end() && *beg != '\0' && *beg == delimiter) beg++;
		end = beg;
	}

	return tokens;
}

DynamicArray<CString> ctokenize(CString str, CString delimiter) {
	DynamicArray<CString> tokens;

	const char* beg = str.beg();
	const char* end = str.beg();

	while (end != str.end() && *end != '\0') {
		while (end != str.end() && *end != '\0' && !char_in_string(*end, delimiter)) end++;
		tokens.push_back(CString(beg, end));
		beg = end;
		while (beg != str.end() && *end != '\0' && char_in_string(*beg, delimiter)) beg++;
		end = beg;
	}

	return tokens;
}