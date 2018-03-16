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
	int64 len = MIN(str_a.count, str_b.count);
	return internal_compare(str_a, str_b, len, ignore_case);
}

bool compare_n(CString str_a, CString str_b, int64 num_chars, bool ignore_case) {
	int64 len = MIN(str_a.count, MIN(str_b.count, num_chars));
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

float32 to_float32(CString str) {
	// Make sure that the string passed into atof is zero-terminated
	char buffer[32] = {};
	memcpy(buffer, str.data, MIN(31, str.count));
	return (float)atof(buffer);
}

float64 to_float64(CString str) {
	// Make sure that the string passed into atof is zero-terminated
	char buffer[32] = {};
	memcpy(buffer, str.data, MIN(31, str.count));
	return atof(buffer);
}

int32 to_int32(CString str) {
	// Make sure that the string passed into atoi is zero-terminated
	char buffer[32] = {};
	memcpy(buffer, str.data, MIN(31, str.count));
	return atoi(buffer);
}

int64 to_int64(CString str) {
	// Make sure that the string passed into atoi is zero-terminated
	char buffer[32] = {};
	memcpy(buffer, str.data, MIN(31, str.count));
	return atoll(buffer);
}

CString trim(CString str) {
	const char* beg = str.data;
	const char* end = str.data + str.count;

	while (beg < end && isspace(*beg)) ++beg;
	while (end > beg && isspace(*(end-1))) --end;

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
    
    while (beg != url.begin() && *beg != '.') beg--;
    
    if (beg == url.begin()) {
        return CString();
    }
    
    return CString(beg, end - beg);
}

DynamicArray<String> tokenize(String str, char delimiter) {
	DynamicArray<String> tokens;

	char* beg = str.beg();
	char* end = str.beg();

	while (end != str.end()) {
		while (end != str.end() && *end != delimiter) end++;
		tokens.push_back(String(beg, end));
		beg = end;
		while (beg != str.end() && *beg == delimiter) beg++;
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
		while (beg != str.end() && *end != '\0' && *beg == delimiter) beg++;
		end = beg;
	}

	return tokens;
}