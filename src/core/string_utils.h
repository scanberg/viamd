#pragma once

#include "types.h"
#include "array.h"
#include "allocator.h"

#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

// A buffer string for wrapping a char buffer[N]
template<int64 Size>
struct StringBuffer {
	static constexpr int64 MAX_LENGTH = Size;
	STATIC_ASSERT(MAX_LENGTH > 1, "Size of StringBuffer must be more than 1");
	char buffer[MAX_LENGTH] = {};
	//int32 length = 0;

	StringBuffer() = default;

	template <int64 N>
	StringBuffer(const char(&cstr)[N]) {
		constexpr auto len = N < MAX_LENGTH ? N : MAX_LENGTH;
		strncpy(buffer, cstr, len);
		buffer[len - 1] = '\0';
	}

	StringBuffer(const char* cstr) {
		int64 len = (int64)strnlen(cstr, MAX_LENGTH);
		strncpy(buffer, cstr, len);
		buffer[len - 1] = '\0';
	}

	StringBuffer(char c) {
		buffer[0] = c;
		buffer[1] = '\0';
	}

	StringBuffer(const StringBuffer& other) {
		memcpy(buffer, other.buffer, MAX_LENGTH);
		buffer[MAX_LENGTH - 1] = '\0';
	}

	template <int64 N>
	StringBuffer(const StringBuffer<N>& other) {
		constexpr auto len = N < MAX_LENGTH ? N : MAX_LENGTH;
		memcpy(buffer, other.buffer, len);
		buffer[len - 1] = '\0';
	}

	StringBuffer(StringBuffer&& other) {
		memcpy(buffer, other.buffer, MAX_LENGTH);
		buffer[MAX_LENGTH - 1] = '\0';
	}

	template <int64 N>
	StringBuffer(StringBuffer<N>&& other) {
		constexpr auto len = N < MAX_LENGTH ? N : MAX_LENGTH;
		memcpy(buffer, other.buffer, len);
		buffer[len - 1] = '\0';
	}

	StringBuffer(const CString& cstr) {
		// @NOTE: MAX_LENGTH - 1 here because we copy from cstring which excludes \0
		auto len = cstr.count < MAX_LENGTH - 1 ? cstr.count : MAX_LENGTH - 1;
		strncpy(buffer, cstr.data, len);
		buffer[len] = '\0';
	}

	StringBuffer& operator =(const StringBuffer& other) {
		if (this != &other) {
			memcpy(buffer, other.buffer, MAX_LENGTH);
			buffer[MAX_LENGTH - 1] = '\0';
		}
		return *this;
	}

	template <int64 N>
	StringBuffer& operator =(const StringBuffer<N>& other) {
		if (this != &other) {
			auto len = other.count < MAX_LENGTH ? other.count : MAX_LENGTH;
			memcpy(buffer, other.buffer, len);
			buffer[len - 1] = '\0';
		}
		return *this;
	}

	StringBuffer& operator =(const CString& cstr) {
		auto len = cstr.count < MAX_LENGTH ? cstr.count : MAX_LENGTH;
		strncpy(buffer, cstr.data, len);
		buffer[len - 1] = '\0';
		return *this;
	}

	StringBuffer& operator =(char c) {
		buffer[0] = c;
		buffer[1] = '\0';
		return *this;
	}

	template<int64 N>
	StringBuffer& operator =(const char(&cstr)[N]) {
		if (buffer != cstr) {
			constexpr auto len = N < MAX_LENGTH ? N : MAX_LENGTH;
			strncpy(buffer, cstr, len);
			buffer[len - 1] = '\0';
		}
		return *this;
	}

    StringBuffer& operator =(const char* cstr) {
        int64 len = (int64)strnlen(cstr, MAX_LENGTH);
        strncpy(buffer, cstr, len);
        buffer[len - 1] = '\0';
    }

    char operator[](int64 i) const {
        ASSERT(i < MAX_LENGTH);
        return buffer[i];
    }

    char& operator[](int64 i) {
        ASSERT(i < MAX_LENGTH);
        return buffer[i];
    }

	operator String() { return String(buffer, MAX_LENGTH); }
	operator CString() const { return CString(buffer, MAX_LENGTH); }
	operator const char*() const { return buffer; }

	int64 size() const { return MAX_LENGTH; }

	const char* begin() const { return buffer; }
	const char* beg() const { return buffer; }
	const char* end() const { return buffer + MAX_LENGTH; }

	char* begin() { return buffer; }
	char* beg() { return buffer; }
	char* end() { return buffer + MAX_LENGTH; }
};

// Comparison of Strings
bool compare(CString str_a, CString str_b, bool ignore_case = false);
bool compare_n(CString str_a, CString str_b, int64 num_chars, bool ignore_case = false);

// Copy String
// Note: will zero terminate the dst String
void copy(String dst, CString src);

// Copy N first characters of src String
// Note: will zero terminate the dst String
void copy_n(String dst, CString src, int64 num_chars);

// Reads the next line from a CString
// Note: Line holds the extracted line, str gets modified
// Note: NO guarantee that the line will be zero terminated. Be careful with printf!
bool extract_line(CString& line, CString& str);

// Copies the next line from a CString
// Note: line holds the buffer which the line will be copied to, str gets modified
// Note: Guaranteed that the line will be zero terminated.
bool copy_line(String& line, CString& str);

// Wrappers around atof
float32 to_float32(CString str);
float64 to_float64(CString str);
inline float to_float(CString str) { return to_float32(str); }

// Wrappers around atoi
int32 to_int32(CString str);
int64 to_int64(CString str);
inline int to_int(CString str) { return to_int32(str); }

// Removes whitespace from begining and end of String
CString trim(CString str);
String trim(String str);

// Reads text file and copies into zero terminated String allocated using supplied allocator
String read_textfile(CString filename, Allocator& alloc = default_alloc);

// Returns directory part from url, ex: func("C:/folder/file.ext") should return "C:/folder/"
CString get_directory(CString url);

// Returns file from url including extension, ex: func("C:/folder/file.ext") should return "file.ext"
CString get_file(CString url);

// Returns file part of url excluding extension, ex: func("C:/folder/file.ext") should return "file"
CString get_file_without_extension(CString url);

// Returns file extension part of url, ex: func("C:/folder/file.ext") should return "ext"
CString get_file_extension(CString url);