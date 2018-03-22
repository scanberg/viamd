#pragma once

#include "types.h"
#include "array.h"

#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

struct CString : Array<const char> {
	CString() {
		data = 0;
		count = 0;
	}

	CString(const char* cstr, int64 length = -1) {
		data = cstr;
		if (length == -1)
			length = strlen(cstr);
		count = length;
	}

	CString(const char* beg, const char* end) {
		data = beg;
		count = end - beg;
		ASSERT(count >= 0);
	}

	template <int64 length>
	CString(const char(&cstr)[length]) {
		data = cstr;
		count = length;
	}

	CString substr(int64 _offset, int64 _count = -1) {
		auto array = sub_array(_offset, _count);
		return { array.data, array.count };
	}

	operator const char*() { return data; }
	operator bool() { return (data != 0 && count != 0); }
};

struct String : Array<char> {
	String() {
		data = 0;
		count = 0;
	}

	String(const String& other) : Array(other.data, other.count) {}

	String(char* cstr, int64 length) {
		data = cstr;
		count = length;
	}

	String(char* beg, char* end) {
		data = beg;
		count = end - beg;
		ASSERT(count >= 0);
	}

	template <int64 length>
	String(char(&cstr)[length]) {
		data = cstr;
		count = length;
	}

	String substr(int64 _offset, int64 _count = -1) {
		auto array = sub_array(_offset, _count);
		return { array.data, array.count };
	}

	operator CString() { return CString(data, count); }
	operator char*() { return data; }
	operator bool() { return (data != 0 && count != 0); }
};

struct DynamicString : DynamicArray<char> {
	DynamicString() = default;

	DynamicString(CString other) : DynamicArray(other) {}

	DynamicString(DynamicString&& other) : DynamicArray(std::forward<DynamicArray<char>>(other)) {};

	DynamicString(char* cstr, int64 length) {
		this->resize(length + 1);
		strncpy(data, cstr, length);
	}

	template <int64 length>
	DynamicString(char(&cstr)[length]) {
		this->resize(length + 1);
		strncpy(data, cstr, length);
	}

	String substr(int64 _offset, int64 _count = -1) {
		auto arr = sub_array(_offset, _count);
		return { arr.data, arr.count };
	}

	operator String() { return String(data, count); }
	operator CString() { return CString(data, count); }
	operator char*() { return data; }
	operator bool() { return (data != 0 && count != 0); }

	DynamicString& operator += (CString other) {
		this->resize(count + other.count + 1);
		strncpy(end(), other.data, other.count);
		return *this;
	}
};

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
		buffer[len] = '\0';
	}

	StringBuffer(const char* cstr) {
		int64 len = (int64)strnlen(cstr, MAX_LENGTH);
		strncpy(buffer, cstr, len);
		buffer[len] = '\0';
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
	operator CString() const { return CString(buffer, strnlen(buffer, MAX_LENGTH)); }
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

template<typename T>
struct ConversionResult {
	T value;
	bool success;

	operator T() {
		return value;
	}
};

// Wrappers around strtof
ConversionResult<float32> to_float32(CString str);
ConversionResult<float64> to_float64(CString str);
inline ConversionResult<float> to_float(CString str) { return to_float32(str); }

// Wrappers around strtol
ConversionResult<int32> to_int32(CString str);
ConversionResult<int64> to_int64(CString str);
inline ConversionResult<int> to_int(CString str) { return to_int32(str); }

// Removes whitespace from begining and end of String
CString trim(CString str);
String trim(String str);

// Reads text file and copies into allocated zero terminated String
String allocate_and_read_textfile(CString filename);

// Returns directory part from url, ex: func("C:/folder/file.ext") should return "C:/folder/"
CString get_directory(CString url);

// Returns file from url including extension, ex: func("C:/folder/file.ext") should return "file.ext"
CString get_file(CString url);

// Returns file part of url excluding extension, ex: func("C:/folder/file.ext") should return "file"
CString get_file_without_extension(CString url);

// Returns file extension part of url, ex: func("C:/folder/file.ext") should return "ext"
CString get_file_extension(CString url);

// Tokenizes a string into shorter strings based on some delimiter 
// @TODO: Implement a CString instead of a char
DynamicArray<String> tokenize(String str, char delimiter = ' ');
DynamicArray<String> tokenize(String str, CString delimiter);
DynamicArray<CString> ctokenize(CString str, char delimiter = ' ');
DynamicArray<CString> ctokenize(CString str, CString delimiter);

