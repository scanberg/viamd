#pragma once

#include "types.h"
#include "array.h"

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

struct CString : Array<const char> {
    CString() = default;

    CString(const char* cstr, int64 length = -1) {
        data = cstr;
        if (length == -1) length = strlen(cstr);
        count = length;
    }

    CString(const char* beg, const char* end) {
        data = beg;
        count = end - beg;
        ASSERT(count >= 0);
    }

    template <int64 length>
    CString(const char (&cstr)[length]) {
        data = cstr;
        count = length;
    }

    CString substr(int64 _offset, int64 _count = -1) {
        auto array = sub_array(_offset, _count);
        return {array.data, array.count};
    }

    const char* cstr() const { return data; }

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
    String(char (&cstr)[length]) {
        data = cstr;
        count = length;
    }

    String substr(int64 _offset, int64 _count = -1) {
        auto array = sub_array(_offset, _count);
        return {array.data, array.count};
    }

    operator CString() { return CString(data, count); }
    operator char*() { return data; }
    operator bool() { return (data != 0 && count != 0); }
};

// A buffer string for wrapping a char buffer[N]
template <int64 Size>
struct StringBuffer {
    static constexpr int64 MAX_LENGTH = Size;
    STATIC_ASSERT(MAX_LENGTH > 1, "Size of StringBuffer must be more than 1");
    char buffer[MAX_LENGTH] = {};
    // int32 length = 0;

    StringBuffer() = default;

    template <int64 N>
    StringBuffer(const char (&cstr)[N]) {
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

    StringBuffer& operator=(const StringBuffer& other) {
        if (this != &other) {
            memcpy(buffer, other.buffer, MAX_LENGTH);
            buffer[MAX_LENGTH - 1] = '\0';
        }
        return *this;
    }

    template <int64 N>
    StringBuffer& operator=(const StringBuffer<N>& other) {
        if (this != &other) {
            auto len = other.count < MAX_LENGTH ? other.count : MAX_LENGTH;
            memcpy(buffer, other.buffer, len);
            buffer[len - 1] = '\0';
        }
        return *this;
    }

    StringBuffer& operator=(const CString& cstr) {
        auto len = cstr.count < (MAX_LENGTH - 1) ? cstr.count : (MAX_LENGTH - 1);
        strncpy(buffer, cstr.data, len);
        buffer[len] = '\0';
        return *this;
    }

    StringBuffer& operator=(char c) {
        buffer[0] = c;
        buffer[1] = '\0';
        return *this;
    }

    template <int64 N>
    StringBuffer& operator=(const char (&cstr)[N]) {
        if (buffer != cstr) {
            constexpr auto len = N < MAX_LENGTH ? N : MAX_LENGTH;
            strncpy(buffer, cstr, len);
            buffer[len - 1] = '\0';
        }
        return *this;
    }

    StringBuffer& operator=(const char* cstr) {
        int64 len = (int64)strnlen(cstr, MAX_LENGTH);
        strncpy(buffer, cstr, len);
        buffer[len - 1] = '\0';
        return *this;
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
    operator bool() const { return buffer[0] != '\0'; }

    int64 size() const { return MAX_LENGTH; }

    const char* cstr() const { return buffer; }
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

String allocate_string(CString str);
String allocate_string(int32 length);
void free_string(String* str);

template <typename T>
struct ConversionResult {
    T value;
    bool success;

    operator T() { return value; }
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

// Creates a relative path from the absolute path specified in the first file to the absolute path of the second file
StringBuffer<256> get_relative_path(CString absolute_from, CString absolute_to);

// Creates an absolute path using an absolute reference path specified in the first file and the relative path specified in the second
StringBuffer<256> get_absolute_path(CString absolute_reference, CString relative_file);

// Converts windows backslashes '\\' to forward slashes '/'
void convert_backslashes(String str);

bool contains_whitespace(CString str);

bool balanced_parentheses(CString str);

CString extract_parentheses(CString str);
CString extract_parentheses_contents(CString str);

const char* find_character(CString str, char c);
bool contains_character(CString str, char c);

// Tokenizes a string into shorter strings based on some delimiter
DynamicArray<String> tokenize(String str, char delimiter = ' ');
DynamicArray<String> tokenize(String str, CString delimiter);
DynamicArray<CString> ctokenize(CString str, char delimiter = ' ');
DynamicArray<CString> ctokenize(CString str, CString delimiter);

// Range extraction functionality
bool is_range(CString arg);
bool extract_range(int* first, int* last, CString arg);
bool extract_ranges(DynamicArray<ivec2>* ranges, Array<const CString> args);

// Temporary string object with managed memory
struct TmpString : CString {
    TmpString(CString str) {
        String tmp = allocate_string(str);
        this->data = tmp.data;
        this->count = tmp.count;
    }
    TmpString(const TmpString& other) = delete;
    TmpString(TmpString&& other) {
        this->data = other.data;
        this->count = other.count;
    }
    ~TmpString() { FREE((void*)this->data); }
};

// This is a hack to generate a zero terminated string from a CString object
// Returns an object with a temporary allocated string which is freed upon its destruction
inline TmpString make_tmp_str(CString str) { return TmpString(str); }

// TODO: Possibly implement a good templated print function in the style of printf as discussed here
// https://stackoverflow.com/questions/17671772/c11-variadic-printf-performance