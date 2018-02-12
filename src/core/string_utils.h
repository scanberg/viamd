#pragma once

#include "types.h"
#include "array.h"
#include "allocator.h"

// Comparison of Strings
bool compare(CString str_a, CString str_b, bool ignore_case = false);
bool compare_n(CString str_a, CString str_b, int64 N, bool ignore_case = false);

// Copy String data
void copy(String dst, CString src);

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
