#pragma once

#include "types.h"
#include "array.h"
#include "allocator.h"

bool compare(CString str_a, CString str_b, bool ignore_case = false);
bool compare_n(CString str_a, CString str_b, int64 N, bool ignore_case = false);

void copy(String dst, CString src);
bool get_line(String& line, CString& str);

float to_float(CString str);
int to_int(CString str);

CString trim(CString str);
String trim(String str);

String read_textfile(CString filename, Allocator& alloc = default_alloc);

CString get_directory(CString url);
CString get_file(CString url);
CString get_file_without_extension(CString url);
CString get_file_extension(CString url);
