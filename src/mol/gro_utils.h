#pragma once

#include <mol/gro.h>
#include <core/common.h>
#include <core/allocator.h>
#include <core/array.h>

struct GroResult {
	bool success = false;
	GroStructure gro = {};

	operator bool() { return success; }
};

GroResult load_gro_from_file(const char* filename, Allocator& alloc = default_alloc);
GroResult parse_gro_from_string(CString string, Allocator& alloc = default_alloc);