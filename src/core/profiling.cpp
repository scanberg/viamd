#include "profiling.h"
#include <core/array.h>
#include <core/string_utils.h>
#include <core/hash.h>
#include <stdlib.h>

namespace profiling {

typedef uint32 Key;

struct Entry {
	Key key;
	float32 ms;
	int32 parent_idx = -1;
	StringBuffer<64> name;
};

static DynamicArray<Entry> entries;
static int32 curr_section = -1;

inline Entry* find_entry(Array<Entry> ent, Key key) {
	for (auto& e : ent) {
		if (key == e.key) return &e;
	}
	return nullptr;
}

void initialize() {}

void shutdown() {}

void push_section(CString sec) {
	uint32 key = hash::crc32(sec);
	auto* e = find_entry(entries, key);
	if (!e) {
		e = &entries.push_back({ key, 0, -1, sec });
	}
}

void pop_section() {}

void finish() {
	entries.clear();
}

void print_log() {}

}  // namespace profiling
