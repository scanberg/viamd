#include "profiling.h"
#include <core/array.h>
#include <core/string_utils.h>
#include <core/hash.h>
#include <stdlib.h>

namespace profiling {

struct Entry {
	uint32 key;
	float32 ms;
	int32 parent_idx = -1;
};

static DynamicArray<StringBuffer<64>> names;
static DynamicArray<Entry> map;
static int32 curr_section = -1;

static Entry* lower_bound(Array<Entry>& data, uint32 key) {
	Entry* first = data.begin();
	Entry* last = data.end();
	size_t count = (size_t)(last - first);
	while (count > 0) {
		size_t count2 = count >> 1;
		Entry* mid = first + count2;
		if (mid->key < key) {
			first = ++mid;
			count -= count2 + 1;
		}
		else {
			count = count2;
		}
	}
	return first;
}

static void sort_map(Array<Entry>& map) {
	qsort(map.data, map.count, sizeof(Entry), [](const void* a, const void* b) -> int {
		if (((Entry*)a)->key < ((Entry*)b)->key) return -1;
		if (((Entry*)a)->key >((Entry*)b)->key) return +1;
		return 0;
	});
}

static Entry* get_or_insert_entry(DynamicArray<Entry>& map, uint32 key) {
	auto ptr = lower_bound(map, key);
	if (ptr == map.end() || ptr->key != key) {
		return map.insert(ptr, { key, 0 });
	}
	return ptr;
}

void profiling::initialize() {}

void profiling::shutdown() {}

void profiling::push_section(CString sec) {
	uint32 key = hash::crc32(sec);
	auto* e = get_or_insert_entry(map, key);
}

void profiling::pop_section() {}

void profiling::finish() {
	map.clear();
}

void profiling::print_log() {}

}  // namespace profiling
