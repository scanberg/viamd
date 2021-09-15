#pragma once

#include <core/md_common.h> // assert
#include <core/md_intrinsics.h>

template <typename Key, typename Type, int Size>
struct Cache {
	struct Node {
		int prev = 0;
		int next = 0;
	};

	Cache() {
		clear();
	}

	void clear() {
		memset(this, 0, sizeof(*this));
		int N = Size - 1;
		for (int i = 0; i < Size; ++i) {
			list[i].prev = (i == 0) ? -1 : i - 1;
			list[i].next = (i == N) ? -1 : i + 1;
		}
		front = 0;
		back  = N;
	}

	int find(Key key) {
		for (int i = 0; i < Size; ++i) {
			if (keys[i] == key) return i;
		}
		return -1;
	}

	int pop_back() {
		int idx = back;
		int prev = list[back].prev;
		if (prev != -1) {
			list[prev].next = -1;
			back = prev;
		}
		else {
			back = -1;
		}
		return idx;
	}

	void pop(int idx) {
		int prev = list[idx].prev;
		int next = list[idx].next;
		if (prev != -1) list[prev].next = next;
		if (next != -1) list[next].prev = prev;
	}

	void push_front(int idx) {
		list[front].prev = idx;
		list[idx].next = front;
		list[idx].prev = -1;
		front = idx;
	}

	void put(Key key, const Type& item) {
		ASSERT(find(key) == -1);
		int idx = pop_back();
		push_front(idx);
		keys[idx] = key;
		data[idx] = item;
	}

	Type* reserve(Key key) {
		ASSERT(find(key) == -1);
		int idx = pop_back();
		push_front(idx);
		keys[idx] = key;
		return &data[idx];
	}

	Type* get(Key key) {
		int idx = find(key);
		if (idx != -1) {
			if (front != idx) {
				if (idx == back)
					pop_back();
				else
					pop(idx);
				push_front(idx);
			}
			return &data[idx];
		}
		return nullptr;
	}

	int front = 0;
	int back = 0;
	Type data[Size] = {0};
	Key  keys[Size] = {0};
	Node list[Size] = {0};
};

template <typename Key, typename Type>
struct LRU_Cache_4 {
	uint16_t ref_matrix = 0x8CE; // Bit matrix with upper triangular bits set
	Type data[4] = {0};
	Key  keys[4] = {0};

	void clear() {
		memset(this, 0, sizeof(*this));
		ref_matrix = 0x8CE;
	}

	int get_lru_idx() {
		// Find first zero nibble
		if ((ref_matrix & 0x000F) == 0) return 0;
		if ((ref_matrix & 0x00F0) == 0) return 1;
		if ((ref_matrix & 0x0F00) == 0) return 2;
		if ((ref_matrix & 0xF000) == 0) return 3;
		ASSERT(false && "Corrupted ref matrix!");
		return -1;
	}

	void mark_lru_idx(int i) {
		ref_matrix = ref_matrix | (0xF << (4 * i));
		ref_matrix = ref_matrix & ~(0x1111 << i);
	}

	int find(Key key) {
		for (int i = 0; i < 4; ++i) {
			if (keys[i] == key) return i;
		}
		return -1;
	}

	Type* get(Key key) {
		int idx = find(key);
		if (idx != -1) {
			mark_lru_idx(idx);
			return &data[idx];
		}
		return nullptr;
	}

	void put(Key key, const Type& item) {
		ASSERT(find(key) == -1);
		int idx = get_lru_idx();
		mark_lru_idx(idx);
		keys[idx] = key;
		data[idx] = item;
	}

	Type* reserve(Key key) {
		ASSERT(find(key) == -1);
		int idx = get_lru_idx();
		mark_lru_idx(idx);
		keys[idx] = key;
		return &data[idx];
	}
};

template <typename Key, typename Type>
struct LRU_Cache_8 {
	uint64_t ref_matrix = 0x0080c0e0f0f8fcfe; // Bit matrix with upper triangular bits set
	Type data[8] = {0};
	Key  keys[8] = {0};

	void clear() {
		memset(this, 0, sizeof(*this));
		ref_matrix = 0x0080c0e0f0f8fcfe;
	}

	int get_lru_idx() {
		int result = (int)find_first_zero_byte64(ref_matrix);
		ASSERT(result < 8);
		return result;
	}

	void mark_lru_idx(int i) {
		ref_matrix = ref_matrix | (0xFFLLU << (8 * i));
		ref_matrix = ref_matrix & ~(0x0101010101010101LLU << i);
	}

	int find(Key key) {
		for (int i = 0; i < 8; ++i) {
			if (keys[i] == key) return i;
		}
		return -1;
	}

	Type* get(Key key) {
		int idx = find(key);
		if (idx != -1) {
			mark_lru_idx(idx);
			return &data[idx];
		}
		return nullptr;
	}

	void set(Key key, const Type val) {
		ASSERT(find(key) == -1);
		int idx = get_lru_idx();
		mark_lru_idx(idx);
		keys[idx] = key;
		data[idx] = val;
	}

	Type* reserve(Key key) {
		ASSERT(find(key) == -1);
		int idx = get_lru_idx();
		mark_lru_idx(idx);
		keys[idx] = key;
		return &data[idx];
	}
};