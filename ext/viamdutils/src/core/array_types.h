#pragma once

#include "types.h"
#include "common.h"

#include <string.h>
#include <type_traits>

/*
  This is an 'array-view' which exposes access to some data which is not owned by the array itself.
  Nothing will be allocated or freed by the constructors and destructors of this object.
*/

template <typename T>
struct Array {
    Array() = default;
    template <size_t N>
    Array(T (&c_arr)[N]) : data(c_arr), count(N) {}
    Array(T* _data, int64 _count) : data(_data), count(_count) {}
    Array(T* _data_beg, T* _data_end) : data(_data_beg), count(_data_end - _data_beg) {}

    Array<T> sub_array(int64 _offset, int64 _count = -1) {
        ASSERT(0 <= _offset);
        ASSERT(_count >= -1);
        if (_count == -1) {
            _count = (this->count - _offset) > 0 ? (this->count - _offset) : 0;
        }
        ASSERT((_offset + _count) <= this->count);
        return {data + _offset, _count};
    }

    Array<const T> sub_array(int64 _offset, int64 _count = -1) const {
        ASSERT(0 <= _offset);
        ASSERT(_count >= -1);
        if (_count == -1) {
            _count = (this->count - _offset) > 0 ? (this->count - _offset) : 0;
        }
        ASSERT((_offset + _count) <= this->count);
        return {data + _offset, _count};
    }

    const T* begin() const { return data; }
    const T* beg() const { return data; }
    const T* end() const { return data + count; }

    T* begin() { return data; }
    T* beg() { return data; }
    T* end() { return data + count; }

    const T& front() const { return data[0]; }
    const T& back() const { return data[count - 1]; }
    T& front() { return data[0]; }
    T& back() { return data[count - 1]; }

    int64 size() const { return count; }
    int64 size_in_bytes() const { return count * sizeof(T); }

    operator bool() const { return data != nullptr && count > 0; }
    const T& operator[](int64 i) const { return data[i]; }
    T& operator[](int64 i) { return data[i]; }

    operator Array<const T>() const { return {data, count}; }

    T* data;
    int64 count;
};

// Light-weight std::vector alternative
// @WARNING: THIS IS NOT A STRAIGHT FORWARD REPLACEMENT TO STD::VECTOR AS CONSTRUCTORS AND DESTRUCTORS ARE NEVER CALLED.
template <typename T>
struct DynamicArray : Array<T> {
    static_assert(std::is_trivially_destructible<T>::value, "DynamicArray only supports trivially destructable data types");

    static constexpr int64 INIT_CAPACITY = 8;
    DynamicArray() : capacity(INIT_CAPACITY) {
        this->data = (T*)CALLOC(capacity, sizeof(T));
        this->count = 0;
    }

    DynamicArray(int64 count) {
        count < 1 ? count = 1 : count;
        capacity = count;
        this->data = (T*)CALLOC(capacity, sizeof(T));
        this->count = capacity;
    }

    DynamicArray(int64 count, T value) : capacity(count > INIT_CAPACITY ? count : INIT_CAPACITY) {
        this->data = (T*)MALLOC(capacity * sizeof(T));
        this->count = count;
        for (int i = 0; i < count; i++) {
            this->data[i] = value;
        }
    }

    DynamicArray(const T* first, const T* last) noexcept {
        this->capacity = last - first;
        this->count = capacity;
        if (this->count > 0) {
            this->data = (T*)MALLOC(capacity * sizeof(T));
            memcpy(this->data, first, this->count * sizeof(T));
        }
    }

    DynamicArray(const Array<const T>& clone_source) : capacity(clone_source.count) {
        this->count = capacity;
        if (this->count > 0) {
            this->data = (T*)MALLOC(capacity * sizeof(T));
            memcpy(this->data, clone_source.data, this->count * sizeof(T));
        }
    }

    DynamicArray(const DynamicArray& other) : capacity(other.count) {
        this->count = capacity;
        if (this->count > 0) {
            this->data = (T*)MALLOC(capacity * sizeof(T));
            memcpy(this->data, other.data, this->count * sizeof(T));
        }
    }

    DynamicArray(DynamicArray&& other) : capacity(other.capacity) {
        this->data = other.data;
        other.data = nullptr;
        other.capacity = 0;
        this->count = other.count;
        other.count = 0;
    }

    ~DynamicArray() {
        if (this->data) {
            FREE(this->data);
        }
        this->count = 0;
    }

    DynamicArray& operator=(const Array<const T>& other) {
        // Is this ok? It probably is since we're only comparing memory adresses...
        if (&other != (const Array<const T>*)this) {
            if (other.count > capacity) {
                reserve(other.count);
            }
            this->count = other.count;
            memcpy(this->data, other.data, this->count * sizeof(T));
        }
        return *this;
    }

    DynamicArray& operator=(const DynamicArray& other) {
        if (&other != this) {
            if (other.count > capacity) {
                reserve(other.count);
            }
            this->count = other.count;
            memcpy(this->data, other.data, this->count * sizeof(T));
        }
        return *this;
    }

    DynamicArray& operator=(DynamicArray&& other) {
        // @NOTE: Is this check needed?
        if (&other != this) {
            if (this->data) {
                FREE(this->data);
            }
            capacity = other.capacity;
            other.capacity = 0;
            this->data = other.data;
            other.data = nullptr;
            this->count = other.count;
            other.count = 0;
        }
        return *this;
    }

    inline int64 _grow_capacity(int64 sz) const {
        int64 new_capacity = capacity ? (capacity + capacity / 2) : INIT_CAPACITY;
        return new_capacity > sz ? new_capacity : sz;
    }

    void append(Array<const T> arr) {
        if (this->count + arr.count >= capacity) {
            reserve(_grow_capacity(this->count + arr.count));
        }
        memcpy(this->end(), arr.data, arr.count * sizeof(T));
        this->count += arr.count;
    }

    void append(Array<T> arr) {
        if (this->count + arr.count >= capacity) {
            reserve(_grow_capacity(this->count + arr.count));
        }
        memcpy(this->end(), arr.data, arr.count * sizeof(T));
        this->count += arr.count;
    }

    T& push_back(const T& item) {
        if (this->count == capacity) {
            reserve(_grow_capacity(this->count + 1));
        }
        this->data[this->count] = item;
        this->count++;
        return this->back();
    }

    T pop_back() {
        ASSERT(this->count > 0);
        this->count--;
        return this->data[this->count];
    }

    void reserve(int64 new_capacity) {
        if (new_capacity < capacity) return;
        T* new_data = (T*)CALLOC(new_capacity, sizeof(T));
        if (this->data) {
            memcpy(new_data, this->data, this->count * sizeof(T));
        }
        FREE(this->data);
        this->data = new_data;
        capacity = new_capacity;
    }

    // Resizes the array to a new size and zeros eventual new slots
    void resize(int64 new_count) {
        if (new_count == this->count) {
            return;
        } else if (new_count < this->count) {
            this->count = new_count;
        } else {
            if (capacity < new_count) {
                reserve(_grow_capacity(new_count));
                // memset(this->data + this->count, 0, (new_count - this->count) * sizeof(T));
            }
            this->count = new_count;
        }
    }

    T* insert(T* it, const T& v) {
        ASSERT(this->beg() <= it && it <= this->end());
        const ptrdiff_t off = it - this->beg();
        if (this->count == capacity) reserve(_grow_capacity(this->count + 1));
        if (off < (int64)this->count) memmove(this->beg() + off + 1, this->beg() + off, ((size_t)this->count - (size_t)off) * sizeof(T));
        this->data[off] = v;
        this->count++;
        return this->beg() + off;
    }

    void remove(T* it, int64 num_items = 1) {
        ASSERT(this->beg() <= it && it < this->end());
        ASSERT(it + num_items <= this->end());
        auto dst = it;
        auto src = it + num_items;
        memmove(dst, src, (this->end() - src) * sizeof(T));
        this->count--;
    }

    void swap_back_and_pop(T* it) {
        ASSERT(this->beg() <= it && it < this->end());
        *it = this->back();
        pop_back();
    }

    void clear() { this->count = 0; }

    void set_mem_to_zero() { memset(this->data, 0, this->count * sizeof(T)); }

private:
    int64 capacity;
};

/*
// @TODO: Implement this (like std::array ish)
// @TODO: Make StringBuffer to use this

template <typename T, int64 Length>
struct StaticArray {
	T data[Length];
};
*/

template <typename T>
Array<T> allocate_array(int64 num_elements) {
    if (num_elements == 0) return {};
    return {(T*)MALLOC(num_elements * sizeof(T)), num_elements};
}

template <typename T>
void free_array(Array<T>* arr) {
    ASSERT(arr);
    if (arr->data) {
        FREE(arr->data);
    }
    arr->data = nullptr;
    arr->count = 0;
}

template <typename T>
void zero_array(Array<T>* arr) {
	ASSERT(arr);
	memset(arr->data, 0, arr->count * sizeof(T));
}

// count = 4
// i = 2
// [0,1,2,3]
template <typename T>
void remove_array_element(Array<T>* arr, int i) {
    ASSERT(arr);
    ASSERT(i < arr->count);
    memmove(arr->data + i, arr->data + (i + 1), arr->count - (i + 1));
}