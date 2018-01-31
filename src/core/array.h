#pragma once

#include <core/types.h>
#include <core/common.h>
#include <core/allocator.h>

template <typename T>
struct Array {
    Array(T* _data = 0, int64 _count = 0) : data(_data), count(_count) {}

    template <size_t N>
    Array(const T (&c_arr)[N]) : data(c_arr), count(N) {}

    Array SubArray(int64 _offset, int64 _count = -1) {
        ASSERT(0 <= _offset);
		ASSERT(_count > 0 || _count == -1);
        if (_count == -1) {
            _count = (this->count - _offset) > 0 ? (this->count - _offset) : 0;
        }
		ASSERT((_offset + _count) <= this->count);
        return {data + _offset, _count};
    }

    T* begin() { return data; }
    T* beg() { return data; }
    T* end() { return data + count; }

    operator bool() const { return data != nullptr && count > 0; }
    const T& operator[](int64 i) const { return data[i]; }
    T& operator[](int64 i) { return data[i]; }

    T* data;
    int64 count;
};

template <typename T>
struct DynamicArray : Array<T> {
    DynamicArray(Allocator& alloc = default_alloc) : capacity(32), allocator(alloc) {
        this->data = (T*)allocator.Alloc(capacity * sizeof(T));
        this->count = 0;
    }

    DynamicArray(const Array<T>& clone_source, Allocator& alloc = default_alloc) : allocator(alloc) {
        capacity = clone_source.count;
        this->count = capacity;
        if (this->count > 0) {
            this->data = (T*)allocator.Alloc(capacity * sizeof(T));
            memcpy(this->data, clone_source.data, this->count * sizeof(T));
        }
    }

    ~DynamicArray() {
        if (this->data) {
            allocator.Free(this->data);
        }
    }

    void PushBack(const T& item) {
        if (this->count >= capacity) {
            // GROW
            Reserve(capacity * 2);
        }
        this->data[this->count] = item;
        this->count++;
    }

    void Reserve(int64 new_capacity) {
        if (new_capacity < capacity) return;
        T* new_data = (T*)allocator.Alloc(new_capacity * sizeof(T));
        if (this->data) {
            memcpy(new_data, this->data, this->count * sizeof(T));
            allocator.Free(this->data);
        }
        this->data = new_data;
        capacity = new_capacity;
    }

    // Resizes the array to a new size and zeros eventual new slots
    void Resize(int64 new_count) {
        ASSERT(new_count > 0);
        if (new_count == this->count)
            return;
        else if (new_count < this->count) {
            this->count = new_count;
            return;
        } else {
            if (capacity < new_count) {
                Reserve(new_count);
                this->count = new_count;
            }
        }
    }

    int64 capacity;
    Allocator& allocator;
};

struct CString : Array<const char> {
    CString() {
        data = 0;
        count = 0;
    }

    CString(const char* cstr, int64 length) {
        data = cstr;
        count = length;
    }

    template <int64 length>
    CString(const char (&cstr)[length]) {
        data = cstr;
        count = length;
    }

    CString substr(int64 _offset, int64 _count = -1) {
        auto array = SubArray(_offset, _count);
        return {array.data, array.count};
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

    template <int64 length>
    String(char (&cstr)[length]) {
        data = cstr;
        count = length;
    }

    String substr(int64 _offset, int64 _count = -1) {
        auto array = SubArray(_offset, _count);
        return {array.data, array.count};
    }

    operator CString() { return CString(data, count); }
    operator const char*() { return data; }
    operator bool() { return (data != 0 && count != 0); }
};

template <typename T>
Array<T> AllocateArray(int64 count, Allocator& alloc = default_alloc) noexcept {
    ASSERT(count > 0);
    return {(T*)alloc.Alloc(sizeof(T) * count), count};
}

template <typename T>
Array<T> AllocateArrayAndZero(int64 count, Allocator& alloc = default_alloc) noexcept {
    Array<T> array = AllocateArray<T>(count, alloc);
    memset(array.data, 0, array.count * sizeof(T));
    return array;
}
