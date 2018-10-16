#pragma once

#include <string.h>

#include "types.h"
#include "array_types.h"

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
    int64 length() const { return count; }

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

    int64 length() const { return count; }

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

    StringBuffer& operator+=(CString txt) {
        int64 offset = (int64)strnlen(buffer, MAX_LENGTH);
        int64 length = MAX_LENGTH - offset;
        if (length > 0) {
            strncpy(buffer + offset, txt, txt.count < length ? txt.count : length);
        }
        buffer[MAX_LENGTH - 1] = '\0';
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
