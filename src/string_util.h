#pragma once

#include <stdint.h>

#include <core/md_str.h>

// StrBuf, Bloated like a true c++ pig in all its 'glory'

namespace cexpr_helper {
    constexpr int64_t cexpr_strlen(const char* cstr) {
        int64_t len = 0;
        while (*cstr != '\0') {
            ++cstr;
            ++len;
        }
        return len;
    }

    constexpr int64_t cexpr_strnlen(const char* cstr, int64_t max_length) {
        int64_t len = 0;
        while (len < max_length && *cstr != '\0') {
            ++cstr;
            ++len;
        }
        return len;
    }

    constexpr int64_t cexpr_strncpy(char* dst, const char* src, int64_t max_length) {
        int64_t i = 0;
        while (i < max_length && src[i] != '\0') {
            dst[i] = src[i];
            ++i;
        }
        return i;
    }

    template <typename T>
    constexpr void cexpr_copy(T* dst, const T* src, int64_t count) {
        int64_t i = 0;
        while (i < count) {
            dst[i] = src[i];
            i++;
        }
    }
}

// A buf string for wrapping a char buf[N]
template <int64_t Size>
struct StrBuf {
    static_assert(Size > 1, "Size of StrBuf must be more than 1");
    char buf[Size] = {};

    constexpr StrBuf() noexcept {};

    template <int64_t N>
    constexpr StrBuf(const char (&cstr)[N]) noexcept  {
        constexpr auto length = N < Size ? N : Size - 1;
        cexpr_helper::cexpr_copy(buf, cstr, length);
        buf[length] = '\0';
    }

    constexpr StrBuf(const char* cstr) noexcept {
        int64_t len = cexpr_helper::cexpr_strnlen(cstr, Size - 1);
        cexpr_helper::cexpr_copy(buf, cstr, len);
        buf[len] = '\0';
    }

    constexpr StrBuf(const char* cstr, int64_t len) noexcept {
        ASSERT(len > 0);
        len = len < Size - 1 ? len : Size - 1;
        cexpr_helper::cexpr_copy(buf, cstr, len);
        buf[len] = '\0';
    }

    constexpr StrBuf(char c) noexcept {
        buf[0] = c;
        buf[1] = '\0';
    }

    constexpr StrBuf(const StrBuf& other) noexcept {
        cexpr_helper::cexpr_copy(buf, other.buf, Size);
        //memcpy(buf, other.buf, MaxSize);
        buf[Size - 1] = '\0';
    }

    template <int64_t N>
    constexpr StrBuf(const StrBuf<N>& other) noexcept {
        constexpr auto len = N < Size ? N : Size;
        cexpr_helper::cexpr_copy(buf, other.buf, len);
        buf[len - 1] = '\0';
    }

    constexpr StrBuf(StrBuf&& other) noexcept {
        cexpr_helper::cexpr_copy(buf, other.buf, Size);
        buf[Size - 1] = '\0';
    }

    template <int64_t N>
    constexpr StrBuf(StrBuf<N>&& other) noexcept {
        constexpr auto len = N < Size ? N : Size;
        cexpr_helper::cexpr_copy(buf, other.buf, len);
        buf[len - 1] = '\0';
    }

    StrBuf(str_t str);

    constexpr StrBuf& operator=(const StrBuf& other) noexcept {
        if (this != &other) {
            cexpr_helper::cexpr_copy(buf, other.buf, Size);
            buf[Size - 1] = '\0';
        }
        return *this;
    }

    template <int64_t N>
    constexpr StrBuf& operator=(const StrBuf<N>& other) noexcept {
        if (this != &other) {
            constexpr auto len = N < (Size - 1) ? N : (Size - 1);
            cexpr_helper::cexpr_copy(buf, other.buf, len);
            buf[len] = '\0';
        }
        return *this;
    }

    StrBuf& operator=(str_t str);

    constexpr StrBuf& operator=(char c) noexcept {
        buf[0] = c;
        buf[1] = '\0';
        return *this;
    }

    template <int64_t N>
    constexpr StrBuf& operator=(const char (&cstr)[N]) noexcept {
        if (&buf[0] != &cstr[0]) {
            constexpr auto len = N < (Size - 1) ? N : (Size - 1);
            cexpr_helper::cexpr_copy(buf, cstr, len);
            buf[len] = '\0';
        }
        return *this;
    }

    constexpr StrBuf& operator=(const char* cstr) noexcept {
        const int64_t len = cexpr_helper::cexpr_strnlen(cstr, Size - 1);
        cexpr_helper::cexpr_copy(buf, cstr, len);
        buf[len] = '\0';
        return *this;
    }

    StrBuf& operator+=(str_t txt);

    constexpr char operator[](int64_t i) const noexcept {
        ASSERT(0 <= i && i < Size);
        return buf[i];
    }

    constexpr char& operator[](int64_t i) noexcept {
        ASSERT(0 <= i && i < Size);
        return buf[i];
    }

    constexpr operator bool() const noexcept { return buf[0] != '\0'; }

    constexpr int64_t capacity() const noexcept { return Size; }
    constexpr int64_t length() const noexcept { return cexpr_helper::cexpr_strnlen(buf, Size); }

    constexpr char* cstr() const noexcept { return (char*)buf; }
    constexpr operator str_t() const noexcept { return {buf, length()}; }

    constexpr const char* begin() const noexcept { return buf; }
    constexpr const char* beg() const noexcept { return buf; }
    constexpr const char* end() const noexcept { return buf + Size; }

    constexpr char* begin() noexcept { return buf; }
    constexpr char* beg() noexcept { return buf; }
    constexpr char* end() noexcept { return buf + Size; }
};

template <int64_t N>
StrBuf<N>::StrBuf(str_t str) {
    // @NOTE: MAX_LENGTH - 1 here because we copy from cstring which excludes \0
    auto len = str.len < (N - 1) ? str.len : (N - 1);
    cexpr_helper::cexpr_copy(buf, str.ptr, len);
    buf[len] = '\0';
}

template <int64_t N>
StrBuf<N>& StrBuf<N>::operator=(str_t str) {
    auto len = str.len < (N - 1) ? str.len : (N - 1);
    cexpr_helper::cexpr_copy(buf, str.ptr, len);
    buf[len] = '\0';
    return *this;
}

template <int64_t N>
StrBuf<N>& StrBuf<N>::operator+=(str_t txt) {
    const int64_t offset = (int64_t)cexpr_helper::cexpr_strnlen((const char*)buf, N);
    const int64_t length = N - offset;
    if (length > 0) {
        cexpr_helper::cexpr_strncpy((char*)buf + offset, txt.ptr, txt.len < length ? txt.len : length);
    }
    buf[N - 1] = '\0';
    return *this;
}
