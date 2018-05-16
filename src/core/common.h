#pragma once

/*
#include <ppk_assert.h>

#define ASSERT PPK_ASSERT
#define ASSERT_WARNING PPK_ASSERT_WARNING
#define ASSERT_DEBUG PPK_ASSERT_DEBUG
#define ASSERT_ERROR PPK_ASSERT_ERROR
#define ASSERT_FATAL PPK_ASSERT_FATAL
#define ASSERT_CUSTOM PPK_ASSERT_CUSTOM
#define ASSERT_USED PPK_ASSERT_USED
#define ASSERT_USED_WARNING PPK_ASSERT_USED_WARNING
#define ASSERT_USED_DEBUG PPK_ASSERT_USED_DEBUG
#define ASSERT_USED_ERROR PPK_ASSERT_USED_ERROR
#define ASSERT_USED_FATAL PPK_ASSERT_USED_FATAL
#define ASSERT_USED_CUSTOM PPK_ASSERT_USED_CUSTOM
#define STATIC_ASSERT PPK_STATIC_ASSERT
*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>

#if !defined(NDEBUG)
inline void _assert(bool cond, const char* file, const char* func, int line, const char* fmt, ...) {
    if (!cond) {
        va_list ap;
        va_start(ap, fmt);
        fprintf(stderr, "ASSERTION FAILED!\n%s:%s:%i\n", file, func, line);
        vfprintf(stderr, fmt, ap);
        va_end(ap);
        assert(false);
    }
}
inline void _assert(bool cond, const char* file, const char* func, int line) { _assert(cond, file, func, line, ""); }

#define ASSERT(cond, ...) _assert(cond, __FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)
#else
#define ASSERT(...) \
    {}
#endif
#define STATIC_ASSERT(cond, msg) static_assert(cond, msg)

#define KILOBYTES(x) (x << 10)
#define MEGABYTES(x) (KILOBYTES(x) << 10)
#define GIGABYTES(x) (MEGABYTES(x) << 10)

#define BIT(x) (1 << x)
#define UNUSED(x) (void)(x)

// @TODO: This has to be changed so that the default version is only used when neither malloc or free is defined
#ifndef MALLOC
#define USE_DEFAULT_MALLOC_AND_FREE
#endif

#ifndef TMP_MALLOC
#define USE_DEFAULT_TMP_MALLOC_AND_FREE
#endif

#if defined USE_DEFAULT_MALLOC_AND_FREE || defined USE_DEFAULT_TMP_MALLOC_AND_FREE
#include <stdlib.h>
#endif

#ifdef USE_DEFAULT_MALLOC_AND_FREE
#define MALLOC(x) malloc(x)
#define REALLOC(x, y) realloc(x, y)
#define CALLOC(x, y) calloc(x, y)
#define FREE(x) free(x)
#endif

#ifdef USE_DEFAULT_TMP_MALLOC_AND_FREE
#define TMP_MALLOC(x) malloc(x)
#define TMP_REALLOC(x, y) realloc(x, y)
#define TMP_CALLOC(x, y) calloc(x, y)
#define TMP_FREE(x) free(x)
#endif
