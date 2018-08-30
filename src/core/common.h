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

#ifndef NDEBUG
inline void _assert(const char* file, const char* func, int line, bool cond, const char* fmt, ...) {
    if (!cond) {
        va_list ap;
        va_start(ap, fmt);
        fprintf(stderr, "ASSERTION FAILED!\n%s:%s:%i\n", file, func, line);
        vfprintf(stderr, fmt, ap);
        va_end(ap);
        assert(false);
    }
}
inline void _assert(const char* file, const char* func, int line, bool cond) { _assert(file, func, line, cond, ""); }

#define ASSERT(...) _assert(__FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)
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

// implementation of 'defer' in c++.
// from here https://pastebin.com/suTkpYp4

#define CONCAT_INTERNAL(x, y) x##y
#define CONCAT(x, y) CONCAT_INTERNAL(x, y)

struct ExitScopeHelp {
    template <typename T>
    struct ExitScope {
        T lambda;
        ExitScope(T lambda) : lambda(lambda) {}
        ~ExitScope() { lambda(); }
        ExitScope& operator=(const ExitScope&) = delete;
    };

    template <typename T>
    ExitScope<T> operator+(T t) {
        return t;
    }
};

#define defer const auto& CONCAT(defer__, __LINE__) = ExitScopeHelp() + [&]()

// //sample code:
// {
// 	defer
// 	{
// 		printf("three\n");
// 	};
// 	printf("one\n");
// 	defer
// 	{
// 		printf("two\n");
// 	};
// }
// prints
// one
// two
// three
