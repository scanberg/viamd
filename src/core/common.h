#pragma once

#include <ppk_assert.h>

#define ASSERT                PPK_ASSERT
#define ASSERT_WARNING        PPK_ASSERT_WARNING
#define ASSERT_DEBUG          PPK_ASSERT_DEBUG
#define ASSERT_ERROR          PPK_ASSERT_ERROR
#define ASSERT_FATAL          PPK_ASSERT_FATAL
#define ASSERT_CUSTOM         PPK_ASSERT_CUSTOM
#define ASSERT_USED           PPK_ASSERT_USED
#define ASSERT_USED_WARNING   PPK_ASSERT_USED_WARNING
#define ASSERT_USED_DEBUG     PPK_ASSERT_USED_DEBUG
#define ASSERT_USED_ERROR     PPK_ASSERT_USED_ERROR
#define ASSERT_USED_FATAL     PPK_ASSERT_USED_FATAL
#define ASSERT_USED_CUSTOM    PPK_ASSERT_USED_CUSTOM
#define STATIC_ASSERT		  PPK_STATIC_ASSERT

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

#ifdef USE_DEFAULT_MALLOC_AND_FREE
#define MALLOC(x) malloc(x)
#define REALLOC(x, y) realloc(x, y)
#define FREE(x) free(x)
#endif

#ifdef USE_DEFAULT_TMP_MALLOC_AND_FREE
#define TMP_MALLOC(x) alloc(x)
#define TMP_REALLOC(x, y) realloc(x, y)
#define TMP_FREE(x) free(x)
#endif