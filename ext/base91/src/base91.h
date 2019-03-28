/*
 * Copyright (c) 2000-2006 Joachim Henke
 *
 * For conditions of distribution and use, see copyright notice in base91.c
 */

#ifndef BASE91_H
#define BASE91_H 1

#include <stddef.h>

struct basE91 {
	unsigned long queue;
	unsigned int nbits;
	int val;
};

void basE91_init(struct basE91 *);

size_t basE91_encode(struct basE91 *, const void *, size_t, void *);

size_t basE91_encode_end(struct basE91 *, void *);

size_t basE91_decode(struct basE91 *, const void *, size_t, void *);

size_t basE91_decode_end(struct basE91 *, void *);

#endif	/* base91.h */
