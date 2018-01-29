#pragma once

#include "types.h"
#include <stdlib.h>
#include <assert.h>

// TODO: Support aligned allocations
// TODO: Encode size in stack-allocator and assert comparison with prev on free

struct Allocator {
    virtual void* Alloc(int64) = 0;
    virtual void Free(void*) = 0;
};

struct StackAllocator : Allocator {
    uint8* base_ptr;
    uint8* curr_ptr;
    uint8* prev_ptr;
    int64 max_size;

    StackAllocator(void* memory, int64 size)
        : base_ptr((uint8*)memory), curr_ptr(base_ptr), max_size(size) {
        assert(base_ptr != 0);
        assert(max_size > 0);
    }

    virtual void* Alloc(int64 size) override {
        assert(curr_ptr + size < base_ptr + max_size);
        prev_ptr = curr_ptr;
        curr_ptr += size;
        return prev_ptr;
    }

    virtual void Free(void* ptr) override {
        assert(curr_ptr <= ptr && ptr < curr_ptr + max_size);
        //assert(ptr == prev_ptr);
        curr_ptr = (uint8*)ptr;
    }
};

struct FrameAllocator : StackAllocator {
    FrameAllocator(void* memory, int64 size) :
        StackAllocator(memory, size) {};

    void Reset() {
        curr_ptr = base_ptr;
    }

    virtual void Free(void*) override {}
};

struct DefaultAllocator : Allocator {
    virtual void* Alloc(int64 size) override {
        return malloc(size);
    }

    virtual void Free(void* ptr) override {
        return free(ptr);
    }
};

static DefaultAllocator default_alloc;
