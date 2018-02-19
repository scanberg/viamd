#pragma once

#include <core/types.h>
#include <core/common.h>
#include <stdlib.h>

// @TODO: Support aligned allocations

struct Allocator {
    virtual void* alloc(int64) = 0;
    virtual void free(void*) = 0;
};

// @TODO: Encode size in stack-allocator and assert comparison with prev on free
// @TODO: Maybe not, let user free a larger chunk of the stack at once.
struct StackAllocator : Allocator {
    uint8* base_ptr;
    uint8* curr_ptr;
    int64 max_size;

    StackAllocator(void* memory, int64 size)
        : base_ptr((uint8*)memory), curr_ptr(base_ptr), max_size(size) {
        ASSERT(base_ptr != 0);
		ASSERT(max_size > 0);
    }

    virtual void* alloc(int64 size) override {
		ASSERT(curr_ptr + size < base_ptr + max_size);
        auto ptr = curr_ptr;
        curr_ptr += size;
        return ptr;
    }

    virtual void free(void* ptr) override {
		ASSERT(base_ptr <= ptr && ptr < base_ptr + max_size);
        curr_ptr = (uint8*)ptr;
    }
};

struct FrameAllocator : StackAllocator {
    FrameAllocator(void* memory, int64 size) :
        StackAllocator(memory, size) {};
	
    void reset() {
        curr_ptr = base_ptr;
    }

    virtual void free(void*) override {}
};

struct DefaultAllocator : Allocator {
    virtual void* alloc(int64 size) override {
        return std::malloc(size);
    }

    virtual void free(void* ptr) override {
        return std::free(ptr);
    }
};

static DefaultAllocator default_alloc;
