#include <core/types.h>
#include <core/array.h>

// @TODO: Implement more hash functions... constexpr is nice!

namespace hash {
uint32 crc32(const char* data, int64 byte_size);

template <typename T>
uint32 crc32(Array<T> arr) {
    return crc32(arr.data, arr.count * sizeof(T));
}
}  // namespace hash