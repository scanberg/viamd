#pragma once

#include <core/md_str.h>

#if VIAMD_TREXIO

#ifdef __cplusplus
extern "C" {
#endif

bool md_trexio_system_init_from_file(struct md_system_t* sys, str_t filename, struct md_allocator_i* alloc);
bool md_trexio_path_is_trexio(str_t filename);

#ifdef __cplusplus
}
#endif

#endif // VIAMD_TREXIO
