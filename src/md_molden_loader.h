/**
 * @file md_molden_loader.h
 * @brief Molden file system loader interface
 * 
 * This header declares the Molden system loader entry points for VIAMD.
 */

#pragma once

#include <core/md_str.h>

struct md_allocator_i;
struct md_system_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Initialize md_system_t from a Molden file
 * 
 * Parses a Molden file (.molden, .mold) and converts it to an md_system_t
 * structure suitable for visualization.
 * 
 * @return true on success, false on failure
 */
bool md_molden_system_init_from_file(struct md_system_t* sys, str_t filename, struct md_allocator_i* alloc);

/**
 * @brief Initialize md_system_t from Molden file content
 */
bool md_molden_system_init_from_str(struct md_system_t* sys, str_t content, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
