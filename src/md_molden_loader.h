/**
 * @file md_molden_loader.h
 * @brief Molden file system loader interface
 * 
 * This header declares the Molden system loader interface for VIAMD.
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

struct md_system_loader_i;

/**
 * @brief Get Molden system loader interface
 * 
 * Returns the system loader interface for Molden files (.molden, .mold).
 * This loader parses Molden files and converts them to md_system_t structures
 * for visualization.
 * 
 * @return Pointer to Molden system loader interface
 */
struct md_system_loader_i* md_molden_system_loader(void);

#ifdef __cplusplus
}
#endif
