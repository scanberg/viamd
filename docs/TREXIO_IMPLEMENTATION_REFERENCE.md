# TREXIO Implementation Reference

This document provides a minimal reference implementation showing the correct allocator usage patterns for the TREXIO loader.

## Correct md_trexio Structure Pattern

Following the VeloxChem pattern from `md_vlx.c`:

```c
// md_trexio.h
#pragma once

#include <core/md_allocator.h>
#include <core/md_str.h>
#include <md_system.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct md_trexio_t md_trexio_t;

// Creation and destruction - MUST pass allocator
md_trexio_t* md_trexio_create(struct md_allocator_i* alloc);
void md_trexio_destroy(md_trexio_t* trexio, struct md_allocator_i* alloc);

// Parsing - MUST pass allocator  
bool md_trexio_parse_file(md_trexio_t* trexio, str_t filename, struct md_allocator_i* alloc);

// Data access
int64_t md_trexio_num_atoms(const md_trexio_t* trexio);
const double* md_trexio_atom_coordinates(const md_trexio_t* trexio);
const int32_t* md_trexio_atomic_numbers(const md_trexio_t* trexio);

// System loader interface
struct md_system_loader_i* md_trexio_system_loader(void);

#ifdef __cplusplus
}
#endif
```

## Correct Implementation Pattern

```c
// md_trexio.c
#include "md_trexio.h"
#include <core/md_log.h>
#include <core/md_array.h>
#include <trexio.h>

struct md_trexio_t {
    int64_t num_atoms;
    double* atom_coordinates;  // Array of 3 * num_atoms
    int32_t* atomic_numbers;   // Array of num_atoms
    char** atom_labels;        // Array of num_atoms
};

// ===== CRITICAL: Proper allocator validation =====
md_trexio_t* md_trexio_create(md_allocator_i* alloc) {
    // STEP 1: Validate allocator
    if (!alloc) {
        MD_LOG_ERROR("NULL allocator passed to md_trexio_create");
        return NULL;
    }
    
    if (!alloc->realloc) {
        MD_LOG_ERROR("Allocator missing realloc function");
        return NULL;
    }
    
    // STEP 2: Allocate structure
    md_trexio_t* trexio = (md_trexio_t*)md_alloc(alloc, sizeof(md_trexio_t));
    if (!trexio) {
        MD_LOG_ERROR("Failed to allocate TREXIO structure");
        return NULL;
    }
    
    // STEP 3: Initialize to zero
    MEMSET(trexio, 0, sizeof(md_trexio_t));
    
    return trexio;
}

void md_trexio_destroy(md_trexio_t* trexio, md_allocator_i* alloc) {
    if (!trexio || !alloc) return;
    
    // Free arrays if allocated
    if (trexio->atom_coordinates) {
        md_free(alloc, trexio->atom_coordinates, sizeof(double) * 3 * trexio->num_atoms);
    }
    
    if (trexio->atomic_numbers) {
        md_free(alloc, trexio->atomic_numbers, sizeof(int32_t) * trexio->num_atoms);
    }
    
    if (trexio->atom_labels) {
        for (int64_t i = 0; i < trexio->num_atoms; i++) {
            if (trexio->atom_labels[i]) {
                // Free individual labels (TREXIO uses malloc)
                free(trexio->atom_labels[i]);
            }
        }
        md_free(alloc, trexio->atom_labels, sizeof(char*) * trexio->num_atoms);
    }
    
    // Free main structure
    md_free(alloc, trexio, sizeof(md_trexio_t));
}

bool md_trexio_parse_file(md_trexio_t* trexio, str_t filename, md_allocator_i* alloc) {
    // ===== CRITICAL: Validate inputs =====
    if (!trexio) {
        MD_LOG_ERROR("NULL TREXIO structure");
        return false;
    }
    
    if (!alloc) {
        MD_LOG_ERROR("NULL allocator - cannot parse TREXIO file");
        return false;
    }
    
    if (!alloc->realloc) {
        MD_LOG_ERROR("Allocator has NULL realloc function");
        return false;
    }
    
    // Convert str_t to C string for TREXIO API
    char path[4096];
    if (filename.len >= sizeof(path)) {
        MD_LOG_ERROR("Filename too long");
        return false;
    }
    MEMCPY(path, filename.ptr, filename.len);
    path[filename.len] = '\0';
    
    // Open TREXIO file
    trexio_exit_code rc;
    trexio_t* file = trexio_open(path, 'r', TREXIO_AUTO, &rc);
    if (!file || rc != TREXIO_SUCCESS) {
        MD_LOG_ERROR("Failed to open TREXIO file: %s", trexio_string_of_error(rc));
        return false;
    }
    
    // Read number of atoms
    int32_t num_atoms_i32;
    rc = trexio_read_nucleus_num(file, &num_atoms_i32);
    if (rc != TREXIO_SUCCESS) {
        MD_LOG_ERROR("Failed to read nucleus number");
        trexio_close(file);
        return false;
    }
    trexio->num_atoms = (int64_t)num_atoms_i32;
    
    // ===== Allocate with validated allocator =====
    
    // Allocate coordinates array
    trexio->atom_coordinates = (double*)md_alloc(alloc, sizeof(double) * 3 * trexio->num_atoms);
    if (!trexio->atom_coordinates) {
        MD_LOG_ERROR("Failed to allocate coordinates array");
        trexio_close(file);
        return false;
    }
    
    // Read coordinates (in Bohr)
    rc = trexio_read_nucleus_coord(file, trexio->atom_coordinates);
    if (rc != TREXIO_SUCCESS) {
        MD_LOG_ERROR("Failed to read nucleus coordinates");
        trexio_close(file);
        return false;
    }
    
    // Allocate atomic numbers array
    trexio->atomic_numbers = (int32_t*)md_alloc(alloc, sizeof(int32_t) * trexio->num_atoms);
    if (!trexio->atomic_numbers) {
        MD_LOG_ERROR("Failed to allocate atomic numbers array");
        trexio_close(file);
        return false;
    }
    
    // Read atomic charges (convert to atomic numbers)
    double* charges = (double*)md_alloc(alloc, sizeof(double) * trexio->num_atoms);
    if (!charges) {
        MD_LOG_ERROR("Failed to allocate temporary charges array");
        trexio_close(file);
        return false;
    }
    
    rc = trexio_read_nucleus_charge(file, charges);
    if (rc == TREXIO_SUCCESS) {
        for (int64_t i = 0; i < trexio->num_atoms; i++) {
            trexio->atomic_numbers[i] = (int32_t)charges[i];
        }
    }
    md_free(alloc, charges, sizeof(double) * trexio->num_atoms);
    
    trexio_close(file);
    return true;
}

// ===== System Loader Interface =====

static bool trexio_sys_init_from_file(
    struct md_system_o* inst, 
    str_t filename, 
    struct md_allocator_i* alloc, 
    uint32_t flags
) {
    (void)flags;
    
    // ===== CRITICAL: Validate allocator FIRST =====
    if (!alloc || !alloc->realloc) {
        MD_LOG_ERROR("Invalid allocator in TREXIO system loader");
        return false;
    }
    
    md_trexio_t* trexio = (md_trexio_t*)inst;
    
    // Parse file with validated allocator
    if (!md_trexio_parse_file(trexio, filename, alloc)) {
        return false;
    }
    
    // Populate md_system_t structure
    md_system_t* sys = (md_system_t*)inst;
    
    // Set atom count
    sys->atom.count = trexio->num_atoms;
    
    // Allocate and populate coordinates
    sys->atom.x = (float*)md_alloc(alloc, sizeof(float) * trexio->num_atoms);
    sys->atom.y = (float*)md_alloc(alloc, sizeof(float) * trexio->num_atoms);
    sys->atom.z = (float*)md_alloc(alloc, sizeof(float) * trexio->num_atoms);
    
    if (!sys->atom.x || !sys->atom.y || !sys->atom.z) {
        MD_LOG_ERROR("Failed to allocate coordinate arrays");
        return false;
    }
    
    // Convert Bohr to Angstrom (CODATA 2018: 1 Bohr = 0.529177210903 Angstrom)
    const double BOHR_TO_ANGSTROM = 0.529177210903;
    for (int64_t i = 0; i < trexio->num_atoms; i++) {
        sys->atom.x[i] = (float)(trexio->atom_coordinates[i * 3 + 0] * BOHR_TO_ANGSTROM);
        sys->atom.y[i] = (float)(trexio->atom_coordinates[i * 3 + 1] * BOHR_TO_ANGSTROM);
        sys->atom.z[i] = (float)(trexio->atom_coordinates[i * 3 + 2] * BOHR_TO_ANGSTROM);
    }
    
    return true;
}

static bool trexio_sys_init_from_str(
    struct md_system_o* inst,
    str_t str,
    struct md_allocator_i* alloc,
    uint32_t flags
) {
    (void)inst;
    (void)str;
    (void)alloc;
    (void)flags;
    
    // TREXIO doesn't support loading from string
    MD_LOG_ERROR("TREXIO loader does not support loading from string");
    return false;
}

static void trexio_sys_free(struct md_system_o* inst, struct md_allocator_i* alloc) {
    if (!inst || !alloc) return;
    
    md_trexio_t* trexio = (md_trexio_t*)inst;
    md_trexio_destroy(trexio, alloc);
}

static md_system_loader_i trexio_system_loader = {
    trexio_sys_init_from_str,
    trexio_sys_init_from_file,
    trexio_sys_free,
};

md_system_loader_i* md_trexio_system_loader(void) {
    return &trexio_system_loader;
}
```

## Correct Test Pattern

```c
// test_trexio.c
#include "utest.h"
#include <md_trexio.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

UTEST(trexio, allocator_validation) {
    // Verify heap allocator is valid
    md_allocator_i* alloc = md_get_heap_allocator();
    ASSERT_NE(alloc, NULL);
    ASSERT_NE(alloc->realloc, NULL);
}

UTEST(trexio, create_and_destroy) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Create TREXIO structure
    md_trexio_t* trexio = md_trexio_create(alloc);
    ASSERT_NE(trexio, NULL);
    
    // Destroy
    md_trexio_destroy(trexio, alloc);
}

UTEST(trexio, parse_water) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Create
    md_trexio_t* trexio = md_trexio_create(alloc);
    ASSERT_NE(trexio, NULL);
    
    // Parse
    bool result = md_trexio_parse_file(
        trexio,
        STR_LIT(MD_UNITTEST_DATA_DIR "/trexio/water.trexio"),
        alloc
    );
    ASSERT_TRUE(result);
    
    // Verify
    EXPECT_EQ(3, md_trexio_num_atoms(trexio));
    
    const double* coords = md_trexio_atom_coordinates(trexio);
    ASSERT_NE(coords, NULL);
    
    const int32_t* atomic_nums = md_trexio_atomic_numbers(trexio);
    ASSERT_NE(atomic_nums, NULL);
    EXPECT_EQ(8, atomic_nums[0]);  // Oxygen
    EXPECT_EQ(1, atomic_nums[1]);  // Hydrogen
    EXPECT_EQ(1, atomic_nums[2]);  // Hydrogen
    
    // Cleanup
    md_trexio_destroy(trexio, alloc);
}

UTEST(trexio, null_allocator_handling) {
    // Test that NULL allocator is properly rejected
    md_trexio_t* trexio = md_trexio_create(NULL);
    EXPECT_EQ(NULL, trexio);  // Should return NULL, not crash
}
```

## Key Takeaways

1. **Always validate allocator before use**:
   ```c
   if (!alloc || !alloc->realloc) {
       MD_LOG_ERROR("Invalid allocator");
       return false;
   }
   ```

2. **Pass allocator to ALL functions** that allocate memory

3. **Follow VeloxChem pattern** for consistency

4. **Check allocation results** before using pointers

5. **Use proper cleanup** in destroy functions

6. **Test with NULL allocator** to ensure graceful handling
