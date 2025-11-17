/**
 * mdlib Allocator Reproducer and Diagnostic Tool
 * 
 * This tool reproduces the md_get_heap_allocator() / md_alloc() call pattern
 * used by the TREXIO loader to investigate allocator crashes.
 * 
 * It emits detailed diagnostic logging including:
 * - Allocator pointer and state
 * - Allocation sizes and return values
 * - Any errors encountered
 * 
 * Usage: ./mdlib_allocator_reproducer [--iterations N]
 */

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_common.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#define DIAGNOSTIC_LOG(fmt, ...) \
    do { \
        printf("[DIAGNOSTIC] " fmt "\n", ##__VA_ARGS__); \
        fflush(stdout); \
    } while(0)

// Mimics the md_trexio_t structure allocation pattern
typedef struct test_data_t {
    md_allocator_i* alloc;
    
    // Nucleus data (similar to TREXIO)
    int64_t nucleus_num;
    double* nucleus_coord;
    double* nucleus_charge;
    char** nucleus_labels;
    
    // Basis set data (similar to TREXIO)
    int64_t basis_shell_num;
    int64_t basis_prim_num;
    int64_t* shell_ang_mom;
    int64_t* shell_factor;
    double* exponent;
    double* coefficient;
    
    // MO data (similar to TREXIO)
    int64_t mo_num;
    double* mo_coeff;
    double* mo_energy;
    double* mo_occupation;
} test_data_t;

// Forward declarations
static test_data_t* create_test_data(md_allocator_i* alloc);
static void free_test_data(test_data_t* data);
static bool simulate_trexio_allocations(test_data_t* data);
static void print_allocator_diagnostics(md_allocator_i* alloc, const char* context);

int main(int argc, char** argv) {
    int iterations = 1;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--iterations") == 0 && i + 1 < argc) {
            iterations = atoi(argv[i + 1]);
            i++;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            printf("Usage: %s [--iterations N]\n", argv[0]);
            printf("  --iterations N   Run N iterations of allocations (default: 1)\n");
            printf("  --help           Show this help message\n");
            return 0;
        }
    }
    
    DIAGNOSTIC_LOG("mdlib Allocator Reproducer Started");
    DIAGNOSTIC_LOG("Iterations: %d", iterations);
    DIAGNOSTIC_LOG("========================================");
    
    for (int iter = 0; iter < iterations; iter++) {
        DIAGNOSTIC_LOG("\n--- Iteration %d/%d ---", iter + 1, iterations);
        
        // Get heap allocator (same as TREXIO loader does)
        DIAGNOSTIC_LOG("Calling md_get_heap_allocator()...");
        md_allocator_i* heap_alloc = md_get_heap_allocator();
        print_allocator_diagnostics(heap_alloc, "After md_get_heap_allocator()");
        
        // Create test data structure
        DIAGNOSTIC_LOG("Creating test data structure...");
        test_data_t* data = create_test_data(heap_alloc);
        if (!data) {
            DIAGNOSTIC_LOG("ERROR: Failed to create test data structure");
            return 1;
        }
        print_allocator_diagnostics(heap_alloc, "After create_test_data()");
        
        // Simulate TREXIO-like allocations
        DIAGNOSTIC_LOG("Simulating TREXIO allocation pattern...");
        bool success = simulate_trexio_allocations(data);
        if (!success) {
            DIAGNOSTIC_LOG("ERROR: Allocation simulation failed");
            free_test_data(data);
            return 1;
        }
        print_allocator_diagnostics(heap_alloc, "After simulate_trexio_allocations()");
        
        // Free everything
        DIAGNOSTIC_LOG("Freeing test data...");
        free_test_data(data);
        print_allocator_diagnostics(heap_alloc, "After free_test_data()");
        
        DIAGNOSTIC_LOG("Iteration %d completed successfully", iter + 1);
    }
    
    DIAGNOSTIC_LOG("\n========================================");
    DIAGNOSTIC_LOG("All iterations completed successfully!");
    DIAGNOSTIC_LOG("No crashes detected in allocator operations.");
    
    return 0;
}

static test_data_t* create_test_data(md_allocator_i* alloc) {
    if (!alloc) {
        DIAGNOSTIC_LOG("ERROR: NULL allocator passed to create_test_data");
        return NULL;
    }
    
    DIAGNOSTIC_LOG("  Allocating test_data_t structure (size: %zu bytes)", sizeof(test_data_t));
    
    // This mimics: md_trexio_t* trexio = (md_trexio_t*)md_alloc(alloc, sizeof(md_trexio_t));
    test_data_t* data = (test_data_t*)md_alloc(alloc, sizeof(test_data_t));
    if (!data) {
        DIAGNOSTIC_LOG("  ERROR: md_alloc returned NULL for test_data_t");
        return NULL;
    }
    
    DIAGNOSTIC_LOG("  test_data_t allocated at address: %p", (void*)data);
    
    // Zero-initialize (like MEMSET in TREXIO)
    MEMSET(data, 0, sizeof(test_data_t));
    data->alloc = alloc;
    
    DIAGNOSTIC_LOG("  test_data_t initialized successfully");
    return data;
}

static bool simulate_trexio_allocations(test_data_t* data) {
    if (!data || !data->alloc) {
        DIAGNOSTIC_LOG("  ERROR: Invalid data or allocator in simulate_trexio_allocations");
        return false;
    }
    
    md_allocator_i* alloc = data->alloc;
    
    // Simulate reading a small molecule (e.g., H2O: 3 atoms)
    data->nucleus_num = 3;
    DIAGNOSTIC_LOG("  Setting nucleus_num = %lld", (long long)data->nucleus_num);
    
    // Allocate coordinates (3 * num_atoms doubles)
    size_t coord_size = 3 * data->nucleus_num * sizeof(double);
    DIAGNOSTIC_LOG("  Allocating nucleus_coord (%zu bytes)...", coord_size);
    data->nucleus_coord = (double*)md_alloc(alloc, coord_size);
    if (!data->nucleus_coord) {
        DIAGNOSTIC_LOG("  ERROR: md_alloc returned NULL for nucleus_coord");
        return false;
    }
    DIAGNOSTIC_LOG("  nucleus_coord allocated at: %p", (void*)data->nucleus_coord);
    
    // Initialize with dummy data
    for (int64_t i = 0; i < 3 * data->nucleus_num; i++) {
        data->nucleus_coord[i] = (double)i * 0.529177;  // Convert like Bohr to Angstrom
    }
    DIAGNOSTIC_LOG("  nucleus_coord initialized with test data");
    
    // Allocate charges
    size_t charge_size = data->nucleus_num * sizeof(double);
    DIAGNOSTIC_LOG("  Allocating nucleus_charge (%zu bytes)...", charge_size);
    data->nucleus_charge = (double*)md_alloc(alloc, charge_size);
    if (!data->nucleus_charge) {
        DIAGNOSTIC_LOG("  ERROR: md_alloc returned NULL for nucleus_charge");
        return false;
    }
    DIAGNOSTIC_LOG("  nucleus_charge allocated at: %p", (void*)data->nucleus_charge);
    
    for (int64_t i = 0; i < data->nucleus_num; i++) {
        data->nucleus_charge[i] = (double)(i + 1);  // Dummy charges
    }
    
    // Allocate labels (array of char pointers)
    size_t labels_size = data->nucleus_num * sizeof(char*);
    DIAGNOSTIC_LOG("  Allocating nucleus_labels array (%zu bytes)...", labels_size);
    data->nucleus_labels = (char**)md_alloc(alloc, labels_size);
    if (!data->nucleus_labels) {
        DIAGNOSTIC_LOG("  ERROR: md_alloc returned NULL for nucleus_labels");
        return false;
    }
    DIAGNOSTIC_LOG("  nucleus_labels array allocated at: %p", (void*)data->nucleus_labels);
    
    // Allocate individual label strings
    for (int64_t i = 0; i < data->nucleus_num; i++) {
        const char* labels[] = {"H", "O", "H"};
        size_t label_len = strlen(labels[i]) + 1;
        DIAGNOSTIC_LOG("  Allocating nucleus_labels[%lld] (%zu bytes)...", (long long)i, label_len);
        data->nucleus_labels[i] = (char*)md_alloc(alloc, label_len);
        if (!data->nucleus_labels[i]) {
            DIAGNOSTIC_LOG("  ERROR: md_alloc returned NULL for nucleus_labels[%lld]", (long long)i);
            return false;
        }
        strcpy(data->nucleus_labels[i], labels[i]);
        DIAGNOSTIC_LOG("  nucleus_labels[%lld] = \"%s\" at %p", (long long)i, 
                      data->nucleus_labels[i], (void*)data->nucleus_labels[i]);
    }
    
    // Simulate basis set allocations (smaller test case)
    data->basis_shell_num = 5;
    data->basis_prim_num = 10;
    DIAGNOSTIC_LOG("  Setting basis_shell_num = %lld, basis_prim_num = %lld", 
                  (long long)data->basis_shell_num, (long long)data->basis_prim_num);
    
    // Allocate shell data
    size_t shell_ang_size = data->basis_shell_num * sizeof(int64_t);
    DIAGNOSTIC_LOG("  Allocating shell_ang_mom (%zu bytes)...", shell_ang_size);
    data->shell_ang_mom = (int64_t*)md_alloc(alloc, shell_ang_size);
    if (!data->shell_ang_mom) {
        DIAGNOSTIC_LOG("  ERROR: md_alloc returned NULL for shell_ang_mom");
        return false;
    }
    DIAGNOSTIC_LOG("  shell_ang_mom allocated at: %p", (void*)data->shell_ang_mom);
    
    // Initialize shell data
    for (int64_t i = 0; i < data->basis_shell_num; i++) {
        data->shell_ang_mom[i] = i % 3;  // 0=s, 1=p, 2=d
    }
    
    // Allocate primitive data
    size_t prim_size = data->basis_prim_num * sizeof(double);
    DIAGNOSTIC_LOG("  Allocating exponent (%zu bytes)...", prim_size);
    data->exponent = (double*)md_alloc(alloc, prim_size);
    if (!data->exponent) {
        DIAGNOSTIC_LOG("  ERROR: md_alloc returned NULL for exponent");
        return false;
    }
    DIAGNOSTIC_LOG("  exponent allocated at: %p", (void*)data->exponent);
    
    for (int64_t i = 0; i < data->basis_prim_num; i++) {
        data->exponent[i] = (double)(i + 1) * 0.5;
    }
    
    // Simulate MO data allocations
    data->mo_num = 7;
    DIAGNOSTIC_LOG("  Setting mo_num = %lld", (long long)data->mo_num);
    
    size_t mo_energy_size = data->mo_num * sizeof(double);
    DIAGNOSTIC_LOG("  Allocating mo_energy (%zu bytes)...", mo_energy_size);
    data->mo_energy = (double*)md_alloc(alloc, mo_energy_size);
    if (!data->mo_energy) {
        DIAGNOSTIC_LOG("  ERROR: md_alloc returned NULL for mo_energy");
        return false;
    }
    DIAGNOSTIC_LOG("  mo_energy allocated at: %p", (void*)data->mo_energy);
    
    for (int64_t i = 0; i < data->mo_num; i++) {
        data->mo_energy[i] = -1.0 * (double)(i + 1);  // Dummy energies
    }
    
    DIAGNOSTIC_LOG("  All allocations completed successfully");
    return true;
}

static void free_test_data(test_data_t* data) {
    if (!data) {
        DIAGNOSTIC_LOG("  Warning: free_test_data called with NULL data");
        return;
    }
    
    if (!data->alloc) {
        DIAGNOSTIC_LOG("  ERROR: data->alloc is NULL in free_test_data");
        return;
    }
    
    md_allocator_i* alloc = data->alloc;
    
    // Free in reverse order of allocation
    DIAGNOSTIC_LOG("  Freeing allocated memory...");
    
    if (data->mo_energy) {
        DIAGNOSTIC_LOG("  Freeing mo_energy (%zu bytes)", data->mo_num * sizeof(double));
        md_free(alloc, data->mo_energy, data->mo_num * sizeof(double));
    }
    
    if (data->exponent) {
        DIAGNOSTIC_LOG("  Freeing exponent (%zu bytes)", data->basis_prim_num * sizeof(double));
        md_free(alloc, data->exponent, data->basis_prim_num * sizeof(double));
    }
    
    if (data->shell_ang_mom) {
        DIAGNOSTIC_LOG("  Freeing shell_ang_mom (%zu bytes)", data->basis_shell_num * sizeof(int64_t));
        md_free(alloc, data->shell_ang_mom, data->basis_shell_num * sizeof(int64_t));
    }
    
    if (data->nucleus_labels) {
        for (int64_t i = 0; i < data->nucleus_num; i++) {
            if (data->nucleus_labels[i]) {
                DIAGNOSTIC_LOG("  Freeing nucleus_labels[%lld]", (long long)i);
                md_free(alloc, data->nucleus_labels[i], strlen(data->nucleus_labels[i]) + 1);
            }
        }
        DIAGNOSTIC_LOG("  Freeing nucleus_labels array (%zu bytes)", data->nucleus_num * sizeof(char*));
        md_free(alloc, data->nucleus_labels, data->nucleus_num * sizeof(char*));
    }
    
    if (data->nucleus_charge) {
        DIAGNOSTIC_LOG("  Freeing nucleus_charge (%zu bytes)", data->nucleus_num * sizeof(double));
        md_free(alloc, data->nucleus_charge, data->nucleus_num * sizeof(double));
    }
    
    if (data->nucleus_coord) {
        DIAGNOSTIC_LOG("  Freeing nucleus_coord (%zu bytes)", 3 * data->nucleus_num * sizeof(double));
        md_free(alloc, data->nucleus_coord, 3 * data->nucleus_num * sizeof(double));
    }
    
    DIAGNOSTIC_LOG("  Freeing test_data_t structure (%zu bytes)", sizeof(test_data_t));
    md_free(alloc, data, sizeof(test_data_t));
    
    DIAGNOSTIC_LOG("  All memory freed successfully");
}

static void print_allocator_diagnostics(md_allocator_i* alloc, const char* context) {
    DIAGNOSTIC_LOG("  [%s]", context);
    
    if (!alloc) {
        DIAGNOSTIC_LOG("    ERROR: Allocator is NULL!");
        return;
    }
    
    DIAGNOSTIC_LOG("    Allocator address: %p", (void*)alloc);
    DIAGNOSTIC_LOG("    Allocator inst:    %p", (void*)alloc->inst);
    DIAGNOSTIC_LOG("    Allocator realloc: %p", (void*)(uintptr_t)alloc->realloc);
    
    // Verify allocator is still functional by doing a small allocation/free
    DIAGNOSTIC_LOG("    Testing allocator with small allocation...");
    void* test_ptr = md_alloc(alloc, 16);
    if (!test_ptr) {
        DIAGNOSTIC_LOG("    ERROR: Test allocation failed!");
        return;
    }
    DIAGNOSTIC_LOG("    Test allocation successful at: %p", test_ptr);
    md_free(alloc, test_ptr, 16);
    DIAGNOSTIC_LOG("    Test free completed");
}
