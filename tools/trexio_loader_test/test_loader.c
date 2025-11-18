#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <core/md_allocator.h>
#include <core/md_str.h>
#include <md_trexio.h>

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <trexio_file>\n", argv[0]);
        return 1;
    }
    
    const char* filename = argv[1];
    printf("Testing TREXIO loader with file: %s\n", filename);
    
    // Get allocator
    md_allocator_i* alloc = md_get_heap_allocator();
    if (!alloc) {
        fprintf(stderr, "ERROR: md_get_heap_allocator() returned NULL\n");
        return 1;
    }
    printf("Allocator obtained: %p\n", (void*)alloc);
    
    // Create TREXIO structure
    printf("Creating TREXIO structure...\n");
    md_trexio_t* trexio = md_trexio_create(alloc);
    if (!trexio) {
        fprintf(stderr, "ERROR: md_trexio_create() failed\n");
        return 1;
    }
    printf("TREXIO structure created: %p\n", (void*)trexio);
    
    // Parse file
    printf("Parsing TREXIO file...\n");
    str_t fname = { filename, strlen(filename) };
    bool success = md_trexio_parse_file(trexio, fname);
    if (!success) {
        fprintf(stderr, "ERROR: md_trexio_parse_file() failed\n");
        md_trexio_destroy(trexio);
        return 1;
    }
    printf("File parsed successfully\n");
    
    // Print basic info
    printf("Number of atoms: %zu\n", md_trexio_number_of_atoms(trexio));
    
    // Clean up
    printf("Cleaning up...\n");
    md_trexio_destroy(trexio);
    printf("Test completed successfully\n");
    
    return 0;
}
