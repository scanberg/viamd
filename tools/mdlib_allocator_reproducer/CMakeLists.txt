cmake_minimum_required(VERSION 3.20)

# mdlib Allocator Reproducer Tool
# A minimal diagnostic harness to investigate allocator crashes

project(mdlib_allocator_reproducer LANGUAGES C)

# This tool depends on mdlib
if(NOT TARGET mdlib)
    message(FATAL_ERROR "mdlib target not found. This tool must be built as part of the viamd project.")
endif()

# Create the reproducer executable
add_executable(mdlib_allocator_reproducer
    main.c
)

# Link against mdlib to get the allocator implementation
target_link_libraries(mdlib_allocator_reproducer PRIVATE mdlib)

# Set C standard
set_property(TARGET mdlib_allocator_reproducer PROPERTY C_STANDARD 11)

# Platform-specific settings
if(MSVC)
    # MSVC-specific flags if needed
    target_compile_options(mdlib_allocator_reproducer PRIVATE /W4)
else()
    # GCC/Clang flags
    target_compile_options(mdlib_allocator_reproducer PRIVATE -Wall -Wextra)
endif()

# Install target (optional - useful for packaging)
install(TARGETS mdlib_allocator_reproducer
    RUNTIME DESTINATION bin
    COMPONENT tools
)

# Add a custom target to run the reproducer easily
add_custom_target(run_allocator_reproducer
    COMMAND mdlib_allocator_reproducer --iterations 100
    DEPENDS mdlib_allocator_reproducer
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running mdlib allocator reproducer with 100 iterations..."
)

message(STATUS "mdlib allocator reproducer tool configured")
message(STATUS "  Build target: mdlib_allocator_reproducer")
message(STATUS "  Run target:   run_allocator_reproducer")
