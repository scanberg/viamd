cmake_minimum_required(VERSION 3.20)
project(trexio_loader_test LANGUAGES C)

if(NOT TARGET mdlib)
    message(FATAL_ERROR "mdlib target not found")
endif()

add_executable(trexio_loader_test test_loader.c)
target_link_libraries(trexio_loader_test PRIVATE mdlib)
set_property(TARGET trexio_loader_test PROPERTY C_STANDARD 11)
