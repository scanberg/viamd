### LIBXDRFILE ###

project(libxdrfile)
#enable_language (C)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile_trr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile_xtc.h
)
                 
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile_trr.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xdrfile_xtc.c
)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
add_library(libxdrfile STATIC ${HEADERS} ${SOURCES})