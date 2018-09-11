add_subdirectory(rsk_platform)


### LIBXDRFILE ###

project(libxdrfile)
#enable_language (C)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/libxdrfile/src/xdrfile.h
    ${CMAKE_CURRENT_SOURCE_DIR}/libxdrfile/src/xdrfile_trr.h
    ${CMAKE_CURRENT_SOURCE_DIR}/libxdrfile/src/xdrfile_xtc.h
)
                 
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/libxdrfile/src/xdrfile.c
    ${CMAKE_CURRENT_SOURCE_DIR}/libxdrfile/src/xdrfile_trr.c
    ${CMAKE_CURRENT_SOURCE_DIR}/libxdrfile/src/xdrfile_xtc.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/libxdrfile/src)
add_library(libxdrfile STATIC ${HEADERS} ${SOURCES})


### TINYEXPR ###

project(tinyexpr)
add_library(tinyexpr STATIC ${CMAKE_CURRENT_SOURCE_DIR}/tinyexpr/tinyexpr.c)
