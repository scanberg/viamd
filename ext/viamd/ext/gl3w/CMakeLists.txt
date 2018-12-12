### GL3W ###
project(gl3w)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/GL/gl3w.h
)
            
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/gl3w.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(gl3w STATIC ${HEADERS} ${SOURCES})