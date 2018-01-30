project(gl3w)
enable_language (C)

find_package(OpenGL REQUIRED)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/gl3w/include/GL/gl3w.h
    ${CMAKE_CURRENT_SOURCE_DIR}/gl3w/include/GL/glcorearb.h
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/gl3w/src/gl3w.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/gl3w/include/)

add_library(gl3w STATIC ${HEADERS} ${SOURCES})

target_link_libraries(gl3w
                      ${OPENGL_LIBRARY})

project(imgui)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/imgui.h
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/imgui.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/imgui_draw.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/imgui_demo.cpp
)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/)

add_library(imgui STATIC ${HEADERS} ${SOURCES})