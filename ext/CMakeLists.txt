### GLFW ###

option(GLFW_BUILD_DOCS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_INSTALL OFF)

add_subdirectory(glfw)





### GL3W ###

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





### IMGUI ###

project(imgui)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/imgui.h
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/imgui.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/imgui_draw.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui/imgui_demo.cpp
)

add_library(imgui STATIC ${HEADERS} ${SOURCES})





### PPK ASSERT ###

project(ppk_assert)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/ppk_assert/ppk_assert.h
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ppk_assert/ppk_assert.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ppk_assert)
add_library(ppk_assert STATIC ${HEADERS} ${SOURCES})



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