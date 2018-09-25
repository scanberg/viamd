### OPENGL ###
find_package(OpenGL REQUIRED)

### GLFW ###
option(GLFW_BUILD_DOCS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_INSTALL OFF)
add_subdirectory(ext/glfw)

### NATIVEFILEDIALOG ###
add_subdirectory(ext/nativefiledialog)

### IMGUI ###
add_subdirectory(ext/imgui)

### RSK PLATFORM ###
project(rsk_platform)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/array_types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/common.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/gl.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/keys.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/log.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/platform.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/string_types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/string_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/core/types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GL/glcorearb.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GL/gl3w.h
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/log.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/platform.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/platform_unix.inl
    ${CMAKE_CURRENT_SOURCE_DIR}/src/platform_linux.inl
    ${CMAKE_CURRENT_SOURCE_DIR}/src/platform_osx.inl
    ${CMAKE_CURRENT_SOURCE_DIR}/src/platform_win32.inl
    ${CMAKE_CURRENT_SOURCE_DIR}/src/string_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gl3w.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/glfw/include/
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/nativefiledialog/src/include/
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/)

add_definitions(-DGLFW_INCLUDE_NONE)
add_library(rsk_platform STATIC ${HEADERS} ${SOURCES})

target_link_libraries(rsk_platform
                      PUBLIC imgui
                      PRIVATE glfw
                      PRIVATE nativefiledialog
                      PRIVATE ${glfw_LIBRARIES}
                      PRIVATE ${OPENGL_LIBRARY})