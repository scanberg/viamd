### GLFW ###
option(GLFW_BUILD_DOCS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_INSTALL OFF)

add_subdirectory(ext/glfw)
add_subdirectory(ext/nativefiledialog)

### OPENGL ###
find_package(OpenGL REQUIRED)

### IMGUI ###
project(imgui)

file(GLOB_RECURSE WIDGET_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/widgets/*.h)
file(GLOB_RECURSE WIDGET_INLINES ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/widgets/*.inl)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/imconfig.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/imgui.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/imgui_internal.h
    ${WIDGET_HEADERS}
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/imgui.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/imgui_draw.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/imgui_demo.cpp
    ${WIDGET_INLINES}
)

add_library(imgui STATIC ${HEADERS} ${SOURCES})

### PLATFORM ###
project(rsk_platform)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/gl3w/include/GL/gl3w.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/gl3w/include/GL/glcorearb.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/array_types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/common.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gl.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/keys.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/platform.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/string_types.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/string_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/types.h 
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/gl3w/src/gl3w.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/platform.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/string_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/log.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/gl3w/include/
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/glfw/include/
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/
                    ${CMAKE_CURRENT_SOURCE_DIR}/ext/nativefiledialog/src/include/
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/)

add_library(${PROJECT_NAME} ${HEADERS} ${SOURCES})

target_link_libraries(${PROJECT_NAME}
                      imgui
                      nativefiledialog
                      ${GLFW_LIBRARIES}
                      ${OPENGL_LIBRARY})