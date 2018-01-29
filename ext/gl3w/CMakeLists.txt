project(gl3w)
enable_language (C)

find_package(OpenGL REQUIRED)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GL/gl3w.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GL/glcorearb.h
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gl3w.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/)

add_library(gl3w STATIC ${HEADERS} ${SOURCES})

target_link_libraries(${PROJECT_NAME}
                      ${OPENGL_LIBRARY})