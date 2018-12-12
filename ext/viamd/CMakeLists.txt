cmake_minimum_required(VERSION 3.0)

project(viamd_lib)

add_subdirectory(ext/glm)
add_subdirectory(ext/gl3w)
add_subdirectory(ext/libxdrfile)
add_subdirectory(ext/tinyexpr)

include_directories(src
                    ext/glm
                    ext/gl3w
                    ext/libxdrfile/src
                    ext/tinyexpr
                    ext/svd3)


file(GLOB CORE_FILES src/core/*)
file(GLOB GFX_FILES src/gfx/*)
file(GLOB MOL_FILES src/mol/*)

file(GLOB PROJECT_CONFIGS CMakeLists.txt
                          Readme.md
                         .gitattributes
                         .gitignore
                         .gitmodules)

source_group("core" FILES ${CORE_FILES})
source_group("gfx" FILES ${GFX_FILES})
source_group("mol" FILES ${MOL_FILES})

add_library(${PROJECT_NAME} STATIC ${CORE_FILES} ${GFX_FILES} ${MOL_FILES})

target_link_libraries(${PROJECT_NAME}
                      gl3w
                      libxdrfile
                      tinyexpr)