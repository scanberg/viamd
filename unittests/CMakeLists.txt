cmake_minimum_required(VERSION 3.0)

project(unittests)

include_directories(../src
                    ../ext/glm
                    ../ext/libxdrfile/src
                    ../ext/tinyexpr
                    ../ext/Eigen
                    ../ext/svd3
                    ../ext/rsk_platform/include
                    ../ext/rsk_platform/ext/imgui)

file(GLOB_RECURSE PROJECT_HEADERS *.h
                                  ../src/*.h)
file(GLOB_RECURSE PROJECT_SOURCES *.cpp
                                  ../src/core/*.cpp
                                  ../src/gfx/*.cpp
                                  ../src/mol/*.cpp)

file(GLOB PROJECT_CONFIGS CMakeLists.txt)

source_group("Headers" FILES ${PROJECT_HEADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})

add_definitions(-DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
add_executable(${PROJECT_NAME} ${PROJECT_HEADERS} ${PROJECT_SOURCES})

target_link_libraries(${PROJECT_NAME}
                      rsk_platform
                      libxdrfile
                      tinyexpr)