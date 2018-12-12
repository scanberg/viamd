cmake_minimum_required(VERSION 3.0)

project(unittests)

include_directories(../ext/viamd/src
                    ../ext/viamd/ext/glm
                    ../ext/viamd/ext/gl3w)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} viamd_lib)