cmake_minimum_required(VERSION 3.0)

project(viamd)

if(APPLE)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
endif()

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /wd4201 /MP")
  #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_subdirectory(ext)

                    #ext/noc
                    #ext/glfw/include
                    #ext/gl3w/include
                    #ext/imgui
include_directories(src
                    ext/glm
                    ext/libxdrfile/src
                    ext/tinyexpr
                    ext/Eigen)

file(GLOB_RECURSE PROJECT_HEADERS src/*.h)
file(GLOB_RECURSE PROJECT_SOURCES src/*.cpp
                                  src/*.inl)

file(GLOB PROJECT_CONFIGS CMakeLists.txt
                          Readme.md
                         .gitattributes
                         .gitignore
                         .gitmodules)

source_group("Headers" FILES ${PROJECT_HEADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})

add_definitions(-DGLFW_INCLUDE_NONE
                -DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")

#add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${PROJECT_HEADERS} ${PROJECT_SOURCES})
add_executable(${PROJECT_NAME} ${PROJECT_HEADERS} ${PROJECT_SOURCES})

target_link_libraries(${PROJECT_NAME}
                      #glfw
                      #gl3w
                      #imgui
                      rsk_platform
                      libxdrfile
                      tinyexpr
                      #${GLFW_LIBRARIES}
                      #${OPENGL_LIBRARY}
                      )

#set_target_properties(${PROJECT_NAME} PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME})