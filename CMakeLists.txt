cmake_minimum_required(VERSION 3.0)

project(viamd)

find_package(OpenGL REQUIRED)

if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif()

if(APPLE)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
endif()

if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /wd4201")
  #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_subdirectory(ext)

include_directories(src
                    ext/glfw/include
                    ext/gl3w/include
                    ext/glm
                    ext/imgui
                    ext/ppk_assert
                    ext/libxdrfile/src
                    ext/nativefiledialog/src/include
                    ext/tinyexpr)

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

add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${PROJECT_HEADERS} ${PROJECT_SOURCES})

target_link_libraries(${PROJECT_NAME}
                      glfw
                      gl3w
                      imgui
                      ppk_assert
                      libxdrfile
                      nativefiledialog
                      tinyexpr
                      ${GLFW_LIBRARIES}
                      ${OPENGL_LIBRARY})

#set_target_properties(${PROJECT_NAME} PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
