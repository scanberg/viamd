cmake_minimum_required(VERSION 3.0)

project(viamd_app)

if(APPLE)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

if(MSVC)
  message(STATUS "MSVC configure was used")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /wd4201 /MP /GR-")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
else()
  message(STATUS "non MSVC configure was used")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types")
endif()

add_subdirectory(ext/glfw)
add_subdirectory(ext/imgui)
add_subdirectory(ext/nativefiledialog)
add_subdirectory(ext/viamd)
add_subdirectory(unittest)

include_directories(src/
                    ext/viamd/src
                    ext/viamd/ext/glm
                    ext/viamd/ext/gl3w
                    ext/Eigen
                    ext/imgui
                    ext/glfw/include
                    ext/nativefiledialog/src/include)

file(GLOB_RECURSE PLATFORM_FILES src/platform/*)

file(GLOB PROJECT_CONFIGS CMakeLists.txt
                          Readme.md
                         .gitattributes
                         .gitignore
                         .gitmodules)

set(APPLICATION_FILES
    src/main.cpp
)

source_group("platform" FILES ${PLATFORM_FILES})
source_group("application" FILES ${APPLICATION_FILES})

add_definitions(-DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
add_definitions(-DDATA_DIR=\"${PROJECT_SOURCE_DIR}/data\")

#add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${PLATFORM_FILES} ${APPLICATION_FILES})
add_executable(${PROJECT_NAME} ${PLATFORM_FILES} ${APPLICATION_FILES})

target_link_libraries(${PROJECT_NAME}
                      glfw
                      imgui
                      nativefiledialog
                      viamd_lib)