cmake_minimum_required(VERSION 3.5)

project(viamd)

if (APPLE)
	set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
endif()

if (MSVC)
	message(STATUS "MSVC configure was used")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /wd4201 /wd4324 /MP /GR-")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
else()
	message(STATUS "non MSVC configure was used")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types")
endif()

option(BUILD_SHARED_LIBS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_BUILD_DOCS OFF)
option(GLFW_INSTALL OFF)
option(GLFW_VULKAN_STATIC OFF)

option(VIAMD_CREATE_MACOSX_BUNDLE "Build a macosx bundle instead of just an executable" OFF)
option(VIAMD_USE_RELATIVE_PATHS "Use relative paths for resource directories" OFF)

add_subdirectory(ext/glfw)
add_subdirectory(ext/imgui)
add_subdirectory(ext/nativefiledialog)
add_subdirectory(ext/tinyexpr)
add_subdirectory(ext/mdutils)
add_subdirectory(unittest)

file(GLOB PLATFORM_FILES src/platform/*)
file(GLOB SRC_FILES src/*)

source_group("platform" FILES ${PLATFORM_FILES})

set(OSX_BUNDLE "")
if (VIAMD_CREATE_MACOSX_BUNDLE)
	set(OSX_BUNDLE "MACOSX_BUNDLE")
endif()

add_executable(viamd ${OSX_BUNDLE} ${SRC_FILES} ${PLATFORM_FILES})

set(RESOURCE_DIR "${PROJECT_SOURCE_DIR}/")
if (VIAMD_USE_RELATIVE_PATHS)
    set(RESOURCE_DIR "")
endif()

target_compile_definitions(viamd PRIVATE
	VIAMD_DATA_DIR=\"${RESOURCE_DIR}data\"
	VIAMD_SHADER_DIR=\"${RESOURCE_DIR}shaders\"
	VIAMD_IMAGE_DIR=\"${RESOURCE_DIR}images\"
	VIAMD_SCREENSHOT_DIR=\"${RESOURCE_DIR}screenshots\")

#target_compile_options(viamd PRIVATE -W4)
target_compile_features(viamd PRIVATE cxx_std_17)

target_include_directories(viamd
#    PUBLIC 
#        $<INSTALL_INTERFACE:src>    
#        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(viamd
    glfw
    imgui
    nativefiledialog
	tinyexpr
    mdutils
)