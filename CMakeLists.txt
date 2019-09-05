cmake_minimum_required(VERSION 3.5)

# add_compile_options_config(<CONFIG> <option> ...)
function(add_compile_options_config CONFIG)
    foreach(opt ${ARGN})
        add_compile_options("$<$<CONFIG:${CONFIG}>:${opt}>")
    endforeach()
endfunction()

set(OSX_BUNDLE "")
if (APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
    if (VIAMD_CREATE_MACOSX_BUNDLE)
        set(OSX_BUNDLE "MACOSX_BUNDLE")
    endif()
endif()

project(viamd)

option(VIAMD_PROFILE_COMPILE_TIMES ON)

if (MSVC)
    message(STATUS "MSVC configuration was used")
    add_compile_options(/MT /W4 /wd4201 /wd4324 /MP /GR- /fp:fast /MP)# /arch:AVX2

if (VIAMD_PROFILE_COMPILE_TIMES)
	add_compile_options(/Bt+ /d1reportTime /d2cgsummary)
    add_link_options(/time+)
endif()
	add_compile_options_config(RELEASE /Oi /GS-)
    add_compile_options_config(DEBUG /MTd)
else()
	message(STATUS "non MSVC configuration was used")
    add_compile_options(
        -Wall -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types -mavx -maes
    )
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types")
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
add_subdirectory(ext/stb)
add_subdirectory(ext/gl3w)
#add_subdirectory(ext/base91)
add_subdirectory(unittest)
add_subdirectory(benchmark)

file(GLOB PLATFORM_FILES src/platform/*)
file(GLOB SRC_FILES src/*)
file(GLOB GFX_FILES src/gfx/*)

source_group("platform" FILES ${PLATFORM_FILES})
source_group("gfx" FILES ${GFX_FILES})

add_executable(viamd ${OSX_BUNDLE} ${SRC_FILES} ${PLATFORM_FILES} ${GFX_FILES})

set(RESOURCE_DIR "${PROJECT_SOURCE_DIR}/")
if (VIAMD_USE_RELATIVE_PATHS)
    set(RESOURCE_DIR "")
endif()

target_compile_definitions(viamd PRIVATE
	VIAMD_DATA_DIR=\"${RESOURCE_DIR}data\"
	VIAMD_SHADER_DIR=\"${RESOURCE_DIR}shaders\"
	VIAMD_IMAGE_DIR=\"${RESOURCE_DIR}images\"
	VIAMD_SCREENSHOT_DIR=\"${RESOURCE_DIR}screenshots\")

target_compile_features(viamd PRIVATE cxx_std_17)

target_include_directories(viamd
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/gl3w
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/Eigen
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/svd3
)

target_link_libraries(viamd
    gl3w
    glfw
    imgui
    nativefiledialog
	tinyexpr
    mdutils
    stb
)