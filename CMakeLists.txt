cmake_minimum_required(VERSION 3.15)

INCLUDE(ext/mdlib/cmake/common.cmake)

project(viamd VERSION 0.1.5 LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDERS OFF)

# GLOBAL OPTIONS
option(BUILD_SHARED_LIBS OFF)
option(USE_MSVC_RUNTIME_LIBRARY_DLL OFF)

if (USE_MSVC_RUNTIME_LIBRARY_DLL)
    set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>DLL)
else()
    set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreaded$<$<CONFIG:Debug>:Debug>)
endif()

# GLFW OPTIONS
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_BUILD_DOCS OFF)
option(GLFW_INSTALL OFF)
option(GLFW_VULKAN_STATIC OFF)

# ENKI OPTIONS
option(ENKITS_BUILD_C_INTERFACE OFF)
option(ENKITS_BUILD_EXAMPLES OFF)
option(ENKITS_BUILD_SHARED OFF)
option(ENKITS_INSTALL OFF)

find_package(OpenGL REQUIRED)

add_subdirectory(ext/mdlib)
add_subdirectory(ext/glfw)
add_subdirectory(ext/imgui)
add_subdirectory(ext/nativefiledialog)
add_subdirectory(ext/stb)
add_subdirectory(ext/enkiTS)
add_subdirectory(ext/ImGuiColorTextEdit)
add_subdirectory(ext/implot)
add_subdirectory(ext/atomic_queue)
add_subdirectory(ext/imgui_notify)

# VIAMD OPTIONS
option(VIAMD_CREATE_MACOSX_BUNDLE "Build a macosx bundle instead of just an executable" OFF)
option(VIAMD_LINK_STDLIB_STATIC "Link against stdlib statically" ${MD_LINK_STDLIB_STATIC})

# Copy many of the fields from mdlib
set(VIAMD_DEFINES ${MD_DEFS})
set(VIAMD_STDLIBS)
set(VIAMD_FLAGS ${MD_FLAGS})
set(VIAMD_FLAGS_DEB ${MD_FLAGS_DEB} "/INCREMENTAL")
set(VIAMD_FLAGS_REL ${MD_FLAGS_REL})

if (MD_LINK_STDLIB_STATIC)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(VIAMD_STDLIBS "-static-libgcc -static-libstdc++")
    endif()
    set(glfw_LIBRARIES ${VIAMD_STDLIBS})
endif()

set(OSX_BUNDLE "")
if (APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
    if (VIAMD_CREATE_MACOSX_BUNDLE)
        set(OSX_BUNDLE "MACOSX_BUNDLE")
    endif()
endif()

if (CMAKE_C_COMPILER_ID STREQUAL "GNU") # GCC
    set(VIAMD_FLAGS ${VIAMD_FLAGS} "-fms-extensions") #Silence pesky warnings of anonymous structs in vector types
endif()

file(GLOB APP_FILES src/application/*)
file(GLOB SRC_FILES src/*)
file(GLOB GFX_FILES src/gfx/*)
file(GLOB_RECURSE SHADER_FILES src/shaders/*)

source_group("app" FILES ${APP_FILES})
source_group("gfx" FILES ${GFX_FILES})
source_group("shaders" FILES ${SHADER_FILES})

add_executable(viamd ${OSX_BUNDLE} ${SRC_FILES} ${APP_FILES} ${GFX_FILES} ${SHADER_FILES})

# Default to linking statically
if (VIAMD_LINK_STDLIB_STATIC)
    set_property(TARGET viamd PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
else()
    set_property(TARGET viamd PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

set(VIAMD_BIN_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
if (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    set(VIAMD_BIN_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>")
endif()

create_copy_resource_dir_target(viamd_copy_shaders  "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders"  "${VIAMD_BIN_DIR}/shaders")
create_copy_resource_dir_target(viamd_copy_datasets "${CMAKE_CURRENT_SOURCE_DIR}/datasets" "${VIAMD_BIN_DIR}/datasets")

add_dependencies(viamd viamd_copy_shaders viamd_copy_datasets)

target_compile_definitions(viamd PRIVATE
	VIAMD_DATASET_DIR=\"datasets\"
	VIAMD_SHADER_DIR=\"shaders\"
	VIAMD_SCREENSHOT_DIR=\"screenshots\")

# We just hijack the warning and compile flags from mdlib
target_compile_options(viamd PRIVATE ${VIAMD_FLAGS} $<$<CONFIG:Debug>:${VIAMD_FLAGS_DEB}> $<$<CONFIG:Release>:${VIAMD_FLAGS_REL}>)

target_compile_features(viamd PRIVATE cxx_std_20)

target_include_directories(viamd
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/gl3w
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/enkiTS/src
)

set_target_properties(viamd PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${VIAMD_BIN_DIR}")

target_link_libraries(viamd
    glfw
    imgui
    nativefiledialog
    mdlib
    stb
    enkiTS
    ImGuiColorTextEdit
    implot
    atomic_queue
    imgui_notify
    ${VIAMD_STDLIBS}
)
