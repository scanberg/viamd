cmake_minimum_required(VERSION 3.10)

INCLUDE(ext/mdlib/cmake/common.cmake)

project(viamd)

set_property(GLOBAL PROPERTY USE_FOLDERS OFF)

# GLOBAL OPTIONS
option(BUILD_SHARED_LIBS OFF)
option(USE_MSVC_RUNTIME_LIBRARY_DLL OFF)

# VIAMD OPTIONS
option(VIAMD_USE_ADDRESS_SANITIZER "Use Address sanitizer" OFF)
option(VIAMD_CREATE_MACOSX_BUNDLE "Build a macosx bundle instead of just an executable" OFF)
option(VIAMD_LINK_STDLIB_STATIC "Link against stdlib statically" ON)

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

set(VIAMD_STDLIBS)
if (MD_LINK_STDLIB_STATIC)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(VIAMD_STDLIBS "-static-libgcc -static-libstdc++")
    endif()
    set(glfw_LIBRARIES ${VIAMD_STDLIBS})
endif()

add_subdirectory(ext/mdlib)
add_subdirectory(ext/glfw)
add_subdirectory(ext/imgui)
add_subdirectory(ext/nativefiledialog)
add_subdirectory(ext/stb)
add_subdirectory(ext/enkiTS)
add_subdirectory(ext/ImGuiColorTextEdit)
add_subdirectory(ext/implot)
add_subdirectory(ext/atomic_queue)

set(OSX_BUNDLE "")
if (APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
    if (VIAMD_CREATE_MACOSX_BUNDLE)
        set(OSX_BUNDLE "MACOSX_BUNDLE")
    endif()
endif()


file(GLOB APP_FILES src/application/*)
file(GLOB SRC_FILES src/*)
file(GLOB GFX_FILES src/gfx/*)

source_group("app" FILES ${APP_FILES})
source_group("gfx" FILES ${GFX_FILES})

add_executable(viamd ${OSX_BUNDLE} ${SRC_FILES} ${APP_FILES} ${GFX_FILES})


set(VIAMD_BIN_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
if (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    set(VIAMD_BIN_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>")
endif()

create_copy_resource_dir_target(viamd_copy_shaders  "${CMAKE_CURRENT_SOURCE_DIR}/shaders"  "${VIAMD_BIN_DIR}/shaders")
create_copy_resource_dir_target(viamd_copy_datasets "${CMAKE_CURRENT_SOURCE_DIR}/datasets" "${VIAMD_BIN_DIR}/datasets")

add_dependencies(viamd viamd_copy_shaders viamd_copy_datasets)

target_compile_definitions(viamd PRIVATE
	VIAMD_DATASET_DIR=\"datasets\"
	VIAMD_SHADER_DIR=\"shaders\"
	VIAMD_SCREENSHOT_DIR=\"screenshots\")

# We just hijack the warning and compile flags from mdlib
target_compile_options(viamd PRIVATE ${MD_WARNING_FLAGS} ${MD_COMMON_COMPILE_FLAGS})
add_compile_options_config(DEBUG ${MD_DEBUG_COMPILE_FLAGS})
add_compile_options_config(RELEASE ${MD_RELEASE_COMPILE_FLAGS})

if (APPLE)
target_compile_features(viamd PRIVATE cxx_std_17)
else()
target_compile_features(viamd PRIVATE cxx_std_20)
endif()

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
    ${VIAMD_STDLIBS}
)
