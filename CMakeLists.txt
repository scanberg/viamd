cmake_minimum_required(VERSION 3.20)

project(viamd VERSION 0.1.15 LANGUAGES CXX)

INCLUDE(ext/mdlib/cmake/common.cmake)

if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

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
option(VIAMD_IMGUI_ENABLE_VIEWPORTS "Enable ImGui Viewports" OFF)
option(VIAMD_IMGUI_ENABLE_DOCKSPACE "Enable ImGui Dockspace in main window" OFF)
option(VIAMD_CREATE_MACOSX_BUNDLE "Build a macosx bundle instead of just an executable" OFF)
option(VIAMD_LINK_STDLIB_STATIC "Link against stdlib statically" ${MD_LINK_STDLIB_STATIC})
set(VIAMD_FRAME_CACHE_SIZE_MB "2048" CACHE STRING "Reserved frame cache size in Megabytes")
set(VIAMD_NUM_WORKER_THREADS "8" CACHE STRING "Number of worker threads (Decrease if you run out of memory during evaluation)")

# Copy many of the fields from mdlib
set(VIAMD_STDLIBS)
set(VIAMD_FLAGS ${MD_FLAGS})
set(VIAMD_FLAGS_DEB ${MD_FLAGS_DEB})
set(VIAMD_FLAGS_REL ${MD_FLAGS_REL})
set(VIAMD_LINK_FLAGS ${MD_LINK_FLAGS})
set(VIAMD_LINK_FLAGS_DEB ${MD_LINK_FLAGS_DEB})
set(VIAMD_LINK_FLAGS_REL ${MD_LINK_FLAGS_REL})

if (MD_LINK_STDLIB_STATIC)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(VIAMD_STDLIBS "-static-libgcc -static-libstdc++")
    endif()
    set(glfw_LIBRARIES ${VIAMD_STDLIBS})
endif()

set(OSX_BUNDLE "")
if (APPLE)
    # Remove once we have proper support for Apple Silicon?
    SET(CMAKE_OSX_ARCHITECTURES "x86_64")
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

set(SHADER_FILES 
    src/shaders/highlight.frag
    src/shaders/scale_hsv.frag
    src/shaders/deferred_shading.frag
    src/shaders/temporal.frag
    src/shaders/tonemap/passthrough.frag
    src/shaders/tonemap/exposure_gamma.frag
    src/shaders/tonemap/uncharted.frag
    src/shaders/tonemap/aces.frag
    src/shaders/tonemap/fast_reversible.frag
    src/shaders/dof/dof_half_res_prepass.frag
    src/shaders/dof/dof.frag
    src/shaders/velocity/blit_velocity.frag
    src/shaders/velocity/blit_tilemax.frag
    src/shaders/velocity/blit_neighbormax.frag
    src/shaders/volume/raycaster.vert
    src/shaders/volume/raycaster.frag
    src/shaders/ssao/ssao.frag
    src/shaders/ssao/blur.frag
)

create_resources("${SHADER_FILES}" "gen/shaders.inl")

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
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(VIAMD_BIN_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>")
endif()

create_copy_resource_dir_target(viamd_copy_datasets "${CMAKE_CURRENT_SOURCE_DIR}/datasets"      "${VIAMD_BIN_DIR}/datasets")

add_dependencies(viamd viamd_copy_datasets)

target_compile_definitions(viamd PRIVATE
	VIAMD_DATASET_DIR=\"datasets\"
	VIAMD_SCREENSHOT_DIR=\"screenshots\"
    VIAMD_NUM_WORKER_THREADS=${VIAMD_NUM_WORKER_THREADS}
    VIAMD_FRAME_CACHE_SIZE=${VIAMD_FRAME_CACHE_SIZE_MB}
    VIAMD_IMGUI_ENABLE_VIEWPORTS=$<BOOL:${VIAMD_IMGUI_ENABLE_VIEWPORTS}>
    VIAMD_IMGUI_ENABLE_DOCKSPACE=$<BOOL:${VIAMD_IMGUI_ENABLE_DOCKSPACE}>
    ${MD_DEFINES}
)

if (WIN32)
    target_sources(viamd PRIVATE icon.rc)
endif()

# We just hijack the warning and compile flags from mdlib
target_compile_options(viamd PRIVATE ${VIAMD_FLAGS} $<$<CONFIG:Debug>:${VIAMD_FLAGS_DEB}> $<$<CONFIG:Release>:${VIAMD_FLAGS_REL}>)

target_compile_features(viamd PRIVATE cxx_std_20)

target_include_directories(viamd
    PRIVATE
        src
        gen
        ext/gl3w
        ext/enkiTS/src
)

set_target_properties(viamd PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${VIAMD_BIN_DIR}")

target_link_options(viamd PRIVATE ${VIAMD_LINK_FLAGS} $<$<CONFIG:Debug>:${VIAMD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${VIAMD_LINK_FLAGS_REL}>)
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
