cmake_minimum_required(VERSION 3.20)

project(viamd VERSION 0.1.15 LANGUAGES CXX)

INCLUDE(ext/mdlib/cmake/common.cmake)

# If configured using single configuration generator, set default build type to Release
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  endif()
endif()

# GLOBAL OPTIONS
set_property(GLOBAL PROPERTY USE_FOLDERS OFF)
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

if (UNIX AND NOT APPLE)
    option(GLFW_BUILD_WAYLAND OFF)
endif()

# ENKI OPTIONS
option(ENKITS_BUILD_C_INTERFACE OFF)
option(ENKITS_BUILD_EXAMPLES OFF)
option(ENKITS_BUILD_SHARED OFF)
option(ENKITS_INSTALL OFF)

# VIAMD OPTIONS
option(VIAMD_ENABLE_DEFAULT_DATASET "Enable Default Dataset" ON)
option(VIAMD_IMGUI_ENABLE_VIEWPORTS "Enable ImGui Viewports" OFF)
option(VIAMD_IMGUI_ENABLE_DOCKSPACE "Enable ImGui Dockspace in main window" OFF)
option(VIAMD_CREATE_MACOSX_BUNDLE "Build a macosx bundle instead of just an executable" OFF)
option(VIAMD_LINK_STDLIB_STATIC "Link against stdlib statically" ON)
option(VIAMD_ENABLE_VELOXCHEM "Enable Veloxchem Module" OFF)
option(VIAMD_ENABLE_TREXIO "Enable TREXIO Module" OFF)
set(VIAMD_FRAME_CACHE_SIZE_MB "2048" CACHE STRING "Reserved frame cache size in Megabytes")
set(VIAMD_NUM_WORKER_THREADS "8" CACHE STRING "Number of worker threads (Decrease if you run out of memory during evaluation)")

# MDLIB OPTIONS
set(MD_LINK_STDLIB_STATIC ${VIAMD_LINK_STDLIB_STATIC} CACHE BOOL "" FORCE)

if (VIAMD_ENABLE_VELOXCHEM)
    set(MD_ENABLE_VLX ON CACHE BOOL "" FORCE)
    set(MD_VLX_BASIS_FOLDER "basis")
endif ()

if (VIAMD_ENABLE_TREXIO)
    set(MD_ENABLE_TREXIO ON CACHE BOOL "" FORCE)
endif ()

if (APPLE)
    # Remove once we have proper support for Apple Silicon?
    # SET(CMAKE_OSX_ARCHITECTURES "x86_64")
    if (VIAMD_CREATE_MACOSX_BUNDLE)
        set(OSX_BUNDLE "MACOSX_BUNDLE")
        set(MACOSX_BUNDLE_ICON_FILE viamd.icns)
        set(MACOSX_ICON ${CMAKE_CURRENT_SOURCE_DIR}/icon/viamd.icns)
        set_source_files_properties(${MACOSX_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    endif()
endif()

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

# Add tools subdirectory if TREXIO is enabled
if(VIAMD_ENABLE_TREXIO)
    add_subdirectory(tools/mdlib_allocator_reproducer)
    add_subdirectory(tools/trexio_loader_test)
endif()

# Copy many of the fields from mdlib
set(VIAMD_STDLIBS)
set(VIAMD_FLAGS ${MD_FLAGS})
set(VIAMD_FLAGS_DEB ${MD_FLAGS_DEB})
set(VIAMD_FLAGS_REL ${MD_FLAGS_REL})
set(VIAMD_LINK_FLAGS ${MD_LINK_FLAGS})
set(VIAMD_LINK_FLAGS_DEB ${MD_LINK_FLAGS_DEB})
set(VIAMD_LINK_FLAGS_REL ${MD_LINK_FLAGS_REL})
set(VIAMD_DEFINES)

if (VIAMD_LINK_STDLIB_STATIC)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(VIAMD_STDLIBS "-static-libgcc -static-libstdc++")
    endif()
    set(glfw_LIBRARIES ${VIAMD_STDLIBS})
endif()

if (CMAKE_C_COMPILER_ID STREQUAL "GNU") # GCC
    set(VIAMD_FLAGS ${VIAMD_FLAGS} "-fms-extensions") #Silence pesky warnings of anonymous structs in vector types
endif()

file(GLOB SRC_FILES src/*.h src/*.cpp src/*.inl)
file(GLOB APP_FILES src/app/*.cpp src/app/*.h)
file(GLOB GFX_FILES src/gfx/*.cpp src/gfx/*.h)

set(COMP_FILES
    src/components/ramachandran/ramachandran.cpp
    src/components/shapespace/shapespace.cpp
    src/components/dataset/dataset.cpp
)

if (VIAMD_ENABLE_VELOXCHEM)
    set(COMP_FILES ${COMP_FILES} src/components/veloxchem/veloxchem.cpp)
endif()

set(SHADER_FILES 
    src/shaders/highlight.frag
    src/shaders/scale_hsv.frag
    src/shaders/compose_deferred.frag
    src/shaders/temporal.frag
    src/shaders/luma.frag
    src/shaders/fxaa.frag
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
    src/shaders/volume/entryexit.vert
    src/shaders/volume/entryexit.frag
    src/shaders/volume/raycaster.frag
    src/shaders/ssao/ssao.frag
    src/shaders/ssao/blur.frag
)

# Bake shaders into a single header file
create_resources("${SHADER_FILES}" "gen/shaders.inl")

source_group("src" FILES ${SRC_FILES})
source_group("src/app" FILES ${APP_FILES})
source_group("src/gfx" FILES ${GFX_FILES})
source_group("shaders" FILES ${SHADER_FILES})
source_group("src/component" FILES ${COMP_FILES})

add_executable(viamd ${OSX_BUNDLE} ${MACOSX_ICON} ${SRC_FILES} ${APP_FILES} ${GFX_FILES} ${SHADER_FILES} ${COMP_FILES})

install(TARGETS viamd DESTINATION bin)

# The question is how this should be treated on Linux...
# For now we simply don't copy the datasets
# Install datasets
#install(FILES datasets/1ALA-500.pdb DESTINATION datasets)

if (UNIX AND NOT APPLE)
    # Linux metadata
    install(FILES viamd.desktop DESTINATION share/applications)
    install(FILES viamd.metainfo.xml DESTINATION share/metainfo)
	install(FILES icon/viamd.png DESTINATION share/icons)
endif()

if (WIN32)
    target_sources(viamd PRIVATE icon.rc)
    if (VIAMD_LINK_STDLIB_STATIC)
        set_property(TARGET viamd PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    else()
        set_property(TARGET viamd PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
    endif()
endif()

set(VIAMD_BIN_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(VIAMD_BIN_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>")
endif()

if (VIAMD_ENABLE_VELOXCHEM)
    create_copy_resource_dir_target(viamd_copy_basis "${CMAKE_CURRENT_SOURCE_DIR}/ext/mdlib/ext/veloxchem/basis"      "${VIAMD_BIN_DIR}/basis")
    add_dependencies(viamd viamd_copy_basis)
endif()

if (VIAMD_ENABLE_DEFAULT_DATASET)
    create_copy_resource_dir_target(viamd_copy_datasets "${CMAKE_CURRENT_SOURCE_DIR}/datasets"      "${VIAMD_BIN_DIR}/datasets")
    add_dependencies(viamd viamd_copy_datasets)
    set(VIAMD_DEFINES ${VIAMD_DEFINES} VIAMD_DATASET_DIR=\"datasets\")
endif()

target_compile_definitions(viamd PRIVATE
    ${VIAMD_DEFINES}
    VIAMD_NUM_WORKER_THREADS=${VIAMD_NUM_WORKER_THREADS}
    VIAMD_FRAME_CACHE_SIZE=${VIAMD_FRAME_CACHE_SIZE_MB}
    VIAMD_IMGUI_ENABLE_VIEWPORTS=$<BOOL:${VIAMD_IMGUI_ENABLE_VIEWPORTS}>
    VIAMD_IMGUI_ENABLE_DOCKSPACE=$<BOOL:${VIAMD_IMGUI_ENABLE_DOCKSPACE}>
    ${MD_DEFINES}
)

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

