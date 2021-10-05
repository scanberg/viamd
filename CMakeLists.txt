
cmake_minimum_required(VERSION 3.10)

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

option(VIAMD_USE_RELATIVE_PATHS ON)
option(VIAMD_PROFILE_COMPILE_TIMES OFF)

if (MSVC)
    message(STATUS "MSVC configuration was used")
    add_compile_options(/MT /W4 /wd4201 /wd4324 /wd26451 /MP /GR- /fp:fast /arch:AVX)# /arch:AVX2

	if (VIAMD_PROFILE_COMPILE_TIMES)	
        add_compile_options(/d2cgsummary)
		#add_compile_options(/Bt+ /d1reportTime /d2cgsummary)
		#add_link_options(/time+)
	endif()
	add_compile_options_config(RELEASE /MT /Oi /GS-)
    add_compile_options_config(DEBUG /MTd /ZI)
else()
	message(STATUS "non MSVC configuration was used")
    add_compile_options(
        -Wall -Wextra -Wno-unused-function -Wno-unused-parameter -Wno-missing-braces -Wno-missing-field-initializers -march=native #-maes
    )
	if (VIAMD_PROFILE_COMPILE_TIMES)
		add_compile_options(-ftime-trace)
	endif()
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-c99-extensions -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types")
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS OFF)

option(BUILD_SHARED_LIBS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_BUILD_DOCS OFF)
option(GLFW_INSTALL OFF)
option(GLFW_VULKAN_STATIC OFF)
option(ENKITS_BUILD_C_INTERFACE OFF)
option(ENKITS_BUILD_EXAMPLES OFF)
option(ENKITS_BUILD_SHARED OFF)
option(ENKITS_INSTALL OFF)
option(VIAMD_CREATE_MACOSX_BUNDLE "Build a macosx bundle instead of just an executable" OFF)
option(VIAMD_USE_RELATIVE_PATHS "Use relative paths for resource directories" OFF)
option(VIAMD_LINK_STDLIB_STATIC "Link against stdlib statically" ON)

if (VIAMD_LINK_STDLIB_STATIC)
    set(VIAMD_STDLIBS "")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(VIAMD_STDLIBS "-static-libgcc -static-libstdc++")
    endif()
    set(glfw_LIBRARIES ${VIAMD_STDLIBS})
endif()


add_subdirectory(ext/glfw)
add_subdirectory(ext/imgui)
add_subdirectory(ext/nativefiledialog)
add_subdirectory(ext/mold)
add_subdirectory(ext/stb)
add_subdirectory(ext/gl3w)
add_subdirectory(ext/enkiTS)
add_subdirectory(ext/ImGuiColorTextEdit)
add_subdirectory(ext/implot)
add_subdirectory(ext/atomic_queue)

file(GLOB APP_FILES src/application/*)
file(GLOB SRC_FILES src/*)
file(GLOB GFX_FILES src/gfx/*)

source_group("app" FILES ${APP_FILES})
source_group("gfx" FILES ${GFX_FILES})

add_executable(viamd ${OSX_BUNDLE} ${SRC_FILES} ${APP_FILES} ${GFX_FILES})

set(RESOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/")
if (VIAMD_USE_RELATIVE_PATHS)
    set(RESOURCE_DIR "")

	# COPY RESOURCE FILES TO EXECUTABLE DIRECTORY TO KEEP RELATIVE PATH
    # CANNOT COPY ENTIRE DIRECTORY -> WILL REMOVE ANY EXISTING DIRECTORY
    file(GLOB_RECURSE COPY_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} shaders/* images/* datasets/*)

    foreach(FILENAME ${COPY_FILES})
    set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}")
    set(DST "${CMAKE_BINARY_DIR}/$<CONFIG>/${FILENAME}")

    add_custom_command(
        TARGET viamd POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
        )
    endforeach(FILENAME)

	# COPY SHADER FILES TO EXECUTABLE DIRECTORY TO KEEP RELATIVE PATH
	#add_custom_command(
	#		TARGET viamd POST_BUILD
	#		COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_BINARY_DIR}/$<CONFIG>/shaders
    #       COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/images ${CMAKE_BINARY_DIR}/$<CONFIG>/images)
endif()



target_compile_definitions(viamd PRIVATE
	VIAMD_DATASET_DIR=\"${RESOURCE_DIR}datasets\"
	VIAMD_SHADER_DIR=\"${RESOURCE_DIR}shaders\"
	VIAMD_IMAGE_DIR=\"${RESOURCE_DIR}images\"
	VIAMD_SCREENSHOT_DIR=\"${RESOURCE_DIR}screenshots\")

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

target_link_libraries(viamd
    gl3w
    glfw
    imgui
    nativefiledialog
    mold
    stb
    enkiTS
    ImGuiColorTextEdit
    implot
    atomic_queue
    ${VIAMD_STDLIBS}
)
