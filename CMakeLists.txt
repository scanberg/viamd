cmake_minimum_required(VERSION 3.5)

project(viamd)

if(APPLE)
	set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

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
add_subdirectory(ext/viamdutils)
add_subdirectory(unittest)

file(GLOB PLATFORM_FILES src/platform/*)

source_group("platform" FILES ${PLATFORM_FILES})

add_definitions(-DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
add_definitions(-DDATA_DIR=\"${PROJECT_SOURCE_DIR}/data\")

#add_executable(viamd MACOSX_BUNDLE src/main.cpp ${PLATFORM_FILES})
add_executable(viamd src/main.cpp src/plot_extended.cpp src/range_slider.cpp ${PLATFORM_FILES})

target_compile_options(viamd PRIVATE -W4)
target_compile_features(viamd PRIVATE cxx_std_17)

target_include_directories(viamd
    PUBLIC 
        $<INSTALL_INTERFACE:src>    
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(viamd
    glfw
    imgui
    nativefiledialog
    viamdutils
)