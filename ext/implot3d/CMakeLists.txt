### IMGUI 3D PLOTTING LIBRARY ###
project(implot3d LANGUAGES CXX)

set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/implot3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implot3d_demo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implot3d_items.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implot3d_meshes.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implot3d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/implot3d_internal.h
)

add_library(implot3d ${SOURCES})

target_include_directories(implot3d
	PUBLIC
		$<INSTALL_INTERFACE:/>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(implot3d PRIVATE cxx_std_11)

target_link_libraries(implot3d PRIVATE imgui ${VIAMD_STDLIBS})
