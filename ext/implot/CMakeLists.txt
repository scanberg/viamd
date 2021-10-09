### IMGUI PLOTTING LIBRARY ###
project(implot LANGUAGES CXX)

set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/implot.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implot_demo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implot_items.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implot.h
    ${CMAKE_CURRENT_SOURCE_DIR}/implot_internal.h
)

add_library(implot ${SOURCES})

target_include_directories(implot
	PUBLIC
		$<INSTALL_INTERFACE:/>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(implot PRIVATE cxx_std_11)

target_link_libraries(implot PRIVATE imgui ${VIAMD_STDLIBS})