### IMGUI ###
project(imgui LANGUAGES CXX)

set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui_widgets.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui_draw.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui_demo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui_tables.cpp
)

add_library(imgui ${SOURCES})

target_include_directories(imgui
	PUBLIC
		$<INSTALL_INTERFACE:/>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(imgui PRIVATE cxx_std_11)

if (MSVC)
    target_compile_options(imgui
	    PRIVATE /W3 /wd4152)
endif()
