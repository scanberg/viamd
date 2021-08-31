### IMGUI COLOR TEXT EDITOR ###
project(ImGuiColorTextEdit LANGUAGES CXX)

set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/TextEditor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TextEditor.h
)

add_library(ImGuiColorTextEdit ${SOURCES})

target_include_directories(ImGuiColorTextEdit
	PUBLIC
		$<INSTALL_INTERFACE:/>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(ImGuiColorTextEdit imgui)