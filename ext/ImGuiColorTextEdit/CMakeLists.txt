### IMGUI COLOR TEXT EDITOR ###
# https://github.com/goossens/ImGuiColorTextEdit at 99fcb0d (2026-09-20), vendored unmodified.
# To update: copy TextEditor.h, TextEditor.cpp, LICENSE and README.md from a newer commit. Everything viamd specific
# (the mdscript language, markers, focus and hover handling) lives in src/script_editor.h/.cpp instead.
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

target_compile_features(ImGuiColorTextEdit PRIVATE cxx_std_20)

target_link_libraries(ImGuiColorTextEdit PRIVATE imgui ${VIAMD_STDLIBS})