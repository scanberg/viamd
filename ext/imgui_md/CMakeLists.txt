### imgui_md: Markdown rendering for Dear ImGui on top of MD4C ###
# Adapted from https://github.com/mekhontsev/imgui_md (MIT), see imgui_md.h for the list of changes.
project(imgui_md LANGUAGES CXX)

add_library(imgui_md STATIC imgui_md.cpp imgui_md.h)

target_include_directories(imgui_md
	PUBLIC
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_compile_features(imgui_md PRIVATE cxx_std_17)
target_link_libraries(imgui_md PUBLIC imgui md4c)
