### TINYEXPR ###
project(libtinyexpr LANGUAGES C)

add_library(tinyexpr ${CMAKE_CURRENT_SOURCE_DIR}/tinyexpr.c)

target_include_directories(tinyexpr
	PUBLIC
		$<INSTALL_INTERFACE:/>
		$<BUILD_INTERFACE:
			${CMAKE_CURRENT_SOURCE_DIR}/;
		>
)