### BASE91 ###
project(libbase91 LANGUAGES C)

add_library(base91 ${CMAKE_CURRENT_SOURCE_DIR}/src/base91.c)

target_include_directories(base91
	PUBLIC
		$<INSTALL_INTERFACE:/>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
)