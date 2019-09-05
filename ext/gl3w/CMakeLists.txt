### GL3W ###
project(gl3w)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(gl3w STATIC ${CMAKE_CURRENT_SOURCE_DIR}/GL/gl3w.h ${CMAKE_CURRENT_SOURCE_DIR}/gl3w.c)

target_include_directories(gl3w
	PUBLIC
		$<INSTALL_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
)

if(MSVC)
	target_compile_options(gl3w
		PRIVATE /W3 /wd4152)
endif()