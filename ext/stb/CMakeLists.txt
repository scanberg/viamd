### stb ###
project(stb LANGUAGES C)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/stb_image.h
)
                  
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stb_image.c
)

add_library(stb ${HEADERS} ${SOURCES})

target_include_directories(stb
	PUBLIC
		$<INSTALL_INTERFACE:/include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
)