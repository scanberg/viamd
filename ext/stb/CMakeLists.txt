### stb ###
project(stb LANGUAGES C)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/stb_image.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/stb_image_write.h
)
                  
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stb_image.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/stb_image_write.c 
)

add_library(stb ${HEADERS} ${SOURCES})

target_include_directories(stb
	PUBLIC
		$<INSTALL_INTERFACE:/include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
)

if(MSVC)
	target_compile_options(stb
		PRIVATE /W3 /wd4100)
endif()