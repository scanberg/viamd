
### NATIVEFILEDIALOG ###
project(nativefiledialoglib LANGUAGES C)
                  
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_common.c
)

IF(UNIX AND NOT APPLE)
    find_package(GTK2 REQUIRED gtk)
    set(SOURCES
        ${SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_gtk.c)
ENDIF()

IF(APPLE)
    set(SOURCES
        ${SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_cocoa.m)
ENDIF()

IF(WIN32)
    set(SOURCES
        ${SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_win.cpp)
ENDIF()

add_library(nativefiledialog ${SOURCES})

target_include_directories(nativefiledialog
	PUBLIC
		$<INSTALL_INTERFACE:/src/include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
)

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/include)

IF(APPLE)
    target_link_libraries(nativefiledialog "-framework AppKit")
ENDIF()
IF (UNIX AND NOT APPLE)
	target_link_libraries(nativefiledialog gtk)
ENDIF()
