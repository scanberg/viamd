
### NATIVEFILEDIALOG ###
project(nativefiledialog)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_common.h
)
                   
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

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/include)
add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

IF(APPLE)
    target_link_libraries(${PROJECT_NAME}
                        "-framework AppKit")
ENDIF()
