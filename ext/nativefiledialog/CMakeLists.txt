### NATIVEFILEDIALOG ###
project(nativefiledialoglib LANGUAGES C)
                  
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_common.c)
set(LIB "")
set(INC "")
set(DEF "")

if(UNIX AND NOT APPLE)
    option(NATIVEFILEDIALOG_USE_GTK "Use GTK for window handling, fallback is Zenity" ON)
    if (NATIVEFILEDIALOG_USE_GTK)
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
	message(${CMAKE_MODULE_PATH})
        #include(cmake/FindGTK3.cmake)
        find_package(GTK3 REQUIRED gtk)
        set(SRC
            ${SRC}
            ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_gtk.c)
        set(INC ${GTK3_INCLUDE_DIRS})
        set(LIB ${GTK3_LIBRARIES})
    else()
        set(SRC
            ${SRC}
            ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_zenity.c)
    endif()
endif()

if(APPLE)
    set(SRC
        ${SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_cocoa.m)
    #set(LIB \"-framework AppKit\")
endif()

if(WIN32)
    set(SRC
        ${SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_win.cpp)
    set(DEF "_CRT_SECURE_NO_WARNINGS")
endif()

add_library(nativefiledialog ${SRC})

target_include_directories(nativefiledialog
	PUBLIC
		$<INSTALL_INTERFACE:/src/include>
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
	PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src
        ${INC}
)

target_compile_definitions(nativefiledialog
    PRIVATE ${DEF})

target_link_libraries(nativefiledialog ${LIB})
