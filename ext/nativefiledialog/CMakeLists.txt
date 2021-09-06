### NATIVEFILEDIALOG ###
project(nativefiledialoglib LANGUAGES C)

cmake_minimum_required(VERSION 3.10)
                  
set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_common.c)
set(LIB "")
set(INC "")
set(DEF "")
set(OPT "")

if(UNIX AND NOT APPLE)
    option(NATIVEFILEDIALOG_USE_GTK "Use GTK for window handling, fallback is Zenity" ON)
    if (NATIVEFILEDIALOG_USE_GTK)
	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
	message(${CMAKE_MODULE_PATH})
        #include(cmake/FindGTK3.cmake)
        #find_package(GTK3 REQUIRED gtk)
        find_package(PkgConfig REQUIRED)
        
        #PKG_CHECK_MODULES(GTK REQUIRED gtk+-3.0)
        #pkg_check_modules(GTK REQUIRED gtk+-3.0)
        find_package(GTK3 3.0 REQUIRED gtk)

        set(SRC
            ${SRC}
            ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_gtk.c)
        set(INC ${GTK3_INCLUDE_DIRS})
        set(LIB ${GTK3_LIBRARIES})
        set(DEF ${GTK3_CFLAGS_OTHER})
    else()
        set(SRC
            ${SRC}
            ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_zenity.c)
    endif()
elseif(APPLE)
    find_library(APPKIT_LIBRARY AppKit)
    set(SRC
        ${SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_cocoa.m)
    set(LIB ${APPKIT_LIBRARY})
else() #windows
    set(SRC
        ${SRC}
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nfd_win.cpp)
		set(DEF "_CRT_SECURE_NO_WARNINGS")
	if(MSVC)
		set(OPT "/W3" "/wd4245")
	endif()
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

target_compile_options(nativefiledialog
	PRIVATE ${OPT})

target_compile_definitions(nativefiledialog
    PRIVATE ${DEF})

target_link_libraries(nativefiledialog
	PRIVATE ${LIB})
