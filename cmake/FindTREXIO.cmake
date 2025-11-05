# FindTREXIO.cmake
# Find the TREXIO library (https://github.com/TREX-CoE/trexio)
#
# This module defines:
#  TREXIO_FOUND - True if TREXIO was found
#  TREXIO_INCLUDE_DIRS - Include directories for TREXIO
#  TREXIO_LIBRARIES - Libraries to link for TREXIO
#  TREXIO_VERSION - Version of TREXIO found

# Try to find TREXIO using pkg-config first
find_package(PkgConfig QUIET)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(PC_TREXIO QUIET trexio)
endif()

# Find the include directory
find_path(TREXIO_INCLUDE_DIR
    NAMES trexio.h
    HINTS
        ${PC_TREXIO_INCLUDEDIR}
        ${PC_TREXIO_INCLUDE_DIRS}
    PATHS
        /usr/include
        /usr/local/include
        $ENV{TREXIO_ROOT}/include
        $ENV{HOME}/.local/include
)

# Find the library
find_library(TREXIO_LIBRARY
    NAMES trexio
    HINTS
        ${PC_TREXIO_LIBDIR}
        ${PC_TREXIO_LIBRARY_DIRS}
    PATHS
        /usr/lib
        /usr/local/lib
        /usr/lib64
        /usr/local/lib64
        $ENV{TREXIO_ROOT}/lib
        $ENV{HOME}/.local/lib
)

# Extract version if found
if(PC_TREXIO_VERSION)
    set(TREXIO_VERSION ${PC_TREXIO_VERSION})
elseif(TREXIO_INCLUDE_DIR)
    # Try to extract version from header file
    if(EXISTS "${TREXIO_INCLUDE_DIR}/trexio.h")
        # Try various version definition patterns
        file(STRINGS "${TREXIO_INCLUDE_DIR}/trexio.h" TREXIO_VERSION_LINES 
             REGEX "#define[ \t]+TREXIO_(PACKAGE_)?VERSION")
        
        # Try semantic versioning first (X.Y.Z)
        string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" TREXIO_VERSION "${TREXIO_VERSION_LINES}")
        
        # If not found, try simpler version patterns (X.Y or just X)
        if(NOT TREXIO_VERSION)
            string(REGEX MATCH "([0-9]+\\.[0-9]+)" TREXIO_VERSION "${TREXIO_VERSION_LINES}")
        endif()
        
        if(NOT TREXIO_VERSION)
            string(REGEX MATCH "([0-9]+)" TREXIO_VERSION "${TREXIO_VERSION_LINES}")
        endif()
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TREXIO
    REQUIRED_VARS TREXIO_LIBRARY TREXIO_INCLUDE_DIR
    VERSION_VAR TREXIO_VERSION
)

if(TREXIO_FOUND)
    set(TREXIO_LIBRARIES ${TREXIO_LIBRARY})
    set(TREXIO_INCLUDE_DIRS ${TREXIO_INCLUDE_DIR})
    
    if(NOT TARGET TREXIO::TREXIO)
        add_library(TREXIO::TREXIO UNKNOWN IMPORTED)
        set_target_properties(TREXIO::TREXIO PROPERTIES
            IMPORTED_LOCATION "${TREXIO_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${TREXIO_INCLUDE_DIR}"
        )
    endif()
endif()

mark_as_advanced(TREXIO_INCLUDE_DIR TREXIO_LIBRARY)
