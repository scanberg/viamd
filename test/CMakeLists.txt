cmake_minimum_required(VERSION 3.15)

project(viamd_test VERSION 0.1 LANGUAGES C CXX)

# viamd is one executable with no library to link against, so the translation units under test are
# compiled into the test binary directly. Only units which are free of ImGui, OpenGL and
# ApplicationState can go here - see README.md for what that means.
#
# The split follows the one the main CMakeLists already makes: the application's own sources are
# listed, because the set of them is part of its structure; components are discovered, because a
# component is meant to be added by dropping in a folder without editing anything central. A test
# build that named components would give that up - and it did, briefly.
set(UNIT_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/serialization_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/event.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/task_system.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/loader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/gfx/camera_utils.cpp
)

# Every .cpp in this folder is a test for the application's own units. Dropping one in is all it
# takes - CONFIGURE_DEPENDS has the build re-check the glob and re-configure when the set changes.
file(GLOB TEST_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# A component declares itself testable by its own layout, and nothing here names one:
#
#   src/components/foo/foo.cpp             the ImGui drawing and the event handling
#   src/components/foo/foo_core.{h,cpp}    the computation, linkable without an application
#   src/components/foo/tests/*.cpp         its tests
#
# The _core suffix is the declaration of separability, and it is the split README.md already asks
# for. A component with no _core unit contributes nothing, which is also what makes a component
# that is mid-split, or half-landed from another branch, inert rather than a broken build: its
# tests have nothing to link against, so they wait for the unit rather than failing the configure.
file(GLOB COMP_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../src/components/*)
foreach(COMP_DIR ${COMP_DIRS})
    if (IS_DIRECTORY ${COMP_DIR})
        file(GLOB COMP_CORE CONFIGURE_DEPENDS ${COMP_DIR}/*_core.cpp)
        if (COMP_CORE)
            file(GLOB COMP_TESTS CONFIGURE_DEPENDS ${COMP_DIR}/tests/*.cpp)
            list(APPEND UNIT_FILES ${COMP_CORE})
            list(APPEND TEST_FILES ${COMP_TESTS})
        endif()
    endif()
endforeach()

# The one component behind a build option, filtered exactly as the main build filters it.
if (NOT VIAMD_ENABLE_QM)
    list(FILTER UNIT_FILES EXCLUDE REGEX "/components/qm/")
    list(FILTER TEST_FILES EXCLUDE REGEX "/components/qm/")
endif()

add_executable(viamd_test ${TEST_FILES} ${UNIT_FILES})

# Same standard as the viamd target. The units under test and the mdlib headers they pull in use
# designated initializers, which gcc and clang accept as an extension in C++17 and MSVC does not.
target_compile_features(viamd_test PRIVATE cxx_std_20)

# The loader decides which reader gets a file partly by looking inside it, so its tests need real
# files. mdlib's fixtures are the ones already in the tree, and they are what the readers under
# that dispatch are themselves tested against.
target_compile_definitions(viamd_test PRIVATE UTEST_USE_OLD_QPC
    VIAMD_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../ext/mdlib/test_data"
    ${VIAMD_DEFINES} ${MD_DEFINES})
target_compile_options(viamd_test PRIVATE ${VIAMD_FLAGS} $<$<CONFIG:Debug>:${VIAMD_FLAGS_DEB}> $<$<CONFIG:Release>:${VIAMD_FLAGS_REL}>)

target_include_directories(viamd_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/components          # a component includes as <foo/foo_core.h>
    ${CMAKE_CURRENT_SOURCE_DIR}/../ext/mdlib/unittest      # utest.h, shared with mdlib's own suite
    ${CMAKE_CURRENT_SOURCE_DIR}/../ext/atomic_queue
    ${CMAKE_CURRENT_SOURCE_DIR}/../ext/enkiTS/src
)

target_link_libraries(viamd_test PRIVATE mdlib enkiTS)

add_test(NAME viamd_test COMMAND viamd_test)
