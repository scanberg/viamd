cmake_minimum_required(VERSION 3.0)

project(viamd)

find_package(OpenGL REQUIRED)

if(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /wd4201")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_subdirectory(ext)

include_directories(src
                    ext/glfw/include
                    ext/gl3w/include
                    ext/glm
                    ext/imgui
                    ext/ppk_assert
                    ext/libxdrfile/src
                    ext/nativefiledialog/src/include)

file(GLOB_RECURSE PROJECT_HEADERS src/*.h)
file(GLOB_RECURSE PROJECT_SOURCES src/*.cpp
                                  src/*.inl)

file(GLOB PROJECT_CONFIGS CMakeLists.txt
                          Readme.md
                         .gitattributes
                         .gitignore
                         .gitmodules)

source_group("Headers" FILES ${PROJECT_HEADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})

add_definitions(-DGLFW_INCLUDE_NONE
                -DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")

add_executable(${PROJECT_NAME} ${PROJECT_HEADERS} ${PROJECT_SOURCES})

target_link_libraries(${PROJECT_NAME}
                      glfw
                      gl3w
                      imgui
                      ppk_assert
                      libxdrfile
                      nativefiledialog
                      ${GLFW_LIBRARIES}
                      ${OPENGL_LIBRARY})

#set_target_properties(${PROJECT_NAME} PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME})

INCLUDE(InstallRequiredSystemLibraries)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "VIAMD")
SET(CPACK_PACKAGE_VENDOR "LIU")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "1")
#SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
#IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
#  SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
#  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
#  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} My Famous Project")
#  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.my-project-home-page.org")
#  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.my-personal-home-page.com")
#  SET(CPACK_NSIS_CONTACT "me@my-personal-home-page.com")
#  SET(CPACK_NSIS_MODIFY_PATH ON)
#ELSE(WIN32 AND NOT UNIX)
#  SET(CPACK_STRIP_FILES "bin/MyExecutable")
#  SET(CPACK_SOURCE_STRIP_FILES "")
#ENDIF(WIN32 AND NOT UNIX)
SET(CPACK_PACKAGE_EXECUTABLES "viamd" "viamd")
INCLUDE(CPack)
