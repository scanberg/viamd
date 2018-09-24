### IMGUI ###
project(imgui)

file(GLOB_RECURSE WIDGET_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/widgets/*.h)
file(GLOB_RECURSE WIDGET_INLINES ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/widgets/*.inl)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/imconfig.h
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui.h
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui_internal.h
    ${WIDGET_HEADERS}
)
                   
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui_draw.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/imgui_demo.cpp
    ${WIDGET_INLINES}
)

add_library(imgui STATIC ${HEADERS} ${SOURCES})