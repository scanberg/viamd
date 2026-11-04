### IMGUIZMO: 3D manipulation gizmos drawn through ImGui ###
# ImGuizmo.cpp and ImGuizmo.h of https://github.com/CedricGuillemet/ImGuizmo at 105e4fde45e4e74aae8973d98750820ef7aca7e2 (MIT, LICENSE)
project(imguizmo LANGUAGES CXX)

add_library(imguizmo
    ${CMAKE_CURRENT_SOURCE_DIR}/ImGuizmo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ImGuizmo.h
)

target_include_directories(imguizmo
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(imguizmo PUBLIC imgui)
