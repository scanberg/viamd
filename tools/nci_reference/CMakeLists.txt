project(viamd_nci_reference LANGUAGES C CXX)

# A headless batch tool: the NCI category references of a VeloxChem calculation, as the NCI panel exports
# them. Like viamd_test it compiles the GUI-free _core units it needs directly.
set(_topology ${CMAKE_CURRENT_SOURCE_DIR}/../../src/components/topology)
add_executable(viamd_nci_reference
    ${CMAKE_CURRENT_SOURCE_DIR}/nci_reference.cpp
    ${_topology}/nci_surface_core.cpp
    ${_topology}/nci_reference_data_core.cpp
    ${_topology}/critical_point_distribution_core.cpp
)
target_compile_features(viamd_nci_reference PRIVATE cxx_std_20)
target_compile_definitions(viamd_nci_reference PRIVATE ${VIAMD_DEFINES} ${MD_DEFINES})
target_compile_options(viamd_nci_reference PRIVATE ${VIAMD_FLAGS} $<$<CONFIG:Debug>:${VIAMD_FLAGS_DEB}> $<$<CONFIG:Release>:${VIAMD_FLAGS_REL}>)
target_include_directories(viamd_nci_reference PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src
    ${CMAKE_CURRENT_SOURCE_DIR}/../../src/components
)
target_link_libraries(viamd_nci_reference PRIVATE mdlib)
