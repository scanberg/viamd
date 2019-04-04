cmake_minimum_required(VERSION 3.5)

project(benchmark)

add_executable(benchmark main.cpp)

target_compile_features(benchmark PRIVATE cxx_std_17)

target_compile_definitions(benchmark PRIVATE
    VIAMD_DATA_DIR=\"${PROJECT_SOURCE_DIR}/../data\")

target_link_libraries(benchmark mdutils)