cmake_minimum_required(VERSION 3.5)

project(unittests)

add_executable(unittests main.cpp)

target_compile_options(unittests PRIVATE -W4)
target_compile_features(unittests PRIVATE cxx_std_17)

target_link_libraries(unittests viamdutils)