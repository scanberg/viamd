cmake_minimum_required(VERSION 3.5)

project(unittests)

add_executable(unittests main.cpp)

target_compile_features(unittests PRIVATE cxx_std_20)
target_compile_definitions(unittests PRIVATE UNITTEST_DATA_DIR=\"${PROJECT_SOURCE_DIR}/data\")
target_link_libraries(unittests mdutils mold)