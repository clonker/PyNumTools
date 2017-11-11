cmake_minimum_required(VERSION 2.8.12)
project(cmake_example)

add_subdirectory(lib/pybind11)

set(spdlog_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include")
set(eigen_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")

pybind11_add_module(pynumtools_binding cpp/binding.cpp)

target_include_directories(pynumtools_binding PUBLIC ${spdlog_include_dir} ${eigen_include_dir})
