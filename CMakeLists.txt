cmake_minimum_required(VERSION 2.8.12)
project(pynumtools)

add_subdirectory(lib/pybind11)

set(spdlog_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include")
set(eigen_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")

set(sources cpp/binding.cpp
        cpp/lgmres.h cpp/lgmres.cpp
        cpp/util.h)

pybind11_add_module(pynumtools_binding ${sources})

target_include_directories(pynumtools_binding PUBLIC ${spdlog_include_dir} ${eigen_include_dir})
