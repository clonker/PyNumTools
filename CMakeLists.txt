cmake_minimum_required(VERSION 2.8.12)
project(pynumtools)

enable_language(C CXX)

find_package(BLAS)
find_package(LAPACK)

add_subdirectory(lib/pybind11)

set(spdlog_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include")
set(eigen_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")

set(sources ${CMAKE_CURRENT_SOURCE_DIR}/cpp/binding.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cpp/lgmres.h ${CMAKE_CURRENT_SOURCE_DIR}/cpp/lgmres.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cpp/util.h
        ${CMAKE_CURRENT_SOURCE_DIR}/cpp/quadrature.h ${CMAKE_CURRENT_SOURCE_DIR}/cpp/quadrature.cpp)

add_subdirectory(pynumtools)
