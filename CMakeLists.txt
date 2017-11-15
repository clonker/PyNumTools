cmake_minimum_required(VERSION 2.8.12)
project(pynumtools)

enable_language(C CXX)

find_package(BLAS)
find_package(LAPACK)

add_subdirectory(lib/pybind11)

set(spdlog_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include")
set(eigen_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen")

set(sources cpp/binding.cpp
        cpp/lgmres.h cpp/lgmres.cpp
        cpp/util.h)

pybind11_add_module(pynumtools_binding ${sources})

target_include_directories(pynumtools_binding PUBLIC ${spdlog_include_dir} ${eigen_include_dir})

target_link_libraries(pynumtools_binding PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

target_compile_definitions(pynumtools_binding PUBLIC EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
