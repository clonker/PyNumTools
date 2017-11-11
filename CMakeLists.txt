cmake_minimum_required(VERSION 2.8.12)
project(cmake_example)

add_subdirectory(lib/pybind11)

set(spdlog_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include")
set(eigen_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen/Eigen")


