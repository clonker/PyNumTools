project(binding)

pybind11_add_module(pynumtools_binding ${sources})

target_include_directories(pynumtools_binding PUBLIC ${spdlog_include_dir} ${eigen_include_dir})

target_link_libraries(pynumtools_binding PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

target_compile_definitions(pynumtools_binding PUBLIC EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
