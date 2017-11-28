project(kmc_binding)

pybind11_add_module(${PROJECT_NAME} ${kmc_sources})

target_include_directories(${PROJECT_NAME} PUBLIC ${spdlog_include_dir})
