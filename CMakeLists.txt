cmake_minimum_required(VERSION 2.8.12)
project(pynumtools_binding)

find_package(Git REQUIRED)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

enable_language(C CXX)
set(CMAKE_CXX_STANDARD 17)

if(SKBUILD)
    # Scikit-Build does not add your site-packages to the search path
    # automatically, so we need to add it _or_ the pybind11 specific directory
    # here.
    execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c
            "import pybind11; print(pybind11.get_cmake_dir())"
            OUTPUT_VARIABLE _tmp_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
    list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
else()
    execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c
            "import skbuild; from pathlib import Path; print(str((Path(skbuild.__file__).parents[0] / 'resources' / 'cmake').absolute()))"
            OUTPUT_VARIABLE _tmp_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT
    )
    list(APPEND CMAKE_MODULE_PATH "${_tmp_dir}")
endif()

set(spdlog_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/spdlog/include")
set(eigen_include_dir "${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen-git-mirror")
find_package(pybind11 REQUIRED)

if(MSVC)
    add_compile_options(/EHsc /bigobj)
endif()

set(sources ${CMAKE_CURRENT_SOURCE_DIR}/cpp/binding.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cpp/lgmres.h ${CMAKE_CURRENT_SOURCE_DIR}/cpp/lgmres.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cpp/util.h
        ${CMAKE_CURRENT_SOURCE_DIR}/cpp/quadrature.h ${CMAKE_CURRENT_SOURCE_DIR}/cpp/quadrature.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/cpp/kmc.cpp)

pybind11_add_module(${CMAKE_PROJECT_NAME} ${sources})

target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC ${spdlog_include_dir} ${eigen_include_dir})
install(TARGETS ${CMAKE_PROJECT_NAME} DESTINATION .)
