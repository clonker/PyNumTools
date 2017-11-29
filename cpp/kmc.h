/**
 * @file kmc.h
 * @brief Functions used in Kinetic Monte Carlo algorithm
 * @author chrisfroe
 * @date 28.11.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace kmc {

using kmc_result_array = py::array_t<std::uint32_t, py::array::c_style>;
using kmc_state_array = py::array_t<std::uint32_t, py::array::c_style>;
using kmc_times_array = py::array_t<double, py::array::c_style>;

void
convert_to_timeseries(kmc_result_array &result, const kmc_times_array &times, const kmc_times_array &times_list,
                      const kmc_state_array &state_list);

}