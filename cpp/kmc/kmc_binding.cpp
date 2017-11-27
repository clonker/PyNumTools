/**
 * @file kmc_binding.cpp
 * @brief Python bindings for functions used in Kinetic Monte Carlo algorithm
 * @author chrisfroe
 * @date 27.11.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using kmc_result_array = py::array_t<std::uint32_t, py::array::c_style>;
using kmc_state_array = py::array_t<std::uint32_t, py::array::c_style>;
using kmc_times_array = py::array_t<double, py::array::c_style>;

static void convert_kmc(kmc_result_array &result, const kmc_times_array &times, const kmc_times_array &times_list, const kmc_state_array &state_list) {
    std::size_t state = 0;

    auto nstates = (std::size_t) state_list.shape()[0];
    auto nframes = (std::size_t) result.shape()[0];
    auto nboxes = (std::size_t) result.shape()[1];
    auto nspecies = (std::size_t) result.shape()[2];

    for(std::size_t ix = 0; ix < nframes; ++ix) {
        auto t = times.at(ix);
        if(t <= times_list.at(state)) {
            for (std::size_t s = 0; s < nspecies; ++s) {
                for (std::size_t b = 0; b < nboxes; ++b) {
                    result.mutable_at(ix, b, s) = state_list.at(state, b, s);
                }
            }
        } else {
            while(state < nstates && t > times_list.at(state)) {
                ++state;
            }
            for (std::size_t s = 0; s < nspecies; ++s) {
                for (std::size_t b = 0; b < nboxes; ++b) {
                    result.mutable_at(ix, b, s) = state_list.at(state, b, s);
                }
            }
        }
    }
}

PYBIND11_MODULE(kmc_binding, m) {
    using namespace py::literals;
    m.def("convert_kmc", &convert_kmc);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

}
