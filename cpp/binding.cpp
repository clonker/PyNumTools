/**
 * << detailed description >>
 *
 * @file binding.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.11.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "lgmres.h"

namespace py = pybind11;

PYBIND11_MODULE(pynumtools_binding, m) {
    using namespace py::literals;

    m.def("lgmres", &pnt::lgmres::sparse::compute,
          "A"_a,
          "b"_a,
          "x0"_a = pnt::Vec(0),
          "tol"_a = 1e-5,
          "maxiter"_a = 1000,
          "M"_a = pnt::Matrix(0, 0),
          "inner_m"_a = 30,
          "outer_k"_a = 3,
          "outer_v"_a = std::vector<std::tuple<pnt::Vec, pnt::Vec>>(),
          "store_outer_Av"_a = true,
          "prepend_outer_Av"_a = false);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

}
