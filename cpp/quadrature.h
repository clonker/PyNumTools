//
// Created by mho on 11/23/17.
//

#ifndef PYNUMTOOLS_QUADRATURE_H
#define PYNUMTOOLS_QUADRATURE_H

#include <pybind11/numpy.h>

namespace py = pybind11;

namespace pnt {

namespace quad {

using array_type = py::array_t<double, py::array::c_style>;

double newton_cotes(const array_type &y, double dx, std::size_t n=1, long axis = -1);

}

}

#endif //PYNUMTOOLS_QUADRATURE_H
