//
// Created by mho on 11/23/17.
//

#include "quadrature.h"

namespace pnt {

namespace quad {

constexpr std::array<double, 2> weights1 {{.5, .5}};
constexpr std::array<double, 3> weights2 {{1./3., 4./3., 1./3.}};
constexpr std::array<double, 4> weights3 {{3./8., 9./8., 9./8., 3./8.}};
constexpr std::array<double, 5> weights4 {{14./45., 64./45., 8./15., 64./45., 14./45.}};
constexpr std::array<double, 6> weights5 {{95./288., 125./96., 125./144., 125./144., 125./96., 95./288.}};


double newton_cotes(const array_type &y, double dx, std::size_t n, long axis) {
    if (n < 1 || n > 5){
        throw std::invalid_argument("the newton-cotes n-point rule is only implemented for 1 <= n <= 5.");
    }
    auto ndim = y.ndim();
    if(axis == -1) {
        axis = ndim -1;
    }
    auto last_dx = dx;
    auto first_dx = dx;
    auto N = y.shape(axis);
    if(N%2 == 0) {

    } else {

    }
}

}

}