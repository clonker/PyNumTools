//
// Created by mho on 11/23/17.
//

#include "quadrature.h"

namespace pnt {

namespace quad {

QuadratureCoefficients<2> q1{1., 2., {{1, 1}}};
QuadratureCoefficients<3> q2{1., 3., {{1, 4, 1}}};
QuadratureCoefficients<4> q3{3., 8., {{1, 3, 3, 1}}};
QuadratureCoefficients<5> q4{2., 45., {{7, 32, 12, 32, 7}}};
QuadratureCoefficients<6> q5{5., 288., {{19, 75, 50, 50, 75, 19}}};
QuadratureCoefficients<7> q6{1., 140., {{41, 216, 27, 272, 27, 216, 41}}};
QuadratureCoefficients<8> q7{7., 17280., {{751, 3577, 1323, 2989, 2989, 1323, 3577, 751}}};
QuadratureCoefficients<9> q8{4., 14175., {{989, 5888, -928, 10496, -4540, 10496, -928, 5888, 989}}};
QuadratureCoefficients<10> q9{9., 89600., {{2857, 15741, 1080, 19344, 5778, 5778, 19344, 1080, 15741, 2857}}};
QuadratureCoefficients<11> q10{5., 299376., {{16067, 106300, -48525, 272400, -260550, 427368, -260550, 272400,
                                              -48525, 106300, 16067}}};
QuadratureCoefficients<12> q11{11., 87091200., {{2171465, 13486539, -3237113, 25226685, -9595542,
                                                 15493566, 15493566, -9595542, 25226685, -3237113,
                                                 13486539, 2171465}}};
QuadratureCoefficients<13> q12{1., 5255250., {{1364651, 9903168, -7587864, 35725120, -51491295,
                                               87516288, -87797136, 87516288, -51491295, 35725120,
                                               -7587864, 9903168, 1364651}}};
QuadratureCoefficients<14> q13{13., 402361344000., {{8181904909, 56280729661, -31268252574,
                                                     156074417954, -151659573325, 206683437987,
                                                     -43111992612, -43111992612, 206683437987,
                                                     -151659573325, 156074417954, -31268252574,
                                                     56280729661, 8181904909}}};
QuadratureCoefficients<15> q14{7., 2501928000., {{90241897, 710986864, -770720657, 3501442784,
                                                  -6625093363, 12630121616, -16802270373, 19534438464,
                                                  -16802270373, 12630121616, -6625093363, 3501442784,
                                                  -770720657, 710986864, 90241897}}};

template<> QuadratureCoefficients<2> weights() { return q1; }
template<> QuadratureCoefficients<3> weights() { return q2; }
template<> QuadratureCoefficients<4> weights() { return q3; }
template<> QuadratureCoefficients<5> weights() { return q4; }
template<> QuadratureCoefficients<6> weights() { return q5; }
template<> QuadratureCoefficients<7> weights() { return q6; }
template<> QuadratureCoefficients<8> weights() { return q7; }
template<> QuadratureCoefficients<9> weights() { return q8; }
template<> QuadratureCoefficients<10> weights() { return q9; }
template<> QuadratureCoefficients<11> weights() { return q10; }
template<> QuadratureCoefficients<12> weights() { return q11; }
template<> QuadratureCoefficients<13> weights() { return q12; }
template<> QuadratureCoefficients<14> weights() { return q13; }
template<> QuadratureCoefficients<15> weights() { return q14; }

double newton_cotes(const array_type &y, double dx, std::size_t n, long axis) {
    if (n < 1 || n > 5) {
        throw std::invalid_argument("the newton-cotes n-point rule is only implemented for 1 <= n <= 5.");
    }
    auto ndim = y.ndim();
    if (axis == -1) {
        axis = ndim - 1;
    }
    auto last_dx = dx;
    auto first_dx = dx;
    auto N = y.shape(axis);
    if (N % 2 == 0) {

    } else {

    }
}

}

}