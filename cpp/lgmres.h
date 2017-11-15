/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file lgmres.h
 * @brief << brief description >>
 * @author clonker
 * @date 11.11.17
 * @copyright GNU Lesser General Public License v3.0
 */
#ifndef PYNUMTOOLS_LGMRES_H
#define PYNUMTOOLS_LGMRES_H

#include "util.h"
#include <Eigen/QR>

namespace pnt {

namespace lgmres {

template<typename VecType>
struct ArnoldiResult {
    // QR
    Matrix Q;
    Matrix R;
    Eigen::HouseholderQR<Matrix> QR;
    // orthogonal projection  coefficients
    Matrix B;
    // columns of matrix V
    std::vector<VecType> vs;
    // columns of matrix Z
    std::vector<VecType> zs;
    // solution to || Hy - e_1||_2 = min!
    VecType y;
};


namespace sparse {

ArnoldiResult<SpVec> arnoldi(const system::MatVec<SpVec> &matvec, const SpVec &v0, std::size_t m,
                             double atol, const system::MatVec<SpVec> &rpsolve,
                             std::vector<std::tuple<Vec, Vec>> &outer_v, bool prependOuterV);

Vec compute(const SpMatrix &A, Eigen::Ref<Vec> b, Eigen::Ref<Vec> x0, double tol, std::size_t maxiter,
            const SpMatrix &M, std::size_t inner_m, std::size_t outer_k,
            std::vector<std::tuple<Vec, Vec>> &outer_v, bool storeOuterAv, bool prependOuterAv);

}
}

/*Vec lgmres(Eigen::Ref<Matrix> A, Eigen::Ref<Vec> b, Vec x0, double tol, std::size_t maxiter,
           Matrix M, std::size_t inner_m, std::size_t outer_k, std::vector<std::tuple<Vec, Vec>> &outer_v,
           bool storeOuterAv);*/

}


#endif //PYNUMTOOLS_LGMRES_H
