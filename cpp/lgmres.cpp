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
 * @file lgmres.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.11.17
 * @copyright GNU Lesser General Public License v3.0
 */
#include "lgmres.h"
#include "util.h"
#include <spdlog/fmt/fmt.h>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <iostream>
#include <spdlog/spdlog.h>

namespace pnt {

using QRMatrix = typename Eigen::HouseholderQR<Matrix>;

constexpr double eps = 2.2204460492503131e-16;

namespace {
auto qr_append(const Matrix &Q, const Matrix &R, const Matrix &u, std::size_t k) {
    // todo this can be further optimized, check
    // https://github.com/scipy/scipy/blob/master/scipy/linalg/_decomp_update.pyx.in
    Matrix A = Q * R;
    A.conservativeResize(Eigen::NoChange, A.cols() + 1);
    for(long ii = A.cols()-1; ii >= static_cast<long>(k)+1; --ii) {
        A.col(ii).swap(A.col(ii-1));
    }

    A.col(k) = u;
    return A.householderQr();
}
}

namespace lgmres {
namespace sparse {

ArnoldiResult<SpVec> arnoldi(const system::MatVec<SpVec> &matvec, const SpVec &v0, std::size_t m, double atol,
                             const system::MatVec<SpVec> &rpsolve, std::vector<std::tuple<Vec, Vec>> &outer_v,
                             bool prependOuterV) {
    auto log = spdlog::get("console");

    auto lpsolve = [](const auto &x) { return x; };

    ArnoldiResult<SpVec> result;

    result.vs = {v0};
    result.zs = {};
    result.y = {};

    m += outer_v.size();

    // orthogonal projection coeffs
    // result.B = Matrix::Zero(0, m);

    // QR
    result.Q = Matrix::Ones(1, 1);
    result.R = Matrix::Zero(1, 0);

    bool breakdown {false};

    // FGMRES Arnoldi process
    std::size_t j = 0;
    for (; j < m; ++j) {
        // L A Z = C B + V H
        Vec z, w;
        bool setW = false;
        if(prependOuterV && j < outer_v.size()) {
            z = std::get<0>(outer_v.at(j));
            w = std::get<1>(outer_v.at(j));
            setW = true;
        }  else if (prependOuterV && j == outer_v.size()) {
            z = rpsolve(v0);
            setW = false;
        } else if( !prependOuterV && j >= m - outer_v.size()) {
            z = std::get<0>(outer_v.at(j - (m - outer_v.size())));
            w = std::get<1>(outer_v.at(j - (m - outer_v.size())));
            setW = true;
        } else {
            z = rpsolve(result.vs.back());
            setW = false;
        }

        if (!setW) {
            w = lpsolve(matvec(z));
        } else {
            // copy!
            w = Vec(w);
        }

        auto wNorm = w.norm();

        // Orthogonalize against V
        Vec hcur = Vec::Zero(j+2);
        std::size_t i = 0;
        for(; i < result.vs.size(); ++i) {
            const auto &v = result.vs.at(i);
            auto alpha = v.dot(w);
            hcur[i] = alpha;
            w -= alpha * v;
        }
        hcur[i+1] = w.norm();
        if(std::abs(hcur[j+1]) > 1e-8) {
            auto alpha = 1./ hcur[j+1];
            w = alpha * w;
        }

        if(!(hcur[j+1] > eps * wNorm)) {
            // w essentially in the span of previous vectors,
            // or we have nans. Bail out after updating the QR solution
            breakdown = true;
        }

        result.vs.push_back(std::move(w));
        result.zs.push_back(std::move(z));

        // arnoldi LSQ problem

        Matrix Q2 = Matrix::Zero(j+2, j+2);
        Q2.topLeftCorner(j+1, j+1) = result.Q;
        Q2(j+1, j+1) = 1;

        Matrix R2 = Matrix::Zero(j+2, j);
        R2.topRows(j+1) = result.R;

        auto QR = qr_append(Q2, R2, hcur, j);

        result.Q = QR.householderQ();
        result.R = QR.matrixQR().triangularView<Eigen::Upper>();

        // residual
        auto res = std::abs(result.Q(0, result.Q.cols()-1));
        if(res < atol || breakdown) {
            break;
        }
    }
    j += 1;

    if(std::isnan(result.R(j,j)) || std::isinf(result.R(j,j))) {
        throw std::runtime_error("nans encountered, bail!");
    }

    {
        // solve lsq problem

        Matrix lsqR = result.R.topLeftCorner(j+1, j+1);
        Matrix lsqQ = result.Q.topLeftCorner(1, j+1).transpose();

        result.y = lsqR.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(lsqQ);
    }

    return result;
}

SpVec compute(const LinearOperator &A, Eigen::Ref<SpVec> b, Eigen::Ref<SpVec> x0, double tol, std::size_t maxiter,
              const LinearOperator &M, std::size_t inner_m, std::size_t outer_k,
              std::vector<std::tuple<Vec, Vec>> &outer_v, bool storeOuterAv, bool prependOuterAv) {
    auto log = spdlog::stdout_color_mt("console");
    log->set_level(spdlog::level::debug);

    SparseSystem system(A, M, x0, b);

    auto matvec = [&system, log](const SpVec &vec) {
        return system.call(vec);
    };
    auto psolve = [&system](const SpVec &vec) {
        return system.prec(vec);
    };

    auto bNorm = b.norm();
    if (bNorm == 0) {
        bNorm = 1;
    }

    SpVec x = system.x0();

    ArnoldiResult<SpVec> arnoldiResult;

    for (std::size_t k_outer = 0; k_outer < maxiter; ++k_outer) {
        SpVec rOuter = matvec(x) - system.b();

        auto rNorm = rOuter.norm();

        // check stopping criterion
        if (rNorm <= tol * bNorm || rNorm <= tol) {
            break;
        }

        // inner lgmres iteration
        SpVec v0 = -psolve(rOuter);
        auto innerRes0 = v0.norm();

        if (innerRes0 == 0) {
            auto rnorm = rOuter.norm();
            throw std::runtime_error(fmt::format("Preconditioner returned a zero vector; |v| ~ {}, |M v| = 0))",
                                                 rnorm));
        }

        v0 /= innerRes0;
        try {
            arnoldiResult = arnoldi(matvec, v0, inner_m, tol * bNorm / rNorm, psolve, outer_v, prependOuterAv);
            arnoldiResult.y = innerRes0*arnoldiResult.y;
        } catch(const std::runtime_error &) {
            // bail
            return x;
        }

        // gmres terminated, eval solution
        Vec dx = arnoldiResult.zs[0] * arnoldiResult.y(0);
        for (std::size_t i = 1; i < std::min(arnoldiResult.zs.size(), static_cast<std::size_t>(arnoldiResult.y.size()));
             ++i) {
            dx += arnoldiResult.zs.at(i) * arnoldiResult.y(i);
        }

        // store augmentation vectors
        auto nx = dx.norm();
        if(nx > 0) {
            if(storeOuterAv) {
                Vec q = arnoldiResult.Q * arnoldiResult.R  * arnoldiResult.y;
                Vec ax = q(0) * arnoldiResult.vs[0];
                for(std::size_t i = 1; i < std::min(arnoldiResult.vs.size(), static_cast<std::size_t>(q.size())); ++i) {
                    ax += q(i) * arnoldiResult.vs.at(i);
                }
                outer_v.emplace_back(std::make_tuple(dx/nx, ax/nx));
            } else {
                outer_v.emplace_back(std::make_tuple(dx/nx, Vec{}));
            }
        }

         // only finite number of augmentation vectors
        while (outer_v.size()>outer_k) {
            outer_v.erase(outer_v.begin());
        }

        x += dx;
    }

    return x;
}

}
}

}
