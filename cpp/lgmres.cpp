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
auto qr_append(Eigen::Ref<Matrix> Q, Eigen::Ref<Matrix> R, Eigen::Ref<Matrix> u)
{   // Returns QR factorization of Q*R with u appended as the last column
    Matrix A = Q * R;
    A.conservativeResize(A.rows(),A.cols()+1);
    A.col(A.cols()-1) = u;
    return A.householderQr();
}
}

SpVec spLgmres(Eigen::Ref<SpMatrix> A, Eigen::Ref<SpVec> b, Eigen::Ref<SpVec> x0, double tol, std::size_t maxiter,
               Eigen::Ref<SpMatrix> M, std::size_t inner_m, std::size_t outer_k,
               std::vector<std::tuple<Vec, Vec>> &outer_v, bool storeOuterAv) {
    auto log = spdlog::stdout_color_mt("console");
    log->set_level(spdlog::level::debug);

    SparseSystem system (A, M, x0, b);

    auto N = A.rows();

    auto matvec = [&system](const SpVec &vec) {
        return system.call(vec);
    };
    auto psolve = [&system](const SpVec &vec) {
        return system.prec(vec);
    };

    auto bNorm = b.norm();
    if(bNorm == 0) {
        bNorm = 1;
    }

    SpVec x = system.x0();

    for(std::size_t k_outer = 0; k_outer < maxiter; ++k_outer) {
        SpVec rOuter = matvec(x) - system.b();

        auto rNorm = rOuter.norm();

        // check stopping criterion
        if(rNorm <= tol * bNorm || rNorm <= tol) {
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

    }

    return x;
}

Vec lgmres(Eigen::Ref<Matrix> A, Eigen::Ref<Vec> b, Vec x0, double tol, std::size_t maxiter,
           Matrix M, std::size_t inner_m, std::size_t outer_k, std::vector<std::tuple<Vec, Vec>> &outer_v,
           bool storeOuterAv) {

    auto console = spdlog::stdout_color_mt("console");
    console->set_level(spdlog::level::debug);

    auto n = b.size();
    Vec x = x0;
    if(x.rows() == 0 || x.cols() == 0) {
        x = Vec::Zero(n);
    }

    if(M.rows() == 0 && M.cols() == 0) {
        std::cout << "using identity as preconditioner" << std::endl;
        M = Matrix::Identity(A.cols(),A.rows());
    }

    auto b_norm = b.norm();
    if(b_norm == 0) {
        b_norm = 1;
    }

    for(std::size_t k_outer = 0; k_outer < maxiter; ++k_outer) {
        Vec r_outer = A*x - b;

        std::cout << "r_outer = (" << r_outer(0) << ", " << r_outer(1) << ", " << r_outer(2) << ")" << std::endl;

        // todo: callback

        // -- check stopping condition
        double r_norm = r_outer.norm();
        if (r_norm <= tol * b_norm || r_norm <= tol) {
            break;
        }

        // -- inner LGMRES iteration
        Vec vs0 = -M*r_outer;
        double inner_res_0 = vs0.norm();

        std::cout << "inner res 0: " << inner_res_0 << std::endl ;

        if(inner_res_0 == 0) {
            auto rnorm = r_outer.norm();
            throw std::runtime_error(fmt::format("Preconditioner returned a zero vector; |v| ~ {}, |M v| = 0))", rnorm));
        }

        vs0 /= inner_res_0;
        std::vector<Vec> vs = {vs0};
        std::vector<Vec> ws;

        Matrix Q = Matrix::Ones(1,1);
        Matrix R = Matrix::Zero(1,0);

        bool breakdown {false};
        std::size_t j = 1;
        for(; j < 1 + inner_m + outer_v.size(); ++j) {
            //     ++ evaluate
            Vec z;
            Vec v_new(0);
            if (j < outer_v.size() + 1) {
                z = std::get<0>(outer_v.at(j-1));
                v_new = std::get<1>(outer_v.at(j-1));
            } else if(j == outer_v.size() + 1) {
                z = vs0;
            } else {
                z = vs.back();
            }

            if(v_new.rows() == 0 || v_new.cols() == 0) {
                v_new = M * (A * z);
            }
            //     ++ orthogonalize
            auto v_new_norm = v_new.norm();

            Vec hcur = Vec::Zero(j+1);
            for(auto it = vs.begin(); it != vs.end(); ++it) {
                const auto ix = std::distance(vs.begin(), it);
                double alpha = it->dot(v_new);
                hcur[ix] = alpha;
                v_new -= alpha * (*it);
            }
            auto v_new_new_norm = v_new.norm();
            hcur[hcur.size()-1] = v_new_new_norm;

            if(std::abs(v_new_new_norm) > 1e-8) {
                v_new /= v_new_new_norm;
            }

            if(!(hcur[hcur.size() - 1] > eps * v_new_norm)) {
                // v_new essentially in the span of previous vectors,
                // or we have nans. Bail out after updating the QR
                // solution.
                breakdown = true;
            }
            vs.push_back(v_new);
            ws.push_back(z);


            // -- GMRES optimization problem

            // Add new column to H=Q*R, padding other columns with zeros

            Matrix Q2 = Matrix::Zero(j+1, j+1);
            Q2.topLeftCorner(j,j) = Q;
            Q2(j,j) = 1.;

            Matrix R2 = Matrix::Zero(j+1, j-1);
            R2.topRows(j) = R;

            QRMatrix QR = qr_append(Q2, R2, hcur);
            Q.conservativeResize(j+1, j+1);
            R.conservativeResize(j+1, j-1);
            Q = QR.householderQ();
            R = QR.matrixQR().triangularView<Eigen::Upper>();

            // Transformed least squares problem
            // || Q R y - inner_res_0 * e_1 ||_2 = min!
            // Since R = [R'; 0], solution is y = inner_res_0 (R')^{-1} (Q^H)[:j,0]

            // Residual is immediately known
            double inner_res = std::abs(Q(0,j)) * inner_res_0;

            // check for termination
            if(inner_res <= tol * inner_res_0 or breakdown) {
                break;
            }
        }

        {
            Vec lsqB = Q.topLeftCorner(1, j).transpose();
            Matrix lsqA = R.topLeftCorner(j, j);
            Vec y = A.jacobiSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(b);
            y *= inner_res_0;

            // todo check if y is finite, otherwise bail

            Vec dx = y(0) * ws[0];
            for (std::size_t i=1; i<std::min(static_cast<std::size_t>(y.size()),
                                             static_cast<std::size_t>(ws.size())); ++i) {
                auto &yi = y(i);
                auto &w = ws[i];
                dx += yi * w;
            }

            double nx = dx.norm();
            if (nx > 0.) {
                if(storeOuterAv) {
                    // todo
                    Vec q = Q * R * y;
                    Vec ax = vs[0]*q[0];
                    console->debug("vs[0] = ({}, {}, {})", vs[0](0), vs[0](1), vs[0](2));
                    console->debug("vs[1] = ({}, {}, {})", vs[1](0), vs[1](1), vs[1](2));
                    for (std::size_t i = 1; i < std::min(static_cast<std::size_t>(vs.size()),
                                                         static_cast<std::size_t>(q.size())); ++i) {
                        ax += q[i] * vs[i];
                    }
                    outer_v.emplace_back(std::make_tuple(dx / nx, ax / nx));
                } else {
                    Vec empty;
                    outer_v.emplace_back(std::make_tuple(dx / nx, empty));
                }
            }

            while (outer_v.size()>outer_k) {
                outer_v.erase(outer_v.begin());
            }

            x += dx;

        }

    }

    return x;
}

}
