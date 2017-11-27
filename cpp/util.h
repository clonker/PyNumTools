//
// Created by mho on 11/13/17.
//

#ifndef PYNUMTOOLS_UTIL_H
#define PYNUMTOOLS_UTIL_H

#if defined(_WIN64) or defined(_WIN32)
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#endif

#include <Eigen/SparseCore>
#include <utility>

namespace pnt {

using Matrix = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;

using SpMatrix = Eigen::SparseMatrix<double>;
using SpVec = Eigen::VectorXd; //Eigen::SparseMatrix<double>;

using LinearOperator = std::function<Vec(const Vec&)>;

namespace system {

/**
 * linear operator
 */
template<typename VecType>
using MatVec = std::function<VecType(const VecType&)>;

using DenseMatVec = MatVec<Vec>;
using SparseMatVec = MatVec<SpVec>;

}

template<typename Mat, typename Vec>
class System {
public:

    System(const Mat &A, const Mat &M, Eigen::Ref<Vec> x0, Eigen::Ref<Vec> b) : _M(M), _x0(x0), _b(b) {
        _A = [A](const Vec &v) {
            return A*v;
        };
        if(M.rows() == 0 || M.cols() == 0) {
            M.conservativeResize(A.rows(), A.cols());
            M.setIdentity();
        }
        _M = [M](const Vec &v) {
            return M*v;
        };
        if(_x0.rows() == 0 || _x0.cols() == 0) {
            _x0.conservativeResize(_b.rows(), _b.cols());
            _x0.setZero();
        }
    }

    System(const LinearOperator &A, const LinearOperator &M, Eigen::Ref<Vec> x0, Eigen::Ref<Vec> b)
            : _A(A), _M(M), _x0(x0), _b(b) {
        if(_x0.rows() == 0 || _x0.cols() == 0) {
            _x0.conservativeResize(_b.rows(), _b.cols());
            _x0.setZero();
        }
    }

    Vec call(const Vec &x) const {
        return _A(x);
    }
    Vec prec(const Vec &x) const {
        return _M(x);
    }

    const Vec &x0() const {
        return _x0;
    }

    const Vec &b() const {
        return _b;
    }

private:
    LinearOperator _A;
    LinearOperator _M;

    Vec _x0;
    Vec _b;
};

using SparseSystem = System<SpMatrix, SpVec>;
using DenseSystem = System<Matrix, Vec>;

}

#endif //PYNUMTOOLS_UTIL_H
