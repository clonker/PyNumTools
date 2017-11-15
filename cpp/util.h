//
// Created by mho on 11/13/17.
//

#ifndef PYNUMTOOLS_UTIL_H
#define PYNUMTOOLS_UTIL_H

#include <Eigen/SparseCore>

namespace pnt {

using Matrix = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;

using SpMatrix = Eigen::SparseMatrix<double>;
using SpVec = Eigen::VectorXd; //Eigen::SparseMatrix<double>;

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

    System(const Mat &A, const Mat &M, Eigen::Ref<Vec> x0, Eigen::Ref<Vec> b) : _A(A), _M(M), _x0(x0), _b(b) {
        if(_M.rows() == 0 && _M.cols() == 0) {
            _M.conservativeResize(_A.rows(), _A.cols());
            _M.setIdentity();
        }
        if(_x0.rows() == 0 && _x0.cols() == 0) {
            _x0.conservativeResize(_b.rows(), _b.cols());
            _x0.setZero();
        }
    }

    Vec call(const Vec &x) const {
        return _A*x;
    }
    Vec prec(const Vec &x) const {
        return _M*x;
    }

    const Mat &A() const {
        return _A;
    }
    Mat &A() {
        return _A;
    }

    const Mat &M() const {
        return _M;
    }
    Mat &M() {
        return _M;
    }

    const Vec &x0() const {
        return _x0;
    }

    const Vec &b() const {
        return _b;
    }

private:
    Mat _A;
    Mat _M;

    Vec _x0;
    Vec _b;
};

using SparseSystem = System<SpMatrix, SpVec>;
using DenseSystem = System<Matrix, Vec>;

}

#endif //PYNUMTOOLS_UTIL_H
