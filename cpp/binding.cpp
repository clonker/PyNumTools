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
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "lgmres.h"

namespace py = pybind11;

PYBIND11_MODULE(pynumtools_binding, m) {
    using namespace py::literals;

    /*m.def("lgmres", &pnt::lgmres::sparse::compute,
          "A"_a,
          "b"_a,
          "x0"_a = pnt::Vec(0),
          "tol"_a = 1e-5,
          "maxiter"_a = 1000,
          "M"_a = id,
          "inner_m"_a = 30,
          "outer_k"_a = 3,
          "outer_v"_a = std::vector<std::tuple<pnt::Vec, pnt::Vec>>(),
          "store_outer_Av"_a = true,
          "prepend_outer_Av"_a = false);*/

    m.def("lgmres", [](const pnt::LinearOperator &A, Eigen::Ref<pnt::SpVec> b, Eigen::Ref<pnt::SpVec> x0, double tol,
                       std::size_t maxiter, const pnt::LinearOperator &M, std::size_t inner_m, std::size_t outer_k,
                       std::vector<std::tuple<pnt::Vec, pnt::Vec>> &outer_v, bool storeOuterAv, bool prependOuterAv) {
        return pnt::lgmres::sparse::compute(A, b, x0, tol, maxiter, M, inner_m, outer_k, outer_v, storeOuterAv,
                                            prependOuterAv);
    });

    m.def("lgmres_mat", [](const pnt::SpMatrix &A, Eigen::Ref<pnt::SpVec> b, Eigen::Ref<pnt::SpVec> x0, double tol,
                       std::size_t maxiter, const pnt::SpMatrix &M, std::size_t inner_m, std::size_t outer_k,
                       std::vector<std::tuple<pnt::Vec, pnt::Vec>> &outer_v, bool storeOuterAv, bool prependOuterAv) {
       return pnt::lgmres::sparse::compute([A](const pnt::SpVec &v) {
           return A*v;
       }, b, x0, tol, maxiter, [M](const pnt::SpVec &v) -> pnt::SpVec {
           if (M.rows() == 0 || M.cols() == 0) {
               return v;
           } else {
               return M*v;
           }
       }, inner_m, outer_k, outer_v, storeOuterAv, prependOuterAv);
    });

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif

}
