import scipy.sparse as _scsp
import numpy as _np
import pynumtools.pynumtools_binding as _binding
import scipy.sparse.linalg as _splin


def lgmres(A, b, x0=None, tol=1e-5, maxiter=1000, M=None, inner_m=30, outer_k=3, outer_v=None,
           store_outer_Av=True, prepend_outer_Av=False):
    _A = A
    if M is None:
        if isinstance(A, _splin.LinearOperator):
            _M = _splin.LinearOperator(A.shape, lambda _x: _x)
        elif callable(A):
            _M = lambda _x: _x
        else:
            _M = _scsp.identity(A.shape[0])
    else:
        _M = M
        if isinstance(A, _splin.LinearOperator):
            # A is a lin op
            if callable(M):
                _M = _splin.LinearOperator(A.shape, lambda _x: M(_x))
            else:
                _M = _splin.LinearOperator(A.shape, lambda _x: M*x)
        elif callable(A):
            # A is a lambda
            if callable(M):
                _M = M
            else:
                _M = lambda v: M*v
        else:
            # A is a matrix
            if isinstance(M, _splin.LinearOperator):
                _A = _splin.LinearOperator(M.shape, lambda _x: A*x)
            elif callable(M):
                _A = lambda v: A*v
            else:
                _A = A

    if callable(_A) and not callable(_M):
        raise ValueError("A and M can only be either both callable or both matrices!")

    if isinstance(_A, _splin.LinearOperator) and not isinstance(_M, _splin.LinearOperator):
        raise ValueError("A and M can only be either both linear operators or both matrices!")

    if x0 is None:
        x0 = _np.zeros_like(b)
    if outer_v is None:
        outer_v = []

    if callable(_A) and callable(_M):
        return _binding.lgmres(_A, b, x0, tol, maxiter, _M, inner_m, outer_k, outer_v, store_outer_Av, prepend_outer_Av)
    elif isinstance(_A, _splin.LinearOperator) and isinstance(_M, _splin.LinearOperator):
        return _binding.lgmres(_A.matvec, b, x0, tol, maxiter, _M.matvec, inner_m, outer_k, outer_v, store_outer_Av,
                               prepend_outer_Av)
    else:
        return _binding.lgmres_mat(_A, b, x0, tol, maxiter, _M, inner_m, outer_k, outer_v, store_outer_Av,
                                   prepend_outer_Av)


if __name__ == '__main__':
    mat = _np.array([[1, 0, 5], [0, 2, 0], [7, 0, 1]], dtype=_np.float64)
    mat = _scsp.csr_matrix(mat)

    op = lambda v: mat * v

    op = _splin.LinearOperator((3, 3), lambda v: mat * v)

    b = _np.array([1, 2, 3], dtype=_np.float64)
    # print(splin.lgmres(mat, b))
    x = lgmres(op, b, tol=1e-12)
    print(x)

    print(_splin.lgmres(mat, b))
