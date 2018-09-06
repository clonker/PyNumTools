import numpy as _np
import pynumtools.finite_differences as _fd
from scipy import sparse as _sparse
from scipy.sparse import linalg as _splin
from pynumtools.lgmres import lgmres as _lgmres


def _cumtrapz_operator(xs):
    """
    Returns matrix representation of the cumulative trapezoidal rule, i.e. \int_0^x f
    :param xs: grid
    :return: (n-1, n)-matrix
    """
    n = len(xs)
    data = _np.zeros(shape=(int(.5 * (n ** 2 + n) - 1),))
    row_data = _np.empty_like(data)
    col_data = _np.empty_like(data)

    current_row = _np.zeros(1)

    offset = 0
    for row in range(n - 1):
        dx = xs[row + 1] - xs[row]
        current_row[-1] += dx
        current_row = _np.append(current_row, dx)
        data[offset:offset + len(current_row)] = current_row
        row_data[offset:offset + len(current_row)] = row
        col_data[offset:offset + len(current_row)] = _np.array([range(len(current_row))])
        offset += len(current_row)
    data *= .5
    assert (len(data) == offset)
    return _sparse.csc_matrix((data, (row_data, col_data)), shape=(n - 1, n))


def tv_derivative(data, xs, u0=None, alpha=10., maxit=1000, linalg_solver_maxit=100., tol=1e-4,
                  verbose=False, solver='lgmres'):
    data = _np.asarray(data, dtype=_np.float64).squeeze()
    xs = _np.asarray(xs, dtype=_np.float64).squeeze()
    n = data.shape[0]
    assert xs.shape[0] == n, "the grid must have the same dimension as data"

    epsilon = 1e-6

    # grid of mid points between xs, extrapolating first and last node:
    #
    #    x--|--x--|---x---|-x-|-x
    #
    midpoints = _np.concatenate((
        [xs[0] - .5 * (xs[1] - xs[0])],
        .5 * (xs[1:] + xs[:-1]),
        [xs[-1] + .5 * (xs[-1] - xs[-2])]
    )).squeeze()
    assert midpoints.shape[0] == n + 1

    diff = _fd.get_fd_matrix_midpoints(midpoints, k=1, window_width=5)
    assert diff.shape[0] == n
    assert diff.shape[1] == n + 1

    diff_t = diff.transpose(copy=True).tocsc()
    assert diff.shape[0] == n
    assert diff.shape[1] == n + 1

    A = _cumtrapz_operator(midpoints)
    AT = A.transpose(copy=True)

    ATA = AT.dot(A)

    if u0 is None:
        u = _np.concatenate(([0], _np.diff(data), [0]))
    else:
        u = u0
    # Aadj_A = lambda v: A_adjoint(A(v))
    Aadj_offset = AT * (data[0] - data)

    E_n = _sparse.dia_matrix((n, n), dtype=xs.dtype)
    midpoints_diff = _np.diff(midpoints)

    for ii in range(1, maxit + 1):
        E_n.setdiag(midpoints_diff * (1. / _np.sqrt(_np.diff(u) ** 2.0 + epsilon)))
        L = diff_t * E_n * diff
        g = ATA.dot(u) + Aadj_offset + alpha * L * u

        # solve linear equation.
        info_i = 0
        if solver == 'lgmres' or solver == 'lgmres_scipy':
            if solver == 'lgmres_scipy':
                s, info_i = _splin.lgmres(A=alpha * L + ATA, b=-g, x0=u, tol=tol, maxiter=linalg_solver_maxit, outer_k=7)
            else:
                s = _lgmres(A=alpha * L + ATA, b=-g, x0=u, tol=tol, maxiter=linalg_solver_maxit)
        elif solver == 'bicgstab':
            [s, info_i] = _splin.bicgstab(A=alpha * L + ATA, b=-g, x0=u, tol=tol, maxiter=linalg_solver_maxit)
        elif solver == 'spsolve':
            s = _splin.spsolve((alpha * L + ATA), -g, use_umfpack=True)
        elif solver == 'np':
            s = _np.linalg.solve((alpha * L + ATA).todense().astype(_np.float64), (-g).astype(_np.float64))

        relative_change = _np.linalg.norm(s[0]) / _np.linalg.norm(u)
        if verbose:
            print('iteration {0:4d}: relative change = {1:.3e}, gradient norm = {2:.3e}'
                  .format(ii, relative_change, _np.linalg.norm(g)))
            if info_i > 0:
                print("WARNING - convergence to tolerance not achieved!")
            elif info_i < 0:
                print("WARNING - illegal input or breakdown")

        # Update current solution
        u = u + s

    return u


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    noise_variance = .08 * .08
    x0 = _np.arange(0, 2.0 * _np.pi, 0.005)
    testf = _np.array([_np.sin(x) for x in x0])
    testf = testf + _np.random.normal(0.0, _np.sqrt(noise_variance), x0.shape)
    true_deriv = [np.cos(x) for x in x0]

    tv_deriv0 = tv_derivative(testf, x0, alpha=.001, maxit=20, linalg_solver_maxit=5e6,
                              verbose=True, solver='spsolve', tol=1e-12)
    tv_deriv0_back = .5 * (tv_deriv0[1:] + tv_deriv0[:-1])
    tv_deriv = tv_derivative(testf, x0, alpha=.01, maxit=20, linalg_solver_maxit=5e6,
                              verbose=True, solver='spsolve', tol=1e-12)
    tv_deriv_back = .5 * (tv_deriv[1:] + tv_deriv[:-1])
    tv_deriv2 = tv_derivative(testf, x0, alpha=.5, maxit=20, linalg_solver_maxit=5e6,
                              verbose=True, solver='spsolve', tol=1e-12)
    tv_deriv2_back = .5 * (tv_deriv2[1:] + tv_deriv2[:-1])

    plt.figure(figsize=(15, 15))
    plt.plot(x0, testf, label='f')
    plt.plot(x0, true_deriv, label='df')

    plt.plot(x0, tv_deriv_back, label='tv, alpha=.01')
    plt.plot(x0, tv_deriv2_back, label='tv, alpha=.5')
    plt.plot(x0, tv_deriv0_back, label='tv, alpha=.001')
