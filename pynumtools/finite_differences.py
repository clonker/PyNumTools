import numpy as _np
import scipy.sparse as _sparse
from pynumtools.util import sliding_window as _sliding_window


def fd_coefficients(x_bar, xs, k=1):
    """
    Calculate finite difference coefficients, so that

    f^(k)(x_bar) \approx f(xs) * fd_coeff(xs, window_around(x_bar), k).

    It is required that len(x) > k.

    See http://epubs.siam.org/doi/10.1137/S0036144596322507.
    :param x_bar: the evaluation point of the derivative
    :param xs: n values of f so that x_1 <= ... <= xbar <= ... <= x_n
    :param k: k-th derivative
    :return: a vector with weights w so that w.dot(f(xs)) \approx f^(k)(x_bar)
    """
    n = len(xs)
    if k >= n:
        raise ValueError("length(xs) = {} has to be larger than k = {}".format(len(xs), k))

    if _np.min(xs) > x_bar or _np.max(xs) < x_bar:
        raise ValueError("the grid xs has to be so that min(xs) <= xbar <= max(xs)")

    # change to m = n - 1 to compute coeffs for all derivatives, then output C
    m = k

    c1 = 1
    c4 = xs[0] - x_bar
    C = _np.zeros(shape=(n, m + 1))
    C[0, 0] = 1
    for i in range(1, n):
        mn = min(i, m)
        c2 = 1
        c5 = c4
        c4 = xs[i] - x_bar
        for j in range(0, i):
            c3 = xs[i] - xs[j]
            c2 = c2 * c3
            if j == i - 1:
                for s in range(mn, 0, -1):
                    C[i, s] = c1 * (s * C[i - 1, s - 1] - c5 * C[i - 1, s]) / c2
                C[i, 0] = -c1 * c5 * C[i - 1, 0] / c2
            for s in range(mn, 0, -1):
                C[j, s] = (c4 * C[j, s] - s * C[j, s - 1]) / c3
            C[j, 0] = c4 * C[j, 0] / c3
        c1 = c2
    c = C[:, -1]
    return c


def get_fd_matrix(xs, k=1, window_width=2):
    """
    Yields a sparse csc matrix D of size (n_nodes)x(n_nodes), such that

    D*fun

    yields an approximation to the k-th derivative.

    :param xs: the grid
    :param k: k-th derivative
    :param window_width: window width - larger means better derivative, usually
    :return: differentiation matrix
    """
    n_nodes = len(xs)
    indices = _np.arange(0, n_nodes, step=1, dtype=int)
    entries_per_row = 2 * window_width + 1
    data = _np.empty(shape=(entries_per_row * n_nodes,), dtype=xs.dtype)
    row_data = _np.empty_like(data)
    col_data = _np.empty_like(data)

    offset = 0
    for row, window in enumerate(_sliding_window(indices, width=window_width, fixed_width=True)):
        window_grid = xs[window]
        window_slice = slice(offset, offset + len(window))
        data[window_slice] = fd_coefficients(xs[row], window_grid, k=k)
        row_data[window_slice] = row
        col_data[window_slice] = window

        offset += len(window)

    return _sparse.csc_matrix((data, (row_data, col_data)), shape=(n_nodes, n_nodes))
