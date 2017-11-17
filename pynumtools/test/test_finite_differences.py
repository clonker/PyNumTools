import unittest
import numpy as np
import matplotlib.pyplot as plt

from pynumtools.finite_differences import fd_coefficients, get_fd_matrix
from pynumtools.util import sliding_window


class TestFiniteDifferences(unittest.TestCase):
    def test_first_derivative(self):
        x0 = np.arange(0, 2.0 * np.pi, 0.05)
        xx = []
        for x in x0:
            if np.random.random() < .3:
                xx.append(x)
        x0 = np.array(xx)

        testf = np.array([np.sin(x) for x in x0])
        true_deriv = [np.cos(x) for x in x0]

        wwidth = 5
        deriv = np.empty_like(x0)
        for ix, (wx, wy) in enumerate(zip(sliding_window(x0, width=wwidth, fixed_width=True),
                                          sliding_window(testf, width=wwidth, fixed_width=True))):
            x = x0[ix]
            coeff = fd_coefficients(x, wx, k=1)
            deriv[ix] = coeff.dot(wy)
        plt.plot(x0, deriv, 'o')
        plt.plot(x0, true_deriv)
        plt.show()
        print(np.array(deriv) - np.array(true_deriv))

    def test_second_derivative(self):
        x0 = np.arange(0, 2.0 * np.pi, 0.05)
        xx = []
        for x in x0:
            if np.random.random() < .3:
                xx.append(x)
        x0 = np.array(xx)

        testf = np.array([np.sin(x) for x in x0])
        true_deriv = [-np.sin(x) for x in x0]

        wwidth = 5
        deriv = np.empty_like(x0)
        for ix, (wx, wy) in enumerate(zip(sliding_window(x0, width=wwidth, fixed_width=True),
                                          sliding_window(testf, width=wwidth, fixed_width=True))):
            x = x0[ix]
            coeff = fd_coefficients(x, wx, k=2)
            deriv[ix] = coeff.dot(wy)
        plt.plot(x0, deriv, 'o')
        plt.plot(x0, true_deriv)
        plt.show()
        print(np.array(deriv) - np.array(true_deriv))

    def test_sliding_window_fixed_width(self):
        n = 0
        seq = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        width = 2
        for ix in sliding_window(seq, width=width, fixed_width=True):
            assert len(ix) == 2 * width + 1
            assert seq[n] >= np.min(ix)
            assert seq[n] <= np.max(ix)
            n += 1

        assert n == len(seq), "need n={} == len(seq)=={}".format(n, len(seq))

    def test_differentiation_matrix(self):
        x0 = np.arange(0, 2.0 * np.pi, 0.05)
        xx = []
        for x in x0:
            if np.random.random() < .3:
                xx.append(x)
        x0 = np.array(xx)

        fun = np.sin(x0)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)

        ax1.set_title('first derivative')
        ax1.plot(x0, get_fd_matrix(x0, k=1, window_width=2)*fun, 'o')
        ax1.plot(x0, np.cos(x0))

        ax2.set_title('second derivative')
        ax2.plot(x0, get_fd_matrix(x0, k=2, window_width=2)*fun, 'o')
        ax2.plot(x0, -np.sin(x0))

        ax3.set_title('third derivative')
        ax3.plot(x0, get_fd_matrix(x0, k=3, window_width=2)*fun, 'o')
        ax3.plot(x0, -np.cos(x0))

        ax4.set_title('fourth derivative')
        ax4.plot(x0, get_fd_matrix(x0, k=4, window_width=4)*fun, 'o')
        ax4.plot(x0, np.sin(x0))

        plt.show()

    def test_sliding_window(self):
        n = 0
        seq = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        width = 2
        for ix in sliding_window(seq, width=width, fixed_width=False):
            assert len(ix) <= 2 * width + 1
            assert seq[n] >= np.min(ix)
            assert seq[n] <= np.max(ix)
            n += 1

        assert n == len(seq), "need n={} == len(seq)=={}".format(n, len(seq))


if __name__ == '__main__':
    unittest.main()
