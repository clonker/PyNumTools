import unittest

import scipy.sparse.linalg as splin
import numpy as np
import pynumtools.pynumtools_binding as binding

class TestLGMRES(unittest.TestCase):

    def test_lgmres_scipy(self):

        mat = np.array([[1,0,0], [0,1,0], [0,0,1]], dtype=np.float64)
        b = np.array([1,2,3], dtype=np.float64)
        splin.lgmres(mat, b)
        x = binding.lgmres(mat, b, np.zeros_like(b))
        print(x)
