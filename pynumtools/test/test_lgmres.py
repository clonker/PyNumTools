import unittest

import numpy as np
import pynumtools.pynumtools_binding as binding

class TestLGMRES(unittest.TestCase):

    # @fixme
    def test_lgmres_scipy(self):
        mat = np.array([[1,0,0], [0,1,0], [0,0,1]], dtype=float)
        b = np.array([1,2,3], dtype=float)
        #x = binding.lgmres(mat, b, np.zeros_like(b))
