import unittest

class TestSanity(unittest.TestCase):

    def test_import(self):
        import pynumtools.pynumtools_binding as binding
        assert binding.__version__ is not None
