import os
import sys
from numpy import get_include as get_numpy_include
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

__version__ = '0.0.1'

def get_eigen_include():
    eigen_inc = os.path.join(os.path.dirname(__file__), 'lib', 'eigen')
    assert os.path.exists(eigen_inc)
    return eigen_inc


def get_pybind_include():
    pybind_inc = os.path.join(os.path.dirname(__file__), 'lib', 'pybind11', 'include')
    assert os.path.exists(pybind_inc)
    return pybind_inc


def get_spdlog_include():
    spdlog_inc = os.path.join(os.path.dirname(__file__), 'lib', 'spdlog', 'include')
    assert os.path.exists(spdlog_inc)
    return spdlog_inc


def get_project_include():
    inc = os.path.join(os.path.dirname(__file__), 'cpp')
    assert os.path.exists(inc)
    return inc


ext_modules = [
    Extension(
        'pynumtools.pynumtools_binding',
        sources=['cpp/binding.cpp', 'cpp/lgmres.cpp'],
        language='c++',
        include_dirs=[
            get_pybind_include(), get_project_include(), get_numpy_include(), get_spdlog_include(), get_eigen_include()
        ],
        extra_compile_args=['-std=c++14', '-O0', '-fvisibility=hidden']
        # -ffast-math
    ),
]


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


setup(
    name='PyNumTools',
    version=__version__,
    author='Moritz Hoffmann',
    author_email='clonker at gmail.com',
    url='https://github.com/clonker/PyNumTools',
    description='estimator skeleton',
    long_description='yes',
    ext_modules=ext_modules,
    packages=find_packages(),
    zip_safe=False,
    install_requires=['numpy']
)
