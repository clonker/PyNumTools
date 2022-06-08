import sys

from setuptools import find_packages

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

cmake_args = []

metadata = \
    dict(
        name="pynumtools",
        version="0.2",
        author='Moritz Hoffmann, Christoph Fröhner',
        author_email='clonker at gmail.com, christoph.froehner at fu-berlin.de',
        url='https://github.com/clonker/PyNumTools',
        description='estimator skeleton',
        long_description='yes',
        zip_safe=False,
        packages=find_packages(where="."),
        package_dir={'pynumtools': 'pynumtools'},
        cmake_install_dir="pynumtools/",
        install_requires=["scipy", "numpy", "tqdm"],
        cmake_args=cmake_args,
        include_package_data=True,
        ext_modules=[]
    )

if __name__ == '__main__':
    setup(**metadata)
