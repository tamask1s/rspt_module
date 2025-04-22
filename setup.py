from setuptools import setup, Extension
import numpy

# pybind11 include dir-t csak build időben töltjük be
try:
    import pybind11
    pybind11_include = pybind11.get_include()
except ImportError:
    pybind11_include = ""

module = Extension(
    'rspt_module',
    sources=['rspt_module/rspt_module.cpp',
             'rspt_module/iir_filter_design.cpp'],
    include_dirs=[numpy.get_include(), pybind11_include],
    language='c++'
)

setup(
    name='rspt_module',
    version='0.1',
    description='2D array squaring functions in C++',
    ext_modules=[module],
    install_requires=['numpy'],
    setup_requires=['pybind11>=2.5.0'],
)
