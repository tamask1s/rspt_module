from distutils.core import setup, Extension
import numpy
import pybind11

module = Extension(
    'rspt_module',
    sources=['rspt_module/rspt_module.cpp'],
    include_dirs=[numpy.get_include(), pybind11.get_include()],
    language='c++'
)

setup(
    name='rspt_module',
    version='0.1',
    description='2D array squaring functions in C++',
    ext_modules=[module]
)
