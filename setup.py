from setuptools import setup, Extension
import pybind11
import numpy

ext_modules = [
    Extension(
        "rspt_module",
        ["rspt_module/rspt_module.cpp"],
        include_dirs=[
            pybind11.get_include(),
            numpy.get_include(),
        ],
        language="c++"
    ),
]

setup(
    name="rspt_module",
    version="0.1",
    author="tamask1s",
    description="2D array squaring functions in C++",
    ext_modules=ext_modules,
    install_requires=["numpy", "pybind11"],
    zip_safe=False,
)
