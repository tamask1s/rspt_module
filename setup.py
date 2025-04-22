from setuptools import setup, Extension
import numpy

ext_modules = [
    Extension(
        "rspt_module",
        ["rspt_module/rspt_module.cpp"],
        include_dirs=[
            numpy.get_include(),
            # pybind11 include dir-t majd a setup_requires biztosítja
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
    install_requires=["numpy"],
    setup_requires=["pybind11>=2.5.0"],
    zip_safe=False,
)
