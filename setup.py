from setuptools import setup, Extension, find_packages
import numpy
import os

# pybind11 include dir-t csak build időben töltjük be
try:
    import pybind11
    pybind11_include = pybind11.get_include()
except ImportError:
    pybind11_include = ""

# Read README for long description
def read_readme():
    try:
        with open('README.md', 'r', encoding='utf-8') as f:
            return f.read()
    except:
        return 'Python module for rspt peak detection algorithm'

module = Extension(
    'rspt_module.rspt_module',
    sources=['rspt_module/rspt_module.cpp',
             'rspt_module/lib_filter/iir_filter_design.cpp',
             'rspt_module/lib_filter/iir_filter.cpp',
             'rspt_module/ecg_analysis.cpp'],
    include_dirs=[numpy.get_include(), pybind11_include],
    language='c++',
    extra_compile_args=['-std=c++17']
)

setup(
    name='rspt_module',
    version='0.1.0',
    description='Python module for rspt peak detection algorithm',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    ext_modules=[module],
    install_requires=['numpy'],
    setup_requires=['pybind11>=2.5.0,<2.11.0', 'numpy'],
    python_requires='>=3.6',
    zip_safe=False,
)
