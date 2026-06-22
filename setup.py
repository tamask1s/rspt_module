from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

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
             'rspt_module/ecg_analysis.cpp',
             'rspt_module/ecg_analysis_c.cpp'],
    include_dirs=[],
    language='c++',
    extra_compile_args=['-std=c++17']
)


class BuildExt(build_ext):
    def finalize_options(self):
        super().finalize_options()

        import numpy
        import pybind11

        for extension in self.extensions:
            extension.include_dirs.append(numpy.get_include())
            extension.include_dirs.append(pybind11.get_include())


setup(
    name='rspt_module',
    version='0.1.0',
    description='Python module for rspt peak detection algorithm',
    long_description=read_readme(),
    long_description_content_type='text/markdown',
    packages=find_packages(),
    ext_modules=[module],
    cmdclass={'build_ext': BuildExt},
    install_requires=['numpy<2'],
    python_requires='>=3.6',
    zip_safe=False,
)
