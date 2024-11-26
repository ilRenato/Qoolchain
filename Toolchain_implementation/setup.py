from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize

setup(
    name='ToolchainVectorVersion',
    ext_modules=cythonize([Extension(name='QuboPreprocessing',
                                     sources=['PreprocessingToolchain.pyx', 'PreQubo.cpp', 'ImplicationNetwork.cpp',
                                              'ResidualNetwork.cpp', 'FlowFunctions.cpp', 'GraphFunctions.cpp'],
                                     language='c++')]),
    cmdclass={'build_ext': build_ext},
)
