import os
import sys
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext


if (float(sys.version[:3])<3):
    sys.stderr.write("ERROR: Python3 required! \n")
    sys.exit(1)

extra_c_args = ["-w","-O3","-ffast-math"]

ext_modules = [Extension("sicer.src.coarsegraining",["sicer/src/coarsegraining.c"],extra_compile_args=extra_c_args)]

setup(
    name='SICER',
    version='0.5.4',
    description='Re-implementation of SICER algorithm. Still in development',
    long_description='Re-implementation of SICER algorithm. Still in development',
    url = 'https://github.com/jeffreyyoo/SICER-2',
    author = 'Jeffrey Yoo',
    author_email = 'jy2ma@virginia.edu',
    license = 'MIT',
    packages=find_packages(),
    scripts=['bin/sicer','bin/sicer_df'],
    setup_requires=['numpy','scipy'],
    keywords = ['ChIP-Seq','SICER'],
    classifiers=["Programming Language :: Python :: 3",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"],
    ext_modules=ext_modules
)
