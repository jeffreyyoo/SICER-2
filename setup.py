import os
import sys
from setuptools import setup, find_packages

setup(
    name='SICER',
    version='0.3.5',
    description='Re-implementation of SICER algorithm. Still in development',
    long_description='Re-implementation of SICER algorithm. Still in development',
    url = 'https://github.com/jeffreyyoo/SICER-2',
    author = 'Jeffrey Yoo',
    author_email = 'jy2ma@virginia.edu',
    license = 'MIT',
    packages=find_packages(),
    scripts=['bin/sicer','bin/sicer_df'],
    install_requires = ['scipy','numpy'],
    keywords = ['ChIP-Seq','SICER'],
    classifiers=["Programming Language :: Python :: 3",
        "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"]
)
