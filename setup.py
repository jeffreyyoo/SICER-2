import os
import sys
from setuptools import setup, find_packages

setup(
    name='SICER',
    version='0.1dev',
    description='Re-implementation of SICER algorithm. Still in development',
    long_description='Re-implementation of SICER algorithm. Still in development',
    url = 'https://github.com/jeffreyyoo/SICER-2'
    author = 'Jeffrey Yoo',
    author_email = 'jy2ma@virginia.edu'
    license = 'MIT',
    packages=find_packages(),
    scripts='bin/sicer'
    install_requires = ['scipy','numpy']
    keywords = ['ChIP-Seq','SICER'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 1 - Alpha",
        "Environment :: Other Environment", "Intended Audience :: De velopers",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ]
)
