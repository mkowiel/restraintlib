#!/usr/bin/env python
import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='RestraintLib',
    version='2021.10.1',
    description='Bond length and angle restraints for DNA and RNA oligonucleotides',
    long_description=read('README.md'),
    author='Marcin Kowiel, Dariusz Brzezinski',
    url='https://github.com/mkowiel/restraintlib',
    packages=[
        'restraintlib',
        'restraintlib.lib',
        'restraintlib.lib.regressors',
        'restraintlib.data',
    ],
    include_package_data=True,  # done via MANIFEST.in under setuptools
    license = "BSD",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    install_requires=[
        'six',
        'cctbx-base>=2020.8',
        'scikit-learn>0.21',
        'pandas',
        'pytest',
    ],
    scripts=[
        'bin/restraintlib',
        'bin/retrain_restraintlib',
    ]
)
