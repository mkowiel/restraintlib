#!/usr/bin/env python
import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='RestraintLib',
    version='2019.11.1',
    description='Bond length and angle restraints for DNA and RNA oligonucleotides',
    long_description=read('README.md'),
    author='Marcin Kowiel, Dariusz Brzezinski',
    url='https://github.com/mkowiel/restraintlib',
    packages=['restraintlib'],
    license = "BSD",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    install_requires=[
        'six',
        #'cctbx>=2020.8',
        'scikit-learn==0.20.3',
        'pandas==0.24.2',
    ]
)