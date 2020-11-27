#!/usr/bin/env python
import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='RestraintLib',
    version='2019.10.5',
    description='Bond length and angle restraints for DNA and RNA oligonucleotides',
    long_description=read('README'),
    author='Marcin Kowiel, Dariusz Brzezinski',
    url='https://github.com/mkowiel/restraintlib',
    packages=['restraintlib'],
    license = "BSD",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 2.7"
    ],
)