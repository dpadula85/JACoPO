#!/usr/bin/env python

import setuptools
from numpy.distutils.core import Extension, setup

setup(
    name="JACoPO",
    version="1.0",
    author="Daniele Padula",
    author_email="dpadula85@yahoo.it",
    description="A python package to compute electrostatic interactions",
    url="https://github.com/dpadula85/JACoPO",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPL License",
        "Operating System :: OS Independent",
    ],
    ext_modules=[ Extension('JACoPO.tdc', ['JACoPO/tdc.f90'],
                  extra_f90_compile_args=['-fopenmp', '-lgomp'],
                  extra_link_args=['-lgomp']) ],
	scripts=['JACoPO/bin/compute_coupling'],
    zip_safe=False
)
