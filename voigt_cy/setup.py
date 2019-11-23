#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 08:54:11 2019

@author: yukiitoh
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

sourcefiles = ['voigt_cy.pyx', 'Faddeeva.c']

extensions = [Extension("voigt_cy", sourcefiles,
                        libraries=["m"],
              extra_compile_args = ["-O3", "-fPIC", "-ffast-math", "-march=native", "-fno-strict-aliasing"],
              )]
# -fopenmp??
#extra_link_args=["-fopenmp"]
setup(
    ext_modules = cythonize(extensions,annotate=True,
                            ),
                            include_dirs=[numpy.get_include()]
    
)