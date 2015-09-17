# -*- coding: utf-8 -*-
"""
:mod:`setup`
============

.. module:: setup
    :platform: Unix, Windows
    :synopsis: The Python Packaging setup file for Empyrean.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-09-17, 12:39

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from setuptools import setup, Extension
import numpy

ellipse_fitting_module = Extension('ellipse_fitter',
                                   sources=['b2ac/ext/src/ellipse_fitter.c',
                                            'b2ac/ext/src/Eigenmethods_double.c',
                                            'b2ac/ext/src/Eigenmethods_float.c',
                                            'b2ac/ext/src/Eigenmethods_int.c',
                                            'b2ac/ext/src/EllipseFit.c',
                                            'b2ac/ext/src/EllipseFit_double.c',
                                            'b2ac/ext/src/EllipseFit_float.c',
                                            'b2ac/ext/src/EllipseFit_int.c',
                                            ],
                                   include_dirs=[numpy.get_include(),
                                                 'b2ac/ext/src/'])

setup(
    name='b2ac',
    version='0.1.0',
    description='Python and C implementations of an ellipse fitting algorithm in double and fixed point precision.',
    author='Henrik Blidh',
    author_email='henrik.blidh@nedomkull.com',
    license='MIT',
    url='https://github.com/hbldh/ellipse-fitting',
    packages=[
        'b2ac',
        'b2ac.fit',
        'b2ac.fit.eigenmethods',
        'b2ac.fit.matrix',
        'b2ac.ext',
    ],
    package_data={'b2ac.ext': [
        'src/*',
        '*.cfg',
        '*.py']},
    install_requires=[
        'numpy>=1.6',
    ],
    dependency_links=[],
    ext_modules=[
        ellipse_fitting_module,
    ],
)
