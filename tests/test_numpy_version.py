#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`test_numpy_version`
==================

.. module:: test_numpy_version
    :synopsis:

.. moduleauthor:: hbldh <henrik.blidh@swedwise.com>

Created on 2015-09-24, 08:06:46

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

import b2ac.fit


class TestNumpySuite():

    def naive_test_1(self):
        points = np.array([[3, 0], [0, 5], [3, 10], [6, 5]])
        e1 = b2ac.fit.fit_improved_B2AC_numpy(points.copy())
        e2 = b2ac.fit.fit_improved_B2AC_double(points.copy())
        assert np.linalg.norm(np.array(e1[0]) - np.array(e2[0])) < 1e-7