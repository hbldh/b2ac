#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`fixed_point` --
======================

.. module:: fixed_point
   :platform: Unix, Windows
   :synopsis: Scaling of matrices for fixed point use.

.. moduleauthor:: hbh <henrik.blidh@nedomkull.com>

Created on 2013-09-09, 16:36

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np


def scale_64bit_matrix(v):
    m, M = np.abs(v).min(), np.abs(v).max()
    if np.log2(M) <= 24:
        return v, 0
    elif np.log2(M) <= 32:
        scale = 8
        return v >> scale, scale
    elif np.log2(M) <= 40:
        scale = 16
        return v >> scale, scale
    elif np.log2(M) <= 48:
        scale = 24
        return v >> scale, scale
    elif np.log2(M) <= 56:
        scale = 32
        return v >> scale, scale
    else:
        scale = 40
        return v >> scale, scale


def scale_64bit_vector(v):
    m, M = np.abs(v).min(), np.abs(v).max()

    scale = 0
    value = M >> 14
    while value != 0:
        value >>= 1
        scale += 1

    return v >> scale, scale


def scale_T_matrix(T_no_det, det_S3):
    m, M = np.abs(T_no_det).min(), np.abs(T_no_det).max()
    if np.log2(M) <= 32:
        scale = 0
    elif np.log2(M) <= 40:
        scale = 8
    elif np.log2(M) <= 48:
        scale = 16
    elif np.log2(M) <= 56:
        scale = 24
    else:
        scale = 32

    det_S3 >>= scale
    T_no_det >>= scale

    return T_no_det, det_S3
