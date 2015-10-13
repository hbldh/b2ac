#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`fit`
==================

.. module:: fit
    :platform: Unix, Windows
    :synopsis:

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-10-13

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.eigenmethods.qr_algorithm import QR_algorithm_shift_Givens_int
from b2ac.eigenmethods.inverse_iteration import inverse_iteration_for_eigenvector_int
from b2ac.matrix import matrix_operations as mo

__all__ = ['fit_improved_B2AC_int']


def fit_improved_B2AC_int(points):
    """Ellipse fitting in Python with improved B2AC algorithm as described in
    this `paper <http://autotrace.sourceforge.net/WSCG98.pdf>`_.

    This version of the fitting uses int64 storage during calculations and performs the
    eigensolver on an integer array.

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: The conic section coefficients array defining the fitted ellipse.
    :rtype: :py:class:`numpy.ndarray`

    """
    e_conds = []
    M, T_no_det, determinant_S3 = _calculate_M_and_T_int64(points)

    e_vals = sorted(QR_algorithm_shift_Givens_int(M)[0])

    a = None
    for ev_ind in [1, 2, 0]:
        # Find the eigenvector that matches this eigenvector.
        eigenvector, e_norm = inverse_iteration_for_eigenvector_int(M, e_vals[ev_ind])
        # See if that eigenvector yields an elliptical solution.
        elliptical_condition = (4 * eigenvector[0] * eigenvector[2]) - (eigenvector[1] ** 2)
        e_conds.append(elliptical_condition)
        if elliptical_condition > 0:
            a = eigenvector
            break

    if a is None:
        raise ArithmeticError("No elliptical solution found.")

    conic_coefficients = np.concatenate((a, np.dot(T_no_det, a) // determinant_S3))
    return conic_coefficients


def _scale_T_matrix(T_no_det, det_S3):
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


def _calculate_M_and_T_int64(points):
    """Part of the B2AC ellipse fitting algorithm, calculating the M and T
     matrices needed.

     This integer implementation also returns the determinant of the
     scatter matrix, which hasn't been applied to the matrix T yet.

     M is exact in integer values, but is truncated towards zero compared
     to the double version.

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: M, T undivided by determinant and the determinant.
    :rtype: tuple

    """
    S = _calculate_scatter_matrix_c(points[:, 0], points[:, 1])
    # Extract the symmetric parts of the S matrix.
    S1 = np.array([S[0, 0], S[0, 1], S[0, 2], S[1, 1], S[1, 2], S[2, 2]], dtype='int64')
    S3 = np.array([S[3, 3], S[3, 4], S[3, 5], S[4, 4], S[4, 5], S[5, 5]], dtype='int64')

    adj_S3, det_S3 = mo.inverse_symmetric_3by3_int64(S3)
    if np.linalg.norm((np.array(adj_S3, 'float') / det_S3) - np.linalg.norm(np.linalg.inv(S[3:, 3:]))) > 0.1:
        print("Norm diff: {0}".format(np.linalg.norm(np.array(adj_S3, 'float') / det_S3) - np.linalg.norm(S[3:,3:])))
    S2 = S[:3, 3:]

    T_no_det = - np.dot(np.array(adj_S3.reshape((3, 3)), 'int64'), np.array(S2.T, 'int64'))
    T_no_det, det_S3 = _scale_T_matrix(T_no_det, det_S3)
    M_term_2 = np.dot(np.array(S2, 'int64'), T_no_det) // det_S3
    M = mo.matrix_add_symmetric(M_term_2, S1)
    M[[0, 2], :] /= 2
    M[1, :] = -M[1, :]

    #T_no_det, scale = eig._scale_matrix(T_no_det)
    #det_S3 = det_S3 >> scale

    return M, T_no_det, det_S3


def _calculate_scatter_matrix_c(x, y):
    """Calculates the upper triangular scatter matrix for the input coordinates.

    :param x: The x coordinates.
    :type x: :py:class:`numpy.ndarray`
    :param y: The y coordinates.
    :type y: :py:class:`numpy.ndarray`
    :return: The upper triangular scatter matrix.
    :rtype: :py:class:`numpy.ndarray`

    """
    S = np.zeros((6, 6), 'int64')

    for i in xrange(len(x)):
        tmp_x2 = x[i] ** 2
        tmp_x3 = tmp_x2 * x[i]
        tmp_y2 = y[i] ** 2
        tmp_y3 = tmp_y2 * y[i]

        S[0, 0] += tmp_x2 * tmp_x2
        S[0, 1] += tmp_x3 * y[i]
        S[0, 2] += tmp_x2 * tmp_y2
        S[0, 3] += tmp_x3
        S[0, 4] += tmp_x2 * y[i]
        S[0, 5] += tmp_x2
        S[1, 2] += tmp_y3 * x[i]
        S[1, 4] += tmp_y2 * x[i]
        S[1, 5] += x[i] * y[i]
        S[2, 2] += tmp_y2 * tmp_y2
        S[2, 4] += tmp_y3
        S[2, 5] += tmp_y2
        S[3, 5] += x[i]
        S[4, 5] += y[i]

    S[5, 5] = len(x)

    # Doubles
    S[1, 1] = S[0, 2]
    S[1, 3] = S[0, 4]
    S[2, 3] = S[1, 4]
    S[3, 3] = S[0, 5]
    S[3, 4] = S[1, 5]
    S[4, 4] = S[2, 5]

    return S