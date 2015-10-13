#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`fit`
==================

.. module:: fit
    :synopsis:

.. moduleauthor:: hbldh <henrik.blidh@swedwise.com>

Created on 2015-09-24, 07:18:22

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

import b2ac.matrix.matrix_operations as mo
import b2ac.eigenmethods.qr_algorithm as qr
import b2ac.eigenmethods.inverse_iteration as inv_iter


def fit_improved_B2AC_double(points):
    """Ellipse fitting in Python with improved B2AC algorithm as described in
    this `paper <http://autotrace.sourceforge.net/WSCG98.pdf>`_.

    This version of the fitting uses float storage during calculations and performs the
    eigensolver on a float array. It only uses `b2ac` package methods for fitting, to
    be as similar to the integer implementation as possible.

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: The conic section array defining the fitted ellipse.
    :rtype: :py:class:`numpy.ndarray`

    """
    e_conds = []
    points = np.array(points, 'float')

    M, T = _calculate_M_and_T_double(points)

    e_vals = sorted(qr.QR_algorithm_shift_Givens_double(M)[0])

    a = None
    for ev_ind in [1, 2, 0]:
        # Find the eigenvector that matches this eigenvector.
        eigenvector = inv_iter.inverse_iteration_for_eigenvector_double(M, e_vals[ev_ind], 5)

        # See if that eigenvector yields an elliptical solution.
        elliptical_condition = (4 * eigenvector[0] * eigenvector[2]) - (eigenvector[1] ** 2)
        e_conds.append(elliptical_condition)
        if elliptical_condition > 0:
            a = eigenvector
            break

    if a is None:
        print("Eigenvalues = {0}".format(e_vals))
        print("Elliptical conditions = {0}".format(e_conds))
        raise ArithmeticError("No elliptical solution found.")

    conic_coefficients = np.concatenate((a, np.dot(T, a)))
    return conic_coefficients


def _calculate_M_and_T_double(points):
    """Part of the B2AC ellipse fitting algorithm, calculating the M and T
     matrices needed.

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: Matrices M and T.
    :rtype: tuple

    """
    S = _calculate_scatter_matrix_py(points[:, 0], points[:, 1])
    S1 = S[:3, :3]
    S3 = np.array([S[3, 3], S[3, 4], S[3, 5], S[4, 4], S[4, 5], S[5, 5]])
    S3_inv = mo.inverse_symmetric_3by3_double(S3).reshape((3, 3))
    S2 = S[:3, 3:]
    T = -np.dot(S3_inv, S2.T)
    M_term_2 = np.dot(S2, T)
    M = S1 + M_term_2
    M[[0, 2], :] = M[[2, 0], :] / 2
    M[1, :] = -M[1, :]

    return M, T


def _calculate_scatter_matrix_py(x, y):
    """Calculates the complete scatter matrix for the input coordinates.

    :param x: The x coordinates.
    :type x: :py:class:`numpy.ndarray`
    :param y: The y coordinates.
    :type y: :py:class:`numpy.ndarray`
    :return: The complete scatter matrix.
    :rtype: :py:class:`numpy.ndarray`

    """
    D = np.ones((len(x), 6), 'int64')
    D[:, 0] = x * x
    D[:, 1] = x * y
    D[:, 2] = y * y
    D[:, 3] = x
    D[:, 4] = y

    return np.dot(D.T, D)


