#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: ellipse_fitting
   :platform: Unix, Windows
   :synopsis: Ellipse fitting algorithms and handling of ellipse information.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2013-05-05, 23:22

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np
from geometry.matrix import inverse_symmetric_3by3_double, inverse_symmetric_3by3_int, add_symmetric_matrix


def fit_B2AC(points):
    """Ellipse fitting in Python with numerically unstable algorithm.

    Described `here <http://research.microsoft.com/pubs/67845/ellipse-pami.pdf>`_.

    N_POLYPOINTS.B. Do not use, since it works with almost singular matrix.

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: The conic section array defining the fitted ellipse.
    :rtype: :py:class:`numpy.ndarray`

    """
    import scipy.linalg as scla

    constraint_matrix = np.zeros((6, 6))
    constraint_matrix[0, 2] = 2
    constraint_matrix[1, 1] = -1
    constraint_matrix[2, 0] = 2

    S = _calculate_scatter_matrix_py(points[:, 0], points[:, 1])

    evals, evect = scla.eig(S, constraint_matrix)
    ind = np.where(evals == (evals[evals > 0].min()))[0][0]
    return evect[:, ind]


def fit_improved_B2AC(points):
    """Ellipse fitting in Python with improved B2AC algorithm as described in
    this `paper <http://autotrace.sourceforge.net/WSCG98.pdf>`_.

    This version of the fitting uses float storage during calculations and performs the
    eigensolver on a float array.

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: The conic section array defining the fitted ellipse.
    :rtype: :py:class:`numpy.ndarray`

    """

    S = _calculate_scatter_matrix_py(points[:, 0], points[:, 1])
    S3_inv = inverse_symmetric_3by3_double(S[3:, 3:]).reshape((3, 3))
    S2 = S[:3, 3:]
    T = -np.dot(S3_inv, S2.T)
    M = S[:3, :3] + np.dot(S2, T)
    M[[0, 2], :] /= 2
    M[1, :] = -M[1, :]

    e_vals, e_vect = np.linalg.eig(M)

    #print("Eigenvalues: " + ", ".join([str(e) for e in e_vals]))
    #print("Eigenvectors: " + ", ".join([str(e) for e in e_vect.flatten()]))

    try:
        elliptical_solution_index = np.where(((4 * e_vect[0, :] * e_vect[2, :]) - ((e_vect[1, :] ** 2))) > 0)[0][0]
    except:
        # No positive eigenvalues. Fit was not ellipse.
        raise ArithmeticError("No elliptical solution found.")

    a = e_vect[:, elliptical_solution_index]
    return np.concatenate((a, np.dot(T, a)))


def fit_improved_B2AC_int(points):
    """Ellipse fitting in Python with improved B2AC algorithm as described in
    this `paper <http://autotrace.sourceforge.net/WSCG98.pdf>`_.

    This version of the fitting uses int64 storage during calculations and performs the
    eigensolver on an integer array.

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: The conic section array defining the fitted ellipse.
    :rtype: :py:class:`numpy.ndarray`

    """
    S = _calculate_scatter_matrix_c(points[:, 0], points[:, 1])
    S1 = np.array([S[0, 0], S[0, 1], S[0, 2], S[1, 1], S[1, 2], S[2, 2]])
    S3 = np.array([S[3, 3], S[3, 4], S[3, 5], S[4, 4], S[4, 5], S[5, 5]])
    adj_S3, det_S3 = inverse_symmetric_3by3_int(S3)
    S2 = S[:3, 3:]
    T_no_det = - np.dot(np.array(adj_S3.reshape((3, 3)), 'int64'), np.array(S2.T, 'int64'))
    M_term2 = np.dot(np.array(S2, 'int64'), T_no_det) // det_S3
    M = add_symmetric_matrix(M_term2, S1)
    M[[0, 2], :] /= 2
    M[1, :] = -M[1, :]

    e_vals, e_vect = np.linalg.eig(M)

    #print("Eigenvalues: " + ", ".join([str(e) for e in e_vals]))
    #print("Eigenvectors: " + ", ".join([str(e) for e in e_vect.flatten()]))

    try:
        elliptical_solution_index = np.where(((4 * e_vect[0, :] * e_vect[2, :]) - ((e_vect[1, :] ** 2))) > 0)[0][0]
    except:
        # No positive eigenvalues. Fit was not ellipse.
        raise ArithmeticError("No elliptical solution found.")
    a = e_vect[:, elliptical_solution_index]
    return np.concatenate((a, np.dot(T_no_det, a) / det_S3))


def conic_to_general(conic_coeffs):
    """

    :param conic_coeffs: The six coefficients defining the ellipse as a conic shape.
    :type conic_coeffs: :py:class:`numpy.ndarray` or tuple
    :return: The general form for the ellipse. Returns tuple :math:`(x_c,\ y_c),\\ (a,\\ b),\\ \\theta`
        that fulfills the equation

        .. math::

            \\frac{((x-x_c)\\cos(\\theta) + (y-y_c)\\sin(\\theta))^2}{a^2} +
            \\frac{((x-x_c)\\sin(\\theta) - (y-y_c)\\sin(\\theta))^2}{b^2} = 1

    :rtype: tuple

    """
    a, b, c, d, e, f = conic_coeffs
    denom = 2 * ((b ** 2) - (a * c))
    x = (c * d - b * e) / denom
    y = (a * e - b * d) / denom
    mu = 1 / ((a * (x ** 2)) + (2 * b * x * y) + (c * (y ** 2)) - f)

    sqrt_expr = np.sqrt(((mu * a - mu * c) ** 2) + (4 * ((mu * b) ** 2)))
    min_axis = 1 / np.sqrt((mu * a + mu * c + sqrt_expr) / 2)
    maj_axis = 1 / np.sqrt((mu * a + mu * c - sqrt_expr) / 2)
    angle = np.arctan2(-2 * b, c - a)

    return (x, y), (maj_axis, min_axis), angle


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


def _calculate_scatter_matrix_c(x, y):
    """Calculates the upper triangular scatter matrix for the input coordinates.

    :param x: The x coordinates.
    :type x: :py:class:`numpy.ndarray`
    :param y: The y coordinates.
    :type y: :py:class:`numpy.ndarray`
    :return: The upper triangular scatter matrix.
    :rtype: :py:class:`numpy.ndarray`

    """
    S = np.zeros((6, 6), 'int32')

    # TODO: Remove mean of x and y coordinates to reduce word length.

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


