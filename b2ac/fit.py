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
import b2ac.conversion as c2gconv

DEBUG = True


def fit_unstable_B2AC(points):
    """Ellipse fitting in Python with numerically unstable algorithm. Requires SciPy to run!

    Described `here <http://research.microsoft.com/pubs/67845/ellipse-pami.pdf>`_.

    N.B. Do not use, since it works with almost singular matrix.

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

    eigenvalues, eigenvalues = scla.eig(S, constraint_matrix)
    ind = np.where(eigenvalues == (eigenvalues[eigenvalues > 0].min()))[0][0]
    return eigenvalues[:, ind]


def fit_improved_B2AC_numpy(points):
    """Ellipse fitting in Python with improved B2AC algorithm as described in
    this `paper <http://autotrace.sourceforge.net/WSCG98.pdf>`_.

    This version of the fitting simply applies NumPy:s methods for calculating
    the conic section, modelled after the Matlab code in the paper:

    .. code-block::

        function a = fit_ellipse(x, y)

        D1 = [x .ˆ 2, x .* y, y .ˆ 2]; % quadratic part of the design matrix
        D2 = [x, y, ones(size(x))]; % linear part of the design matrix
        S1 = D1’ * D1; % quadratic part of the scatter matrix
        S2 = D1’ * D2; % combined part of the scatter matrix
        S3 = D2’ * D2; % linear part of the scatter matrix
        T = - inv(S3) * S2’; % for getting a2 from a1
        M = S1 + S2 * T; % reduced scatter matrix
        M = [M(3, :) ./ 2; - M(2, :); M(1, :) ./ 2]; % premultiply by inv(C1)
        [evec, eval] = eig(M); % solve eigensystem
        cond = 4 * evec(1, :) .* evec(3, :) - evec(2, :) .ˆ 2; % evaluate a’Ca
        a1 = evec(:, find(cond > 0)); % eigenvector for min. pos. eigenvalue
        a = [a1; T * a1]; % ellipse coefficients

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: The conic section array defining the fitted ellipse.
    :rtype: :py:class:`numpy.ndarray`

    """
    points = np.array(points, 'float')
    points, x_mean, y_mean = _remove_mean_values(points)

    x = points[:, 0]
    y = points[:, 1]
    D1 = np.vstack([x ** 2, x * y, y ** 2]).T
    D2 = np.vstack([x, y, np.ones((len(x), ), dtype=x.dtype)]).T
    S1 = D1.T.dot(D1)
    S2 = D1.T.dot(D2)
    S3 = D2.T.dot(D2)
    T = -np.linalg.inv(S3).dot(S2.T)
    M = S1 + S2.dot(T)
    #M2 = np.array([M[2, :] / 2, -M[1, :], M[0, :] / 2])
    inv_mat = np.array([[0, 0, 0.5], [0, -1, 0], [0.5, 0, 0]], 'float')
    M = inv_mat.dot(M)
    eigenvalues, eigenvectors = np.linalg.eig(M)
    cond = (4 * eigenvectors[:, 0] * eigenvectors[:, 2]) - (eigenvectors[:, 1] ** 2)
    I = np.where(cond > 0)[0]
    ev_ind = I[np.argmin(cond[I])]
    a1 = eigenvectors[:, ev_ind]
    if DEBUG:
        print("Numpy Version, Elliptical solution = {0}: {1}".format(eigenvalues[ev_ind], a1))

    conic_coefficients = np.concatenate((a1, np.dot(T, a1)))
    rotated_euclidian_coefficients = list(c2gconv.conic_to_general_1(conic_coefficients))
    rotated_euclidian_coefficients[0] = (rotated_euclidian_coefficients[0][0] + x_mean,
                                         rotated_euclidian_coefficients[0][1] + y_mean)
    return rotated_euclidian_coefficients


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
    points, x_mean, y_mean = _remove_mean_values(points)

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

    if DEBUG:
        print("Double precision, Elliptical solution = {0}: {1}".format(e_vals[ev_ind], a))

    conic_coefficients = np.concatenate((a, np.dot(T, a)))
    rotated_euclidian_coefficients = list(c2gconv.conic_to_general_1(conic_coefficients))
    rotated_euclidian_coefficients[0] = (rotated_euclidian_coefficients[0][0] + x_mean,
                                         rotated_euclidian_coefficients[0][1] + y_mean)
    return rotated_euclidian_coefficients


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
    points, x_mean, y_mean = _remove_mean_values(points)

    M, T_no_det, determinant_S3 = _calculate_M_and_T_int64(points)

    e_vals = sorted(qr.QR_algorithm_shift_Givens_int(M)[0])

    a = None
    for ev_ind in [1, 2, 0]:
        # Find the eigenvector that matches this eigenvector.
        eigenvector, e_norm = inv_iter.inverse_iteration_for_eigenvector_int(M, e_vals[ev_ind])
        # See if that eigenvector yields an elliptical solution.
        elliptical_condition = (4 * eigenvector[0] * eigenvector[2]) - (eigenvector[1] ** 2)
        e_conds.append(elliptical_condition)
        if elliptical_condition > 0:
            a = eigenvector
            break

    if a is None:
        raise ArithmeticError("No elliptical solution found.")

    if DEBUG:
        print("Integer precision, Elliptical solution = {0}: {1}".format(e_vals[ev_ind], a))

    conic_coefficients = np.concatenate((a, np.dot(T_no_det, a) // determinant_S3))
    rotated_euclidian_coefficients = list(c2gconv.conic_to_general_int(conic_coefficients, True))
    rotated_euclidian_coefficients[0] = (rotated_euclidian_coefficients[0][0] + x_mean,
                                         rotated_euclidian_coefficients[0][1] + y_mean)
    return rotated_euclidian_coefficients


def _remove_mean_values(points):
    x_mean = int(points[:, 0].mean())
    y_mean = int(points[:, 1].mean())
    return points - (x_mean, y_mean), x_mean, y_mean


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
