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

from b2ac.compat import *
from b2ac.matrix.matrix_ref import inverse_symmetric_3by3_double, inverse_symmetric_3by3_int, add_symmetric_matrix


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
    x = points[:, 0]
    y = points[:, 1]

    D1 = np.vstack([x ** 2, x * y, y ** 2]).T
    D2 = np.vstack([x, y, np.ones((len(x), ), dtype=x.dtype)]).T
    S1 = D1.T.dot(D1)
    S2 = D1.T.dot(D2)
    S3 = D2.T.dot(D2)
    T = -np.linalg.inv(S3).dot(S2.T)
    M = S1 + S2.dot(T)
    M = np.array([M[2, :] / 2, -M[1, :], M[0, :] / 2])
    eval, evec = np.linalg.eig(M)
    cond = (4 * evec[:, 0] * evec[:, 2]) - (evec[:, 1] ** 2)
    I = np.where(cond > 0)[0]
    a1 = evec[:, I[np.argmin(cond[I])]]
    return np.concatenate([a1, T.dot(a1)])


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
    points = np.array(points, 'float')
    S = _calculate_scatter_matrix_py(points[:, 0], points[:, 1])
    S3 = S[3:, 3:]
    S3 = np.array([S3[0, 0], S3[0, 1], S3[0, 2], S3[1, 1], S3[1, 2], S3[2, 2]])
    S3_inv = inverse_symmetric_3by3_double(S3).reshape((3, 3))
    S2 = S[:3, 3:]
    T = -np.dot(S3_inv, S2.T)
    M = S[:3, :3] + np.dot(S2, T)
    inv_mat = np.array([[0, 0, 0.5], [0, -1, 0], [0.5, 0, 0]], 'float')
    M = inv_mat.dot(M)

    e_vals, e_vect = np.linalg.eig(M)

    try:
        elliptical_solution_index = np.where(((4 * e_vect[0, :] * e_vect[2, :]) - ((e_vect[1, :] ** 2))) > 0)[0][0]
    except:
        # No positive eigenvalues. Fit was not ellipse.
        raise ArithmeticError("No elliptical solution found.")

    a = e_vect[:, elliptical_solution_index]
    if a[0] < 0:
        a = -a
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

    try:
        elliptical_solution_index = np.where(((4 * e_vect[0, :] * e_vect[2, :]) - ((e_vect[1, :] ** 2))) > 0)[0][0]
    except:
        # No positive eigenvalues. Fit was not ellipse.
        raise ArithmeticError("No elliptical solution found.")
    a = e_vect[:, elliptical_solution_index]
    return np.concatenate((a, np.dot(T_no_det, a) / det_S3))


def _calculate_scatter_matrix_py(x, y):
    """Calculates the complete scatter matrix for the input coordinates.

    :param x: The x coordinates.
    :type x: :py:class:`numpy.ndarray`
    :param y: The y coordinates.
    :type y: :py:class:`numpy.ndarray`
    :return: The complete scatter matrix.
    :rtype: :py:class:`numpy.ndarray`

    """
    D = np.ones((len(x), 6), dtype=x.dtype)
    D[:, 0] = x * x
    D[:, 1] = x * y
    D[:, 2] = y * y
    D[:, 3] = x
    D[:, 4] = y

    return D.T.dot(D)


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


