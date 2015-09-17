#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: ellipse_fitting
   :platform: Unix, Windows
   :synopsis: Ellipse fitting algorithms and handling of ellipse information.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2013-09-09, 16:32

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

import b2ac.matrix.matrix_operations as mo
import b2ac.matrix.matrix_algorithms as ma
import b2ac.eigenmethods.qr_algorithm as qr
import b2ac.eigenmethods.inverse_iteration as inv_iter

DEBUG = True


def fit_B2AC(points):
    """Ellipse fitting in Python with numerically unstable algorithm.

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


def fit_improved_B2AC_double(points):
    """Ellipse fitting in Python with improved B2AC algorithm as described in
    this `paper <http://autotrace.sourceforge.net/WSCG98.pdf>`_.

    This version of the fitting uses float storage during calculations and performs the
    eigensolver on a float array.

    :param points: The [Nx2] array of points to fit ellipse to.
    :type points: :py:class:`numpy.ndarray`
    :return: The conic section array defining the fitted ellipse.
    :rtype: :py:class:`numpy.ndarray`

    """
    e_conds = []
    points, x_mean, y_mean = _remove_mean_values(points)

    M, T = _calculate_M_and_T_double(points)

    e_vals = sorted(qr.QR_algorithm_shift_Givens_double(M)[0])

    a = None
    for ev_ind in [1, 2, 0]:
        # Find the eigenvector that matches this eigenvector.
        eigenvector = inv_iter.inverse_iteration_for_eigenvector_double(M, e_vals[ev_ind], 5)

        # See if that eigenvector yields an elliptical solution.
        elliptical_condition = (4 * eigenvector[0] * eigenvector[2]) - ((eigenvector[1] ** 2))
        e_conds.append(elliptical_condition)
        if elliptical_condition > 0:
            a = eigenvector
            break

    if a is None:
        print("Eigenvalues = {0}".format(e_vals))
        print("Elliptical conditions = {0}".format(e_conds))
        raise ArithmeticError("No elliptical solution found.")

    if DEBUG:
        print("Elliptical solution = {0}: {1}".format(e_vals[ev_ind], a))

    conic_coefficients = np.concatenate((a, np.dot(T, a)))
    rotated_euclidian_coefficients = list(_general_to_rotated(conic_coefficients))
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
        elliptical_condition = (4 * eigenvector[0] * eigenvector[2]) - ((eigenvector[1] ** 2))
        e_conds.append(elliptical_condition)
        if elliptical_condition > 0:
            a = eigenvector
            break

    if a is None:
        raise ArithmeticError("No elliptical solution found.")

    if DEBUG:
        print("Elliptical solution = {0}: {1}".format(e_vals[ev_ind], a))

    conic_coefficients = np.concatenate((a, np.dot(T_no_det, a) // determinant_S3))
    rotated_euclidian_coefficients = list(_general_to_rotated_int(conic_coefficients, True))
    rotated_euclidian_coefficients[0] = (rotated_euclidian_coefficients[0][0] + x_mean,
                                         rotated_euclidian_coefficients[0][1] + y_mean)
    return rotated_euclidian_coefficients


def _general_to_rotated(conic_coeffs):
    """Transform from conic section format to general format.

    :param conic_coeffs: The six coefficients defining the ellipse as a conic shape in
     :math:`ax^2 + bxy + cy^2 + dx + ey + f = 0`.
    :type conic_coeffs: :py:class:`numpy.ndarray` or tuple
    :return: The general form for the ellipse. Returns tuple :math:`(x_c,\ y_c),\\ (a,\\ b),\\ \\theta`
        that fulfills the equation

        .. math::

            \\frac{((x-x_c)\\cos(\\theta) + (y-y_c)\\sin(\\theta))^2}{a^2} +
            \\frac{((x-x_c)\\sin(\\theta) - (y-y_c)\\sin(\\theta))^2}{b^2} = 1

    :rtype: tuple

    """
    a, b, c, d, e, f = conic_coeffs

    angle = np.arctan2(b, a - c) / 2

    # Obtaining angels using Givens rotations.
    # cos_2t, sin_2t = ma.Givens_rotation_double(b, a-c)
    # cos_2t, sin_2t = -sin_2t, cos_2t
    #
    # if cos_2t < 0:
    #     sin_t = np.sqrt((1 - cos_2t)/2)
    #     cos_t = sin_2t / (2 * sin_t)
    # else:
    #     cos_t = np.sqrt((cos_2t + 1)/2)
    #     sin_t = sin_2t / (2 * cos_t)
    cos_t = np.cos(angle)
    sin_t = np.sin(angle)

    a_prime = a * (cos_t ** 2) + b * cos_t * sin_t + c * (sin_t ** 2)
    c_prime = a * (sin_t ** 2) - b * cos_t * sin_t + c * (cos_t ** 2)
    d_prime = (d * cos_t) + (e * sin_t)
    e_prime = (-(d * sin_t)) + (e * cos_t)
    f_prime = f

    x_prime = (-d_prime) / (2 * a_prime)
    y_prime = (-e_prime) / (2 * c_prime)
    x = (x_prime * cos_t) - (y_prime * sin_t)
    y = (x_prime * sin_t) + (y_prime * cos_t)

    maj_axis = np.sqrt(((-4 * f_prime * a_prime * c_prime) +
                        (c_prime * (d_prime ** 2)) +
                        (a_prime * (e_prime ** 2))) / (4 * a_prime * (c_prime ** 2)))
    min_axis = np.sqrt(((-4 * f_prime * a_prime * c_prime) +
                        (c_prime * (d_prime ** 2)) +
                        (a_prime * (e_prime ** 2))) / (4 * (a_prime ** 2) * c_prime))

    if a_prime > c_prime:
        x_axis = min_axis
        y_axis = maj_axis
    else:
        x_axis = maj_axis
        y_axis = min_axis

    general_coeffs = (x, y), (x_axis, y_axis), angle

    if DEBUG:
        print("Cos: {0}, Sin: {1}".format(cos_t, sin_t))
        print("Conic form    = {0:.4f}x^2 + {1:.4f}xy + "
              "{2:.4f}y^2 + {3:.4f}x + {4:.4f}y + {5:.4f}".format(*conic_coeffs))
        print("Conic form 2  = {0:.4f}x^2 + {1:.4f}xy + "
            "{2:.4f}y^2 + {3:.4f}x + {4:.4f}y + {5:.4f}".format(a_prime, 0, c_prime,
                                                                d_prime, e_prime, f_prime))
        print("Elliptical form: Center = {0}, Radii = {1}, Angle = {2}\n".format(*general_coeffs))

    return general_coeffs


def _general_to_rotated_int(conic_coeffs, return_float=False):
    """Transform from conic section format to general format, in integer precision.

    :param conic_coeffs: The six coefficients defining the ellipse as a conic shape in
     :math:`ax^2 + bxy + cy^2 + dx + ey + f = 0`.
    :type conic_coeffs: :py:class:`numpy.ndarray` or tuple
    :return: The general form for the ellipse. Returns tuple :math:`(x_c,\ y_c),\\ (a,\\ b),\\ \\theta`
        that fulfills the equation

        .. math::

            \\frac{((x-x_c)\\cos(\\theta) + (y-y_c)\\sin(\\theta))^2}{a^2} +
            \\frac{((x-x_c)\\sin(\\theta) - (y-y_c)\\sin(\\theta))^2}{b^2} = 1

    :rtype: tuple

    """
    # These conic coefficients have a scaling e_norm = sqrt(a**2 + b**2 + c**2) from
    # the eigenvector not being normalized before turning it into conic coefficients.
    a, b, c, d, e, f = conic_coeffs
    angle = np.arctan2(b, (a - c)) / 2

    if b == 0:
        if (a-c) < 0:
            unity = 1
            cos_t = 0
            sin_t = 1
        else:
            unity = 1
            cos_t = 1
            sin_t = 0
    elif (a-c) == 0:
        if b < 0:
            unity = 50
            cos_t = 5
            sin_t = -5
        else:
            unity = 50
            cos_t = 5
            sin_t = 5
    else:
        # unity here becomes a value with scaling e_norm.
        cos_2t, sin_2t, unity = ma.Givens_rotation_int(b, (a-c))

        cos_2t, sin_2t = -sin_2t, cos_2t

        if cos_2t < 0:
            sin_t = ma.sqrt_int64(int((unity - cos_2t) >> 1), False)  # sin_t is scaled up sqrt(unity)
            cos_t = sin_2t // (2 * sin_t)                       # cos_t also becomes sqrt(unity) scaled
            # Now that we have cos and sin values we establish a new scaling with
            # unity value obtained through the sin^2(x) + cos^2(x) = 1 trigonometric identity.
            unity = sin_t ** 2 + cos_t ** 2

        else:
            cos_t = ma.sqrt_int64(int((cos_2t + unity) >> 1), False)  # cos_t is scaled up sqrt(unity)
            sin_t = sin_2t // (2 * cos_t)                       # sin_t also becomes sqrt(unity) scaled
            # Now that we have cos and sin values we establish a new scaling with
            # unity value obtained through the sin^2(x) + cos^2(x) = 1 trigonometric identity.
            unity = sin_t ** 2 + cos_t ** 2

    a_prime = a * (cos_t ** 2) + b * cos_t * sin_t + c * (sin_t ** 2)  # Is unity + e_norm scaled
    c_prime = a * (sin_t ** 2) - b * cos_t * sin_t + c * (cos_t ** 2)  # Is unity + e_norm scaled
    d_prime = (d * cos_t) + (e * sin_t)  # Is sqrt(unity) + e_norm scaled
    e_prime = (-(d * sin_t)) + (e * cos_t)  # Is sqrt(unity) + e_norm scaled
    f_prime = f  # Is e_norm scaled only.

    x_prime_num = (-d_prime)
    x_prime_denom = (2 * a_prime)

    y_prime_num = (-e_prime)
    y_prime_denom = (2 * c_prime)

    if return_float:
        # At sub-pixel precision, treat values as floats.
        x = (((x_prime_num * cos_t) / x_prime_denom) - ((y_prime_num * sin_t) / y_prime_denom))
        y = (((x_prime_num * sin_t) / x_prime_denom) + ((y_prime_num * cos_t) / y_prime_denom))
    else:
        # At pixel precision, perform integer division.
        x = int(((x_prime_num * cos_t) // x_prime_denom) -
                ((y_prime_num * sin_t) // y_prime_denom))
        y = int(((x_prime_num * sin_t) // x_prime_denom) +
                ((y_prime_num * cos_t) // y_prime_denom))

    sqrt_unity = ma.sqrt_int64(int(unity))

    a_prime = a_prime // (unity)
    c_prime = c_prime // (unity)
    d_prime = d_prime // (sqrt_unity)
    e_prime = e_prime // (sqrt_unity)

    if return_float:
        # At sub-pixel precision, use float divison and square root.
        numerator = ((-4 * f_prime * a_prime * c_prime) +
                     (c_prime * (d_prime ** 2)) +
                     (a_prime * (e_prime ** 2)))
        denominator_major = (4 * a_prime * (c_prime ** 2))
        tmp_axis = ma.sqrt_int64((numerator // denominator_major), True)
        maj_axis = tmp_axis[0] + tmp_axis[1]
        denominator_minor = (4 * (a_prime ** 2) * c_prime)
        tmp_axis = ma.sqrt_int64(numerator // denominator_minor, True)
        min_axis = tmp_axis[0] + tmp_axis[1]
    else:
        # At pixel precision, use integer division and perform
        # fixed point square root.
        numerator = ((-4 * f_prime * a_prime * c_prime) +
                     (c_prime * (d_prime ** 2)) +
                     (a_prime * (e_prime ** 2))) * unity
        denominator_major = (4 * a_prime * (c_prime ** 2))
        maj_axis = ma.sqrt_int64(int(numerator // denominator_major))
        denominator_minor = (4 * (a_prime ** 2) * c_prime)
        min_axis = ma.sqrt_int64(int(numerator // denominator_minor))

    if a_prime > c_prime:
        x_axis = min_axis
        y_axis = maj_axis
    else:
        x_axis = maj_axis
        y_axis = min_axis

    general_coeffs = (x, y), (x_axis, y_axis), angle

    if DEBUG:
        print("Cos: {0} => {1}, Sin: {2} => {3}".format(cos_t, cos_t / np.sqrt(unity), sin_t, sin_t / np.sqrt(unity)))
        print("Conic form    = {0:.4f}x^2 + {1:.4f}xy + "
              "{2:.4f}y^2 + {3:.4f}x + {4:.4f}y + {5:.4f}".format(*conic_coeffs))
        print("Conic form 2  = {0:.4f}x^2 + {1:.4f}xy + "
            "{2:.4f}y^2 + {3:.4f}x + {4:.4f}y + {5:.4f}".format(
            a_prime * unity, 0, c_prime * unity, d_prime * np.sqrt(unity), e_prime * np.sqrt(unity), f_prime))
        print("Elliptical form: Center = {0}, Radii = {1}, Angle = {2}\n".format(*general_coeffs))

    return general_coeffs


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
    S3 = np.array([S[3, 3], S[3, 4], S[3, 5], S[4, 4], S[4, 5], S[5, 5]])
    S3_inv = mo.inverse_symmetric_3by3_double(S3).reshape((3, 3))
    S2 = S[:3, 3:]
    T = -np.dot(S3_inv, S2.T)
    M_term_2 = np.dot(S2, T)
    M = S[:3, :3] + M_term_2
    M[[0, 2], :] /= 2
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
