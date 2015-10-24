#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`conversion`
==================

.. module:: conversion
    :platform: Unix, Windows
    :synopsis:

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-10-07

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

import b2ac.matrix.matrix_algorithms as ma


def conic_to_general_reference(conic_coeffs):
    """Transform from conic section format to general format.

    :param conic_coeffs: The six coefficients defining the ellipse as a conic shape.
    :type conic_coeffs: :py:class:`numpy.ndarray` or tuple
    :param verbose: If debug printout is desired.
    :type verbose: bool
    :return:
    :rtype: tuple

    """
    a, b, c, d, e, f = conic_coeffs

    angle = np.arctan2(b, a - c) / 2

    cos_theta = np.cos(angle)  # np.sqrt((1 + np.cos(2*angle)) / 2)
    sin_theta = np.sin(angle)  # np.sqrt((1 - np.cos(2*angle)) / 2)

    a_prime = a * (cos_theta ** 2) + (b * cos_theta * sin_theta) + c * (sin_theta ** 2)
    c_prime = a * (sin_theta ** 2) - (b * cos_theta * sin_theta) + c * (cos_theta ** 2)
    d_prime = (d * cos_theta) + (e * sin_theta)
    e_prime = (-(d * sin_theta)) + (e * cos_theta)
    f_prime = f

    x_prime = (-d_prime) / (2 * a_prime)
    y_prime = (-e_prime) / (2 * c_prime)
    major_axis = np.sqrt(
        ((-4 * f_prime * a_prime * c_prime) + (c_prime * (d_prime ** 2)) + (a_prime * (e_prime ** 2))) /
        (4 * a_prime * (c_prime ** 2)))
    minor_axis = np.sqrt(
        ((-4 * f_prime * a_prime * c_prime) + (c_prime * (d_prime ** 2)) + (a_prime * (e_prime ** 2))) /
        (4 * (a_prime ** 2) * c_prime))

    if a_prime > c_prime:
        angle += np.pi / 2

    x = x_prime * cos_theta - y_prime * sin_theta
    y = x_prime * sin_theta + y_prime * cos_theta

    return [x, y], [major_axis, minor_axis], angle


def conic_to_general_1(conic_coeffs):
    """Transform from conic section format to general format.

    Adopted from http://math.stackexchange.com/a/423272

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

    # Obtaining angles using Givens rotations.
    # TODO: Evaluate this method of calculating angles.
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
    if np.isnan(maj_axis):
        maj_axis = 0.0
    if np.isnan(min_axis):
        min_axis = 0.0

    if a_prime > c_prime:
        x_axis = min_axis
        y_axis = maj_axis
    else:
        x_axis = maj_axis
        y_axis = min_axis

    general_coeffs = [x, y], [x_axis, y_axis], angle

    return general_coeffs


def conic_to_general_2(conic_coeffs):
    """Transform from conic section format to general format.

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

    return [x, y], [maj_axis, min_axis], angle


def conic_to_general_int(conic_coeffs, return_float=False, verbose=False):
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

    a_prime //= unity
    c_prime //= unity
    d_prime //= sqrt_unity
    e_prime //= sqrt_unity

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

    general_coeffs = [x, y], [x_axis, y_axis], angle

    if verbose:
        print("Cos: {0} => {1}, Sin: {2} => {3}".format(cos_t, cos_t / np.sqrt(unity), sin_t, sin_t / np.sqrt(unity)))
        print("Conic form    = {0:.4f}x^2 + {1:.4f}xy + "
              "{2:.4f}y^2 + {3:.4f}x + {4:.4f}y + {5:.4f}".format(*conic_coeffs))
        print("Conic form 2  = {0:.4f}x^2 + {1:.4f}xy + "
              "{2:.4f}y^2 + {3:.4f}x + {4:.4f}y + {5:.4f}".format(
            a_prime * unity, 0, c_prime * unity, d_prime * np.sqrt(unity), e_prime * np.sqrt(unity), f_prime))
        print("Elliptical form: Center = {0}, Radii = {1}, Angle = {2}\n".format(*general_coeffs))

    return general_coeffs
