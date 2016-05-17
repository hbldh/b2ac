#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: matrix
   :platform: Unix, Windows
   :synopsis: Operations on matrices.

.. moduleauthor:: hbldh <henrik.blidh@nedomull.com>

Created on 2013-05-15, 10:45

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.compat import *


def inverse_symmetric_3by3_double(M):
    """C style inverse of a symmetric, flattened 3 by 3 matrix.

    Returns the adjoint matrix and the determinant, s.t.

    .. math::

        M^{-1} = \\frac{1}{\\det(M)} \\cdot \\text{adj}(M).


    For integer matrices, this then returns exact results by avoiding
    to apply the division.

    :param M: The matrix to find inverse to. Assumes array with shape (6,).
    :type M: :py:class:`numpy.ndarray`
    :return: The inverse matrix, flattened.
    :rtype: :py:class:`numpy.ndarray`

    """

    determinant = 0
    inverse = np.zeros((9,), dtype='float')

    # First row of adjunct matrix
    inverse[0] = (M[3] * M[5] - (M[4] ** 2))  # Det #0
    inverse[1] = -(M[1] * M[5] - M[4] * M[2])  # Det #1
    inverse[2] = (M[1] * M[4] - M[3] * M[2])  # Det #2

    # Second row of adjunct matrix
    inverse[3] = inverse[1]
    inverse[4] = (M[0] * M[5] - (M[2] ** 2))
    inverse[5] = -(M[0] * M[4] - M[1] * M[2])

    # Third row of adjunct matrix
    inverse[6] = inverse[2]
    inverse[7] = inverse[5]
    inverse[8] = (M[0] * M[3] - (M[1] ** 2))

    determinant += M[0] * inverse[0]
    determinant += M[1] * inverse[1]  # Using addition since minus is integrated in adjunct matrix.
    determinant += M[2] * inverse[2]

    for i in xrange(len(inverse)):
        inverse[i] /= determinant

    return inverse


def inverse_symmetric_3by3_int(M):
    """C style inverse of a symmetric, flattened 3 by 3 matrix.

    Returns the adjoint matrix and the determinant, s.t.

    .. math::

        M^{-1} = \\frac{1}{\\det(M)} \\cdot \\text{adj}(M).


    For integer matrices, this then returns exact results by avoiding
    to apply the division.

    :param M: The matrix to find inverse to. Assumes array with shape (6,).
    :type M: :py:class:`numpy.ndarray`
    :return: The adjoint flattened matrix and the determinant.
    :rtype: tuple

    """

    determinant = 0
    adj_M = np.zeros((9,), dtype='int32')

    # First row of adjunct matrix
    adj_M[0] = (M[3] * M[5] - (M[4] ** 2))  # Det #0
    adj_M[1] = -(M[1] * M[5] - M[4] * M[2])  # Det #1
    adj_M[2] = (M[1] * M[4] - M[3] * M[2])  # Det #2

    # Second row of adjunct matrix
    adj_M[3] = adj_M[1]
    adj_M[4] = (M[0] * M[5] - (M[2] ** 2))
    adj_M[5] = -(M[0] * M[4] - M[1] * M[2])

    # Third row of adjunct matrix
    adj_M[6] = adj_M[2]
    adj_M[7] = adj_M[5]
    adj_M[8] = (M[0] * M[3] - (M[1] ** 2))

    determinant += np.int64(M[0]) * np.int64(adj_M[0])
    determinant += np.int64(M[1]) * np.int64(adj_M[1])  # Using addition since minus is integrated in adjunct matrix.
    determinant += np.int64(M[2]) * np.int64(adj_M[2])

    return adj_M, determinant


def inverse_3by3_int(M):
    """C style inverse of a flattened 3 by 3 matrix.

    Returns the adjoint matrix and the determinant, s.t.

    .. math::

        M^{-1} = \\frac{1}{\\det(M)} \\cdot \\text{adj}(M).


    For integer matrices, this then returns exact results by avoiding
    to apply the division.

    :param M: The matrix to find inverse to.
    :type M: :py:class:`numpy.ndarray`
    :return: The adjoint flattened matrix and the determinant.
    :rtype: tuple

    """
    if len(M.shape) > 1:
        M = M.flatten()

    determinant = 0
    adj_M = np.zeros((9,), 'int')

    # First row of adjunct matrix
    adj_M[0] = (M[4] * M[8] - M[7] * M[5])  # Det #0
    adj_M[1] = -(M[1] * M[8] - M[7] * M[2])
    adj_M[2] = (M[1] * M[5] - M[4] * M[2])

    # Second row of adjunct matrix
    adj_M[3] = -(M[3] * M[8] - M[6] * M[5])  # Det #1
    adj_M[4] = (M[0] * M[8] - M[6] * M[2])
    adj_M[5] = -(M[0] * M[5] - M[3] * M[2])

    # Third row of adjunct matrix
    adj_M[6] = (M[3] * M[7] - M[6] * M[4])  # Det #2
    adj_M[7] = -(M[0] * M[7] - M[6] * M[1])
    adj_M[8] = (M[0] * M[4] - M[3] * M[1])

    determinant += np.int64(M[0]) * np.int64(adj_M[0])
    determinant += np.int64(M[1]) * np.int64(adj_M[3])  # Using addition since minus is integrated in adjunct matrix.
    determinant += np.int64(M[2]) * np.int64(adj_M[6])

    return adj_M, determinant


def inverse_3by3_double(M):
    """C style inverse of a flattened 3 by 3 matrix.

    .. math::

        M^{-1} = \\frac{1}{\\det(M)} \\cdot \\text{adj}(M).


    For integer matrices, this then returns exact results by avoiding
    to apply the division.

    :param M: The matrix to find inverse to.
    :type M: :py:class:`numpy.ndarray`
    :return: The inverse matrix.
    :rtype: :py:class:`numpy.ndarray`

    """
    if len(M.shape) > 1:
        M = M.flatten()

    determinant = 0
    adj_M = np.zeros((9,), 'int')

    # First row of adjunct matrix
    adj_M[0] = (M[4] * M[8] - M[7] * M[5])  # Det #0
    adj_M[1] = -(M[1] * M[8] - M[7] * M[2])
    adj_M[2] = (M[1] * M[5] - M[4] * M[2])

    # Second row of adjunct matrix
    adj_M[3] = -(M[3] * M[8] - M[6] * M[5])  # Det #1
    adj_M[4] = (M[0] * M[8] - M[6] * M[2])
    adj_M[5] = -(M[0] * M[5] - M[3] * M[2])

    # Third row of adjunct matrix
    adj_M[6] = (M[3] * M[7] - M[6] * M[4])  # Det #2
    adj_M[7] = -(M[0] * M[7] - M[6] * M[1])
    adj_M[8] = (M[0] * M[4] - M[3] * M[1])

    determinant += M[0] * adj_M[0]
    determinant += M[1] * adj_M[3]  # Using addition since minus is integrated in adjunct matrix.
    determinant += M[2] * adj_M[6]

    return adj_M / determinant


def add_symmetric_matrix(M, M_sym):
    """Add a regular matrix and a symmetric one.

    :param M: A [3x3] matrix to add with symmetric matrix.
    :type M: :py:class:`numpy.ndarray`
    :param M_sym: A [6x1] array to add with M.
    :type M_sym: :py:class:`numpy.ndarray`
    :return: The sum of the two matrices.
    :rtype: :py:class:`numpy.ndarray`

    """
    M[0, 0] += M_sym[0]
    M[0, 1] += M_sym[1]
    M[1, 0] += M_sym[1]
    M[0, 2] += M_sym[2]
    M[2, 0] += M_sym[2]

    M[1, 1] += M_sym[3]
    M[1, 2] += M_sym[4]
    M[2, 1] += M_sym[4]

    M[2, 2] += M_sym[5]

    return M
