#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: matrix_operations
   :platform: Unix, Windows
   :synopsis: Operations on matrices in different arithmetic.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2013-09-09, 16:32

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np


def inverse_symmetric_3by3_double(M):
    """C style inverse of a symmetric, flattened 3 by 3 matrix, returning full matrix in
    floating-point, also in flattened format.

    :param M: The matrix to find inverse to. Assumes array with shape (6,).
    :type M: :py:class:`numpy.ndarray`
    :return: The inverse matrix, flattened.
    :rtype: :py:class:`numpy.ndarray`

    """

    determinant = 0
    adj_M = np.zeros((9,), dtype='float')

    # First row of adjugate matrix
    adj_M[0] = (M[3] * M[5] - (M[4] ** 2))  # Det #0
    adj_M[1] = -(M[1] * M[5] - M[4] * M[2])  # Det #1
    adj_M[2] = (M[1] * M[4] - M[3] * M[2])  # Det #2

    # Second row of adjugate matrix
    adj_M[3] = adj_M[1]
    adj_M[4] = (M[0] * M[5] - (M[2] ** 2))
    adj_M[5] = -(M[0] * M[4] - M[1] * M[2])

    # Third row of adjugate matrix
    adj_M[6] = adj_M[2]
    adj_M[7] = adj_M[5]
    adj_M[8] = (M[0] * M[3] - (M[1] ** 2))

    determinant += M[0] * adj_M[0]
    determinant += M[1] * adj_M[1]  # Using addition since minus is integrated in adjugate matrix.
    determinant += M[2] * adj_M[2]

    return adj_M / determinant


def inverse_symmetric_3by3_int64(M):
    """C style inverse of a symmetric, flattened 3 by 3 matrix, represented
    by the six unique values

    .. code-block:: python

        np.array([S[0, 0], S[0, 1], S[0, 2], S[1, 1], S[1, 2], S[2, 2]])

    Returns the flattened full adjugate matrix and the determinant, s.t.

    .. math::

        M^{-1} = \\frac{1}{\\det(M)} \\cdot \\text{adj}(M).


    For integer matrices, this then returns exact results by avoiding
    to apply the division.

    :param M: The matrix to find inverse to. Assumes array with shape (6,).
    :type M: :py:class:`numpy.ndarray`
    :return: The adjugate flattened matrix and the determinant.
    :rtype: tuple

    """

    determinant = 0
    adj_M = np.zeros((9,), dtype='int64')

    # First row of adjugate matrix
    adj_M[0] = (M[3] * M[5] - (M[4] ** 2))  # Det #0
    adj_M[1] = -(M[1] * M[5] - M[4] * M[2])  # Det #1
    adj_M[2] = (M[1] * M[4] - M[3] * M[2])  # Det #2

    # Second row of adjugate matrix
    adj_M[3] = adj_M[1]
    adj_M[4] = (M[0] * M[5] - (M[2] ** 2))
    adj_M[5] = -(M[0] * M[4] - M[1] * M[2])

    # Third row of adjugate matrix
    adj_M[6] = adj_M[2]
    adj_M[7] = adj_M[5]
    adj_M[8] = (M[0] * M[3] - (M[1] ** 2))

    determinant += np.int64(M[0]) * np.int64(adj_M[0])
    determinant += np.int64(M[1]) * np.int64(adj_M[1])  # Using addition since minus is integrated in adjugate matrix.
    determinant += np.int64(M[2]) * np.int64(adj_M[2])

    return adj_M, determinant


def inverse_3by3_double(M):
    """C style inverse of a flattened 3 by 3 matrix, returning full matrix in
    floating-point, also in flattened format.

    :param M: The matrix to find inverse to.
    :type M: :py:class:`numpy.ndarray`
    :return: The inverse matrix.
    :rtype: :py:class:`numpy.ndarray`

    """
    if len(M.shape) > 1:
        M = M.flatten()

    M = np.array(M, 'float')

    determinant = 0.
    adj_M = np.zeros((9,), 'float')

    # First row of adjugate matrix
    adj_M[0] = (M[4] * M[8] - M[7] * M[5])  # Det #0
    adj_M[1] = -(M[1] * M[8] - M[7] * M[2])
    adj_M[2] = (M[1] * M[5] - M[4] * M[2])

    # Second row of adjugate matrix
    adj_M[3] = -(M[3] * M[8] - M[6] * M[5])  # Det #1
    adj_M[4] = (M[0] * M[8] - M[6] * M[2])
    adj_M[5] = -(M[0] * M[5] - M[3] * M[2])

    # Third row of adjugate matrix
    adj_M[6] = (M[3] * M[7] - M[6] * M[4])  # Det #2
    adj_M[7] = -(M[0] * M[7] - M[6] * M[1])
    adj_M[8] = (M[0] * M[4] - M[3] * M[1])

    determinant += M[0] * adj_M[0]
    determinant += M[1] * adj_M[3]  # Using addition since minus is integrated in adjugate matrix.
    determinant += M[2] * adj_M[6]

    return (adj_M / determinant)


def inverse_3by3_int64(M, return_determinant=True):
    """C style inverse of a flattened 3 by 3 matrix.

    Returns the full adjugate matrix and the determinant, s.t.

    .. math::

        M^{-1} = \\frac{1}{\\det(M)} \\cdot \\text{adj}(M).


    For integer matrices, this then returns exact results by avoiding
    to apply the division.

    :param M: The matrix to find inverse to.
    :type M: :py:class:`numpy.ndarray`
    :return: The adjugate flattened matrix and the determinant.
    :rtype: tuple

    """
    if len(M.shape) > 1:
        M = M.flatten()

    determinant = np.int64(0)
    adj_M = np.zeros((9,), 'int64')

    # First row of adjugate matrix
    adj_M[0] = (M[4] * M[8] - M[7] * M[5])  # Det #0
    adj_M[1] = -(M[1] * M[8] - M[7] * M[2])
    adj_M[2] = (M[1] * M[5] - M[4] * M[2])

    # Second row of adjugate matrix
    adj_M[3] = -(M[3] * M[8] - M[6] * M[5])  # Det #1
    adj_M[4] = (M[0] * M[8] - M[6] * M[2])
    adj_M[5] = -(M[0] * M[5] - M[3] * M[2])

    # Third row of adjugate matrix
    adj_M[6] = (M[3] * M[7] - M[6] * M[4])  # Det #2
    adj_M[7] = -(M[0] * M[7] - M[6] * M[1])
    adj_M[8] = (M[0] * M[4] - M[3] * M[1])

    if return_determinant:
        if ((np.log2(np.abs(M[0])) + np.log2(np.abs(adj_M[0]))) > 63 or
                (np.log2(np.abs(M[1])) + np.log2(np.abs(adj_M[1]))) > 63 or
                (np.log2(np.abs(M[2])) + np.log2(np.abs(adj_M[6]))) > 63):
            print("inverse_3by3_int64: Overflow in determinant calculation!")
            determinant += int(M[0]) * int(adj_M[0])
            determinant += int(M[1]) * int(adj_M[3])  # Using addition since minus is integrated in adjugate matrix.
            determinant += int(M[2]) * int(adj_M[6])
        else:
            determinant += np.int64(M[0]) * np.int64(adj_M[0])
            determinant += np.int64(M[1]) * np.int64(adj_M[3])  # Using addition since minus is integrated in adjugate matrix.
            determinant += np.int64(M[2]) * np.int64(adj_M[6])
        return adj_M, determinant
    else:
        return adj_M


def matrix_add_symmetric(M, M_sym):
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
