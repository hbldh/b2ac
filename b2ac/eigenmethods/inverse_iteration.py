#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`inverse_iteration` -- Eigenvector finding
===============================================

.. module:: inverse_iteration
   :platform: Unix, Windows
   :synopsis: Inverse iteration to find eigenvectors.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2013-09-09, 16:34

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.compat import *
import b2ac.matrix.matrix_operations as mo
import b2ac.matrix.fixed_point as fp


def inverse_iteration_for_eigenvector_double(A, eigenvalue, n_iterations=1):
    """Performs a series of inverse iteration steps with a known
    eigenvalue to produce its eigenvector.

    :param A: The 3x3 matrix to which the eigenvalue belongs.
    :type A: :py:class:`numpy.ndarray`
    :param eigenvalue: One eigenvalue of the matrix A.
    :type eigenvalue: float
    :param n_iterations: Number of iterations to perform the multiplication
     with the inverse. For a accurate eigenvalue, one iteration is enough
     for a ~1e-6 correct eigenvector. More than five is usually unnecessary.
    :type n_iterations: int
    :return: The eigenvector of this matrix and eigenvalue combination.
    :rtype: :py:class:`numpy.ndarray`

    """
    A = np.array(A, 'float')
    # Subtract the eigenvalue from the diagonal entries of the matrix.
    # N_POLYPOINTS.B. Also slightly perturb the eigenvalue so the matrix will
    # not be so close to singular!
    for k in xrange(A.shape[0]):
        A[k, k] -= eigenvalue + 0.001
    # Obtain the inverse of the matrix.
    A_inv = mo.inverse_3by3_double(A).reshape((3, 3))
    # Instantiate the eigenvector to iterate with.
    eigenvector = np.ones((A.shape[0], ), 'float')
    eigenvector /= np.linalg.norm(eigenvector)
    # Perform the desired number of iterations.
    for k in xrange(n_iterations):
        eigenvector = np.dot(A_inv, eigenvector)
        eigenvector /= np.linalg.norm(eigenvector)

    if np.any(np.isnan(eigenvector)) or np.any(np.isinf(eigenvector)):
        print("Nan and/or Infs in eigenvector!")

    if (eigenvector[0] < 0) and (eigenvector[2] < 0):
        eigenvector = -eigenvector
    return eigenvector


def inverse_iteration_for_eigenvector_int(A, eigenvalue):
    """Performs a series of inverse iteration steps with a known
    eigenvalue to produce its eigenvector.

    :param A: The 3x3 matrix to which the eigenvalue belongs.
    :type A: :py:class:`numpy.ndarray`
    :param eigenvalue: One approximate eigenvalue of the matrix A.
    :type eigenvalue: int
    :return: The eigenvector of this matrix and eigenvalue combination.
    :rtype: :py:class:`numpy.ndarray`

    """
    A = np.array(A, 'int64')

    # Subtract the eigenvalue from the diagonal entries of the matrix.
    for k in xrange(A.shape[0]):
        A[k, k] -= eigenvalue
    A, scale = fp.scale_64bit_matrix(A)

    # Obtain the inverse of the matrix.
    adj_A = mo.inverse_3by3_int64(A.flatten(), False)
    eigenvector = adj_A.reshape((3, 3)).sum(1)
    eigenvector, scale = fp.scale_64bit_vector(eigenvector)

    e_norm = int(np.sqrt((eigenvector ** 2).sum()))
    if (eigenvector[0] < 0) and (eigenvector[2] < 0):
        eigenvector = -eigenvector
    return eigenvector, e_norm
