#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: matrix_algorithms
   :platform: Unix, Windows
   :synopsis: Algorithms operating on matrices e.g. factorisation.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2013-09-09, 16:32

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.compat import *


def QR_factorisation_Householder_double(A):
    """Perform QR factorisation in double floating-point precision.

    :param A: The matrix to factorise.
    :type A: :py:class:`numpy.ndarray`
    :returns: The matrix Q and the matrix R.
    :rtype: tuple

    """
    A = np.array(A, 'float')

    n, m = A.shape
    V = np.zeros_like(A, 'float')
    for k in xrange(n):
        V[k:, k] = A[k:, k].copy()
        V[k, k] += np.sign(V[k, k]) * np.linalg.norm(V[k:, k], 2)
        V[k:, k] /= np.linalg.norm(V[k:, k], 2)
        A[k:, k:] -= 2 * np.outer(V[k:, k], np.dot(V[k:, k], A[k:, k:]))
    R = np.triu(A[:n, :n])

    Q = np.eye(m, n)
    for k in xrange((n - 1), -1, -1):
        Q[k:, k:] -= np.dot((2 * (np.outer(V[k:, k], V[k:, k]))), Q[k:, k:])
    return Q, R


def QR_factorisation_Givens_double(A):

    n, m = A.shape
    R = np.array(A, dtype='float')
    Q = np.eye(n)
    for i in xrange(m - 1):
        for j in xrange(n - 1, i, -1):
            G = Givens_rotation_matrix_double(R[j - 1, i], R[j, i])
            R[(j - 1):(j + 1), :] = np.dot(G, R[(j - 1):(j + 1), :])
            Q[(j - 1):(j + 1), :] = np.dot(G, Q[(j - 1):(j + 1), :])
    return Q.T, np.triu(R)


def QR_factorisation_int64(A):
    raise NotImplementedError("Never implemented this...")


def convert_to_Hessenberg_Givens_double(A):
    n, m = A.shape
    A = np.array(A, 'float')
    for i in xrange(m):
        for k in xrange(m - 1, i + 1, -1):
            c, s = Givens_rotation_double(A[k - 1, i], A[k, i])
            for j in xrange(m):
                tau_1 = A[k-1, j]
                tau_2 = A[k, j]
                A[k-1, j] = ((tau_1 * c) - (tau_2 * s))
                A[k, j] = ((tau_1 * s) + (tau_2 * c))
            for j in xrange(n):
                tau_1 = A[j, k-1]
                tau_2 = A[j, k]
                A[j, k-1] = ((tau_1 * c) - (tau_2 * s))
                A[j, k] = ((tau_1 * s) + (tau_2 * c))
    return np.triu(A, -1)


def convert_to_Hessenberg_Givens_int(A):
    m, n = A.shape
    A = np.array(A, 'int64')
    for i in xrange(m):
        for k in xrange(m - 1, i + 1, -1):
            c_n, s_n, denominator = Givens_rotation_int(A[k - 1, i], A[k, i])
            for j in xrange(m):
                tau_1 = A[k-1, j]
                tau_2 = A[k, j]
                A[k-1, j] = ((tau_1 * c_n) - (tau_2 * s_n)) // denominator
                A[k, j] = ((tau_1 * s_n) + (tau_2 * c_n)) // denominator
            for j in xrange(n):
                tau_1 = A[j, k-1]
                tau_2 = A[j, k]
                A[j, k-1] = ((tau_1 * c_n) - (tau_2 * s_n)) // denominator
                A[j, k] = ((tau_1 * s_n) + (tau_2 * c_n)) // denominator

    return np.triu(A, -1)


def convert_to_Hessenberg_double(A):
    """ Tridiagonalize a square matrix to upper Hessenberg form using Householder reflections.

    Costs 10/3 * n^3 + O(n^2) operations, where n is the size of the matrix.

    .. code-block:: matlab

        function A = mytridiag(A)

        [m,n] = size(A);
        if (m ~= n)
            error('Input matrix is not square.')
        end

        for k = 1:(m - 2)
            vk = A((k+1):m,k);
            vk(1) = vk(1) + sign(vk(1)) * norm(vk);
            vk = vk / norm(vk);
            A((k+1):m,k:m) = A((k+1):m,k:m) - ...
                                2 * vk * (vk' * A((k+1):m,k:m));
            A(1:m,(k+1):m) = A(1:m,(k+1):m) - ...
                               2 * (A(1:m,(k+1):m) * vk) * vk';
        end

        end

    :param A: The matrix to convert.
    :type A: :py:class:`numpy.ndarray`
    :return: The Hessenberg matrix.
    :rtype: :py:class:`numpy.ndarray`

    """
    m, n = A.shape
    A = np.array(A, 'float')

    for k in xrange(m - 1):
        vk = A[(k + 1):, k].copy()
        vk[0] += np.sign(vk[0]) * np.linalg.norm(vk, 2)
        vk /= np.linalg.norm(vk, 2)
        A[(k + 1):, k:] -= 2 * np.outer(vk, np.dot(vk, A[(k + 1):, k:]))
        A[:, (k + 1):] -= 2 * np.outer(np.dot(A[:, (k + 1):], vk), vk)

    return np.triu(A, -1)


def convert_to_Hessenberg_symmetric_double(A):
    """Tridiagonalize a symmetric square matric to Hessenberg form using Householder reflections.

    Costs 4/3 * n^3 + O(n^2) operations.

    .. code-block:: matlab

        function A = mytridiagturbo(A)

        [m,n] = size(A);
        if (m ~= n)
            error('Input matrix is not square.')
        end

        for k = 1:(m - 2)
            vk = A((k+1):m,k);
            beta = sign(vk(1))*norm(vk);
            vk(1) = vk(1) + beta;
            vk = vk / norm(vk);

            % Calculate A_hat for current iteration.
            u = -2*(A((k+1:m),(k+1:m)) * vk);
            d = -vk' * u;
            w = u + d * vk;
            wv = w * vk';
            A((k+1:m),(k+1:m)) = A((k+1:m),(k+1:m)) + wv + wv';

            % Adjust tridiagonality.
            A((k+1):m,k) = [-beta; zeros(m-k-1,1)];
            A(k,(k+1):m) = A((k+1):m,k)';
         end

    :param A: The matrix to convert.
    :type A: :py:class:`numpy.ndarray`
    :return: The Hessenberg matrix.
    :rtype: :py:class:`numpy.ndarray`

    """

    m, n = A.shape
    A = np.array(A, 'float')

    for k in xrange(m - 1):
        vk = A[(k + 1):, k].copy()
        beta = np.sign(vk[0]) * np.linalg.norm(vk)
        vk[0] += beta
        vk = vk / np.linalg.norm(vk)

        # Calculate A_hat for current iteration.
        u = -2 * np.dot(A[(k + 1):, (k + 1):], vk)
        d = np.dot(-vk, u)
        w = u + d * vk
        wv = np.outer(w, vk)
        A[(k + 1):, (k + 1):] += wv + wv.T

        # Adjust tridiagonality.
        t = np.zeros((m - (k + 1), ), dtype='float')
        t[0] = -beta
        A[(k + 1):, k] = t
        A[k, (k + 1):] = t
    return A


def Givens_rotation_matrix_double(a, b):
    c, s = Givens_rotation_double(a, b)
    return np.array([[c, -s], [s, c]])


def Givens_rotation_double(a, b):
    # Working with actual trigonometric functions
    # angle = np.arctan2(-a, b)
    # c = np.cos(angle)
    # s = np.sin(angle)

    # Using naive definitions
    # root = np.sqrt(a ** 2 + b ** 2)
    # c = a / root
    # s = -b / root

    # Using Matrix Computations solution
    if b == 0:
        c = 1.0
        s = 0.0
    else:
        if np.abs(b) > np.abs(a):
            tau = - a / b
            s = 1 / (np.sqrt(1 + tau ** 2))
            c = s * tau
        else:
            tau = - b / a
            c = 1 / (np.sqrt(1 + tau ** 2))
            s = c * tau

    return c, s


def Givens_rotation_int(a, b):
    # Using naive definitions
    if np.abs(a) > 0 and np.log2(np.abs(a)) > 32:
        print("Too large value a = {0} for squaring and rooting!".format(a))
    if np.abs(b) > 0 and np.log2(np.abs(b)) > 32:
        print("Too large value b = {0} for squaring and rooting!".format(b))
    root = sqrt_int64(a ** 2 + b ** 2)
    return a, -b, root


def sqrt_int64(x, return_fraction=False):

    BITSPERLONG = 64
    BPL_HALF = BITSPERLONG // 2

    def top2bits(x):
        return (x & (3 << (BITSPERLONG - 2))) >> (BITSPERLONG - 2)

    a = 0
    r = 0

    for i in xrange(BITSPERLONG):

        r = (r << 2) + top2bits(int(x))
        x <<= 2
        a <<= 1
        e = (a << 1) + 1
        if r >= e:
            r -= e
            a += 1
    if return_fraction:
        return a >> BPL_HALF, (a & (2**32 - 1)) / (2**32 - 1)
    else:
        return a >> BPL_HALF
