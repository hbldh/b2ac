#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`eigenmethods` -- Eigenvalue and Eigenvector methods
=========================================================

.. module:: eigenmethods
   :platform: Unix, Windows
   :synopsis: A plethora of eigensystem solvers.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2013-09-09, 16:46

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

import b2ac.matrix.matrix_algorithms as ma
import b2ac.matrix.fixed_point as fp

QR_ALGORITHM_TOLERANCE = 1e-10  #sys.float_info.epsilon


def QR_algorithm(A):
    """The straight QR algorithm for finding eigenvalues.

    N.B. Has very slow convergence. Use shifted versions instead.

    .. code-block:: matlab

        function [myeig tmvec] = QRalgorithm(A)

        %------ Phase 1 ------
        T = mytridiag(A);

        %------ Phase 2 ------
        m = length(A);
        myeig = zeros(m,1);
        % Size of tmvec is not known in advance.
        % Estimating an initial size.
        tmvec = zeros(m^2,1);
        counter = 0;

        while (m > 1)
            counter = counter + 1;
            [Q,R] = myQR(T);
            T = R*Q;
            tmvec(counter) = abs(T(m,m-1));
            if (tmvec(counter) < 1e-8)
                myeig(m) = T(m,m);
                m = m - 1;
                T = T(1:m,1:m);
            end
        end

        myeig(1) = T;
        tmvec = tmvec(1:counter);

        end

    :param A: A square matrix to find eigenvalues in.
    :type A: :py:class:`numpy.ndarray`
    :return: The eigenvalues.
    :rtype: list

    """

    # First, tridiagonalize the input matrix to a Hessenberg matrix.
    T = ma.convert_to_Hessenberg_double(A.copy())

    m, n = T.shape
    if n != m:
        raise np.linalg.LinAlgError("Array must be square.")

    convergence_measure = []
    eigenvalues = np.zeros((n, ), dtype='float')
    n -= 1

    while n > 0:
        # Perform QR decomposition.
        Q, R = ma.QR_factorisation_Householder_double(T)
        # Multiply R and Q.
        T = np.dot(R, Q)
        # Add convergence information and extract eigenvalue if close enough.
        convergence_measure.append(np.abs(T[n, n - 1]))
        if convergence_measure[-1] < QR_ALGORITHM_TOLERANCE:
            eigenvalues[n] = T[n, n]
            T = T[:n, :n]
            n -= 1

    eigenvalues[0] = T
    return eigenvalues, convergence_measure


def QR_algorithm_shift(A):
    """The QR algorithm with largest value shift for finding eigenvalues.

    Using Householder reflection QR decomposition for finding RQ.

    .. code-block:: matlab

        function [myeig tmvec] = QRshift(A)

        %------ Phase 1 ------
        T = mytridiag(A);

        %------ Phase 2 ------
        m = length(A);
        myeig = zeros(m,1);
        % Size of tmvec is not known in advance.
        % Estimating an initial size.
        tmvec = zeros(m*8,1);
        mu = 0;
        counter = 0;

        while (m > 1)
            counter = counter + 1;
            muMat = diag(mu * ones(m,1));
            [Q,R] = myQR(T - muMat);
            T = R*Q + muMat;
            tmvec(counter) = abs(T(m,m-1));
            mu = T(m,m)
            if (tmvec(counter) < 1e-8)
                myeig(m) = T(m,m);
                m = m - 1;
                T = T(1:m,1:m);
            end
        end

        myeig(1) = T;
        tmvec = tmvec(1:counter);

        end

    :param A:
    :type A:
    :return: The eigenvalues.
    :rtype: list

    """

    # First, tridiagonalize the input matrix to a Hessenberg matrix.
    T = ma.convert_to_Hessenberg_double(A.copy())

    m, n = T.shape
    if n != m:
        raise np.linalg.LinAlgError("Array must be square.")

    convergence_measure = []
    eigenvalues = np.zeros((n, ), dtype='float')
    n -= 1

    while n > 0:
        # Obtain the shift from the lower right corner of the matrix.
        mu_matrix = np.eye(T.shape[0]) * T[n, n]
        # Perform QR decomposition on the shifted matrix.
        Q, R = ma.QR_factorisation_Householder_double(T - mu_matrix)
        # Multiply R and Q and shift the matrix back.
        T = np.dot(R, Q) + mu_matrix
        # Add convergence information and extract eigenvalue if close enough.
        convergence_measure.append(np.abs(T[n, n - 1]))
        if convergence_measure[-1] < QR_ALGORITHM_TOLERANCE:
            eigenvalues[n] = T[n, n]
            T = T[:n, :n]
            n -= 1

    eigenvalues[0] = T
    return eigenvalues, convergence_measure


def QR_algorithm_shift_Givens_double(A):
    """The QR algorithm with largest value shift for finding eigenvalues.

    Using Givens rotations for finding RQ.

    :param A: The square matrix to find eigenvalues of.
    :type A: :py:class:`numpy.ndarray`
    :return: The eigenvalues.
    :rtype: list

    """

    # First, tridiagonalize the input matrix to a Hessenberg matrix.
    T = ma.convert_to_Hessenberg_Givens_double(A.copy())

    m, n = T.shape
    if n != m:
        raise np.linalg.LinAlgError("Array must be square.")

    convergence_measure = []
    eigenvalues = np.zeros((m, ), dtype='float')
    m -= 1

    while m > 0:
        # Obtain the shift from the lower right corner of the matrix.
        mu_matrix = np.eye(T.shape[0]) * T[m, m]
        # Perform Givens QR step (which returns RQ) on the shifted
        # matrix and then shift it back.
        T = Givens_QR_step_double(T - mu_matrix) + mu_matrix
        # Add convergence information and extract eigenvalue if close enough.
        convergence_measure.append(np.abs(T[m, m - 1]))
        if convergence_measure[-1] < QR_ALGORITHM_TOLERANCE:
            eigenvalues[m] = T[m, m]
            T = T[:m, :m]
            m -= 1

    eigenvalues[0] = T
    return eigenvalues, convergence_measure


def QR_algorithm_shift_Givens_int(A):
    """The QR algorithm with largest value shift for finding eigenvalues, in integer precision.

    Using Givens rotations for finding RQ.

    :param A: The square matrix to find eigenvalues of.
    :type A: :py:class:`numpy.ndarray`
    :return: The eigenvalues.
    :rtype: list

    """

    A, scale = fp.scale_64bit_matrix(A.copy())

    # First, tridiagonalize the input matrix to a Hessenberg matrix.
    T = ma.convert_to_Hessenberg_Givens_int(A)

    m, n = T.shape
    if n != m:
        raise np.linalg.LinAlgError("Array must be square.")

    convergence_measure = []
    eigenvalues = np.zeros((m, ), dtype='int64')
    m -= 1

    while m > 0:
        # Obtain the shift from the lower right corner of the matrix.
        mu_matrix = np.eye(T.shape[0], dtype='int64') * T[m, m]
        # Perform Givens QR step (which returns RQ) on the shifted
        # matrix and then shift it back.
        T = Givens_QR_step_int(T - mu_matrix) + mu_matrix
        # Add convergence information and extract eigenvalue if close enough.
        convergence_measure.append(np.abs(T[m, m - 1]))
        if convergence_measure[-1] <= 1:
            eigenvalues[m] = T[m, m]
            T = T[:m, :m]
            m -= 1

    eigenvalues[0] = T
    return eigenvalues << scale, convergence_measure


def QR_algorithm_Wilkinson_shift(A):
    """The QR algorithm with Wilkinson shift for finding eigenvalues.

    .. code-block:: matlab

        function [myeig tmvec] = QRwilkinson(A)

        %------ Phase 1 ------
        T = mytridiag(A);

        %------ Phase 2 ------
        m = length(A);
        myeig = zeros(m,1);
        % Size of tmvec is not known in advance.
        % Estimating an initial size.
        tmvec = zeros(m*8,1);
        mu = 0;
        counter = 0;

        while (m > 1)
            counter = counter + 1;
            muMat = diag(mu * ones(m,1));
            [Q,R] = myQR(T - muMat);
            T = R*Q + muMat;
            tmvec(counter) = abs(T(m,m-1));
            delta = (T(m-1,m-1) - T(m,m)) / 2;
            mu = T(m,m) - sign(delta) * T(m,m-1) / ...
                    (abs(delta) + norm([delta T(m,m-1)]));
            if (tmvec(counter) < 1e-8)
                myeig(m) = T(m,m);
                m = m - 1;
                T = T(1:m,1:m);
            end
        end

        myeig(1) = T;
        tmvec = tmvec(1:counter);

        end

    :param A: The matrix to find eigenvalues in.
    :type A: :py:class:`numpy.ndarray`
    :return: The eigenvalues.
    :rtype: list

    """

    # First, tridiagonalize the input matrix to a Hessenberg matrix.
    T = ma.convert_to_Hessenberg_double(A.copy())

    m, n = T.shape
    convergence_measure = []
    eigenvalues = np.zeros((m, ), dtype='float')
    m -= 1
    mu = 0

    while m > 0:
        mu_matrix = np.eye(T.shape[0]) * mu
        Q, R = ma.QR_factorisation_Householder_double(T - mu_matrix)
        T = np.dot(R, Q) + mu_matrix
        convergence_measure.append(np.abs(T[m, m - 1]))
        delta = (T[m - 1, m - 1] - T[m, m]) / 2
        mu = T[m, m] - np.sign(delta) * T[m, m - 1] / \
             (np.abs(delta) + np.linalg.norm([delta, T[m, m - 1]]))

        if convergence_measure[-1] < QR_ALGORITHM_TOLERANCE:
            eigenvalues[m] = T[m, m]
            T = T[:m, :m]
            m -= 1

    eigenvalues[0] = T
    return eigenvalues, convergence_measure


def Givens_QR_step_double(A):
    """Performs a Givens rotation QR algorithm step in double precision.

    :param A: Square matrix in Hessenberg form.
    :type A: :py:class:`numpy.ndarray`
    :return: Square matrix in Hessenberg form.
    :rtype: :py:class:`numpy.ndarray`

    """
    _, n = A.shape

    G_storage = []
    for k in xrange(n - 1):
        c, s = ma.Givens_rotation_double(A[k, k], A[k + 1, k])
        G_storage.append((c, s))
        for j in xrange(n):
            tau_1 = A[k, j]
            tau_2 = A[k + 1, j]
            A[k, j] = ((tau_1 * c) - (tau_2 * s))
            A[k + 1, j] = ((tau_1 * s) + (tau_2 * c))
    for k in xrange(n - 1):
        c, s = G_storage.pop(0)
        for j in xrange(n):
            tau_1 = A[j, k]
            tau_2 = A[j, k + 1]
            A[j, k] = ((tau_1 * c) - (tau_2 * s))
            A[j, k + 1] = ((tau_1 * s) + (tau_2 * c))

    return np.triu(A, -1)


def Givens_QR_step_int(A):
    """Performs a Givens rotation QR algorithm step in integer precision.

    :param A: Square matrix in Hessenberg form.
    :type A: :py:class:`numpy.ndarray`
    :return: Square matrix in Hessenberg form.
    :rtype: :py:class:`numpy.ndarray`

    """
    _, n = A.shape

    G_storage = []
    for k in xrange(n - 1):
        c_n, s_n, denominator = ma.Givens_rotation_int(A[k, k], A[k + 1, k])
        G_storage.append((c_n, s_n, denominator))
        for j in xrange(n):
            tau_1 = A[k, j]
            tau_2 = A[k + 1, j]
            A[k, j] = ((tau_1 * c_n) - (tau_2 * s_n)) // denominator
            A[k + 1, j] = ((tau_1 * s_n) + (tau_2 * c_n)) // denominator
    for k in xrange(n - 1):
        c_n, s_n, denominator = G_storage.pop(0)
        for j in xrange(n):
            tau_1 = A[j, k]
            tau_2 = A[j, k + 1]
            A[j, k] = ((tau_1 * c_n) - (tau_2 * s_n)) // denominator
            A[j, k + 1] = ((tau_1 * s_n) + (tau_2 * c_n)) // denominator
    return A


