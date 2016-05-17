#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`unstable`
==================

.. module:: unstable
    :platform: Unix, Windows
    :synopsis:

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-10-13

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.compat import *
from b2ac.fit.double import _calculate_scatter_matrix_double


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

    S = _calculate_scatter_matrix_double(points[:, 0], points[:, 1])

    eigenvalues, eigenvalues = scla.eig(S, constraint_matrix)
    ind = np.where(eigenvalues == (eigenvalues[eigenvalues > 0].min()))[0][0]
    return eigenvalues[:, ind]
