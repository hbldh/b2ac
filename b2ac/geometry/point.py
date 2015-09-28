#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`point`
==================

.. module:: point
    :platform: Unix, Windows
    :synopsis: 

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-09-24, 23:51

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.geometry.shape import B2ACGeometricShape


class B2ACPoint(B2ACGeometricShape):
    """A class for representing points.

    :param point: An array of two values: x and y cooedinate for point.
    :type point: list, tuple or :py:class:`numpy.ndarray`

    Example Input:

    .. code-block:: python

        {u'algorithm': u'Is1FpgaGlintPupilDetector',
         u'name': u'glint',
         u'points': [{u'x': 453, u'y': 215}],
         u'type': u'point'}

    """

    def __init__(self, point):
        """Constructor for B2ACEllipse"""
        super(B2ACPoint, self).__init__()
        if len(point) != 2:
            raise ValueError("Only 2D-points supported.")
        self.point = np.array(point, 'float')

    def __str__(self):
        return "Point: ({0:.2f}, {1:.2f})".format(*self.point)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if isinstance(self, other):
            return np.testing.assert_allclose(self.point, other.point)
        return False

    @property
    def mpl_patch_arguments(self):
        """Returns Matplotlib Patch arguments. Can then be plotted by with following
        snippet:

        .. code-block:: python

            import matplotlib.pyplot as plt
            from matplotlib.patches import Point
            from b2ac.geometry.point import B2ACPoint

            p = B2ACPoint((3.0, 5.0))
            plt.gca().add_patch(Point(p.mpl_patch_arguments, color='g'))
            plt.show()

        """
        return self.point,

    def get_center_point(self):
        """Returns the B2ACPoints value."""
        return self.point

    def get_area(self):
        """Returns the area covered of the B2ACPoint."""
        return 0.
