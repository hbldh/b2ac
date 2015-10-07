#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: ellipse
    :platform: Unix, Windows
    :synopsis: Ellipse class.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2013-09-09, 16:32

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.geometry.shape import B2ACGeometricShape


class B2ACEllipse(B2ACGeometricShape):
    """An class for representing elliptical shapes."""

    def __init__(self, center, radii, rotation_angle):
        """Constructor for B2ACEllipse"""

        super(B2ACEllipse, self).__init__()
        self.center_point = np.array(center, 'float')
        self.radii = np.array(radii, 'float')
        self.rotation_angle = float(rotation_angle)

    def __str__(self):
        return "Ellipse: Center = ({0:.2f}, {1:.2f}), Radii = ({2:.2f}, {3:.2f}), Angle = {4:.4f}".format(
            self.center_point[0], self.center_point[1], self.radii[0], self.radii[1], self.rotation_angle)

    def __repr__(self):
        return str(self)

    @property
    def mpl_patch_arguments(self):
        """Returns Matplotlib patch arguments. Can then be plotted by with following
        snippet:

        .. code-block:: python

            import matplotlib.pyplot as plt
            from matplotlib.patches import Ellipse
            from b2ac.ellipse import B2ACEllipse

            e = B2ACEllipse((3.0, 5.0), (2.0, 6.0, 1.57))
            plt.gca().add_patch(Ellipse(
                e.mpl_patch_arguments, fill=False, edgecolor='g'))
            plt.show()

        """
        return self.center_point, self.radii[0] * 2, self.radii[1] * 2, np.rad2deg(self.rotation_angle)

    def get_center_point(self):
        """Returns a center of weight for the object."""
        return self.center_point

    def get_area(self):
        """Returns the area covered of the B2ACEllipse."""
        return np.pi * self.radii[0] * self.radii[1]

    def polygonize(self, n=73):
        """Gets a approximate polygon array representing the ellipse.

        Note that the last point is the same as the first point, creating a closed
        polygon.

        :param n: The number of points to generate. Default is 73 (one vertex every 5 degrees).
        :type n: int
        :return: An [n  x 2] numpy array, describing the boundary vertices of
                 the polygonized ellipse.
        :rtype: :py:class:`numpy.ndarray`

        """
        t = np.linspace(0, 2 * np.pi, num=n, endpoint=True)
        out = np.zeros((len(t), 2), dtype='float')
        out[:, 0] = (self.center_point[0] +
                     self.radii[0] * np.cos(t) * np.cos(self.rotation_angle) -
                     self.radii[1] * np.sin(t) * np.sin(self.rotation_angle))
        out[:, 1] = (self.center_point[1] +
                     self.radii[0] * np.cos(t) * np.sin(self.rotation_angle) +
                     self.radii[1] * np.sin(t) * np.cos(self.rotation_angle))
        return out


