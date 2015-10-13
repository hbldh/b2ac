#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`polygon`
==================

.. module:: polygon
    :platform: Unix, Windows
    :synopsis: 

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-09-24, 23:50

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import warnings
import numpy as np

from b2ac.geometry.shape import B2ACGeometricShape


class B2ACPolygon(B2ACGeometricShape):
    """A class for representing polygons."""

    def __init__(self, points):
        """Constructor for B2ACPolygon"""
        super(B2ACPolygon, self).__init__()

        self.polygon_points = np.array(points, 'float')
        if self.polygon_points.shape[1] != 2:
            raise ValueError("Polygon must be entered as a [n x 2] array, i.e. a 2D polygon.")

    def __str__(self):
        return "Polygon with {0} vertices:\n{1}".format(
            self.polygon_points.shape[0], self.polygon_points)

    def __repr__(self):
        return str(self)

    @property
    def mpl_patch_arguments(self):
        """Returns Matplotlib patch arguments. Can then be plotted by with following
        snippet:

        .. code-block:: python

            import matplotlib.pyplot as plt
            from matplotlib.patches import Polygon
            from b2ac.geometry.polygon import B2ACPolygon

            p = B2ACPolygon([(2.1, 1.2), (5.2, 2.4), (1.0, 0.2)])
            plt.gca().add_patch(Polygon(
                p.mpl_patch_arguments, closed=True, fill=False, edgecolor='r')
            plt.show()

        """
        return self.polygon_points,

    def get_center_point(self, use_centroid=True):
        """Returns a center of weight for the object.

        :param use_centroid: Uses a centroid finding method instead of pure mean of vertices.
        :type use_centroid: bool

        """
        if use_centroid:
            with warnings.catch_warnings(record=False) as w:
                # Cause all warnings to never be triggered.
                warnings.simplefilter("ignore")

                pnt_array = self.get_closed_polygon()

                A = self._area_help_function()
                D = (pnt_array[:-1, 0] * pnt_array[1:, 1] -
                     pnt_array[1:, 0] * pnt_array[:-1, 1])

                c_x = ((pnt_array[:-1, 0] + pnt_array[1:, 0]) * D).sum() / (6 * A)
                c_y = ((pnt_array[:-1, 1] + pnt_array[1:, 1]) * D).sum() / (6 * A)

                if np.isnan(c_x) or np.isinf(c_x) or np.isnan(c_y) or np.isinf(c_y):
                    # If centroid calculations fails (e.g. due to zero-valued area) then use the
                    # mean of the vertices as center point instead.
                    return np.mean(self.get_open_polygon(), 0)
                else:
                    return np.array([c_x, c_y])
        else:
            return np.mean(self.get_open_polygon(), 0)

    def get_area(self):
        """Returns the area covered by the B2ACPolygon.

        :return: The area of the polygon.
        :rtype: float

        """
        # Abs to handle counter-clockwise ordering.
        return np.abs(self._area_help_function())

    @property
    def is_clockwise_ordered(self):
        """Property for checking if the polygon points are ordered clockwise.

        :return: Clockwise ordering status.
        :rtype: bool

        """
        return self._area_help_function() > 0

    @property
    def is_counter_clockwise_ordered(self):
        """Property for checking if the polygon points are ordered counter-clockwise.

        :return: Counter-clockwise ordering status.
        :rtype: bool

        """
        return self._area_help_function() < 0

    def get_bounding_box(self):
        """Get a four point Polygon describing the bounding box of the current Polygon.

        :return: A new polygon, describing this one's bounding box.
        :rtype: :py:class:`b2ac.geometry.polygon.B2ACPolygon`

        """
        mins = self.polygon_points.min(axis=0)
        maxs = self.polygon_points.max(axis=0)
        bb_pnts = np.zeros((4, 2), dtype=self.polygon_points.dtype)
        bb_pnts[0, :] = mins
        bb_pnts[1, :] = [mins[0], maxs[1]]
        bb_pnts[2, :] = maxs
        bb_pnts[3, :] = [maxs[0], mins[1]]
        out = B2ACPolygon(bb_pnts)
        return out

    def _area_help_function(self):
        """Performing the actual area calculation.

        :return: The area of the polygon, negative if counter-clockwise orientation of polygon points.
        :rtype: float

        """
        # If the polygon is not closed, append first point to the end of polygon to ensure that.
        pnt_array = self.get_closed_polygon()

        return (pnt_array[:-1, 0] * pnt_array[1:, 1] -
                pnt_array[1:, 0] * pnt_array[:-1, 1]).sum() / 2

    def get_closed_polygon(self):
        """Appends the first point to the end of point array, in order to "close" the polygon."""
        if not self.is_closed:
            return np.concatenate([self.polygon_points, [self.polygon_points[0, :]]])
        else:
            return self.polygon_points

    def get_open_polygon(self):
        """Removes the last point if it is close enough to the first, in order to "open" the polygon."""
        if self.is_closed:
            return self.polygon_points[:-1, :]
        else:
            return self.polygon_points

    @property
    def is_closed(self):
        return np.linalg.norm(self.polygon_points[0, :] - self.polygon_points[-1, :]) < 1e-10
