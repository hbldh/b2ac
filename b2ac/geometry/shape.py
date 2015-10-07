#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`shape`
==================

.. module:: shape
    :platform: Unix, Windows
    :synopsis: 

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-09-24, 23:51

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

#from b2ac.geometry.overlap.overlap_functions import overlap
#from b2ac.geometry.distance.distance_functions import distance


class B2ACGeometricShape(object):
    """Parent class for geometric shapes."""

    def __init__(self):
        """Constructor for parent class."""
        pass

    def __eq__(self, other):
        """Equality operator."""
        # TODO: Evaluate better equality measure...
        self.distance(other) < 0.01

    def __lt__(self, other):
        # TODO: Do these have any reasonable meaning?
        raise NotImplementedError()

    def __le__(self, other):
        # TODO: Do these have any reasonable meaning?
        raise NotImplementedError()

    def __gt__(self, other):
        # TODO: Do these have any reasonable meaning?
        raise NotImplementedError()

    def __ge__(self, other):
        # TODO: Do these have any reasonable meaning?
        raise NotImplementedError()

    @property
    def mpl_patch_arguments(self):
        """Returns Matplotlib Patch arguments."""
        raise NotImplementedError()

    def get_center_point(self):
        """Returns a center of weight for the object."""
        raise NotImplementedError()

    def get_area(self):
        """Returns the area covered of the B2ACGeometricShape object."""
        raise NotImplementedError()

    def distance(self, other):
        """Calculates distance between B2ACGeometricShapes."""
        #return distance(self, other)
        raise NotImplementedError()

    def overlap(self, other):
        """Returns ratio of overlap between this B2ACGeometricShape and another.

        :param other: The B2ACGeometricShape subclass to check overlap with.
        :type other: :py:class:`b2ac.geometry.shape.B2ACGeometricShape`
        :return: Value in [0, 1] describing how much of ``self`` that is enclosed in ``other``.
        :rtype: float

        """
        #return overlap(self, other)
        raise NotImplementedError()

    def intersection_union_ratio(self, other, N=73):
        """Returns ratio between the two shapes intersection and their union.

        :param other: The B2ACGeometricShape subclass to check intersection and union ratio with.
        :type other:
        :param N: The number of points to generate in case of ellipse to polygon conversion is needed.
         Default is 73 (one vertex every 5 degrees).
        :type N: int
        :return: Value in [0, 1].
        :rtype: float

        """
        raise NotImplementedError()
        # if isinstance(self, B2ACPoint) or isinstance(other, B2ACPoint):
        #     # Points cannot be used with this function.
        #     return 0.
        # if isinstance(self, B2ACEllipse):
        #     self_polygon = self.convert_to_AlgoTrackPolygon(N_POLYPOINTS)
        #     if isinstance(other, B2ACEllipse):
        #         other_polygon = other.convert_to_AlgoTrackPolygon(N_POLYPOINTS)
        #     elif isinstance(other, B2ACPolygon):
        #         other_polygon = other
        #     else:
        #         raise RuntimeError("Invalid comparison object.")
        # elif isinstance(self, B2ACPolygon):
        #     self_polygon = self
        #     if isinstance(other, B2ACEllipse):
        #         other_polygon = other.convert_to_AlgoTrackPolygon(N_POLYPOINTS)
        #     elif isinstance(other, B2ACPolygon):
        #         other_polygon = other
        #     else:
        #         raise RuntimeError("Invalid comparison object.")
        # else:
        #     # Cannot happen...
        #     raise RuntimeError("Invalid comparison object.")
        #
        # intersection_polygon = of.polygon_intersection(self_polygon, other_polygon)
        # if intersection_polygon is None:
        #     # No intersection. Return 0.
        #     return 0.
        #
        # union_polygon = of.polygon_union(self_polygon, other_polygon)
        # return intersection_polygon.get_area() / union_polygon.get_area()


