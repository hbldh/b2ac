#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`distance_functions`
==================

.. module:: distance_functions
    :platform: Unix, Windows
    :synopsis: 

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-09-28, 23:25

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.geometry.point import B2ACPoint
from b2ac.geometry.ellipse import B2ACEllipse
from b2ac.geometry.polygon import B2ACPolygon


def distance(this, other):
    if isinstance(this, B2ACPoint):
        if isinstance(other, B2ACPoint):
            # Distance defined as distance between the Point objects.
            return np.linalg.norm(this.point - other.point)
        elif isinstance(other, B2ACEllipse):
            # Distance defined as distance from this point to center point of ellipse.
            return np.linalg.norm(this.point - other.center_point)
        elif isinstance(other, B2ACPolygon):
            # Distance defined as distance from this point to center point of polygon.
            return np.linalg.norm(this.point - np.mean(other.polygon_points, 0))
            # TODO: Better to evaluate distance from any polygon node point?
        else:
            raise ValueError("Cannot compare B2ACPoint to {0}.".format(type(other)))
    elif isinstance(this, B2ACEllipse):
        if isinstance(other, B2ACPoint):
            # Distance defined as distance from this point to center point of ellipse.
            return np.linalg.norm(self.center_point - other.point)
        elif isinstance(other, B2ACEllipse):
            # Distance defined as distance between center points of ellipses.
            return np.linalg.norm(self.center_point - other.center_point)
        elif isinstance(other, B2ACPolygon):
            # Distance defined as distance from ellipse center point and center of polygon.
            return np.linalg.norm(self.center_point - np.mean(other.polygon_points, 0))
            # return np.min(np.sqrt(((other.polygon_points - self.center_point) ** 2).sum(1)))
        else:
            raise ValueError("Cannot compare B2ACEllipse to {0}.".format(type(other)))
    elif isinstance(this, B2ACPolygon):
        if isinstance(other, B2ACPoint):
            # Distance defined as distance from the other point to center point of polygon.
            return np.linalg.norm(other.point - np.mean(self.polygon_points, 0))

            # TODO: Redefine to minimal distance from center point or any polygon node?
            # p1 = np.min(np.sqrt(((self.polygon_points - other.point) ** 2).sum(1)))
            # p2 = np.linalg.norm(np.mean(self.polygon_points, 0) - other.point)
            # return np.min([p1, p2])

        elif isinstance(other, B2ACEllipse):
            # Distance defined as distance from the center point of the other ellipse
            # to center point of polygon.
            return np.linalg.norm(other.center_point - np.mean(self.polygon_points, 0))

            # TODO: Redefine to minimal distance from center point or any polygon node?
            # p1 = np.min(np.sqrt(((self.polygon_points - other.center_point) ** 2).sum(1)))
            # p2 = np.linalg.norm(np.mean(self.polygon_points, 0) - other.center_point)
            # return np.min([p1, p2])

        elif isinstance(other, B2ACPolygon):
            # Distance defined as distances between center points of polygons.
            return np.linalg.norm(np.mean(self.polygon_points, 0) - np.mean(other.polygon_points, 0))

            # TODO: Redefine to minimal distance from center points or any polygon nodes?
            # p1 = np.min([np.min(np.sqrt(((v - other.polygon_points) ** 2).sum(1))) for v in self.polygon_points])
            # p2 = np.min(np.sqrt(((self.polygon_points - np.mean(other.polygon_points, 0)) ** 2).sum(1)))
            # p3 = np.min(np.sqrt(((other.polygon_points - np.mean(self.polygon_points, 0)) ** 2).sum(1)))
            # return np.min([p1, p2, p3])
        else:
            raise ValueError("Cannot compare B2ACPolygon to {0}.".format(type(other)))
    else:
        raise ValueError("Cannot call this method with a {0}.".format(type(this)))
