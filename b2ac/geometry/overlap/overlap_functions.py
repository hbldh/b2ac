#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: overlap_functions
   :platform: Unix, Windows
   :synopsis: Overlapping functions for the different B2ACGeometricShape objects.

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2013-04-15, 11:37

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

from b2ac.geometry.point import B2ACPoint
from b2ac.geometry.ellipse import B2ACEllipse
from b2ac.geometry.polygon import B2ACPolygon

N_POLYPOINTS = 73


def overlap(this, other):
    if isinstance(this, B2ACPoint):
        if isinstance(other, B2ACPoint):
            # Point to point overlap is really non-defined. Using parent class
            # __eq__ to test this.
            return float(this == other)
        elif isinstance(other, B2ACEllipse):
            return overlap_point_ellipse(this, other)
        elif isinstance(other, B2ACPolygon):
            return overlap_point_polygon(this, other)
    elif isinstance(this, B2ACEllipse):
        if isinstance(other, B2ACPoint):
            return overlap_point_ellipse(other, this)
        elif isinstance(other, B2ACEllipse):
            return overlap_ellipse_ellipse(this, other)
        elif isinstance(other, B2ACPolygon):
            return overlap_ellipse_polygon(this, other)
    elif isinstance(this, B2ACPolygon):
        if isinstance(other, B2ACPoint):
            return overlap_point_polygon(other, this)
        elif isinstance(other, B2ACEllipse):
            return overlap_polygon_ellipse(this, other)
        elif isinstance(other, B2ACPolygon):
            return overlap_polygon_polygon(this, other)


def overlap_point_ellipse(point, ellipse):
    """Test by definition of ellipse.

    Implementation reference:
    `Point in ellipse <http://www.maa.org/joma/volume8/kalman/general.html>`_

    :param point: The point.
    :type point: :py:class:`b2ac.geometry.point.B2ACPoint`
    :param ellipse: The ellipse.
    :type ellipse: :py:class:`b2ac.geometry.ellipse.B2ACEllipse`
    :return: 1.0 if point is inside ellipse, 0.0 otherwise.
    :rtype: float

    """
    tr_point = point.point - ellipse.center_point
    x_coeff = ((np.cos(ellipse.rotation_angle) / ellipse.radii[0]) ** 2 +
               (np.sin(ellipse.rotation_angle) / ellipse.radii[1]) ** 2)
    xy_coeff = (2 * np.cos(ellipse.rotation_angle) * np.sin(ellipse.rotation_angle) *
                ((1 / ellipse.radii[0]) + (1 / ellipse.radii[1])))
    y_coeff = ((np.sin(ellipse.rotation_angle) / ellipse.radii[0]) ** 2 +
               (np.cos(ellipse.rotation_angle) / ellipse.radii[1]) ** 2)
    return float((x_coeff * (tr_point[0] ** 2) + xy_coeff * (tr_point[0] * tr_point[1]) +
                  y_coeff * (tr_point[1])) <= 1)


def overlap_point_polygon(point, polygon):
    """Performs simple ray casting to test if the point is inside the polygon.

    .. note::

        Can only handle convex polygons! Additionally, if point is on polygon boundary,
        this solution might treat it as outside the polygon...

    Method reference: `Ray casting <http://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm>`_

    Line intersection equation:

    .. math::

        \\left[ \\begin{array}{c} x_{i} \\\\ y_{i} \\end{array} \\right] + d \\cdot
        \\left[ \\begin{array}{c} x_{i+1}-x_{i} \\\\ y_{i+1}-y_{i} \\end{array} \\right] =
        \\left[ \\begin{array}{c} \hat{x} \\\\ \hat{y} \\end{array} \\right] + t \\cdot
        \\left[ \\begin{array}{c} 1 \\\\ 0 \\end{array} \\right].

    From this we can solve for :math:`d` and :math:`t` and estimate where
    intersection occurs:

    .. math::

        d = \\frac{\hat{y} - y_{i}}{y_{i+1} - y_{i}}

    .. math::

        t = x_{i} - \\hat{x} + \\frac{\hat{y} - y_{i}}{y_{i+1} - y_{i}}(x_{i+1} - x_{i})

    If :math:`d` is smaller than 0, then intersection does not happen at along the polygon
    edge, likewise if :math:`d` is greater than the distance between the two vertices.
    If it is in between these two, then intersection occurs if :math:`t` is positive,
    i.e. it happens in the positive x-axis direction. If edge is parallel, then no crossing
    occurs either.

    :param point: The point.
    :type point: :py:class:`b2ac.geometry.point.B2ACPoint`
    :param polygon: The polygon.
    :type polygon: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :return: 1.0 if point is inside polygon, 0.0 otherwise.
    :rtype: float

    """
    # Add first point again to the array of polygon points, to simplify loop.
    # Check whether the first and the last point are equal, i.e. the polygon is closed
    # properly. Otherwise, add first point to the end of the array to close the polygon.
    if np.linalg.norm(polygon.polygon_points[-1, :] - polygon.polygon_points[0, :]) > 1e-10:
        pnt_arr = np.concatenate((polygon.polygon_points, [polygon.polygon_points[0, :]]))
    else:
        pnt_arr = polygon.polygon_points

    # Start ray tracing.
    nbr_of_crossings = 0
    for i in xrange(pnt_arr.shape[0] - 1):
        # Estimate distance along edge direction that line from point along x-axis
        # intersects this polygon edge.
        denominator = pnt_arr[i + 1, 1] - pnt_arr[i, 1]
        if denominator == 0:
            # Parallel lines.
            continue
        d = (point.point[1] - pnt_arr[i, 1]) / denominator
        # The equality means that if intersection occurs at first polygon vertex, it counts
        # as a crossing. The omission of it at the "far end" of the edge means that it will
        # not count as a crossing at that end.
        if (d >= 0.0) and (d < 1.0):
            # Have to check if it is in the negative x-direction the crossing occurs;
            # if so we can neglect this.
            t = pnt_arr[i, 0] - point.point[0] + d * (pnt_arr[i + 1, 0] - pnt_arr[i, 0])
            # Greater than 0 means that this point counts as "being on the right hand
            # side" of the edge.
            if t > 0.0:
                nbr_of_crossings += 1

    # If crossings of edges only occurs once, then the point must lie inside the
    # convex polygon. If 0 or 2 times, it is outside.
    if nbr_of_crossings == 1:
        return 1.0
    else:
        return 0.0


def overlap_ellipse_ellipse(e1, e2):
    """Calculates the intersection of two ellipses, returning the ratio of how
    much of ``e1`` that is contained in ``e2``.

    Current implementation:
        Converts both ellipses to polygons and uses :py:meth:`overlap_polygon_polygon`.

    .. note::

        Overlap with itthis should return 1.0, but is currently suffering from some
        numerical instabilities that returns unstable results!

    Implementation proposals:

    1. http://www.geometrictools.com/Documentation/IntersectionOfEllipses.pdf
    2. http://www.geometrictools.com/Documentation/AreaIntersectingEllipses.pdf

    :param e1: The first ellipse.
    :type e1: :py:class:`b2ac.geometry.ellipse.B2ACEllipse`
    :param e1: The second ellipse.
    :type e2: :py:class:`b2ac.geometry.ellipse.B2ACEllipse`
    :return: Value in [0, 1] describing how much of ``e1``
     that is contained in ``e2``.
    :rtype: float

    """
    return overlap_polygon_polygon(B2ACPolygon(e1.polygonize(N_POLYPOINTS)),
                                   B2ACPolygon(e2.polygonize(N_POLYPOINTS)))


def overlap_ellipse_polygon(ellipse, polygon):
    """Get ratio of how much of ``ellipse`` that is contained in ``polygon``.

    Current implementation:
    Convert ellipse to polygon and use :py:meth:`overlap_polygon_polygon`

    :param ellipse: The ellipse.
    :type ellipse: :py:class:`b2ac.geometry.ellipse.B2ACEllipse`
    :param polygon: The polygon.
    :type polygon: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :return: Value in [0, 1] describing how much of ``ellipse``
     that is contained in ``polygon``.
    :rtype: float

    """

    return overlap_polygon_polygon(B2ACPolygon(ellipse.polygonize(N_POLYPOINTS)), polygon)


def overlap_polygon_ellipse(polygon, ellipse):
    """Get ratio of how much of ``polygon`` that is contained in ``ellipse``.

    Current implementation:
    Convert ellipse to polygon and use :py:meth:`overlap_polygon_polygon`

    :param polygon: The polygon.
    :type polygon: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :param ellipse: The ellipse.
    :type ellipse: :py:class:`b2ac.geometry.ellipse.B2ACEllipse`
    :return: Value in [0, 1] describing how much of ``polygon``
     that is contained in ``ellipse``.
    :rtype: float

    """
    return overlap_polygon_polygon(polygon, B2ACPolygon(ellipse.polygonize(N_POLYPOINTS)))


def overlap_polygon_polygon(p1, p2):
    """Get ratio of how much of ``p1`` that is contained in ``p2``.

    .. note::

        Can only handle convex polygons since it uses Sutherland-Hodgman algorithm!

    Method reference:
    `Sutherland-Hodgman <http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm>`_

    :param p1: The first polygon.
    :type p1: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :param p2: The second polygon
    :type p2: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :return: Value in [0, 1] describing how much of ``p1``
     that is contained in ``p2``.
    :rtype: float

    """
    # First, check if the polygons are equal. If they are, the intersection calculation
    # might encounter divide-by-zero errors and is furthermore unnecessary to run...
    if np.allclose(p1.polygon_points.shape, p2.polygon_points.shape) and \
            np.allclose(p1.polygon_points, p2.polygon_points):
        return 1.0

    # Return the ratio of size of the intersection and size of the subject polygon.
    intersection_polygon = polygon_intersection(p1, p2)
    if intersection_polygon is None:
        return 0.
    else:
        return intersection_polygon.get_area() / p1.get_area()


def polygon_intersection(p1, p2):
    """Returns the polygon representing the intersection of two other.

    .. note::

        Can only handle convex polygons since it uses Sutherland-Hodgman algorithm!

    Method reference:
    `Sutherland-Hodgman <http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm>`_

    :param p1: The first polygon.
    :type p1: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :param p2: The second polygon
    :type p2: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :return: The intersection polygon or None if no intersection exists.
    :rtype: :py:class:`b2ac.geometry.polygon.B2ACPolygon` or :py:class:`None`.

    """
    # First remove the "closing points" of the polygon if such a vertex is present.
    p1_pnts = p1.get_open_polygon()
    p2_pnts = p2.get_open_polygon()

    # Perform clipping of polygons and create a new polygon instance from that.
    # TODO: Write function for ordering points clockwise in standard Cartesian plane.
    clipped_polygon_points = sutherland_hodgman_polygon_clipping(p1_pnts, p2_pnts)
    if clipped_polygon_points is not None and len(clipped_polygon_points):
        return p1.__class__(clipped_polygon_points)
    else:
        return None


def sutherland_hodgman_polygon_clipping(subject_polygon, clip_polygon):
    """Sutherland-Hodgman polygon clipping.

    .. note

        This algorithm works in regular Cartesian plane, not in inverted y-axis image plane,
        so make sure that polygons sent in are ordered clockwise in regular Cartesian sense!

    Method reference:
    `Sutherland-Hodgman <http://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm>`_

    Reference code found at `Rosettacode
    <http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Python>`_

    :param subject_polygon: A [n x 2] array of points representing the non-closed polygon to reduce.
    :type subject_polygon: :py:class:`numpy.ndarray`
    :param clip_polygon: A [m x 2] array of points representing the non-closed polygon to clip with.
    :type clip_polygon: :py:class:`numpy.ndarray`
    :return: A [r x 2] array of points representing the intersection polygon or :py:class:`None`
     if no intersection is present.
    :rtype: :py:class:`numpy.ndarray` or :py:class:`None`

    """
    TOLERANCE = 1e-14

    def inside(p):
        # This ``inside`` function assumes y-axis pointing upwards. If one would
        # like to rewrite this function to work with clockwise ordered coordinates
        # in the image style, then reverse the comparison from ``>`` to ``<``.
        return (cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0])

    def compute_intersection():
        dc = [cp1[0] - cp2[0], cp1[1] - cp2[1]]
        dp = [s[0] - e[0], s[1] - e[1]]
        n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
        n2 = s[0] * e[1] - s[1] * e[0]
        denominator = (dc[0] * dp[1] - dc[1] * dp[0])
        if np.abs(denominator) < TOLERANCE:
            # Lines were parallel.
            return None
        n3 = 1.0 / denominator
        return [(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3]

    output_list = list(subject_polygon)
    cp1 = clip_polygon[-1]

    for clip_vertex in clip_polygon:
        cp2 = clip_vertex
        input_list = output_list
        if not input_list:
            return None
        output_list = []
        s = input_list[-1]

        for subject_vertex in input_list:
            e = subject_vertex
            if inside(e):
                if not inside(s):
                    intersection = compute_intersection()
                    if intersection is not None:
                        output_list.append(intersection)
                output_list.append(e)
            elif inside(s):
                intersection = compute_intersection()
                if intersection is not None:
                    output_list.append(intersection)
            s = e
        cp1 = cp2

    # TODO: Verify that points are clockwise sorted here.
    pnts_out = []
    while len(output_list):
        pnt = output_list.pop(0)
        if not any([np.all(np.equal(pnt, unique_pnt)) for unique_pnt in pnts_out]):
            pnts_out.append(pnt)

    return np.array(pnts_out)


def polygon_union(p1, p2, use_graham_scan=False):
    """Calculates the union of two polygons and returns this union polygon.

    Implementation proposals:

    1. `Polygon union <http://stackoverflow.com/questions/6844462/polygon-union-without-holes>`_
    2. Find all intersection points (Sutherland-Hodgman) to get a point cloud and
       then run the `Graham scan <http://en.wikipedia.org/wiki/Graham_scan>`_ algorithm
       on that set of points.

    :param p1: The first polygon.
    :type p1: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :param p2: The second polygon
    :type p2: :py:class:`b2ac.geometry.polygon.B2ACPolygon`
    :param use_graham_scan: Boolean that can be used for selecting the slower Graham Scan convex hull algorithm
     instead of Quickhull.
    :type use_graham_scan: bool
    :return: The union polygon or :py:class:`None` if the polygons are not connected.
    :rtype: :py:class:`b2ac.geometry.polygon.B2ACPolygon` or :py:class:`None`

    """
    p_intersection = polygon_intersection(p1, p2)
    if p_intersection is None:
        return None
    points = np.concatenate((p1.polygon_points,
                             p2.polygon_points,
                             p_intersection.polygon_points))
    if use_graham_scan:
        return p1.__class__(graham_scan(points))
    else:
        return p1.__class__(quickhull(points))


def graham_scan(points):
    """Calculates the convex hull of an arbitrary 2D point cloud.

    Method reference:
    `Graham scan <http://en.wikipedia.org/wiki/Graham_scan>`_

    Code adapted from:
    `Google Code <https://mycodeplayground.googlecode.com/files/graham_scan.py>`_

    :param points: A [n x 2] array of points from which to estimate the convex hull.
    :type points: :py:class:`numpy.ndarray`
    :return: A [m x 2] array defining the convex hull polygon of the points sent in.
    :rtype: :py:class:`numpy.ndarray`

    """
    def angle_cmp(pivot):
        """Receive a coordinate as the pivot and return a
        function for comparing angles formed by another
        two coordinates around the pivot.
        """
        def _angle_cmp(c1, c2):
            v1 = c1[0] - pivot[0], c1[1] - pivot[1]
            v2 = c2[0] - pivot[0], c2[1] - pivot[1]
            cp = np.cross(v1, v2)
            if cp < 0:
                return 1
            elif cp == 0:
                return 0
            else:
                return -1
        return _angle_cmp

    def turning(c1, c2, c3):
        """Determine which way does c1 -> c2 -> c3 turns."""
        v1 = c2[0] - c1[0], c2[1] - c1[1]
        v2 = c3[0] - c2[0], c3[1] - c2[1]
        cp = np.cross(v1, v2)
        if cp < 0:
            return 'RIGHT'
        elif cp == 0:
            return 'STRAIGHT'
        else:
            return 'LEFT'

    def point_cmp(p1, p2):
        """Compares 2D points with regard to y coordinate value first, then x."""
        cmp_val = cmp(p1[1], p2[1])
        if cmp_val == 0:
            return cmp(p1[0], p2[0])
        else:
            return cmp_val

    num = len(points)
    if num < 3:
        raise Exception('Too few coordinates sent in.')

    # sort the coords according to y
    points = sorted(points.tolist(), cmp=point_cmp)

    # select the leftmost coord as the pivot
    pivot = points[0]
    coords = points[1:]

    # for remaining coords, sort them by polar angle
    # in counterclockwise order around pivot
    coords.sort(angle_cmp(pivot))

    # push the first three coords in a stack
    stack = [pivot, coords[0], coords[1]]

    # for the rest of the coords, while the angle formed by
    # the coord of the next-to-top of the stack, coord of
    # top of stack and the next coord makes a nonleft turn,
    # pop the stack
    # also, push the next coord into the stack at each loop
    for i in range(2, num - 1):
        while len(stack) >= 2 and \
              turning(stack[-2], stack[-1], coords[i]) != 'LEFT':
            stack = stack[:-1]
        stack.append(coords[i])

    return np.array(stack)


def quickhull(sample):
    """Calculates the convex hull of an arbitrary 2D point cloud.

    This is a pure Python version of the Quick Hull algorithm.
    It's based on the version of ``literateprograms``, but fixes some
    old-style Numeric function calls.

    This version works with numpy version > 1.2.1

    References:

    * `Literateprograms <http://en.literateprograms.org/Quickhull_(Python,_arrays)>`_
    * `Wikipedia <http://en.wikipedia.org/wiki/QuickHull>`_

    Code adapted from:

    `<http://members.home.nl/wim.h.bakker/python/quickhull2d.py>`_

    :param sample: Points to which the convex hull is desired to be found.
    :type sample: :py:class:`numpy.ndarray`
    :return: The convex hull of the points.
    :rtype: :py:class:`numpy.ndarray`

    """

    def calculate_convex_hull(sample):
        link = lambda a, b: np.concatenate((a, b[1:]))
        edge = lambda a, b: np.concatenate(([a], [b]))

        def dome(sample, base):
            h, t = base
            dists = np.dot(sample-h, np.dot(((0, -1), (1, 0)), (t - h)))
            outer = np.repeat(sample, dists > 0, axis=0)

            if len(outer):
                pivot = sample[np.argmax(dists)]
                return link(dome(outer, edge(h, pivot)),
                            dome(outer, edge(pivot, t)))
            else:
                return base

        if len(sample) > 2:
            axis = sample[:, 0]
            base = np.take(sample, [np.argmin(axis), np.argmax(axis)], axis=0)
            return link(dome(sample, base),
                        dome(sample, base[::-1]))
        else:
            return sample

    # Perform a reversal of points here to get points ordered clockwise instead of
    # counter clockwise that the QuickHull above returns.
    return calculate_convex_hull(sample)[::-1, :]


