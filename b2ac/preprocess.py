#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`preprocess`
==================

.. module:: preprocess
    :platform: Unix, Windows
    :synopsis:

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-10-13

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import


def remove_mean_values(points):
    """Calculate integer mean values of the 2D array and removing it from the data.

    :param points: The points defining the ellipse.
    :type points:
    :return: Points, mean X integer coordinate, mean Y integer coordinate
    :rtype: tuple

    """
    x_mean = int(points[:, 0].mean())
    y_mean = int(points[:, 1].mean())
    return points - (x_mean, y_mean), x_mean, y_mean