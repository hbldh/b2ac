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
