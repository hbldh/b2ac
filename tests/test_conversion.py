#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`test_conversion`
==================

.. module:: test_conversion
    :platform: Unix, Windows
    :synopsis:

.. moduleauthor:: hbldh <henrik.blidh@nedomkull.com>

Created on 2015-10-02

"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy as np

import b2ac.fit.reference as efr
import b2ac.conversion as c2gconv
from b2ac.geometry.ellipse import B2ACEllipse
from b2ac.geometry.overlap.overlap_functions import overlap


class TestConversion(object):

    def setup(self):
        self.e = B2ACEllipse(center=(50.0, 75.0), radii=(50.0, 20.0), rotation_angle=0.707)
        self.points = np.array(self.e.polygonize(), 'int32')

    def test_double_conversions(self):
        output = efr.fit_improved_B2AC(self.points.copy())
        general_1 = c2gconv.conic_to_general_1(output.copy())
        general_2 = c2gconv.conic_to_general_2(output.copy())
        general_3 = c2gconv.conic_to_general_reference(output.copy())
        e_1 = B2ACEllipse(*general_1)
        e_2 = B2ACEllipse(*general_2)
        e_3 = B2ACEllipse(*general_3)

        # Test correct center point.
        assert np.linalg.norm(self.e.center_point - e_1.center_point) < 1
        # assert np.linalg.norm(self.e.center_point - e_2.center_point) < 1
        assert np.linalg.norm(self.e.center_point - e_3.center_point) < 1

        # Test correct radii.
        assert np.linalg.norm(max(self.e.radii) - max(e_1.radii)) < 0.25
        assert np.linalg.norm(min(self.e.radii) - min(e_1.radii)) < 0.25
        # assert np.linalg.norm(max(self.e.radii) - max(e_2.radii)) < 0.25
        # assert np.linalg.norm(min(self.e.radii) - min(e_2.radii)) < 0.25
        assert np.linalg.norm(max(self.e.radii) - max(e_3.radii)) < 0.25
        assert np.linalg.norm(min(self.e.radii) - min(e_3.radii)) < 0.25

        # Test overlap
        assert overlap(self.e, e_1) > 0.98
        assert overlap(e_1, self.e) > 0.98
        # assert overlap(self.e, e_2) > 0.98
        # assert overlap(e_2, self.e) > 0.98
        assert overlap(self.e, e_3) > 0.98
        assert overlap(e_3, self.e) > 0.98
