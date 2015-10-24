#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`test_int`
==================

.. module:: test_int
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
import b2ac.ext.ellipse_fitter as fitext

from b2ac.fit import fit_improved_B2AC_int
import b2ac.conversion as c2gconv
from b2ac.preprocess import remove_mean_values
from b2ac.geometry.ellipse import B2ACEllipse
from b2ac.geometry.overlap.overlap_functions import overlap


class TestIntegerImplementations(object):

    def setup(self):
        self.e = B2ACEllipse(center=(50.0, 75.0), radii=(50.0, 20.0), rotation_angle=0.707)
        self.points = np.array(self.e.polygonize(), 'int32')

    def test_fit_int_version_1(self):
        output = fit_improved_B2AC_int(self.points.copy())
        ellipse_data = c2gconv.conic_to_general_int(output, return_float=True, verbose=True)
        e_fitted = B2ACEllipse(*ellipse_data)
        assert np.linalg.norm(self.e.center_point - e_fitted.center_point) < 1
        assert np.linalg.norm(max(self.e.radii) - max(e_fitted.radii)) < 1
        assert np.linalg.norm(min(self.e.radii) - min(e_fitted.radii)) < 1
        assert overlap(self.e, e_fitted) > 0.95
        assert overlap(e_fitted, self.e) > 0.95

    def test_fit_int_version_2(self):
        points, x_mean, y_mean = remove_mean_values(self.points.copy())
        output = fit_improved_B2AC_int(points)
        ellipse_data = c2gconv.conic_to_general_int(output, return_float=True, verbose=True)
        e_fitted = B2ACEllipse(*ellipse_data)
        e_fitted.center_point += (x_mean, y_mean)
        print(e_fitted)
        assert np.linalg.norm(self.e.center_point - e_fitted.center_point) < 1.0
        assert np.linalg.norm(max(self.e.radii) - max(e_fitted.radii)) < 0.1
        assert np.linalg.norm(min(self.e.radii) - min(e_fitted.radii)) < 0.1
        assert overlap(self.e, e_fitted) > 0.98
        assert overlap(e_fitted, self.e) > 0.98

    def test_fit_ext_int_version(self):
        output = fitext.fit_ellipse_int(self.points.copy())
        e_fitted = B2ACEllipse(*output)
        assert np.linalg.norm(self.e.center_point - e_fitted.center_point) < 1
        assert np.linalg.norm(max(self.e.radii) - max(e_fitted.radii)) < 0.25
        assert np.linalg.norm(min(self.e.radii) - min(e_fitted.radii)) < 0.25
        assert overlap(self.e, e_fitted) > 0.98
        assert overlap(e_fitted, self.e) > 0.98

    def test_py_and_ext_similarity(self):
        output = fitext.fit_ellipse_int(self.points.copy())
        e_ext = B2ACEllipse(*output)
        points, x_mean, y_mean = remove_mean_values(self.points.copy())
        output = fit_improved_B2AC_int(points)
        ellipse_data = c2gconv.conic_to_general_int(output, return_float=True, verbose=True)
        e_py = B2ACEllipse(*ellipse_data)
        e_py.center_point += (x_mean, y_mean)
        assert overlap(e_ext, e_py) > 0.99
        assert overlap(e_py, e_ext) > 0.99
