#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
:mod:`test_ellipse_fit_reference`
==================

.. module:: test_ellipse_fit_reference
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

import b2ac.fit.reference as efr
from b2ac.fit import fit_improved_B2AC_double
import b2ac.conversion as c2gconv
from b2ac.geometry.ellipse import B2ACEllipse
from b2ac.geometry.overlap.overlap_functions import overlap


# Bengt Sändh FInn Zetterholm


class TestReferenceImplementation(object):

    def setup(self):
        self.e = B2ACEllipse(center=(50.0, 75.0), radii=(50.0, 20.0), rotation_angle=0.707)
        self.points = np.array(self.e.polygonize(), 'int32')

#    def test_fit_numpy_version(self):
#        # Fails
#        output = fit.fit_improved_B2AC_numpy(self.points.copy())
#        e_fitted = B2ACEllipse(*output)
#        assert np.linalg.norm(self.e.center_point - e_fitted.center_point) < 1
#        assert np.linalg.norm(max(self.e.radii) - max(e_fitted.radii)) < 0.25
#        assert np.linalg.norm(min(self.e.radii) - min(e_fitted.radii)) < 0.25
#        assert overlap(self.e, e_fitted) > 0.98
#        assert overlap(e_fitted, self.e) > 0.98

    def test_fit_double_version(self):
        output = fit_improved_B2AC_double(self.points.copy())
        output = c2gconv.conic_to_general_1(output)
        e_fitted = B2ACEllipse(*output)
        assert np.linalg.norm(self.e.center_point - e_fitted.center_point) < 1
        assert np.linalg.norm(max(self.e.radii) - max(e_fitted.radii)) < 0.25
        assert np.linalg.norm(min(self.e.radii) - min(e_fitted.radii)) < 0.25
        assert overlap(self.e, e_fitted) > 0.98
        assert overlap(e_fitted, self.e) > 0.98

    def test_fit_ext_double_version(self):
        output = fitext.fit_ellipse_double(self.points.copy())
        e_fitted = B2ACEllipse(*output)
        assert np.linalg.norm(self.e.center_point - e_fitted.center_point) < 1
        assert np.linalg.norm(max(self.e.radii) - max(e_fitted.radii)) < 0.25
        assert np.linalg.norm(min(self.e.radii) - min(e_fitted.radii)) < 0.25
        assert overlap(self.e, e_fitted) > 0.98
        assert overlap(e_fitted, self.e) > 0.98

    def test_fit_ext_float_version(self):
        output = fitext.fit_ellipse_float(self.points.copy())
        e_fitted = B2ACEllipse(*output)
        assert np.linalg.norm(self.e.center_point - e_fitted.center_point) < 1
        assert np.linalg.norm(max(self.e.radii) - max(e_fitted.radii)) < 0.25
        assert np.linalg.norm(min(self.e.radii) - min(e_fitted.radii)) < 0.25
        assert overlap(self.e, e_fitted) > 0.98
        assert overlap(e_fitted, self.e) > 0.98

    def test_fit_ref(self):
        points = self.points.copy()
        output = efr.fit_improved_B2AC(points)
        general_form = c2gconv.conic_to_general_1(output)
        e_fitted = B2ACEllipse(*general_form)
        assert np.linalg.norm(self.e.center_point - e_fitted.center_point) < 1
        assert np.linalg.norm(max(self.e.radii) - max(e_fitted.radii)) < 0.25
        assert np.linalg.norm(min(self.e.radii) - min(e_fitted.radii)) < 0.25
        assert overlap(self.e, e_fitted) > 0.98
        assert overlap(e_fitted, self.e) > 0.98
