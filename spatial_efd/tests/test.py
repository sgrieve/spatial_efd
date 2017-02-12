import spatial_efd
import matplotlib
import os
import os.path as path
import numpy as np
import shapefile as shp
from unittest import TestCase


class TestEFD(TestCase):
    def test_area(self):
        a = spatial_efd.ContourArea([0, 10, 10, 0], [0, 0, 10, 10])
        self.assertEqual(a, 100)

    def test_centroid(self):
        c = spatial_efd.ContourCentroid([0, 10, 10, 0], [0, 0, 10, 10])
        self.assertTupleEqual(c, (5, 5))

    def test_nyquist(self):
        n = spatial_efd.Nyquist(np.array([[0, 0], [10, 0], [10, 10], [10, 0]]))
        self.assertEqual(n, 2)

    def test_plot_init(self):
        a = spatial_efd.InitPlot()
        self.assertTrue(isinstance(a, matplotlib.axes.Axes))

    def test_rotate_contour(self):
        x, y = spatial_efd.RotateContour([0, 10, 10, 0], [0, 0, 10, 10], 30.,
                                         (5, 5))
        self.assertAlmostEqual(x, [6.8301270189221928, -1.8301270189221928,
                                   3.1698729810778077, 11.830127018922193],
                               places=7)
        self.assertAlmostEqual(y, [11.830127018922193, 6.8301270189221928,
                                   -1.8301270189221928, 3.1698729810778077],
                               places=7)

    def test_norm_contour(self):
        x, y, c = spatial_efd.NormContour([0., 10., 10., 0.],
                                          [0., 0., 10., 10.], (5., 5.))
        self.assertTupleEqual(c, (0.5, 0.5))
        self.assertListEqual(x, [0.0, 1.0, 1.0, 0.0])
        self.assertListEqual(y, [0.0, 0.0, 1.0, 1.0])

    def test_get_bbox_dimensions(self):
        xw, yw, xmin, ymin = spatial_efd.getBBoxDimensions([0, 10, 10, 0],
                                                           [0, 0, 10, 10])
        self.assertEqual(xw, 10)
        self.assertEqual(yw, 10)
        self.assertEqual(xmin, 0)
        self.assertEqual(ymin, 0)

    def test_load_geometry(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)

        self.assertTrue(isinstance(s[0], shp._ShapeRecord))
