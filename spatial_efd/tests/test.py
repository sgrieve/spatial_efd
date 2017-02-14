#!/usr/bin/env python
import spatial_efd
import os
import os.path as path
import numpy as np
import numpy.testing as ntest
import shapefile as shp
from unittest import TestCase
import matplotlib
matplotlib.use('Agg')


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

    def test_process_geometry(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, cont, c = spatial_efd.ProcessGeometry(s[1])
        self.assertTupleEqual(c, (0.4729141652616648, 0.22570629971140485))
        self.assertAlmostEqual(cont[25][0], 0.41450777)
        self.assertAlmostEqual(cont[25][1], 0.43523316)
        self.assertListEqual(x[:10], [280587.0, 280598.0, 280598.0, 280599.0,
                                      280599.0, 280600.0, 280600.0, 280601.0,
                                      280601.0, 280602.0])
        self.assertListEqual(y[:10], [3882424.0, 3882424.0, 3882423.0,
                                      3882423.0, 3882422.0, 3882422.0,
                                      3882421.0, 3882421.0, 3882420.0,
                                      3882420.0])

    def test_calculate_efd(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _, _ = spatial_efd.ProcessGeometry(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        ntest.assert_almost_equal(coeffs[6].tolist(),
                                  [-0.18756333, -0.08402254, 0.04527809,
                                   0.27131757])

    def test_inverse_transform(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _, _ = spatial_efd.ProcessGeometry(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs)

        self.assertListEqual(a[:5].tolist(), [-44.62007085419053,
                                              -44.50580647372712,
                                              -44.298299327087825,
                                              -44.000212529943006,
                                              -43.61533943655958])
        self.assertListEqual(b[:5].tolist(), [58.029733557517176,
                                              58.103977415636,
                                              58.07079350363471,
                                              57.92885887469977,
                                              57.67852100429082])

    def test_average_coefficients(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)

        coeffsList = []

        for i in xrange(3):
            x, y, _, _ = spatial_efd.ProcessGeometry(s[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList, 10)

        ntest.assert_almost_equal(avg[6].tolist(),
                                  [0.14906656481505043, 1.1743839538604195,
                                  -0.2589142773979526, -0.2106045300687506])

    def test_average_sd(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)

        coeffsList = []

        for i in xrange(3):
            x, y, _, _ = spatial_efd.ProcessGeometry(s[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList, 10)
        sd = spatial_efd.AverageSD(coeffsList, avg, 10)

        ntest.assert_almost_equal(sd[3].tolist(),
                                  [14.525962062245984, 13.07688969920094,
                                  2.8185995453824284, 4.387278049883412])

    def test_fourier_power(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, contour, _ = spatial_efd.ProcessGeometry(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 500)
        n = spatial_efd.FourierPower(coeffs, contour)
        self.assertEqual(n, 19)

    def test_normalize_efd(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _, _ = spatial_efd.ProcessGeometry(s[0])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        coeffs = spatial_efd.normalize_efd(coeffs)

        ntest.assert_almost_equal(coeffs[9].tolist(),
                                  [-0.004300377673482293,
                                   0.00884561305918755,
                                  -0.013450240117972431,
                                  -0.0029657314108907686])

    def test_calculate_dc_coefficients(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _, _ = spatial_efd.ProcessGeometry(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        dc = spatial_efd.calculate_dc_coefficients(x, y)
        self.assertTupleEqual(dc, (280011.35930735932, 3882261.8852813854))

    def test_plotting(self):
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        figpath = path.join(path_, 'Test')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, contour, _ = spatial_efd.ProcessGeometry(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs)

        ax = spatial_efd.InitPlot()
        spatial_efd.PlotEllipse(ax, a, b, color='k', width=1.)
        spatial_efd.PlotContour(ax, contour, color='b', width=1.)
        spatial_efd.SavePlot(ax, 5, figpath, 'png')
        os.remove('{0}_5.png'.format(figpath))
        self.assertTrue(isinstance(ax, matplotlib.axes.Axes))
