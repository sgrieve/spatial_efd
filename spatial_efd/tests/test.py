#!/usr/bin/env python
import spatial_efd
import os
import os.path as path
import numpy as np
import numpy.testing as ntest
import shapefile as shp
from unittest import TestCase
import matplotlib



class TestEFD(TestCase):
    def test_area(self):
        a = spatial_efd.ContourArea([0, 10, 10, 0], [0, 0, 10, 10])
        self.assertEqual(a, 100)

    def test_centroid(self):
        c = spatial_efd.ContourCentroid([0, 10, 10, 0], [0, 0, 10, 10])
        self.assertTupleEqual(c, (5, 5))

    def test_centroid_closed(self):
        c = spatial_efd.ContourCentroid([0, 10, 10, 0, 0], [0, 0, 10, 10, 0])
        self.assertTupleEqual(c, (5, 5))

    def test_area_closed(self):
        a = spatial_efd.ContourArea([0, 10, 10, 0, 0], [0, 0, 10, 10, 0])
        self.assertEqual(a, 100)

    def test_close_contour_closed_input(self):
        X, Y = spatial_efd.CloseContour([0, 10, 10, 0, 0], [0, 0, 10, 10, 0])
        self.assertEqual(X[0], X[-1])
        self.assertEqual(Y[0], Y[-1])

    def test_close_contour_open_input(self):
        X, Y = spatial_efd.CloseContour([0, 10, 10, 0], [0, 0, 10, 10])
        self.assertEqual(X[0], X[-1])
        self.assertEqual(Y[0], Y[-1])

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

    def test_rotate_point(self):
        rx, ry = spatial_efd.rotatePoint((3., 2.), (1., 1.), 73.)
        self.assertAlmostEqual(rx, 0.628438653482)
        self.assertAlmostEqual(ry, 3.20498121665)

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
        self.assertListEqual(x[:10], [0.29533678756476683, 0.35233160621761656,
                                      0.35233160621761656, 0.35751295336787564,
                                      0.35751295336787564, 0.3626943005181347,
                                      0.3626943005181347, 0.36787564766839376,
                                      0.36787564766839376, 0.37305699481865284])
        self.assertListEqual(y[:10], [0.49740932642487046, 0.49740932642487046,
                                      0.49222797927461137, 0.49222797927461137,
                                      0.48704663212435234, 0.48704663212435234,
                                      0.48186528497409326, 0.48186528497409326,
                                      0.47668393782383417, 0.47668393782383417])

    def test_calculate_efd(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _, _ = spatial_efd.ProcessGeometry(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        ntest.assert_almost_equal(coeffs[6].tolist(),
                                  [-0.00134937648, -0.000604478718,
                                   0.0003257416778, 0.001951924972])

    def test_inverse_transform(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _, _ = spatial_efd.ProcessGeometry(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs)

        self.assertListEqual(a[:5].tolist(), [-0.32100770398698036,
                                              -0.3201856580843658,
                                              -0.318692800914299,
                                              -0.3165482915823218,
                                              -0.31377942040690143])
        self.assertListEqual(b[:5].tolist(), [0.4174800975360954,
                                              0.4180142260117704,
                                              0.417775492831905,
                                              0.416754380393524,
                                              0.41495338852007846])

    def test_average_coefficients(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)

        coeffsList = []

        for i in xrange(3):
            x, y, _, _ = spatial_efd.ProcessGeometry(s[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList)

        ntest.assert_almost_equal(avg[6].tolist(),
                                  [0.00049541617818, 0.00515338138093,
                                  -0.0005087032263, 9.7046992097e-05])

    def test_average_sd(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)

        coeffsList = []

        for i in xrange(3):
            x, y, _, _ = spatial_efd.ProcessGeometry(s[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList)
        sd = spatial_efd.AverageSD(coeffsList, avg)

        ntest.assert_almost_equal(sd[3].tolist(),
                                  [0.000381631249123, 0.00018247277186,
                                   4.6821200993e-05, 9.3013816155e-05])

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
        self.assertTupleEqual(dc, (0.34071444143386936, 0.56752000996605101))

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
