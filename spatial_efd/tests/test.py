import spatial_efd
import matplotlib
import os
import os.path as path
import numpy as np
import numpy.testing as ntest
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

    def test_process_geometry(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, cont, c = spatial_efd.ProcessGeometry(s[1])
        self.assertTupleEqual(c, (0.4729141652616648, 0.22570629971140485))
        self.assertAlmostEqual(cont[25][0], 0.43523316)
        self.assertAlmostEqual(cont[25][1], 0.41450777)
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
        pass

    def test_average_sd(self):
        pass

    def test_fourier_power(self):
        pass

    def test_normalize_efd(self):
        pass

    def test_calculate_dc_coefficients(self):
        pass
