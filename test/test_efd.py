#!/usr/bin/env python
from spatial_efd import spatial_efd
import matplotlib
import os
import os.path as path
import numpy as np
import numpy.testing as ntest
import shapefile as shp
import pytest


@pytest.fixture
def open_square():
    return [0, 10, 10, 0], [0, 0, 10, 10]

@pytest.fixture
def closed_square():
    return [0, 10, 10, 0, 0], [0, 0, 10, 10, 0]

class TestEFD():
    def test_area(self, open_square):
        a = spatial_efd.ContourArea(*open_square)
        assert a == 100

    def test_centroid(self, open_square):
        c = spatial_efd.ContourCentroid(*open_square)
        assert c == (5, 5)

    def test_centroid_closed(self, closed_square):
        c = spatial_efd.ContourCentroid(*closed_square)
        assert c == (5, 5)

    def test_area_closed(self, closed_square):
        a = spatial_efd.ContourArea(*closed_square)
        assert a == 100

    def test_close_contour_closed_input(self):
        X, Y = spatial_efd.CloseContour([0, 10, 10, 0, 0], [0, 0, 10, 10, 0])
        assert X[0] == X[-1]
        assert Y[0] == Y[-1]

    def test_close_contour_open_input(self):
        X, Y = spatial_efd.CloseContour([0, 10, 10, 0], [0, 0, 10, 10])
        assert X[0] == X[-1]
        assert Y[0] == Y[-1]

    def test_nyquist(self):
        n = spatial_efd.Nyquist([0, 10, 10, 0, 0])
        assert n == 2

    def test_plot_init(self):
        a = spatial_efd.InitPlot()
        assert isinstance(a, matplotlib.axes.Axes)

    def test_rotate_contour(self):
        x, y = spatial_efd.RotateContour([0, 10, 10, 0], [0, 0, 10, 10], 30.,
                                         (5, 5))
        ntest.assert_almost_equal(x, [3.1698729810778059, 11.830127018922193,
                                   6.8301270189221945, -1.8301270189221928])
        ntest.assert_almost_equal(y, [-1.8301270189221928, 3.1698729810778059,
                                   11.830127018922193, 6.8301270189221945])

    def test_rotate_point(self):
        rx, ry = spatial_efd.rotatePoint((3., 2.), (1., 1.), 73.)
        assert pytest.approx(rx) == 0.628438653482
        assert pytest.approx(ry) == 3.20498121665

    def test_norm_contour(self):
        x, y, c = spatial_efd.NormContour([0., 10., 10., 0.],
                                          [0., 0., 10., 10.], (5., 5.))
        assert pytest.approx(c) == (0.5, 0.5)
        assert pytest.approx(x) == [0.0, 1.0, 1.0, 0.0]
        assert pytest.approx(y) == [0.0, 0.0, 1.0, 1.0]

    def test_get_bbox_dimensions(self):
        xw, yw, xmin, ymin = spatial_efd.getBBoxDimensions([0, 10, 10, 0],
                                                           [0, 0, 10, 10])
        assert xw == 10
        assert yw == 10
        assert xmin == 0
        assert ymin == 0

    def test_load_geometry(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        assert isinstance(s[0], shp._ShapeRecord)

    def test_process_geometry(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, c = spatial_efd.ProcessGeometry(s[1])
        ntest.assert_almost_equal(c, (280621.2724338955, 3882371.5613158443))
        ntest.assert_almost_equal(x[:10], [280587.0, 280598.0, 280598.0, 280599.0,
                                      280599.0, 280600.0, 280600.0, 280601.0,
                                      280601.0, 280602.0])
        ntest.assert_almost_equal(y[:10], [3882424.0, 3882424.0, 3882423.0,
                                      3882423.0, 3882422.0, 3882422.0,
                                      3882421.0, 3882421.0, 3882420.0,
                                      3882420.0])

    def test_process_geometry_norm(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, c = spatial_efd.ProcessGeometryNorm(s[1])
        ntest.assert_almost_equal(c, (0.4729141652616648, 0.22570629971140485))
        ntest.assert_almost_equal(x[:10], [0.29533678756476683, 0.35233160621761656,
                                      0.35233160621761656, 0.35751295336787564,
                                      0.35751295336787564, 0.3626943005181347,
                                      0.3626943005181347, 0.36787564766839376,
                                      0.36787564766839376, 0.37305699481865284])
        ntest.assert_almost_equal(y[:10], [0.49740932642487046, 0.49740932642487046,
                                      0.49222797927461137, 0.49222797927461137,
                                      0.48704663212435234, 0.48704663212435234,
                                      0.48186528497409326, 0.48186528497409326,
                                      0.47668393782383417, 0.47668393782383417])

    def test_calculate_efd(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometryNorm(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        ntest.assert_almost_equal(coeffs[6],
                                  [-0.00134937648, -0.000604478718,
                                   0.0003257416778, 0.001951924972])

    def test_inverse_transform(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometryNorm(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs)

        ntest.assert_almost_equal(a[:5], [-0.32100770398698036,
                                          -0.3201856580843658,
                                          -0.318692800914299,
                                          -0.3165482915823218,
                                          -0.31377942040690143])
        ntest.assert_almost_equal(b[:5], [0.4174800975360954,
                                          0.4180142260117704,
                                          0.417775492831905,
                                          0.416754380393524,
                                          0.41495338852007846])

    def test_inverse_transform_locus(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometryNorm(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs, locus=(0.5, 0.9))

        ntest.assert_almost_equal(a[:5], [0.1789922960130197,
                                          0.17981434191563422,
                                          0.181307199085701,
                                          0.18345170841767827,
                                          0.1862205795930986]
                                          )
        ntest.assert_almost_equal(b[:5], [1.3174800975360959,
                                          1.3180142260117702,
                                          1.317775492831905,
                                          1.3167543803935242,
                                          1.3149533885200781]
                                          )

    def test_average_coefficients(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)

        coeffsList = []

        for i in xrange(3):
            x, y, _ = spatial_efd.ProcessGeometryNorm(s[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList)

        ntest.assert_almost_equal(avg[6],
                                  [0.00049541617818, 0.00515338138093,
                                  -0.0005087032263, 9.7046992097e-05])

    def test_average_sd(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)

        coeffsList = []

        for i in xrange(3):
            x, y, _ = spatial_efd.ProcessGeometryNorm(s[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList)
        sd = spatial_efd.AverageSD(coeffsList, avg)

        ntest.assert_almost_equal(sd[3],
                                  [0.000381631249123, 0.00018247277186,
                                   4.6821200993e-05, 9.3013816155e-05])

    def test_fourier_power(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometryNorm(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 500)
        n = spatial_efd.FourierPower(coeffs, x)
        assert n == 19

    def test_normalize_efd(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometryNorm(s[0])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        coeffs, rotation = spatial_efd.normalize_efd(coeffs)

        ntest.assert_almost_equal(coeffs[9].tolist(),
                                  [-0.004300377673482293,
                                   0.00884561305918755,
                                  -0.013450240117972431,
                                  -0.0029657314108907686])
        assert pytest.approx(rotation) == 14.5510829786

    def test_calculate_dc_coefficients(self):
        filepath = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(filepath, 'example_data.shp')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometryNorm(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        dc = spatial_efd.calculate_dc_coefficients(x, y)
        assert pytest.approx(dc) == (0.34071444143386936, 0.56752000996605101)

    def test_plotting_savefig(self):
        matplotlib.pyplot.clf()
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        figpath = path.join(path_, 'Test')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometryNorm(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs)

        ax = spatial_efd.InitPlot()
        spatial_efd.PlotEllipse(ax, a, b, color='k', width=1.)
        spatial_efd.SavePlot(ax, 5, figpath, 'png')
        assert path.isfile('{0}_5.png'.format(figpath))
        os.remove('{0}_5.png'.format(figpath))

    def test_plotting(self):
        matplotlib.pyplot.clf()
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        figpath = path.join(path_, 'Test')

        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometryNorm(s[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs)

        ax = spatial_efd.InitPlot()
        spatial_efd.PlotEllipse(ax, a, b, color='k', width=1.)
        spatial_efd.SavePlot(ax, 5, figpath, 'png')
        assert path.isfile('{0}_5.png'.format(figpath))
        os.remove('{0}_5.png'.format(figpath))

    def test_plot_comparison(self):
        matplotlib.pyplot.clf()
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        figpath = path.join(path_, 'Test')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometry(s[0])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        ax = spatial_efd.InitPlot()
        spatial_efd.plotComparison(ax, coeffs, 10, x, y)
        spatial_efd.SavePlot(ax, 10, figpath, 'png')
        assert path.isfile('{0}_10.png'.format(figpath))
        os.remove('{0}_10.png'.format(figpath))

    def test_plot_comparison_norm(self):
        matplotlib.pyplot.clf()
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        figpath = path.join(path_, 'Test')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometry(s[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        coeffs, rotation = spatial_efd.normalize_efd(coeffs,
                                                     size_invariant=False)
        ax = spatial_efd.InitPlot()
        spatial_efd.plotComparison(ax, coeffs, 7, x, y, rotation=rotation)
        spatial_efd.SavePlot(ax, 7, figpath, 'png')
        assert path.isfile('{0}_7.png'.format(figpath))
        os.remove('{0}_7.png'.format(figpath))

    def test_plot_comparison_norm_size_invariant(self):
        matplotlib.pyplot.clf()
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        figpath = path.join(path_, 'Test')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometry(s[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        coeffs, rotation = spatial_efd.normalize_efd(coeffs,
                                                     size_invariant=True)
        ax = spatial_efd.InitPlot()
        spatial_efd.plotComparison(ax, coeffs, 7, x, y, rotation=rotation)
        spatial_efd.SavePlot(ax, 8, figpath, 'png')
        assert path.isfile('{0}_8.png'.format(figpath))
        os.remove('{0}_8.png'.format(figpath))

    def test_write_geometry(self):
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        figpath = path.join(path_, 'Test')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometry(s[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)
        assert isinstance(shape, shp.Writer)

    def test_generate_shapefile(self):
        shape = spatial_efd.generateShapefile()
        assert isinstance(shape, shp.Writer)

    def test_save_shapefile(self):
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        outpath = path.join(path_, 'OutShape')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometry(s[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)
        spatial_efd.saveShapefile(outpath, shape, prj=None)
        assert path.isfile('{0}.shp'.format(outpath))
        os.remove('{0}.shp'.format(outpath))
        os.remove('{0}.dbf'.format(outpath))
        os.remove('{0}.shx'.format(outpath))

    def test_save_shapefile_prj(self):
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        prjpath = path.join(path_, 'example_data.prj')
        outpath = path.join(path_, 'OutShape')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometry(s[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)
        spatial_efd.saveShapefile(outpath, shape, prj=prjpath)
        assert path.isfile('{0}.shp'.format(outpath))
        assert path.isfile(prjpath)
        os.remove('{0}.shp'.format(outpath))
        os.remove('{0}.dbf'.format(outpath))
        os.remove('{0}.shx'.format(outpath))
        os.remove('{0}.prj'.format(outpath))

    def test_save_shapefile_prj_missing(self):
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        prjpath = path.join(path_, 'example_data_missing.prj')
        outpath = path.join(path_, 'OutShape')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometry(s[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)
        spatial_efd.saveShapefile(outpath, shape, prj=prjpath)
        assert path.isfile('{0}.shp'.format(outpath))
        os.remove('{0}.shp'.format(outpath))
        os.remove('{0}.dbf'.format(outpath))
        os.remove('{0}.shx'.format(outpath))

    def test_save_shapefile_prj_wrong(self):
        path_ = path.realpath(path.join(os.getcwd(), path.dirname(__file__)))
        filepath = path.join(path_, 'example_data.shp')
        prjpath = filepath
        outpath = path.join(path_, 'OutShape')
        s = spatial_efd.LoadGeometries(filepath)
        x, y, _ = spatial_efd.ProcessGeometry(s[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)
        spatial_efd.saveShapefile(outpath, shape, prj=prjpath)
        assert path.isfile('{0}.shp'.format(outpath))
        os.remove('{0}.shp'.format(outpath))
        os.remove('{0}.dbf'.format(outpath))
        os.remove('{0}.shx'.format(outpath))
