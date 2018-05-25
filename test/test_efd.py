#!/usr/bin/env python
from spatial_efd import spatial_efd
import matplotlib
import os
import numpy as np
import numpy.testing as ntest
import shapefile as shp
import pytest
import warnings
import re


@pytest.fixture
def open_square():
    return [0, 10, 10, 0], [0, 0, 10, 10]

@pytest.fixture
def closed_square():
    return [0, 10, 10, 0, 0], [0, 0, 10, 10, 0]

@pytest.fixture
def shp_paths():
    filepath = os.path.realpath(os.path.join(os.getcwd(),
                                             os.path.dirname(__file__)))
    shppath = os.path.join(filepath, 'example_data.shp')
    prjpath = os.path.join(filepath, 'example_data.prj')
    return shppath, prjpath

@pytest.fixture
def example_shp():
    return spatial_efd.LoadGeometries(shp_paths()[0])

@pytest.fixture
def warn_wrong_prj():
    return 'The file supplied is not a prj file. No .prj file will be written'

@pytest.fixture
def warn_missing_prj():
    return 'The .prj file supplied does not exist. No .prj file will be written'

def clean_warning(message):
    '''
    Helper function to format warning messages so they can be used for tests
    '''
    return re.sub(r'\(.*?\)\s', '', str(message))


class TestEFD():
    @pytest.mark.parametrize('shape, area',
                             [(open_square(), 100), (closed_square(), 100)])
    def test_area(self, shape, area):
        a = spatial_efd.ContourArea(*shape)
        assert a == area

    @pytest.mark.parametrize('shape, centre',
                             [(open_square(), (5, 5)),
                              (closed_square(), (5, 5))])
    def test_centroid(self, shape, centre):
        c = spatial_efd.ContourCentroid(*shape)
        assert c == (5, 5)

    @pytest.mark.parametrize('shape', [open_square(), closed_square()])
    def test_close_contour(self, shape):
        X, Y = spatial_efd.CloseContour(*shape)
        assert X[0] == X[-1]
        assert Y[0] == Y[-1]

    def test_nyquist(self, closed_square):
        n = spatial_efd.Nyquist(closed_square[0])
        assert n == 2

    def test_plot_init(self):
        a = spatial_efd.InitPlot()
        assert isinstance(a, matplotlib.axes.Axes)

    def test_rotate_contour(self, open_square):
        x, y = spatial_efd.RotateContour(*open_square, rotation=30.,
                                         centroid=(5, 5))
        ntest.assert_almost_equal(x, [3.1698729810778059, 11.830127018922193,
                                      6.8301270189221945, -1.8301270189221928])
        ntest.assert_almost_equal(y, [-1.8301270189221928, 3.1698729810778059,
                                      11.830127018922193, 6.8301270189221945])

    def test_rotate_point(self):
        rx, ry = spatial_efd.rotatePoint((3., 2.), (1., 1.), 73.)
        assert pytest.approx(rx) == 0.628438653482
        assert pytest.approx(ry) == 3.20498121665

    def test_norm_contour(self, open_square):
        x, y, c = spatial_efd.NormContour(*open_square, rawCentroid=(5., 5.))
        assert pytest.approx(c) == (0.5, 0.5)
        assert pytest.approx(x) == [0.0, 1.0, 1.0, 0.0]
        assert pytest.approx(y) == [0.0, 0.0, 1.0, 1.0]

    def test_get_bbox_dimensions(self, open_square):
        xw, yw, xmin, ymin = spatial_efd.getBBoxDimensions(*open_square)
        assert xw == 10
        assert yw == 10
        assert xmin == 0
        assert ymin == 0

    def test_load_geometry(self, example_shp):
        assert isinstance(example_shp[0], shp._ShapeRecord)

    def test_process_geometry(self, example_shp):
        x, y, c = spatial_efd.ProcessGeometry(example_shp[1])
        ntest.assert_almost_equal(c, (280621.2724338955, 3882371.5613158443))
        ntest.assert_almost_equal(x[:10], [280587.0, 280598.0, 280598.0, 280599.0,
                                      280599.0, 280600.0, 280600.0, 280601.0,
                                      280601.0, 280602.0])
        ntest.assert_almost_equal(y[:10], [3882424.0, 3882424.0, 3882423.0,
                                      3882423.0, 3882422.0, 3882422.0,
                                      3882421.0, 3882421.0, 3882420.0,
                                      3882420.0])

    def test_process_geometry_norm(self, example_shp):
        x, y, c = spatial_efd.ProcessGeometryNorm(example_shp[1])
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

    def test_calculate_efd(self, example_shp):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        ntest.assert_almost_equal(coeffs[6],
                                  [-0.00134937648, -0.000604478718,
                                   0.0003257416778, 0.001951924972])

    def test_inverse_transform(self, example_shp):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
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

    def test_inverse_transform_locus(self, example_shp):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
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

    def test_average_coefficients(self, example_shp):
        coeffsList = []

        for i in xrange(3):
            x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList)

        ntest.assert_almost_equal(avg[6],
                                  [0.00049541617818, 0.00515338138093,
                                  -0.0005087032263, 9.7046992097e-05])

    def test_average_sd(self, example_shp):
        coeffsList = []

        for i in xrange(3):
            x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList)
        sd = spatial_efd.AverageSD(coeffsList, avg)

        ntest.assert_almost_equal(sd[3],
                                  [0.000381631249123, 0.00018247277186,
                                   4.6821200993e-05, 9.3013816155e-05])

    def test_fourier_power(self, example_shp):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 500)
        n = spatial_efd.FourierPower(coeffs, x)
        assert n == 19

    def test_normalize_efd(self, example_shp):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[0])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        coeffs, rotation = spatial_efd.normalize_efd(coeffs)

        ntest.assert_almost_equal(coeffs[9],
                                  [-0.004300377673482293,
                                   0.00884561305918755,
                                  -0.013450240117972431,
                                  -0.0029657314108907686])
        assert pytest.approx(rotation) == 14.5510829786

    def test_calculate_dc_coefficients(self, example_shp):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        dc = spatial_efd.calculate_dc_coefficients(x, y)
        assert pytest.approx(dc) == (0.34071444143386936, 0.56752000996605101)

    def test_plotting_savefig(self, example_shp, tmpdir):
        matplotlib.pyplot.clf()
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs)
        ax = spatial_efd.InitPlot()
        spatial_efd.PlotEllipse(ax, a, b, color='k', width=1.)
        spatial_efd.SavePlot(ax, 5, tmpdir, 'png')
        assert os.path.isfile('{0}_5.png'.format(tmpdir))

    def test_plot_comparison(self, example_shp, tmpdir):
        matplotlib.pyplot.clf()
        x, y, _ = spatial_efd.ProcessGeometry(example_shp[0])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        ax = spatial_efd.InitPlot()
        spatial_efd.plotComparison(ax, coeffs, 10, x, y)
        spatial_efd.SavePlot(ax, 10, tmpdir, 'png')
        assert os.path.isfile('{0}_10.png'.format(tmpdir))

    def test_plot_comparison_norm(self, example_shp, tmpdir):
        x, y, _ = spatial_efd.ProcessGeometry(example_shp[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        coeffs, rotation = spatial_efd.normalize_efd(coeffs,
                                                     size_invariant=False)
        ax = spatial_efd.InitPlot()
        spatial_efd.plotComparison(ax, coeffs, 7, x, y, rotation=rotation)
        spatial_efd.SavePlot(ax, 7, tmpdir, 'png')
        assert os.path.isfile('{0}_7.png'.format(tmpdir))

    def test_plot_comparison_norm_size_invariant(self, example_shp, tmpdir):
        matplotlib.pyplot.clf()
        x, y, _ = spatial_efd.ProcessGeometry(example_shp[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        coeffs, rotation = spatial_efd.normalize_efd(coeffs,
                                                     size_invariant=True)
        ax = spatial_efd.InitPlot()
        spatial_efd.plotComparison(ax, coeffs, 7, x, y, rotation=rotation)
        spatial_efd.SavePlot(ax, 8, tmpdir, 'png')
        assert os.path.isfile('{0}_8.png'.format(tmpdir))

    def test_write_geometry(self, example_shp):
        x, y, _ = spatial_efd.ProcessGeometry(example_shp[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)
        assert isinstance(shape, shp.Writer)

    def test_generate_shapefile(self):
        shape = spatial_efd.generateShapefile()
        assert isinstance(shape, shp.Writer)

    def test_save_shapefile(self, example_shp, tmpdir):
        x, y, _ = spatial_efd.ProcessGeometry(example_shp[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)
        spatial_efd.saveShapefile(tmpdir.strpath, shape, prj=None)
        assert os.path.isfile('{}.shp'.format(tmpdir))

    def test_save_shapefile_prj(self, example_shp, tmpdir, shp_paths):
        x, y, _ = spatial_efd.ProcessGeometry(example_shp[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)
        spatial_efd.saveShapefile(tmpdir.strpath, shape, prj=shp_paths[1])
        assert os.path.isfile('{}.shp'.format(tmpdir))
        assert os.path.isfile('{}.prj'.format(tmpdir))

    def test_save_shapefile_prj_missing(self, example_shp, tmpdir, warn_missing_prj):
        x, y, _ = spatial_efd.ProcessGeometry(example_shp[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)

        with warnings.catch_warnings(record=True) as w:
            spatial_efd.saveShapefile(tmpdir.strpath, shape, prj='missing.prj')

            assert os.path.isfile('{0}.shp'.format(tmpdir))
            assert not os.path.isfile('{0}.prj'.format(tmpdir))
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)
            assert clean_warning(w[0].message) == warn_missing_prj

    def test_save_shapefile_prj_wrong(self, example_shp, tmpdir, shp_paths, warn_wrong_prj):
        x, y, _ = spatial_efd.ProcessGeometry(example_shp[1])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        shape = spatial_efd.generateShapefile()
        shape = spatial_efd.writeGeometry(coeffs, x, y, 4, shape, 1)

        with warnings.catch_warnings(record=True) as w:
            spatial_efd.saveShapefile(tmpdir.strpath, shape, prj=shp_paths[0])

            assert os.path.isfile('{0}.shp'.format(tmpdir))
            assert not os.path.isfile('{0}.prj'.format(tmpdir))
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)
            assert clean_warning(w[0].message) == warn_wrong_prj
