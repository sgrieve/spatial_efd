#!/usr/bin/env python
import os
import re
import json
import pytest
import warnings
import matplotlib
import shapefile as shp
import numpy.testing as ntest
from spatial_efd import spatial_efd


@pytest.fixture
def expected():
    filepath = os.path.realpath(os.path.join(os.getcwd(),
                                             os.path.dirname(__file__)))
    with open(os.path.join(filepath, 'fixtures/expected.json')) as f:
        return json.loads(f.read())


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
    shppath = os.path.join(filepath, 'fixtures/example_data.shp')
    prjpath = os.path.join(filepath, 'fixtures/example_data.prj')
    return shppath, prjpath


@pytest.fixture
def example_shp():
    return spatial_efd.LoadGeometries(shp_paths()[0])


@pytest.fixture
def warn_wrong_prj():
    return 'The file supplied is not a prj file. No .prj file will be written'


@pytest.fixture
def warn_missing_prj():
    return ('The .prj file supplied does not exist.'
            ' No .prj file will be written')


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
        assert c == centre

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

    def test_rotate_contour(self, open_square, expected):
        x, y = spatial_efd.RotateContour(*open_square, rotation=30.,
                                         centroid=(5, 5))

        ntest.assert_almost_equal(x, expected['rotate_contour']['x'])
        ntest.assert_almost_equal(y, expected['rotate_contour']['y'])

    def test_rotate_point(self, expected):
        rx, ry = spatial_efd.rotatePoint((3., 2.), (1., 1.), 73.)
        assert pytest.approx(rx) == expected['rotate_point']['rx']
        assert pytest.approx(ry) == expected['rotate_point']['ry']

    def test_norm_contour(self, open_square, expected):
        x, y, c = spatial_efd.NormContour(*open_square, rawCentroid=(5., 5.))
        assert pytest.approx(c) == expected['norm_contour']['c']
        assert pytest.approx(x) == expected['norm_contour']['x']
        assert pytest.approx(y) == expected['norm_contour']['y']

    def test_get_bbox_dimensions(self, open_square, expected):
        xw, yw, xmin, ymin = spatial_efd.getBBoxDimensions(*open_square)
        assert xw == expected['get_bbox_dimensions']['xw']
        assert yw == expected['get_bbox_dimensions']['yw']
        assert xmin == expected['get_bbox_dimensions']['xmin']
        assert ymin == expected['get_bbox_dimensions']['ymin']

    def test_load_geometry(self, example_shp):
        assert isinstance(example_shp[0], shp._ShapeRecord)

    def test_process_geometry(self, example_shp, expected):
        x, y, c = spatial_efd.ProcessGeometry(example_shp[1])
        ntest.assert_almost_equal(c, expected['process_geometry']['c'])
        ntest.assert_almost_equal(x[:10], expected['process_geometry']['x'])
        ntest.assert_almost_equal(y[:10], expected['process_geometry']['y'])

    def test_process_geometry_norm(self, example_shp, expected):
        x, y, c = spatial_efd.ProcessGeometryNorm(example_shp[1])
        ntest.assert_almost_equal(c, expected['process_geometry_norm']['c'])
        ntest.assert_almost_equal(x[:10],
                                  expected['process_geometry_norm']['x'])
        ntest.assert_almost_equal(y[:10],
                                  expected['process_geometry_norm']['y'])

    def test_calculate_efd(self, example_shp, expected):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        ntest.assert_almost_equal(coeffs[6],
                                  expected['calculate_efd']['coeffs'])

    def test_inverse_transform(self, example_shp, expected):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs)

        ntest.assert_almost_equal(a[:5], expected['inverse_transform']['a'])
        ntest.assert_almost_equal(b[:5], expected['inverse_transform']['b'])

    def test_inverse_transform_locus(self, example_shp, expected):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        a, b = spatial_efd.inverse_transform(coeffs, locus=(0.5, 0.9))

        ntest.assert_almost_equal(a[:5],
                                  expected['inverse_transform_locus']['a'])
        ntest.assert_almost_equal(b[:5],
                                  expected['inverse_transform_locus']['b'])

    def test_average_coefficients(self, example_shp, expected):
        coeffsList = []

        for i in range(3):
            x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList)
        ntest.assert_almost_equal(avg[6],
                                  expected['average_coefficients']['avg'])

    def test_average_sd(self, example_shp, expected):
        coeffsList = []

        for i in range(3):
            x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[i])
            coeffsList.append(spatial_efd.CalculateEFD(x, y, 10))

        avg = spatial_efd.AverageCoefficients(coeffsList)
        sd = spatial_efd.AverageSD(coeffsList, avg)
        ntest.assert_almost_equal(sd[3], expected['average_sd']['sd'])

    def test_fourier_power(self, example_shp):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        coeffs = spatial_efd.CalculateEFD(x, y, 500)
        n = spatial_efd.FourierPower(coeffs, x)
        assert n == 19

    def test_normalize_efd(self, example_shp, expected):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[0])
        coeffs = spatial_efd.CalculateEFD(x, y, 10)
        coeffs, rotation = spatial_efd.normalize_efd(coeffs)

        ntest.assert_almost_equal(coeffs[9],
                                  expected['normalize_efd']['coeffs'])
        assert pytest.approx(rotation) == expected['normalize_efd']['rotation']

    def test_calculate_dc_coefficients(self, example_shp, expected):
        x, y, _ = spatial_efd.ProcessGeometryNorm(example_shp[2])
        dc = spatial_efd.calculate_dc_coefficients(x, y)
        assert pytest.approx(dc) == expected['calculate_dc_coefficients']['dc']

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

    def test_save_shapefile_prj_missing(self, example_shp, tmpdir,
                                        warn_missing_prj):
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

    def test_save_shapefile_prj_wrong(self, example_shp, tmpdir, shp_paths,
                                      warn_wrong_prj):
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
