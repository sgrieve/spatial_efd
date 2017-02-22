Spatial Elliptical Fourier Descriptors
=======================================

.. image:: https://travis-ci.org/sgrieve/spatial_efd.svg?branch=master
    :target: https://travis-ci.org/sgrieve/spatial_efd

.. image:: https://img.shields.io/codecov/c/github/sgrieve/spatial_efd.svg
    :target: https://codecov.io/github/sgrieve/spatial_efd

.. image:: https://requires.io/github/sgrieve/spatial_efd/requirements.svg?branch=master
     :target: https://requires.io/github/sgrieve/spatial_efd/requirements/?branch=master

.. image:: https://readthedocs.org/projects/spatial-efd/badge/?version=latest
     :target: http://spatial-efd.readthedocs.io/en/latest/?badge=latest

.. image:: https://img.shields.io/badge/License-MIT-green.svg
    :target: https://opensource.org/licenses/MIT


A pure python implementation of the elliptical Fourier analysis method described by Kuhl and Giardina (1982). This package is designed to allow the rapid analysis of spatial data stored as ESRI shapefiles, handling all of the geometric conversions. The resulting data can be written back to shapefiles to allow analysis with other spatial data or can be plotted using matplotlib.

The code is built upon the `pyefd module <https://github.com/hbldh/pyefd>`_ and it is hoped that this package will allow more geoscientists to apply this technique to analyze spatial data using the elliptical Fourier descriptor technique as there is no longer a data conversion barrier to entry. This package is also more feature rich than previous implementations, providing calculations of Fourier power and spatial averaging of collections of ellipses.

.. figure:: docs/figure_1.png
    :width: 600px
    :align: center
    :alt: spatial_efd example
    :figclass: align-center

    Examples of Fourier ellipses (black) being fitted to a shapefile outline (red), for increasing numbers of harmonics.

Features
--------

- Built-in geometry processing, just pass in a shapefile and get results quickly!
- Fourier coefficient average and standard devation calculation
- Handles spatial input data through the pyshp library
- Compute an appropriate number of harmonics for a given polygon
- Basic plotting for analysis and debugging through matplotlib
- Write Fourier ellipses as shapefiles

Installation
------------

Install ``spatial_efd`` by running::

  $ pip install spatial_efd

Dependencies
------------

This package is developed on Linux for Python 2.7 and requires ``matplotlib``, ``numpy``, ``nose`` and ``pyshp``. These packages will all install automatically if ``spatial_efd`` is installed using ``pip``.

Dependencies can be tracked by visiting `requires.io <https://requires.io/github/sgrieve/spatial_efd/requirements/?branch=master>`_

Tests
----------

A range of unit tests are included in the /spatial/tests/ directory. These can
be run using nose or directly from setup.py::

  $ python setup.py test
  $ nosetests


Many of these tests make use of the ``example_data.shp`` file which is a shapefile containing six polygons taken from a real dataset of landslide source areas.

Usage
----------

Here are examples of how to use the code.

Contribute
----------

.. image:: https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat
    :target: https://codecov.io/github/sgrieve/spatial_efd/issues

I welcome contributions to the code, head to the issue tracker on github to get involved!

- `Issue Tracker <github.com/sgrieve/spatial_efd/issues>`_
- `Source Code <github.com/sgrieve/spatial_efd>`_

Support
-------

If you find any bugs, have any questions or would like to see a feature in a new version, drop me a line:

- Twitter: `@GIStuart <https://www.twitter.com/GIStuart>`_
- Email: s.grieve@ed.ac.uk

License
-------

The project is licensed under the MIT license.

References
-----------

Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed contour. Computer graphics and image processing, 18(3), 236-258.
