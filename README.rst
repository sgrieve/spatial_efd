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


A pure python implementation of the elliptical Fourier analysis method described by Kuhl and Giardina (1982). This package is designed to allow the rapid analysis of spatial data stored as ESRI shapefiles, handling all of the geometric conversions. The code is built upon the pyefd module and it is hoped that this package will allow more geoscientists to apply this technique to analyze spatial data using the elliptical Fourier descriptor technique.

[Image of a drainage basin fitted fitted by the first 2,4,10,20 harmonics of its Fourier series]

Features
--------

- Built-in geometry processing, just pass in a shapefile and get results quickly!
- Fourier coefficient average and standard devation calculation
- Handles spatial input data through the pyshp library
- Compute an appropriate number of harmonics for a given polygon
- Basic plotting for analysis and debugging through matplotlib

Installation
------------

Install spatial_efd by running::

  $ pip install spatial_efd

Dependencies
------------

This package is developed on linux for python 2.7 and requires matplotlib, numpy and pyshp. These packages will all install automatically if spatial_efd is installed using pip.


Tests
----------

A range of unit tests are included in the /spatial/tests/ directory. These can
be run using nose or directly from setup.py::

  $ python setup.py test
  $ nosetests


Usage
----------

Here are examples of how to use the code.

Contribute
----------

.. image:: https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat
    :target: https://codecov.io/github/sgrieve/spatial_efd/issues

I welcome contributions to the code, head to the issue tracker on github to get involved!

- Issue Tracker: github.com/sgrieve/spatial_efd/issues
- Source Code: github.com/sgrieve/spatial_efd

Support
-------

If you find any bugs, have any questions or would like to see a feature in a new version, drop me a line:

- Twitter: @GIStuart
- Email: s.grieve@ed.ac.uk

License
-------

The project is licensed under the MIT license.

References
-----------

Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed contour. Computer graphics and image processing, 18(3), 236-258.


API
----

:ref:`Click here <API-ref>` for the module level documentation.
