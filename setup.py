#!/usr/bin/env python
import string
from setuptools import setup


def readme():
    '''
    Turn relative figure paths to urls
    '''
    url = ('https://raw.githubusercontent.com/sgrieve/'
           'spatial_efd/master/_static/figure_')
    lines = []
    with open('README.rst') as f:
        for l in f.readlines():

            lines.append(string.replace(l, '_static/figure_', url))

    return string.join(lines)


setup(name='spatial_efd',
      version='1.0.1',
      description='Spatial elliptical fourier analysis',
      url='http://github.com/sgrieve',
      long_description=readme(),
      keywords='GIS elliptical fourier analysis shapefile',
      classifiers=['Development Status :: 4 - Beta',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 2.7',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: GIS'],
      author='Stuart WD Grieve',
      author_email='s.grieve@ed.ac.uk',
      license='MIT',
      packages=['spatial_efd'],
      install_requires=['matplotlib>=2.0.0', 'numpy>=1.12.0',
                        'pyshp>=1.2.10'],
      include_package_data=True,
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose>=1.3.7'])
