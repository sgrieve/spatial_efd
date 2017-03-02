#!/usr/bin/env python
from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='spatial_efd',
      version='1.0.4',
      description='Spatial elliptical fourier analysis',
      url='http://github.com/sgrieve/spatial-efd',
      long_description=readme(),
      keywords='GIS elliptical fourier analysis shapefile',
      classifiers=['Development Status :: 5 - Production/Stable',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 2.7',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: GIS'],
      author='Stuart WD Grieve',
      author_email='s.grieve@qmul.ac.uk',
      license='MIT',
      packages=['spatial_efd'],
      install_requires=['matplotlib>=2.0.0', 'numpy>=1.12.0',
                        'pyshp>=1.2.10'],
      include_package_data=True,
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose>=1.3.7'])
