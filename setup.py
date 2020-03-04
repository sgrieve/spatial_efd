#!/usr/bin/env python
from setuptools import setup


def readme():

    # Need to replace the local image paths with the external links so the
    # readme works on pypi
    lut = {
           '_static/figure_1.png': 'https://raw.githubusercontent.com/sgrieve/spatial_efd/master/_static/figure_1.png',
           '_static/figure_2.png': 'https://raw.githubusercontent.com/sgrieve/spatial_efd/master/_static/figure_2.png',
           '_static/figure_3.png': 'https://raw.githubusercontent.com/sgrieve/spatial_efd/master/_static/figure_3.png',
           '_static/figure_4.png': 'https://raw.githubusercontent.com/sgrieve/spatial_efd/master/_static/figure_4.png'
          }

    with open('README.rst') as f:
        text = f.read()

        for k, v in lut.items():
            text = text.replace(k, v)

        return text


setup(name='spatial_efd',
      version='1.2.1',
      description='Spatial elliptical fourier analysis',
      url='http://github.com/sgrieve/spatial_efd',
      long_description=readme(),
      keywords='GIS elliptical fourier analysis shapefile',
      classifiers=['Development Status :: 5 - Production/Stable',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: GIS',
                   'Operating System :: OS Independent'],
      author='Stuart WD Grieve',
      author_email='stuart@swdg.io',
      license='MIT',
      packages=['spatial_efd'],
      setup_requires=['pytest-runner'],
      install_requires=['matplotlib>=2.0.0', 'numpy>=1.12.0',
                        'pyshp>=1.2.10', 'future'],
      include_package_data=True,
      zip_safe=False,
      test_suite='pytest-runner',
      tests_require=['pytest'])
