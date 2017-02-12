from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='spatial_efd',
      version='0.1',
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
      install_requires=['matplotlib', 'numpy', 'pyshp'],
      include_package_data=True,
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose'])
