#!/usr/bin/env python

from setuptools import setup
from setuptools.command.build_ext import build_ext as _build_ext

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='parody_py',
      version='0.0.6',
      description='Python package to process data from PARODY-JA4.3 dynamo simulations.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/jnywong/parody_py',
      author='Jenny Wong',
      author_email='jenny.wong@univ-grenoble-alpes.fr',
      license='MIT',
      packages=['paropy'],
      package_data={'paropy': ['scripts/*.py']},
      setup_requires=["numpy"],
      install_requires=['pytest','cartopy'],
      )
