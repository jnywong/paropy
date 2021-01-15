#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='paropy',
      version='0.0.2',
      description='Python package to process data from PARODY-JA4.3 dynamo simulations.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/jnywong/paropy',
      author='Jenny Wong',
      author_email='jenny.wong@univ-grenoble-alpes.fr',
      license='MIT',
      packages=['paropy'],
      package_data={'paropy': ['scripts/*.py']},
      setup_requires=["numpy"],
      install_requires=['pytest','cartopy'],
      )
