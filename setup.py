#!/usr/bin/env python
from setuptools import setup

entry_points = {'console_scripts': 'acisfp_check = acisfp_check.acisfp_check:main'}

setup(name='acisfp_check',
      packages=["acisfp_check"],
      use_scm_version=True,
      setup_requires=['setuptools_scm', 'setuptools_scm_git_archive'],
      description='ACIS Thermal Model for FPTEMP',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/acisops/acisfp_check',
      include_package_data=True,
      entry_points=entry_points,
      )
