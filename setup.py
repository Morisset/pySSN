#!/usr/bin/env python

# to be used with python setup.py dist

from distutils.core import setup
import pyssn

setup(name='pySSN', 
      version=pyssn.__version__,
      description='Python Spectral Synthesis for Nebulae',
      author='Christophe Morisset, Daniel Pequignot',
      author_email='chris.morisset@gmail.com',
      url='',
      py_modules=[],
      packages=['pyssn','pyssn.core','pyssn.utils','pyssn.qt'],
      package_data={'pyssn':['data/*']}
     )

