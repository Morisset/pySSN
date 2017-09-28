#!/usr/bin/env python

# to be used with python setup.py dist

#from distutils.core import setup
from os import path
from setuptools import setup, find_packages
from pyssn.version import __version__

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(name='pySSN', 
      version=__version__,
      description='Python Spectral Synthesis for Nebulae',
      long_description=long_description,
      author='Christophe Morisset, Daniel Pequignot',
      author_email='chris.morisset@gmail.com',
      url='https://github.com/Morisset/pySSN',
      py_modules=[],
      packages=['pyssn','pyssn.core','pyssn.utils','pyssn.qt',
		'pyssn.phyat_lists'],
      package_data={'pyssn':['data/*']},
      entry_points={'console_scripts': ['pySSN = pyssn.pySSN_exec:main',
					'pySSN_generate_col = pyssn.phyat_lists/generate_phyat_list:generate-col']},
     )


