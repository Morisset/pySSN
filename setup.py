#!/usr/bin/env python


#from distutils.core import setup
from os import path
from setuptools import setup, find_packages
from pyssn.version import __version__

# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
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
                'pyssn.phyat_lists', 'pyssn.fortran'],
      package_data={'pyssn':['data/*', 'pyneb_data/*'],
                    'pyssn.fortran':['XSSN_Phyat.f', 'outputcond.dat',
                                     'liste_phyat_rec.dat',
                                     'ions_rec.dat',
                                     'res/*',
                                     'data_lab/*'],
                    'pyssn.phyat_lists':['asplund_2009.dat', 'phy_cond.dat',
                                         'liste_phyat_others.dat',
                                         'ionfracs/*']
                    },
      entry_points={'console_scripts': ['pySSN = pyssn.pySSN_exec:main',
                            'pySSN_link2data = pyssn.fortran.link_data:link2data',
                            'pySSN_compile = pyssn.fortran.compileit:compile_XSSN',
#                            'pySSN_write_phycond = pyssn.phyat_lists.entries:print_phy_cond',
#                            'pySSN_write_ionfrac = pyssn.phyat_lists.entries:print_ionfrac',
#                            'pySSN_write_asplund = pyssn.phyat_lists.entries:print_asplund',
#                            'pySSN_write_outputcond = pyssn.phyat_lists.entries:print_outputcond',
                            'pySSN_write_files = pyssn.phyat_lists.entries:print_files',
                            'pySSN_makelists = pyssn.phyat_lists.generate_phyat_list:make_all_lists']},
     )


