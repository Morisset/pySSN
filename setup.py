#!/usr/bin/env python

# to be used with python setup.py dist

from distutils.core import setup
from pyssn.version import __version__

setup(name='pySSN', 
      version=__version__,
      description='Python Spectral Synthesis for Nebulae',
      author='Christophe Morisset, Daniel Pequignot',
      author_email='chris.morisset@gmail.com',
      url='',
      py_modules=[],
      packages=['pyssn','pyssn.core','pyssn.utils','pyssn.qt'],
      package_data={'pyssn':['data/*']},
      entry_points={'console_scripts': ['pySSN = pyssn.qt.pyssn_qt:main']},
     )


"""
What to put in pySSN executable:
#!/usr/bin/env python
if __name__ == '__main__':
    import sys
    from pyssn.qt.pyssn_qt import main

    sys.exit(main())

"""
