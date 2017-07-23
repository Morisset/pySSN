# pySSN


Installation
=======

Easily installed form github:

`pip install -U git+https://github.com/Morisset/pySSN.git`

Dependencies
=========

pySSN requieres numpy, matplotlib, scipy

It also requires PyNeb if you want to have reddening correction. https://github.com/Morisset/PyNeb_devel

pySSN is better running with Qt4.

Qt5 is still not supported

Usage
====

Once installed, you can call it from command line:

`pySSN [-f s6302_n_c_init.py]`

The init file is optional. If not given, the program will first ask you for one.

You can also call pySSN from within a python session:

`from pyssn.qt.pyssn_qt import main_loc`

`sp = main_loc('./s6302_n_c_init.py')`

You can interact with sp.


If you don't have Qt4 installed, you can still call pySSN with an init file (mandatory in this case):

`pySSN -f s6302_n_c_init.py`

Acknowledgements
================

This project is partly supported by grants DGAPA/PAPIIT-107215 and CONACyT-CB2015-254132.
