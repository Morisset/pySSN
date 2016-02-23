# pySSN

pySSN is better running with Qt4. Once installed, you can call it from command line:

`pySSN [-f s6302_n_c_init.py]`

The init file is optional. If not given, the program will first ask you for one.

You can also call pySSN from within a python session:

`from pyssn.qt.pyssn_qt import main_loc`
`sp = main_loc('./s6302_n_c_init.py')`

You can interact with sp.


If you don't have Qt4 installed, you can still call pySSN with an init file (mandatory in this case):

`pySSN -f s6302_n_c_init.py` !!! This one is still in develepment.

And you can also call it from python session:

`from pyssn.core.spectrum import main_loc`

`sp = main_loc('./s6302_n_c_init.py')`

You will have the possibility to interact with sp.
