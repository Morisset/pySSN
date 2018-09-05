PYTHON SETUP
======

If you want to run pySSN in a separate python installation 
(because you will need python 2.7 and QT4 and you do not want to change your actual 
python environment), you can just create a new environment using:

`conda create -n env_pySSN python=2.7 numpy scipy matplotlib pyqt=4`

Then swich to this envirponment:

`source activate env_pySSN`

Install pyneb:

`pip install pyneb` 

INSTALLATION OF PYSSN (DO IT ONLY ONCE)
======


1. Install the last version (development one) of PyNeb by `pip install -U git+https://github.com/Morisset/PyNeb_devel.git`. You will need the Chianti database installed (`http://www.chiantidatabase.org/download/CHIANTI_8.0.6_data.tar.gz`) and the environment variable XUVTOP pointing to it.
1. Install the distribution by `pip install -U git+https://github.com/Morisset/pySSN.git`
1. Download the data file (107Mo) from http://132.248.3.66/~morisset/data.tgz
1. Uncompress it somewhere, e.g. /user/myname/pySSN_data/
1. Link this directory with the main program: `pySSN_link2data /user/myname/pySSN_data/data`
1. Compile the fortran program: `pySSN_compile --compiler=gfortran`

DEFINING SOME CONDITIONS (IN EACH WORKING DIRECTORY)
======
From a working directory:

1. Define the physical conditions (Te, Ne) by IP ranges in the `phy_cond.dat`. You can obtain an example by `pySSN_write_files -P`
1. Define the ionic fraction. Examples can be obtained by `pySSN_write_files -I --Teff=50 --logU=-2 --B=R`, where Teff is in 25, 50, 75, 100, 150, 300, logU in -1, -2, -3 and B (Bounded) in R or M60 (matter-bounded 60%). This creates for example a file named `50_-2_R_ionfrac.dat`.
1. Define the elemental abundances. You can have an example with `pySSN_write_files -A` or `pySSN_write_files -D`. It creates a file named `abunds.dat` from Asplund 2009 data and Asplund data depleted following Jenkins 2014 respectively.
1. Define the output conditions. An example is obtained with `pySSN_write_files -O ` which generates a file named `outputcond.dat`.
1. All these commands can be done simultaneously: `pySSN_write_files -P -I -D -O --Teff=50 --logU=-2 --B=R`

GENERATE THE `liste_phyat` AND `liste_model` FILES
====

`pySSN_phyat --abund_file=abunds.dat --ion_frac_file=50_-2_R_ionfrac.dat --phy_cond_file=phy_cond.dat --outputcond_file=outputcond.dat --phyat_file=liste_phyat_test1.dat [--pynebdatafiles=pynebatomicdata.dat] [--ion_only='H_I, C_II']`

The optional pynebdatafiles argument is a file containing the list of atomic data files to be used instead of the default ones.

`pySSN_model --abund_file=abunds.dat --ion_frac_file=50_-2_R_ionfrac.dat --phyat_file=liste_phyat_test1.dat --model_file=liste_model_test1.dat --norm_hbeta=10000 --ion_frac_min=0.0000000001`


RUN SYNTHESIS
=====


Define for example test1_init.py:

	spectr_obs = None
	limit_sp = (4000, 7000)
	lambda_pix = 0.1 
	fic_modele = 'liste_model_test1.dat' #
	fic_cosmetik = ''
	do_cosmetik = False
	phyat_file = 'liste_phyat_test1.dat'

From LINUX: 	

`pySSN -f test1_init.py`


