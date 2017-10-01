1. Install the distribution by: `pip install -U git+https://github.com/Morisset/pySSN.git`
1. Download the data file (107Mo) from http://132.248.1.102/~morisset/data.tgz
1. Uncompress it somewhere, e.g. /user/myname/pySSN_data/
1. Link this directory with the main program: pySSN_link2data /user/myname/pySSN_data/data
1. Compile the fortran program: pySSN_compile --compiler=gfortran

Pour faire une synthese:

1. donner les conditions physiques de la nébuleuse dans le fichier phy_cond.dat du style:
    `IP 1000000   CR 1e4    1e2`

1. determiner quel template va etre utilisé pour les fractions ioniques, voir README_ionfrac.dat.
    Par exemple 1789409_ionfrac.dat

1. determiner quelles abondances vont etre utilisées, par exemple: asplund_2009.dat

1. generer le fichier ions_rec.dat pour XSSN, depuis IPYTHON:

    `import pyssn`
    `pyssn.make_ionrec_file(abund_file='asplund_2009.dat', ion_frac_file='1789409_ionfrac.dat')`

1. Verifier que outputcond.dat est adapté au probleme (seuils de detection des raies en fonction du domaine observé).

1. Faire tourner le programme fortran XSSN_Phyat.exe obtenu avec 
    On obtient un ficher liste_phyat_rec.dat

1. faire tourner depuis IPYTHON:

    `%run generate_phyat_list.py`

    On obtient un fichier liste_phyat_coll.dat

1. ajouter tous les fichiers liste_phyat pour en former un seul, depuis IPYTHON:

    `from pyssn.phyat_lists.manage_phyat_list import merge_files`
    `merge_files(('../fortran/liste_phyat_rec.dat', 'liste_phyat_coll.dat', 'liste_phyat_others.dat'), 'liste_phyat_test1.dat')`

    On obtient un fichier liste_phyat_test1.dat

1. Generer un liste_model:

    `from pyssn.phyat_lists.manage_phyat_list import phyat2model`
    `phyat2model('liste_phyat_test1.dat', 'liste_model1.dat', norm_hbeta=1e4, ion_frac_file='1789409_ionfrac.dat', abund_file='asplund_2009.dat', ion_frac_min=1e-4)`

1. READY for the synthesis!!!

1. Par exemple, mettre dans test1_init.dat:
		#-------------------------------------------------------------------
		#   Initialisation
		#-------------------------------------------------------------------
		log_level = 3
		spectr_obs = None
		limit_sp = (4000, 7000)
		lambda_pix = 0.1 
		fic_modele = 'liste_model1.dat' #
		fic_cosmetik = ''
		do_cosmetik = False
		phyat_file = 'liste_phyat_test1.dat'
		plot_ax2 = False
		#----------------------------END -----------------------------------

1. Depuis LINUX: 
	pySSN -f test1_init.py


