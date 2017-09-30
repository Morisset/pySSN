Pour faire une synthese:

# donner les conditions physiques de la nébuleuse dans le fichier phy_cond.dat du style:
IP 1000000   CR 1e4    1e2

# determiner quel template va etre utilisé pour les fractions ioniques, voir README_ionfrac.dat.
Par exemple 1789409_ionfrac.dat

# determiner quelles abondances vont etre utilisées, par exemple: asplund_2009.dat

# generer le fichier ions_rec.dat pour XSSN, depuis IPYTHON:

import pyssn
pyssn.make_ionrec_file(abund_file='asplund_2009.dat', ion_frac_file='1789409_ionfrac.dat')

# Verifier que outputcond.dat est adapté au probleme (seuils de detection des raies en fonction du domaine observé).

# Compiler le programme fortran:

import pyssn
pyssn.compile_XSSN()


# Faire tourner le programme fortran XSSN_Phyat.exe obtenu avec 
On obtient un ficher liste_phyat_rec.dat

# faire tourner depuis IPYTHON:

%run generate_phyat_list.py
On obtient un fichier 
liste_phyat_coll.dat

# ajouter tous les fichiers liste_phyat pour en former un seul, depuis IPYTHON:

from pyssn.phyat_lists.manage_phyat_list import merge_files
merge_files(('../fortran/liste_phyat_rec.dat', 'liste_phyat_coll.dat', 'liste_phyat_others.dat'), 'liste_phyat_test1.dat')
On obtient un fichier liste_phyat_test1.dat

# Generer un liste_model:

from pyssn.phyat_lists.manage_phyat_list import phyat2model
phyat2model('liste_phyat_test1.dat', 'liste_model1.dat', norm_hbeta=1e4, ion_frac_file='1789409_ionfrac.dat', abund_file='asplund_2009.dat', ion_frac_min=1e-4)

# READY for the synthesis!!!

Par exemple, mettre dans test1_init.dat:
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

Depuis LINUX: 
pySSN -f test1_init.py

