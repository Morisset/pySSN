from pyssn import make_phyat_list
from pyssn.phyat_lists.manage_phyat_list import get_extra_atoms, get_atom_str

# Set to None (atoms = None) to generate all the available atoms
atoms = ['C4', 'C3', 'C2', 'O3', 'O2', 'Ne3', 'Ne2', 'Fe3']
#atoms = None

# Name of the output
filename = 'liste3.dat'
extra_file = 'phyat_list_DP_01.dat'

ref_lines_dic = {'Fe7': ((4, 2),),
		 'Fe6': ((2, 1), (5, 1),),
		 'Fe5': ((3, 2), (7, 5),),
		 'Fe4': ((6, 1),),
		 'Fe3': ((2, 1), (12, 1),),
		 }

NLevels_dic = {'S1': 8,
	       'Ni2': 17,
	       'Fe1': 9,
	       'Ni3': 9,
	       'Ni4': 17,
	       'Fe3': 25,
	       'Fe4': 26,
	       'Fe5': 20,
	       'Fe6': 19,
	       'Fe7': 9,
	       'Fe8': 2
	       }

up_lev_rule_dic = {'Ni2': 'split 3',
		   'Fe1': 'split 5',
		   'Ni3': 'split 3',
		   'Ni4': 'split 4',
		   'Fe3': 'split 5',
		   'Fe4': 'all',
		   'Fe5': 'split 5',
		   'Fe6': 'split 4',
		   'Fe7': 'split 3'
		   }

# Transitions set to 0.0
Aij_zero_dic = {'C3': ((2,1),)
		}

# Temperatures and densities for IP below the given value
tem_den_dic = {0.:   (1e4, 1e3),
	       13.6: (1e4, 1e3),
	       24.0: (1e4, 1e3),
	       1e6:  (1e4, 1e3)
	       }

if extra_file is not None and atoms is not None:
	atoms.extend(get_extra_atoms(extra_file, uniq=True))
	atoms = sorted(atoms, key=get_atom_str)

# Here we run the script the produce the phyat_list
make_phyat_list(filename, atoms=atoms, cut=1e-4, E_cut=20, cut_inter=1e-5, 
	       verbose=False, notry=False, NLevels=50, 
	       ref_lines_dic=ref_lines_dic,
	       NLevels_dic=NLevels_dic,
	       up_lev_rule_dic=up_lev_rule_dic,
	       Aij_zero_dic=Aij_zero_dic,
	       tem_den_dic = tem_den_dic,
	       extra_file=extra_file)
