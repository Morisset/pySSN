from pyssn import make_phyat_list
import pyneb as pn

# Set to None (atoms = None) to generate all the available atoms
atoms = None
#atoms = ['C4', 'C3', 'C2', 'O3', 'O2', 'Ne3', 'Ne2', 'Fe3']
atoms = ['Fe14']

# File containing additional data in the phyat_list format
extra_file = None
# extra_file = 'phyat_list_DP_nouv.dat'
#extra_file = 'phyat_list_bidon.dat'

# Name of the output
filename = 'listeFe14.dat'
#filename = 'listeCIIIt1.5_1.dat'

ref_lines_dic = {'Fe7': ((4, 2),),
		 'Fe6': ((2, 1), (5, 1),),
		 'Fe5': ((3, 2), (7, 5),),
		 'Fe4': ((6, 1),),
		 'Fe3': ((2, 1), (12, 1),),
		 'Fe2': ((2, 1), (14, 6))
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
	       'Fe8': 2,
	       'Fe14': 10
	       }

up_lev_rule_dic = {'Ni2': 'split 3',
		   'Fe1': 'split 5',
		   'Ni3': 'split 3',
		   'Ni4': 'split 4',
		   'Fe2': 'split 9',
		   'Fe3': 'split 5',
		   'Fe4': 'all',
		   'Fe5': 'split 5',
		   'Fe6': 'split 4',
		   'Fe7': 'split 3',
		   'Ca2': (2, 3 , (4, 5), 6)
		   }

# Transitions set to 0.0
Aij_zero_dic = {'C3': ((2,1),)
		}

Del_ion = ('Ne7', 'Na8', 'Mg9', 'Al10', 'Ar7', 'Ca9', 'Ca14', 'Ti16', 'Cr18')

# Temperatures and densities for IP below the given value
tem_den_dic = {0.:   (1e4, 1e3),
	       13.6: (1e4, 1e3),
	       24.0: (1e4, 1e3),
	       1e6:  (1e4, 1e3)
	       }
tem_den_dic = {1e6:  (4e4, 1e8)}

# Here we run the script the produce the phyat_list
make_phyat_list(filename, atoms=atoms, cut=1e-4, E_cut=20,
	       verbose=False, notry=True, NLevels=50, 
	       ref_lines_dic=ref_lines_dic,
	       NLevels_dic=NLevels_dic,
	       up_lev_rule_dic=up_lev_rule_dic,
	       Aij_zero_dic=Aij_zero_dic,
	       Del_ion = Del_ion,
	       tem_den_dic = tem_den_dic,
	       extra_file=extra_file)
