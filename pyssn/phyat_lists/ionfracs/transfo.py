'''
Created on 2 oct. 2017

@author: christophemorisset
'''

import os
import numpy as np

data = np.genfromtxt('all_ionfracs.dat', dtype=None, names=True)
for d in data:
    os.rename('{}{}{}_ionfrac.dat'.format(d['Teff'], d['logU'], d['MB']), '{}_{}_{}_ionfrac.dat'.format(d['Teff'], d['logU'], d['MB']))
