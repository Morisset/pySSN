'''
Created on 21 juin 2016

@author: christophemorisset
'''
import numpy as np
import pyneb as pn

def print_liste_phyat(atom, tem, den, cut=1e-3, ij_ref = None):
    emis = atom.getEmissivity(tem, den)
    if ij_ref is None:
        emis_ref = np.max(emis)
    else:
        emis_ref = emis[ij_ref[0]-1, ij_ref[1]-1]
    wls = atom.wave_Ang
    NLevels = emis.shape[0]
    for i in np.arange(NLevels):
        emis_i = emis[i,:]
        j_emis_ref_loc = np.argmax(emis_i)
        print('{:.2f} {:.4f} {} {}'.format(wls[i,j_emis_ref_loc], emis[i, j_emis_ref_loc]/emis_ref, i+1, j_emis_ref_loc+1)) 
        for j in np.arange(i):
            if emis[i,j] > cut * emis_ref:
                print('{:.2f} {:.4f} {} {}'.format(wls[i,j], emis[i,j]/emis[i, j_emis_ref_loc], i+1, j+1)) 
    