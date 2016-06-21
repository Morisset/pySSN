'''
Created on 21 juin 2016

@author: christophemorisset
'''
import numpy as np
import pyneb as pn
from pyneb.utils.physics import Z

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
        print('{}{:02d}000000{:02d}{:02d} {:10.2f} {:8.4f}'.format(Z[atom.elem], atom.spec, i+1, j_emis_ref_loc+1,wls[i,j_emis_ref_loc], emis[i, j_emis_ref_loc]/emis_ref)) 
        for j in np.arange(i):
            if emis[i,j] > cut * emis_ref:
                
                print('{}{:02d}000000{:02d}{:02d} {:10.2f} {:8.4f}'.format(Z[atom.elem], atom.spec, i+1, j+1, wls[i,j], emis[i,j]/emis[i, j_emis_ref_loc])) 
    