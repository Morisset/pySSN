'''
Created on 21 juin 2016

@author: christophemorisset
'''
import numpy as np
import pyneb as pn
from pyneb.utils.physics import Z
from pyneb.utils.misc import int_to_roman

def print_liste_phyat(atom, tem, den, cut=1e-3, ij_ref = None):
    emis = atom.getEmissivity(tem, den)
    if ij_ref is None:
        emis_ref = np.max(emis)
    else:
        emis_ref = emis[ij_ref[0]-1, ij_ref[1]-1]
    wls = atom.wave_Ang
    NLevels = emis.shape[0]
    str_print = ' {:02d}{:02d}00000{:02d}{:02d} {:2s}_{:5s} {:11.3f} 0.000{:10.3e}  1.000  {:02d}{:02d}00000{:02d}{:02d}   1   1.00 {} {}'
    for i in np.arange(NLevels):
        emis_i = emis[i,:]
        j_emis_ref_loc = np.argmax(emis_i)
        print(str_print.format(Z[atom.elem], 
                               atom.spec, 
                               i+1, 
                               0, 
                               atom.elem,
                               int_to_roman(atom.spec), 
                               1.0, 
                               emis[i, j_emis_ref_loc]/emis_ref,
                               Z[atom.elem], 
                               atom.spec, 
                               0, 0,
                               i+1, 
                               j_emis_ref_loc+1))
        for j in np.arange(i):
            if emis[i,j] > cut * emis_ref:
                print(str_print.format(Z[atom.elem], 
                                       atom.spec, 
                                       i+1, 
                                       j+1, 
                                       atom.elem, 
                                       int_to_roman(atom.spec), 
                                       wls[i,j], 
                                       emis[i,j]/emis[i, j_emis_ref_loc], Z[atom.elem], 
                                       atom.spec, 
                                       i+1, 0, 
                                       i+1, 
                                       j+1))


