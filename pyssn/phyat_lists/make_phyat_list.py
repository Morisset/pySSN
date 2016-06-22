'''
Created on 21 juin 2016

@author: christophemorisset
'''
import numpy as np
import pyneb as pn
from pyneb.utils.physics import Z
from pyneb.utils.misc import int_to_roman

def print_liste_phyat(atom, tem, den, cut=1e-3, ij_ref = None, filename=None):
    
    if filename is None:
        def myprint(s, f):
            print(s)
    else:
        def myprint(s, f):
            f.write(s)
            
    emis = atom.getEmissivity(tem, den)
    if ij_ref is None:
        emis_ref = np.max(emis)
    else:
        emis_ref = emis[ij_ref[0]-1, ij_ref[1]-1]
    wls = atom.wave_Ang
    NLevels = emis.shape[0]
    str_print = '{}{:02d}{:02d}00000{:02d}{:02d} {:>2s}_{:<5s} {:11.3f} 0.000{:10.3e}  1.000  {:02d}{:02d}00000{:02d}{:02d}   1   1.00 {}'
    if filename is None:
        f = None
    else:
        str_print += '\n'
        if type(filename) is str:
            f = open(filename, 'w')
        elif type(filename) is file:
            f = filename
    
    myprint(str_print.format('9', Z[atom.elem], atom.spec, 0, 0, atom.elem, int_to_roman(atom.spec), 1.0, 1.0, 0, 0, 9, 99, ''), f)
    
    for i in np.arange(NLevels-1)+1:
        emis_i = emis[i,:]
        j_emis_ref_loc = np.argmax(emis_i)
        print_it = False
        for j in np.arange(i):
            if emis[i,j] > cut * emis_ref:
                print_it = True
        if print_it:
            com = '{:.1f}+'.format(wls[i,j_emis_ref_loc])
            myprint(str_print.format(' ', Z[atom.elem], atom.spec, i+1, 0, atom.elem, int_to_roman(atom.spec), 
                                     1.0, emis[i, j_emis_ref_loc]/emis_ref, Z[atom.elem], atom.spec, 0, 0, com), f)
            for j in np.arange(i):
                if emis[i,j] > cut * emis_ref:
                    com = '{} {}'.format(i+1, j+1)
                    myprint(str_print.format(' ', Z[atom.elem], atom.spec, i+1, j+1, atom.elem, int_to_roman(atom.spec), 
                                             wls[i,j], emis[i,j]/emis[i, j_emis_ref_loc], Z[atom.elem], atom.spec, i+1, 0, com), f)
    

def make_all(tem, den, filename=None):
    aa = pn.atomicData.getAllAtoms()
    f = open(filename, 'w')
    for a in aa:
        try:
            atom = pn.Atom(atom=a)
            print_liste_phyat(atom, tem, den, filename=f)
        except:
            pass
    