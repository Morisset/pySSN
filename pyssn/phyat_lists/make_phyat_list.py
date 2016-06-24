import numpy as np
import pyneb as pn
from pyneb.utils.physics import Z, IP
from pyneb.utils.misc import int_to_roman

pn.atomicData.setDataFile('fe_ii_atom.chianti')
pn.atomicData.setDataFile('fe_ii_coll.chianti')

# ToDo : faire les ions suivant en reprenant les rapports connus et les energies adaptees.

def print_liste_phyat(atom, tem, den, cut=1e-3, ij_ref = None, filename=None, norm_by='sum'):
    
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
    str_print = '{}{:02d}{:02d}000{:03d}{:03d} {:>2s}_{:<5s} {:11.3f} 0.000{:10.3e}  1.000  {:02d}{:02d}000{:03d}{:03d}   1   1.00 {}'
    if filename is None:
        f = None
    else:
        str_print += '\n'
        if type(filename) is str:
            f = open(filename, 'w')
        elif type(filename) is file:
            f = filename
    com = '{:.0f}K 10**{:.0f}cm-3'.format(tem, np.log10(den))
    myprint(str_print.format('9', Z[atom.elem], atom.spec, 0, 0, atom.elem, int_to_roman(atom.spec), 1.0, 1.0, 0, 0, 9, 99, com), f)
    
    for i in np.arange(NLevels-1)+1:
        emis_i = emis[i,:]
        j_emis_ref_loc = np.argmax(emis_i)
        print_it = False
        sum = 0
        for j in np.arange(i):
            if emis[i,j] > cut * emis_ref:
                sum += emis[i,j]
                print_it = True
        if norm_by == 'sum':
            norm = sum
        else:
            norm = emis[i, j_emis_ref_loc]
        if print_it:
            com = '{:.1f}+'.format(wls[i,j_emis_ref_loc])
            myprint(str_print.format(' ', Z[atom.elem], atom.spec, i+1, 0, atom.elem, int_to_roman(atom.spec), 
                                     1.0, emis[i, j_emis_ref_loc]/emis_ref, Z[atom.elem], atom.spec, 0, 0, com), f)
            for j in np.arange(i):
                if emis[i,j] > cut * emis_ref:
                    com = '{} {}'.format(i+1, j+1)
                    myprint(str_print.format(' ', Z[atom.elem], atom.spec, i+1, j+1, atom.elem, int_to_roman(atom.spec), 
                                             wls[i,j], emis[i,j]/norm, Z[atom.elem], atom.spec, i+1, 0, com), f)
    

def get_tem_den(IP):
    if IP < 14:
        return 1e4, 1e3
    elif IP < 25:
        return 1.2e4, 3e3
    else:
        return 1.3e4, 5e3

def make_all(tem1=None, den1=None, cut=1e-3, filename=None, norm_by='sum'):
    aa = pn.atomicData.getAllAtoms()
    f = open(filename, 'w')
    for a in aa:
        try:
            print('trying {}'.format(a))
            atom = pn.Atom(atom=a)
            atom.NLevels = 30
            if tem1 is None:
                tem, den = get_tem_den(IP[atom.atom])
            else:
                tem = tem1
                den = den1
            print_liste_phyat(atom, tem, den, cut=cut, filename=f)
            print('{} done'.format(a))
        except:
            pass
    f.close()
