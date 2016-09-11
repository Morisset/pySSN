import numpy as np
import pyneb as pn
import matplotlib.pyplot as plt
from pyneb.utils.physics import Z, IP, gsFromAtom
from pyneb.utils.misc import int_to_roman, parseAtom
import os
from ..core.spectrum import read_data

pn.atomicData.setDataFile('fe_ii_atom.chianti')
pn.atomicData.setDataFile('fe_ii_coll.chianti')
#pn.atomicData.setDataFile('o_iii_atom.chianti')
#pn.atomicData.setDataFile('o_iii_coll.chianti')

# ToDo : faire les ions suivant en reprenant les rapports connus et les energies adaptees.
# Il faut verifier que l'on commence a ecrire uniquement quand on est sur d'avoir qqchose a ecrire.

def print_liste_phyat(atom, tem, den, cut=1e-3, ij_ref = None, filename=None, up_lev_rule=None, help_file=None, NLevels=None):
    """
    atom: a string e.g. 'O3' or a pn.Atom object
    up_lev_rule: 'each': each upper level are used to define a different master line
            'all': all the transitions are grouped into a single master line
            list, e.g. [2, 3, [4, 5]]: levels 2, 3 and 4-5 together are used to define a master line. 
            None: default hard coded value is used, see up_levs below
    """
    
    up_levs = {'s1': [[2,3]],  
               's2': 'each',
               'p1': [2, [3, 4, 5], [6,7]], 
               'p2': 'each',
               'p3': [2, 3, [4, 5]],
               'p4': 'each',
               'p5': 'each'} # The others are 'all' by default (i.e. only one master line, with all the lines depending on it

    if filename is None:
        def myprint(s, f):
            print(s)
    else:
        def myprint(s, f):
            f.write(s)

    str_print = '{}{:02d}{:02d}{:02d}0{:03d}{:03d} {:<8s} {:11.3f} 0.000{:10.3e}  1.000  {:02d}{:02d}{:1d}00{:03d}{:03d}   1   1.00 {}'  
    if filename is None:
        f = None
    else:
        str_print += '\n'
        if type(filename) is str:
            f = open(filename, 'w')
        elif type(filename) is file:
            f = filename
    
    if type(atom) is str:
        atom = pn.Atom(atom=atom, NLevels=NLevels) 
    emis = atom.getEmissivity(tem, den)
    emis_max_tot = np.max(emis)
    wls = atom.wave_Ang
    NLevels = emis.shape[0]
    gs = gsFromAtom(atom.atom)
    
    if up_lev_rule is None:
        if gs in up_levs:
            up_lev_rule = up_levs[gs]
        else:
            up_lev_rule = 'all'
    
    if up_lev_rule == 'each':
        up_lev_list = range(2, NLevels+1)
    elif up_lev_rule == 'all':
        up_lev_list = [range(2, NLevels+1)]
    else:
        up_lev_list = up_lev_rule
        
    #print('up_lev_rule: {}, up_lev_list: {}'.format(up_lev_rule, up_lev_list))
    NLevels_max = 0
    print_any = False
    to_print = []
    ion_name = '{:>2s}_{:<5s}'.format(atom.elem, int_to_roman(atom.spec)).strip()
    for up_lev in up_lev_list:
        if type(up_lev) is int:
            up_lev = [up_lev]
            
        emis_max = 0
        i_emis_ref_loc = None
        for i0 in up_lev:
            i = i0-1
            for j in np.arange(i):
                if emis[i,j] > emis_max:
                    i_emis_ref_loc = i
                    j_emis_ref_loc = j
                    emis_max = emis[i,j]

        if i_emis_ref_loc is not None:
            com_ref = '{:.0f}K 10**{:.0f}cm-3 {:.1f}'.format(tem, np.log10(den), wls[i_emis_ref_loc,j_emis_ref_loc])
            com_ref = '{} {} {:.1f}'.format(i_emis_ref_loc+1, j_emis_ref_loc+1, wls[i_emis_ref_loc,j_emis_ref_loc])
            ref_str=str_print.format('9', Z[atom.elem], atom.spec,i_emis_ref_loc+1 , 0, 0, ion_name, 1.0, 1.0, 0, 0, 0, 0, 999, com_ref)
            print_ref = True 
            
            for i0 in up_lev:
                i = i0-1            
                print_it = False
                for j in np.arange(i):
                    if emis[i,j] > cut * emis_max_tot and wls[i,j] > 912:
                        print_it = True
                        print_any = True
                        NLevels_max = i0+1
                if print_it:                
                    for j in np.arange(i):
                        if emis[i,j] > cut * emis_max_tot and wls[i,j] > 912:
                            com = '{} {}'.format(i+1, j+1)
                            ion_name = '{:>2s}_{:<5s}'.format(atom.elem, int_to_roman(atom.spec)).strip()
                            if print_ref:
                                to_print.append(ref_str)
                                print_ref = False
                            to_print.append(str_print.format(' ', Z[atom.elem], atom.spec, i+1, i+1, j+1, ion_name, 
                                                     wls[i,j], emis[i,j]/emis_max, Z[atom.elem], atom.spec, i_emis_ref_loc+1, 0, 0, com))
    

    if print_any:
        for str_ in to_print:
            myprint(str_, f)

    if help_file is not None and print_any:
        if type(help_file) is str:
            hf = open(help_file, 'w')
        else:
            hf = help_file
#         hf.write(' ------------------------- \n'.format())
        hf.write('{} \n'.format(atom.atom))
        hf.write('Temperature = {} K, Density = {} cm-3 \n'.format(tem, den))
                
        hf.write('Level           ' + ' '.join(('{:4}'.format(l) for l in np.arange(NLevels)[1:NLevels_max])) + '\n')
        hf.write('Critical density' + ' '.join(('{:4.1f}'.format(l) for l in np.log10(atom.getCritDensity(tem))[1:NLevels_max])) + '\n')
        hf.write(' ------------------------- \n'.format())

    if type(help_file) is str:
        hf.close()

def get_tem_den(IP):
    if IP < 14:
        return 1e4, 1e3
    elif IP < 25:
        return 1.e4, 1e3
    else:
        return 1.e4, 1e3

def make_all(tem1=None, den1=None, cut=1e-4, filename=None, help_file=None):
    atoms = get_atoms_by_conf()
    f = open(filename, 'w')
    if help_file is not None:
        hf = open(help_file, 'w')
    else:
        hf = None
        
    for a in atoms:
        try:
            print('trying {}'.format(a))
            atom = pn.Atom(atom=a, NLevels=30)
            
            if tem1 is None:
                try:
                    tem, den = get_tem_den(IP[atom.atom])
                except:
                    tem = 1e4
                    den = 1e3
            else:
                tem = tem1
                den = den1
            print_liste_phyat(atom, tem, den, cut=cut, filename=f, help_file=hf)
            print('{} done'.format(a))
        except:
            pass
    
    f.close()


def phyat2model(phyat_file, model_file):
    """  

    """
    
    str_print = '{}{:02d}{:02d}000{:03d}{:03d} {:<8s} {:11.3f} 0.000{:10.3e}  1.000  {:02d}{:02d}000{:03d}{:03d}   1   1.00 {}'

    list_phyat = read_data(phyat_file)
    with open(model_file, 'w') as f:
        for line in list_phyat:
            if line['ref'] == 999:
                f.write('{0:14d} {1[id]:8s}      1.000 0.000 1.000e+04  1.000  0000000000000   1   1.00 {1[comment]:15s}'.format(line['num']-90000000000000, line))

                        
def make_conf_plots():
    
    at_list = [pn.Atom(atom=at) for at in ('C4', 'C3', 'C2', 'O3', 'O2', 'Ne3', 'Ne2')]
    f, axes = plt.subplots(4,2, figsize=(15, 30))
    for at, ax in zip(at_list, axes.ravel()):
        at.plotGrotrian(ax=ax)
    f.savefig('conf_all.pdf')

def get_atoms_by_name():
    atoms = pn.atomicData.getAllAtoms()
    atoms.remove('3He2')
    Chianti_path = os.environ['XUVTOP']
    masterlist = '{0}/masterlist/masterlist.ions'.format(Chianti_path)
    with open(masterlist, 'r') as f:
        lines = f.readlines()
    chianti_ions = [l[0:7].strip() for l in lines]
    chianti_ions = [i.capitalize().replace('_','') for i in chianti_ions if i[-1] != 'd']
    for cion in chianti_ions:
        if cion not in atoms and cion not in ('H1', 'He1', 'He2'):
            try:
                atfiles = pn.atomicData.getAllAvailableFiles(cion)
                for atfile in atfiles:
                    pn.atomicData.setDataFile(atfile)
                atoms.append(cion)
                #print('Adding {} from Chianti'.format(cion))
            except:
                pass      
    return sorted(atoms)

def get_atoms_by_Z():
    atoms = get_atoms_by_name()
    return sorted(atoms, key=lambda k:Z[parseAtom(k)[0]])
    
def get_atoms_by_conf():
    atoms = get_atoms_by_Z()
    res = []
    gss = {}
    for atom in atoms:
        gss[atom] = gsFromAtom(atom)
    for gs in ('s1', 's2', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'unknown'):        
        for atom in atoms:
            if gss[atom] == gs:
                res.append(atom)
    return res
                
        
"""
To make a new liste_phyat, from an ipython session:

from pyssn.phyat_lists.make_phyat_list import make_all
make_all(tem1=1e4, den1=1e3, cut=1e-4, filename='liste_phyat_4_3.dat', help_file='HF_4_3.dat')

"""