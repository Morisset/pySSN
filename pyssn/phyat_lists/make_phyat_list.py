import numpy as np
import pyneb as pn
import matplotlib.pyplot as plt
from pyneb.utils.physics import Z, IP, gsFromAtom
from pyneb.utils.misc import int_to_roman, parseAtom
import os
import time
from ..core.spectrum import read_data

#pn.atomicData.setDataFile('o_iii_atom.chianti')

#pn.atomicData.setDataFile('o_iii_atom.chianti')
#pn.atomicData.setDataFile('o_iii_coll.chianti')

# ToDo : faire les ions suivant en reprenant les rapports connus et les energies adaptees.

def print_liste_phyat(atom, tem, den, cut=1e-3, cut_inter=1e-5, ij_ref = None, filename=None, up_lev_rule=None,
                      NLevels=None, E_cut=20, log_file=None, verbose=False):
    """
    atom: a string e.g. 'O3' or a pn.Atom object
    tem in K
    den in cm-3
    cut: relative intensity (to the master one) of a line to be printed
    cut_inter: largest dynamic between ref lines.
    ij_ref: None: is computed. Otherwise must be of the form e.g. (2, 4)
    filename: where to output the result. May be None, a string or a file object  
    up_lev_rule: 'each': each upper level are used to define a different master line
            'all': all the transitions are grouped into a single master line
            list, e.g. [2, 3, [4, 5]]: levels 2, 3 and 4-5 together are used to define a master line. 
            None: default hard coded value is used, see up_levs below
    NLevels : max number of levels for the atom
    E_cut: max energy in eV for the energy level to give a line
    log_file: for hte log
    """
    
    up_levs = {'s1': [[2, 3], 4, 5, 6],  
               's2': 'each',
               'p1': [2, [3, 4, 5], [6, 7]], 
               'p2': 'each',
               'p3': [2, 3, [4, 5]],
               'p4': 'each',
               'p5': 'each'} # The others are 'all' by default (i.e. only one master line, with all the lines depending on it

    if filename is None:
        def myprint(s):
            print(s)
    else:
        def myprint(s):
            filename.write(s)
    if log_file is None:
        def logprint(s):
            print(s)
    else:
        def logprint(s):
            log_file.write(s)
        
    str_print = '{}{:02d}{:02d}{:01d}{:01d}0{:03d}{:03d} {:<8s} {:11.3f} 0.000{:10.3e}  1.000  {:02d}{:02d}{:01d}00{:03d}{:03d}   1   1.00 {}'  
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
    if atom.NLevels == 0:
        logprint('Atom without data. '.format(atom.atom))
        return None        
    energies = atom.getEnergy(unit='eV')
    N_energies = (energies < E_cut).sum()
    if NLevels is not None:
        this_NLevels = np.min((NLevels,atom.NLevels, N_energies))
    else:
        this_NLevels = np.min((atom.NLevels, N_energies))
        
    if this_NLevels <=1:
        logprint('Atom without interesting lines. '.format(atom.atom))
        return None
    if verbose:
        print('this_NLevels={}'.format(this_NLevels))
    #print('Doing {} with NLevels={}'.format(atom.atom, this_NLevels))
    atom = pn.Atom(atom=atom.atom, NLevels=this_NLevels)
    logprint('{} levels. '.format(this_NLevels))
    emis = atom.getEmissivity(tem, den)
    emis_max_tot = np.max(emis)
    wls = atom.wave_Ang[0:this_NLevels, 0:this_NLevels]
    wls_mask = wls < 911.0
    emis[wls_mask] = 0.0
    gs = gsFromAtom(atom.atom)
    
    if up_lev_rule is None:
        if gs in up_levs:
            up_lev_rule = up_levs[gs]
        else:
            up_lev_rule = 'all'
    
    if up_lev_rule == 'each':
        up_lev_list = range(2, this_NLevels+1)
    elif up_lev_rule == 'all':
        up_lev_list = [range(2, this_NLevels+1)]
    else:
        up_lev_list = up_lev_rule
        
    #print('up_lev_rule: {}, up_lev_list: {}'.format(up_lev_rule, up_lev_list))
    NLevels_max = 0
    print_any = False
    to_print = []
    ion_name = '{:>2s}_{:<5s}'.format(atom.elem, int_to_roman(atom.spec)).strip()
    N_lines = 0
    for up_lev in up_lev_list:
        if type(up_lev) is int:
            up_lev = [up_lev]
        emis_ref = 0
        i_emis_ref_loc = None
        j_emis_ref_loc = None
        if ij_ref is None:
            for i in up_lev:
                if i < this_NLevels+1:
                    for j in 1+np.arange(i):                        
                        if emis[i-1,j-1] > emis_ref:
                            i_emis_ref_loc = i
                            j_emis_ref_loc = j
                            emis_ref = emis[i-1,j-1]
                        if verbose:
                            print('{} {} {}'.format(i,j,emis[i-1, j-1]))
            if emis_ref < cut_inter * np.max(emis):
                i_emis_ref_loc = None
        else:
            i_emis_ref_loc, j_emis_ref_loc = ij_ref
            emis_ref = emis[i_emis_ref_loc, j_emis_ref_loc]
        if i_emis_ref_loc is not None:
            com_ref = '{} {} {:.1f}'.format(i_emis_ref_loc, j_emis_ref_loc, wls[i_emis_ref_loc-1,j_emis_ref_loc-1])
            if i_emis_ref_loc > 9:
                print('!!! Reference line for {} corresponds to a level {} > 9'.format(atom.atom, i_emis_ref_loc))
                logprint('Ref line level is {}. '.format(i_emis_ref_loc))
                ref_str=str_print.format('9', Z[atom.elem], atom.spec,i_emis_ref_loc/10 , i_emis_ref_loc%10, 0, 0, ion_name, 1.0, 1.0, 0, 0, 0, 0, 999, com_ref)
            else:      
                ref_str=str_print.format('9', Z[atom.elem], atom.spec,i_emis_ref_loc , 0, 0, 0, ion_name, 1.0, 1.0, 0, 0, 0, 0, 999, com_ref)
            print_ref = True 
            
            for i in up_lev:
                print_it = False
                for j in 1+np.arange(i):
                    if emis[i-1,j-1] > cut * emis_max_tot and wls[i-1,j-1] > 912:
                        print_it = True
                        print_any = True
                        if i > NLevels_max:
                            NLevels_max = i
                if print_it:                
                    for j in 1+np.arange(i):
                        if emis[i-1,j-1] > cut * emis_max_tot and wls[i-1,j-1] > 912:
                            com = '{} {}'.format(i, j)
                            ion_name = '{:>2s}_{:<5s}'.format(atom.elem, int_to_roman(atom.spec)).strip()
                            if print_ref:
                                to_print.append(ref_str)
                                print_ref = False
                            i_to_print = np.min((i, 9))
                            if i_emis_ref_loc > 9:
                                to_print.append(str_print.format(' ', Z[atom.elem], atom.spec, i_to_print, 0, i, j, ion_name, 
                                                     wls[i-1,j-1], emis[i-1,j-1]/emis_ref, Z[atom.elem], atom.spec, i_emis_ref_loc/10, i_emis_ref_loc%10, 0, com))
                            else:
                                to_print.append(str_print.format(' ', Z[atom.elem], atom.spec, i_to_print, 0, i, j, ion_name, 
                                                     wls[i-1,j-1], emis[i-1,j-1]/emis_ref, Z[atom.elem], atom.spec, i_emis_ref_loc, 0, 0, com))                                
                            N_lines += 1
    

    if print_any:
        myprint('# {} - Temp. = {} K, Dens. = {} cm-3 \n'.format(atom.atom, tem, den))
        #myprint('# Level           ' + ' '.join(('{:4}'.format(l) for l in 2+np.arange(NLevels_max-1))) + '\n')
        #myprint('# log crit. dens. ' + ' '.join(('{:4.1f}'.format(l) for l in np.log10(atom.getCritDensity(tem))[1:NLevels_max])) + '\n')
        for str_ in to_print:
            myprint(str_)
        logprint('{} lines printed. '.format(N_lines))
    else:
        logprint('No print. ')
        if len(up_lev_list) == 0:
            logprint('No up lev. ')
        if i_emis_ref_loc is None:
            logprint('No ref loc. ')

def get_tem_den(IP=-1):
    if IP < 0:
        return 1e4, 1e3
    elif IP < 13.6:
        return 1e4, 1e3
    elif IP < 24:
        return 1.e4, 1e3
    else:
        return 1.e4, 1e3

def make_all(tem1=None, den1=None, cut=1e-4, filename=None, E_cut=20, cut_inter=1e-5, verbose=False):
    f = open(filename, 'w')
    f.write("# liste_phyat automatically generated on {} \n".format(time.ctime()))
    log_file = open('log.dat', 'w')
    log_file.write("# log_file automatically generated on {} \n".format(time.ctime()))
    atoms = get_atoms_by_conf()
    printed_confs = []
    #atoms = ['O3', 'Fe2']
    for a in atoms:
        print(a)
        conf = gsFromAtom(a)
        log_file.write('{}, conf={}, '.format(a, conf))
        if conf not in printed_confs:
            f.write('# Starting Conf {}\n'.format(conf))
            printed_confs.append(conf)
        try:
            atom = pn.Atom(atom=a, NLevels=50)
            if atom.NLevels > 0:
                do_it = True
            else:
                do_it = False
        except:
            log_file.write('NIST missing')
            do_it = False
        if do_it:
            if tem1 is None:
                tem, den = get_tem_den(atom.IP)
            else:
                tem = tem1
                den = den1
            try:
                print_liste_phyat(atom, tem, den, cut=cut, filename=f, E_cut=E_cut, log_file=log_file, 
                                  cut_inter=cut_inter, verbose=verbose)
                log_file.write('Done.')
            except:
                log_file.write('plp error.')
        log_file.write('\n')
    f.close()
    log_file.close()

def phyat2model(phyat_file, model_file):
    """  
    generate a model_file file from a phyat_file, setting all the master lines to I=1.0
    """
    
    list_phyat = read_data(phyat_file)
    with open(model_file, 'w') as f:
        for line in list_phyat:
            if line['ref'] == 999:
                new_num = line['num']-90000000000000
                f.write('{0:14d} {1[id]:9s}      1.000 0.000 1.000e+04  1.000  0000000000000   1   1.00 {1[comment]:>15s}'.format(new_num, line))

                        
def make_conf_plots():
    
    at_list = [pn.Atom(atom=at) for at in ('C4', 'C3', 'C2', 'O3', 'O2', 'Ne3', 'Ne2')]
    f, axes = plt.subplots(4,2, figsize=(15, 30))
    for at, ax in zip(at_list, axes.ravel()):
        at.plotGrotrian(ax=ax)
    f.savefig('conf_all.pdf')

def get_atoms_by_name():
    pn.atomicData.addAllChianti()
    atoms = pn.atomicData.getAllAtoms()
    atoms.remove('3He2')
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
    for gs in ('s1', 's2', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'd8', 'd7', 'd6', 'd5', 'd4', 'd3', 'd2', 'd1', 'unknown'):        
        for atom in atoms:
            if gss[atom] == gs:
                res.append(atom)
    return res
                
        
"""
To make a new liste_phyat, from an ipython session:

import pyssn
pyssn.make_all(tem1=1e4, den1=1e3, cut=1e-4, filename='liste_phyat_4_3.dat')

"""
