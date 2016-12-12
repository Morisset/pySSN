import numpy as np
import pyneb as pn
import matplotlib.pyplot as plt
from pyneb.utils.physics import Z, IP, gsFromAtom
from pyneb.utils.misc import int_to_roman, parseAtom
import os
import time
from pyssn.utils.misc import split_atom, read_data

pn.atomicData.addAllChianti()


def unique(array, orderby='first'):
    if len(array) == 0:
        return array
    array = np.asarray(array)
    order = array.argsort(kind='mergesort')
    array = array[order]
    diff = array[1:] != array[:-1]
    if orderby == 'first':
        diff = np.concatenate([[True], diff])
    elif orderby == 'last':
        diff = np.concatenate([diff, [True]])
    else:
        raise ValueError
    uniq = array[diff]
    index = order[diff]
    return uniq[index.argsort()]

def print_phyat_list(atom, tem, den, cut=1e-3, cut_inter=1e-5, ij_ref = None, filename=None, up_lev_rule=None,
                      NLevels=None, E_cut=20, log_file=None, help_file= None, verbose=False, Aij_zero = None):
    """
    atom: a string e.g. 'O3' or a pn.Atom object
    tem in K
    den in cm-3
    cut: relative intensity (to the master one) for a line to be printed
    cut_inter [deprecated]: largest dynamic between ref lines.
    ij_ref: None: is computed. Otherwise must be of the form e.g. (4, 2) or ((4, 2), (5, 3))
    filename: where to output the result. May be None, a string or a file object  
    up_lev_rule: 'each': each upper level are used to define a different master line
            'all': all the transitions are grouped into a single master line
            list, e.g. [2, 3, [4, 5]]: levels 2, 3 and 4-5 together are used to define a master line. 
            None: default hard coded value is used, see up_levs below
    NLevels : max number of levels for the atom
    E_cut: max energy in eV for the energy level to give a line
    log_file: for the log
    help_file: for the comments on the Te, Ne for each ion.
    """
    
    up_levs = {'s1': [[2, 3], 4, 5, 6],  
               's2': 'each',
               'p1': [2, [3, 4, 5], [6, 7]], 
               'p2': 'each',
               'p3': [2, 3, [4, 5]],
               'p4': 'each',
               'p5': 'each'} # The others are 'all' by default (i.e. only one master line, with all the lines depending on it

    ij_refs = {'s1': ((3, 1),),
               'p1': ((3, 1), (7, 1)),
               'p3': ((5, 1),)}

    str_print = '{}{:02d}{:02d}{:01d}{:01d}0{:03d}{:03d} {:<8s} {:11.3f} 0.000{:10.3e}  1.000  {:02d}{:02d}{:01d}00{:03d}{:03d}   1   1.00 {}'  
    if filename is None:
        f = None
    else:
        str_print += '\n'
        if type(filename) is str:
            f = open(filename, 'w')
        elif type(filename) is file:
            f = filename

    if filename is None:
        def myprint(s):
            print(s)
    else:
        def myprint(s):
            f.write(s)

    if log_file is None:
        def logprint(s):
            print(s)
    elif type(log_file) is str:
        lf = open(log_file, 'w')
        def logprint(s):
            lf.write(s)
    elif type(log_file) is file:
        lf = log_file            
        def logprint(s):
            lf.write(s)

    if help_file is None:
        def hprint(s):
            print(s)
    elif type(help_file) is str:
        hf = open(log_file, 'w')
        def hprint(s):
            hf.write(s)
    elif type(help_file) is file:
        hf = help_file            
        def hprint(s):
            hf.write(s)
        
    
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
        if N_energies <= 1:
            logprint('Atom without level below {} eV'.format(E_cut))
        else:
            logprint('Atom without energy levels')
        return None
    #print('Doing {} with NLevels={}'.format(atom.atom, this_NLevels))
    atom = pn.Atom(atom=atom.atom, NLevels=this_NLevels)
    if Aij_zero  is not None:
        for ij in Aij_zero:
            atom._A[ij[0]-1, ij[1]-1] = 0.0
    try:
        NIST_gsconf = atom.NIST[0][0].split('.')[-1]
        if NIST_gsconf[-1] in ('0123456789'):
            NIST_gsconf = NIST_gsconf[-2:]
        else:
            NIST_gsconf = NIST_gsconf[-1] + '1'
    except:
        NIST_gsconf = 'unknown'
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
    elif 'split' in up_lev_rule:
        N_split = int(up_lev_rule.split()[-1])
        up_lev_list = [range(2, N_split+1), range(N_split+1, this_NLevels+1)]
    else:
        up_lev_list = up_lev_rule
        
    if ij_ref is not None:
        if type(ij_ref[0]) is not tuple:
            ij_ref = (ij_ref,)
    else:
        if gs in ij_refs:
            ij_ref = ij_refs[gs]
    
    #print('up_lev_rule: {}, up_lev_list: {}'.format(up_lev_rule, up_lev_list))
    NLevels_max = 0
    print_any = False
    to_print = []
    ion_name = '{:>2s}_{:<5s}'.format(atom.elem, int_to_roman(atom.spec)).strip()
    N_lines = 0
    for i_up_lev, up_lev in enumerate(up_lev_list):
        if type(up_lev) is int:
            up_lev = [up_lev]
        emis_ref = 0
        i_emis_ref_loc = None
        j_emis_ref_loc = None
        
        for i in up_lev:
            if i < this_NLevels+1:
                for j in 1+np.arange(i):
                    if ij_ref is not None and (i, j) in ij_ref:
                            i_emis_ref_loc = i
                            j_emis_ref_loc = j
                            emis_ref = emis[i-1,j-1]
        if i_emis_ref_loc is None:
            for i in up_lev:
                if i < this_NLevels+1:
                    for j in 1+np.arange(i):
                        if emis[i-1,j-1] > emis_ref:
                            i_emis_ref_loc = i
                            j_emis_ref_loc = j
                            emis_ref = emis[i-1,j-1]
        """
        if emis_ref < cut_inter * np.max(emis):
            if verbose:
                print('reset emis ref loc {} < {} * {} (emis_ref < cut_inter * np.max(emis))'.format(emis_ref, cut_inter, np.max(emis)))
            i_emis_ref_loc = None
        """
        if i_emis_ref_loc is not None:
            com_ref = '{} {} {:.1f}'.format(i_emis_ref_loc, j_emis_ref_loc, wls[i_emis_ref_loc-1,j_emis_ref_loc-1])
            if ("split" in up_lev_rule) or ("all" in up_lev_rule):
                i_emis_ref_loc = i_up_lev + 1
            if i_emis_ref_loc > 9:
                logprint('Ref line level is {}. '.format(i_emis_ref_loc))
                ref_str=str_print.format('9', Z[atom.elem], atom.spec,i_emis_ref_loc/10 , i_emis_ref_loc%10, 0, 0, ion_name, 1.0, 1.0, 0, 0, 0, 0, 999, com_ref)
            else:      
                ref_str=str_print.format('9', Z[atom.elem], atom.spec,i_emis_ref_loc , 0, 0, 0, ion_name, 1.0, 1.0, 0, 0, 0, 0, 999, com_ref)
            print_ref = True 
            
            for i in up_lev:
                print_it = False
                if i > this_NLevels:
                    break
                for j in 1+np.arange(i):
                    if emis[i-1,j-1] > cut * emis_ref and wls[i-1,j-1] > 912:
                        print_it = True
                        print_any = True
                        if i > NLevels_max:
                            NLevels_max = i
                if print_it:                
                    for j in 1+np.arange(i):
                        if emis[i-1,j-1] > cut * emis_ref and wls[i-1,j-1] > 912:
                            com = '{} {}'.format(i, j)
                            ion_name = '{:>2s}_{:<5s}'.format(atom.elem, int_to_roman(atom.spec)).strip()
                            if print_ref:
                                to_print.append(ref_str)
                                print_ref = False
                            i_to_print = np.min((i, 9))
                            i_to_print = i_emis_ref_loc
                            if i_emis_ref_loc > 9:
                                to_print.append(str_print.format(' ', Z[atom.elem], atom.spec, i_to_print, 0, i, j, ion_name, 
                                                     wls[i-1,j-1], emis[i-1,j-1]/emis_ref, Z[atom.elem], atom.spec, i_emis_ref_loc/10, i_emis_ref_loc%10, 0, com))
                            else:
                                to_print.append(str_print.format(' ', Z[atom.elem], atom.spec, i_to_print, 0, i, j, ion_name, 
                                                     wls[i-1,j-1], emis[i-1,j-1]/emis_ref, Z[atom.elem], atom.spec, i_emis_ref_loc, 0, 0, com))                                
                            N_lines += 1
    

    if print_any:
        hprint('# {} - Temp. = {} K, Dens. = {} cm-3 \n'.format(atom.atom, tem, den))
        for str_ in to_print:
            myprint(str_)
        logprint('{} lines printed. '.format(N_lines))
    else:
        logprint('No print. ')
        if len(up_lev_list) == 0:
            logprint('No up lev. ')
        if i_emis_ref_loc is None:
            logprint('No ref loc. ')
    
    if type(help_file) is str:
        hf.close()
    if type(log_file) is str:
        lf.close()

def make_phyat_list(filename, tem1=None, den1=None, cut=1e-4, E_cut=20, cut_inter=1e-5, 
             verbose=False, notry=False, NLevels=50, atoms=None, 
             ref_lines_dic=None, NLevels_dic=None, up_lev_rule_dic=None, Aij_zero_dic=None,
             tem_den_dic = None, extra_file=None):
    
    """
    filename: output file name
    tem1: temperature in K
    den1: density in cm-3
    cut: relative intensity (to the master one) for a line to be printed
    E_cut: max energy in eV for the energy level to give a line
    cut_inter: largest dynamic between ref lines.
    NLevels: max number of levels for the atom
    verbose: set to True of more info needed
    notry: set to True to force calling print_phyat_list outside of a try (debig mode)
    atoms: a list of strings, e.g. ['C4', 'C3', 'C2', 'O3', 'O2', 'Ne3', 'Ne2', 'Fe3']. If None, all atoms are done
    ref_lines_dic: e.g. {'Fe7': ((4, 2),),
                         'Fe6': ((2, 1), (5, 1),),
                         'Fe5': ((3, 2), (7, 5),),
                         'Fe4': ((6, 1),),
                         'Fe3': ((2, 1), (12, 1),)
                         }
    NLevels_dic: e.g. {'S1': 8,
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
    up_lev_rule_dic: e.g. {'Ni2': 'split 3',
                           'Fe1': 'split 5',
                           'Ni3': 'split 3',
                           'Ni4': 'split 4',
                           'Fe3': 'split 5',
                           'Fe4': 'all',
                           'Fe5': 'split 5',
                           'Fe6': 'split 4',
                           'Fe7': 'split 3'
                           }
    Aij_zero_dic: e.g. {'C3': ((2,1),)
                        }
    tem_den_dic: e.g. {0.:   (1e4, 1e3),
                       13.6: (1e4, 1e3),
                       24.0: (1e4, 1e3),
                       1e6:  (1e4, 1e3)
                       }
    extra_file: a file of pySSN data format containing data to include.
    """
    
    
    if ref_lines_dic is None:
        ref_lines_dic = {'Fe7': ((4, 2),),
                         'Fe6': ((2, 1), (5, 1),),
                         'Fe5': ((3, 2), (7, 5),),
                         'Fe4': ((6, 1),),
                         'Fe3': ((2, 1), (12, 1),)
                         }
    if NLevels_dic is None:
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
        
    if up_lev_rule_dic is None:
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
        
    if Aij_zero_dic is None:
        Aij_zero_dic = {'C3': ((2,1),)
                        }
    
    if tem_den_dic is None:
        tem_den_dic = {0.:   (1e4, 1e3),
                       13.6: (1e4, 1e3),
                       24.0: (1e4, 1e3),
                       1e6:  (1e4, 1e3)
                       }
    
    def get_tem_den(IP):
        for k in sorted(tem_den_dic.keys()):
            if IP < k: 
                return tem_den_dic[k]
    
    f = open(filename, 'w')
    f.write("# liste_phyat automatically generated on {} \n".format(time.ctime()))
    log_file = open('log.dat', 'w')
    log_file.write("# log_file automatically generated on {} \n".format(time.ctime()))
    help_file = open('help.dat', 'w')
    help_file.write("# help_file automatically generated on {} \n".format(time.ctime()))

    # The following is very badly programmed, need to make it better...
    if extra_file is not None:
        extra_atoms = np.array(get_extra_atoms(extra_file=extra_file))
        with open(extra_file, 'r') as fextra:
            extra_data = fextra.readlines()
        extra_data = np.array(extra_data)
    else:
        extra_data = []
        extra_atoms = []
    
    if atoms is None:
        atoms = get_atoms_by_conf(extra_file=extra_file)
    else:
        atoms.extend(get_extra_atoms(extra_file, uniq=True))
        atoms = get_atoms_by_conf(atoms=atoms)
        
    atoms = unique(atoms)
    
    printed_confs = []
    #atoms = []
    for a in atoms:
        print(a)
        if a in extra_atoms:
            n_lines = 0
            for line in extra_data[extra_atoms == a]:
                if line[0] != '9':
                    n_lines += 1
                f.write(line)
            log_file.write('{}, {} lines from file {}\n'.format(a, n_lines, extra_file))
        else:
            if a in ref_lines_dic:
                ref_lines = ref_lines_dic[a]
            else:
                ref_lines = None
            if a in NLevels_dic:
                this_NLevels = NLevels_dic[a]
            else:
                this_NLevels = NLevels
            if a in up_lev_rule_dic:
                up_lev_rule = up_lev_rule_dic[a]
            else:
                up_lev_rule = None
            conf = gsFromAtom(a)
            log_file.write('{}, conf={}, '.format(a, conf))
            if conf not in printed_confs:
                printed_confs.append(conf)
            try:
                atom = pn.Atom(atom=a, NLevels=this_NLevels)
                if atom.NLevels > 0:
                    do_it = True
                else:
                    do_it = False
            except:
                log_file.write('NIST missing. \n')
                do_it = False
            if do_it:
                if tem1 is None:
                    tem, den = get_tem_den(atom.IP_up)
                else:
                    tem = tem1
                    den = den1
                if a in Aij_zero_dic:
                    Aij_zero = Aij_zero_dic[a]
                else:
                    Aij_zero = None
                if notry:
                    print_phyat_list(atom, tem, den, cut=cut, filename=f, E_cut=E_cut, log_file=log_file, 
                                      cut_inter=cut_inter, verbose=verbose, help_file=help_file, ij_ref=ref_lines,
                                      up_lev_rule=up_lev_rule, Aij_zero=Aij_zero)                
                else:
                    try:
                        print_phyat_list(atom, tem, den, cut=cut, filename=f, E_cut=E_cut, log_file=log_file, 
                                      cut_inter=cut_inter, verbose=verbose, help_file=help_file, ij_ref=ref_lines,
                                      up_lev_rule=up_lev_rule, Aij_zero=Aij_zero)    
                    except:
                        log_file.write('plp error.')
                log_file.write('\n')
                
    f.write('\n')
    f.close()
    log_file.write('\n')
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

def remove_iso(atom):
    res = ''
    elem_done = False
    for char in atom:
        if not elem_done:
            if char not in '1234567890':
                res += char
                elem_done = True
        else:
            res += char
    return res

def get_atom_str(atom):
    """
    return atom with leading 0.
    eg: get_atom_str('Mg5') -> Mg005
    """
    ato, ion = parseAtom(atom)
    return('{}{:03d}'.format(remove_iso(ato), int(ion)))

def get_extra_atoms(extra_file=None,uniq=False):
    
    if extra_file is None:
        return []
    extra_data = read_data(extra_file)
    atoms=[]
    for ID in extra_data.id:
        IDs = split_atom(ID)
        if IDs[0] is None:
            atoms.append('{0[1]}{0[2]}'.format(IDs))
        else:
            atoms.append('{0[0]}{0[1]}{0[2]}'.format(IDs))
    if uniq:
        res = unique(atoms)
    else:
        res = atoms
    return res

def get_atoms_by_name(extra_file=None):
    """
    'phyat_list_DP_01.dat'
    """
    atoms = pn.atomicData.getAllAtoms()
    atoms.extend(get_extra_atoms(extra_file, uniq=True))
    atoms.remove('3He2')
    return sorted(atoms, key=get_atom_str)

def get_atoms_by_Z(extra_file=None, atoms=None):
    if atoms is None:
        atoms = get_atoms_by_name(extra_file=extra_file)
    return sorted(atoms, key=lambda k:Z[parseAtom(remove_iso(k))[0]])

ZmI = lambda atom: Z[parseAtom(remove_iso(atom))[0]] - int(parseAtom(remove_iso(atom))[1])

def get_atoms_by_ZmI(extra_file=None, atoms=None):
    atoms = get_atoms_by_Z(extra_file=extra_file, atoms=atoms)
    return sorted(atoms, key=ZmI)
    
def get_atoms_by_conf(extra_file=None, atoms=None):
    if atoms is None:
        atoms = get_atoms_by_ZmI(extra_file=extra_file)
    res = []
    gss = {}
    for atom in atoms:
        gss[remove_iso(atom)] = gsFromAtom(remove_iso(atom))
    for gs in ('s1', 's2', 
               'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 
               'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10', 
               'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'f11', 'f12', 'f13', 'f14', 
               'unknown'):        
        for atom in atoms:
            if gss[remove_iso(atom)] == gs:
                res.append(atom)
    return res
                
        
"""
To make a new liste_phyat, from an ipython session:

import pyssn
pyssn.make_phyat_list('liste_phyat_4_3.dat', tem1=1e4, den1=1e3, cut=1e-4)

"""
