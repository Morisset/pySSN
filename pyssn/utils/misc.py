'''
Created on 16/01/2014

@author: morisset
'''
import os
import sys
import re
import argparse
from scipy import interpolate
import numpy as np
import pyssn
from pyneb.utils.physics import vactoair

def execution_path(filename):
    return os.path.join(os.path.dirname(sys._getframe(1).f_code.co_filename), filename)

def change_size(tab, fact):
    if fact == 1:
        return tab.copy()
    size = len(tab)
    old_tab_pix = np.linspace(0, size-1, size)
    new_tab_pix = np.linspace(0, size-1, size*fact)
    interp = interpolate.interp1d(old_tab_pix, tab)
    tab_out = interp(new_tab_pix)
    return tab_out

def rebin(tab, fact):
    if fact == 1:
        return tab.copy()
    if len(tab) % fact != 0:
        pyssn.log_.error('Dimension of tab ({0}) is not a multiple of fact ({1})'.format(len(tab), fact), 
                         calling = 'pyssn.misc.rebin')
        return None
    return tab.reshape(len(tab)/fact, fact).sum(1) / fact
    
def convol(array, kernel, method='numpy'):
    
    if method == 'same':
        return array.copy()
    elif method == 'numpy':
        if len(array) >= len(kernel):
            result = np.convolve(array, kernel, mode='same')
        else: 
            result = np.convolve(array, kernel, mode='valid')
        return result
    elif method == 'local':
        return None
    else:
        return None

def convolgauss(spectrum, w, lambda_0, fwhm):
    """
    Convolution with a Gaussian
    """
    
    pix_0 = np.argmin(np.abs(w - lambda_0))
    if pix_0 == 0:
        pix_0 = 1
    if pix_0 == len(w)-1:
        pix_0 = len(w)-2
        
    lam_pix = np.abs(w[pix_0-1] - w[pix_0+1]) / 2.
    fwhm_pix = fwhm / lam_pix
    sig = fwhm_pix / 2.35482 # sqrt(2.*alog(2.))
    nres = 21
    while True:
        nres = nres * 2 + 1
        wkernel = np.arange(nres) - nres / 2
        kernel = np.exp(-(wkernel / (np.sqrt(2.) * sig))**2)
        if (np.min(kernel) < 1e-9) or ((nres * 2 - 3) > len(w)):
            break
    kernel = kernel/kernel.sum()
    cspectrum = convol(spectrum,kernel)
    
    return cspectrum
    
def is_absorb(raie):
    """
    if raie is a numpy recarray, it reports a table of indexes where the condition is filled
    """
    i_min = 9500000000000
    i_max = 9600000000000
    if type(raie) == np.core.records.recarray:
        index_abs = np.where((raie['num'] >= i_min) &
                             (raie['num'] < i_max) &
                             (raie['i_rel'] < 0.0))[0]
        return index_abs
    if (raie['num'] >= i_min) & (raie['num'] < i_max) and (raie['i_rel'] < 0.):
        return True
    else:
        return False

def no_red_corr(raie):
    """
    if raie is a numpy recarray, it reports a table of indexes where the condition is filled
    """
    
    i_min = 9500000000000
    i_max = 9900000000000
    if type(raie) is np.core.records.recarray:
        index_abs = np.where((raie['num'] >= i_min) &
                             (raie['num'] < i_max))[0]
        return index_abs
    if (raie['num'] >= i_min) & (raie['num'] < i_max):
        return True
    else:
        return False


def read_HITRAN_stick(file_):
    """
#MI WN,cm-1     S,cm/mol  A,s-1    Vair Vself El,cm-1  TexpAshift  GQNup          GQNlow         LQNup          LQNlow         Ierr  Iref         SWup   SWlow  

 7157027.590000 3.556e-20 3.147e+03.04000.038    0.00000.000.000000       B     19       X      0                R  1R  0     d22020032 7 0 2 0 0     6.0    1.0
    
    """
    data = np.genfromtxt(file_, 
                         dtype=['a2', 'a1','float', 'float','float', 'float','float'], 
                         delimiter=[2, 1, 12, 10, 10, 5, 5], 
                         skip_header=5, 
                         usecols=(0, 1, 2, 3, 4, 5, 6), 
                         names=('ID', 'ISO','nu','I','A', 'V', 'Vself'))
    return data

def vactoair(wl):
    mask = ( wl >= 2000. ) & ( wl < 20000. )
    fact = 1. + 2.735182e-4 + 131.4182/wl**2 + 2.76249e8/(wl**4.)
    return wl / (mask * fact + (1 - mask))

def make_abs_from_hitran(file_in, file_out, ID_ref, ID_start, label, fac_tau, cut=1e-4, wl_min=1200, wl_max=14000):
    """
    generate a liste_phyat formated file from an hitran data file.
    USAGE:
        make_abs_from_hitran('O2_296_1.stick.gz', 'liste_O2.dat', 9510000000000, 9510000000001, 'At_O2', 1e23, cut=0.0001, wl_min=3200, wl_max=18000)
    """

    data = read_HITRAN_stick(file_in)
    wls = vactoair(1e8 / data['nu'])
    taus = fac_tau *data['I']
    tau_max = 0.
    ID = ID_start
    count=0
    f = open(file_out, 'w')
    f.write('{:14d} {:8s}         1.0 0.000 1.000e+00  1.000            999   1   1.00\n'.format(ID_ref+90000000000000, label))
    for wl, tau in zip(wls, taus):
        if (tau > cut) & (wl > wl_min ) & (wl < wl_max) :
            f.write('{:14d} {:8s}  {:10.3f} 0.000 {:9.3e}  1.000  {:13d}   1   1.00\n'.format(ID, label, wl, tau, ID_ref))
            count += 1
            if tau > tau_max:
                tau_max = tau
        ID += 1
    f.close()
    print('Wrote {} lines in file {} with max(tau) = {}.'.format(count, file_out, tau_max))
      
def gauss(w, I, w_shift, width): 
    return I * np.exp(-((w + w_shift)/(width))**2)

def lorentz(w, I, w_shift, width): 
    return I / (1. + ((w + w_shift)/width)**2)

def carre(w, I, w_shift, width): 
    prof = np.zeros_like(w)
    prof[np.abs(w + w_shift) < (width)] = I
    return prof
    
def rebin_tapas(file_ = None, cut=1e-5):
    
    d = np.genfromtxt(file_, skip_header=26, dtype=None)
    grad = np.ones_like(d[:,0])
    grad[1::] = np.abs(d[1::,1] - d[0:-1,1])
    mask = grad > cut
    res = d[mask,:]
    fout = file_.split('.')[0] + '_rebin.' + file_.split('.')[1]
    with open(fout, 'w') as f:
        for r in res:
            f.write('{} {}\n'.format(r[0]*10., r[1]))

def clean_label(str_):
    return(re.sub("_", " ", str_))


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="init file .")
    parser.add_argument("-V", "--version", action="version", version=pyssn.__version__,
                        help="Display version information and exit.")
    return parser

    