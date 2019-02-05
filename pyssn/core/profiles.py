'''
Created on 28 janv. 2014

@author: christophemorisset
'''

import numpy as np
from ..utils.physics import CST
from ..utils.misc import convolgauss


def profil_emis(w, raie, lambda_shift=0.):

    """
    raie = {'lambda' : 6200,
            'l_shift': 0.,
            'vitesse' : 5.,
            'profile' : 1,
            'num' : 101000001010}
    """

    lambda_0 = raie['lambda'] + raie['l_shift'] + lambda_shift
    w_norm = w - lambda_0
    
    profil = np.zeros_like(w)
    largeur = raie['vitesse'] * lambda_0 / CST.CLIGHT * 1e5
    sp_type =  raie['profile']
    charge = (raie['num'] % 100000000000)/1000000000
    masse =  2 * (raie['num'] - raie['num'] % 100000000000)/100000000000 
    if (masse == 2 and (raie['num'] - raie['num'] % 101000000000)/100000000 == 1010) : 
        masse = 1
    gauss = lambda I, w_shift, width: I * np.exp(-((w_norm + w_shift*largeur)/(width * largeur))**2)

    if sp_type == 1:
        profil = gauss(1., 0., 1.)  
        T4 = 1.
    else:
        profil = gauss(1., 0., 1.)  
        T4 = 1.
    if T4 > 0.0:
        fwhm_therm = 21.4721 * np.sqrt(T4 / masse) * lambda_0 / 299792.50 #km/s
        profil = convolgauss(profil, w, lambda_0, fwhm_therm)
    
    profil[~np.isfinite(profil)] = 0.0

    
    return profil

def translate_intr_prof(prof):    
    
    keys = prof.keys()
    old_keys = [ 'largeur', 'B_1l', 'B_2l', 'B_3l', 'B_4l', 'B_1r', 'B_2r', 'B_3r', 'B_4r', 'decroiss_1', 'decroiss_2', 'decroiss_3', 'decroiss_4']
    new_keys = [ 'width', 'Bb_1', 'Bb_2', 'Bb_3', 'Bb_4', 'Br_1', 'Br_2', 'Br_3', 'Br_4', 'beta_1', 'beta_2', 'beta_3', 'beta_4' ]
    key_dict = dict(zip(old_keys, new_keys))
    for key in old_keys:
        if key in keys:
            prof[key_dict[key]] = prof[key]
    return prof

def profil_instr(filter_size, prof, lambda_pix):
    
    prof = translate_intr_prof(prof)
    w_norm = np.arange(filter_size)-filter_size/2
    w_norm_abs = np.abs(w_norm)
    profil = np.zeros(filter_size)
    prof_add = lambda dec, alpha, Bl, Br: (Bl*winf + Br*wsup)*np.exp(-(w_norm_abs/dec*lambda_pix)**alpha)

    if prof['width'] > 0:
        profil = np.exp(-(w_norm/prof['width']*lambda_pix)**2)
    else:
        profil[filter_size/2] = 1.0
        profil[w_norm_abs <= abs(prof['width']/lambda_pix)] = 1.0
    
    winf = np.zeros(filter_size)
    wsup = np.zeros(filter_size)
    
    winf[w_norm <= 0.0] = 1.0
    wsup[w_norm > 0.0] = 1.0

    indexes = sorted([indexed_key.replace('alpha','') for indexed_key in prof.keys() if 'alpha' in indexed_key])
    for i in indexes:
        profil += prof_add(prof['beta'+i], prof['alpha'+i], prof['Bb'+i], prof['Br'+i])
    return profil
