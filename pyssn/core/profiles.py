'''
Created on 28 janv. 2014

@author: christophemorisset
'''

import numpy as np
from pyssn.utils.physics import CST
from pyssn.utils.misc import convolgauss


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
    

def profil_instr(filter_size, prof, lambda_pix):
    
    w_norm = np.arange(filter_size)-filter_size/2
    w_norm_abs = np.abs(w_norm)
    profil = np.zeros(filter_size)
    prof_add = lambda dec, alpha, Bl, Br: (Bl*winf + Br*wsup)*np.exp(-(w_norm_abs/dec*lambda_pix)**alpha)
    
    if prof['largeur'] > 0:
        profil = np.exp(-(w_norm/prof['largeur']*lambda_pix)**2)
    else:
        profil[filter_size/2] = 1.0
        profil[w_norm_abs <= abs(prof['largeur']/lambda_pix)] = 1.0
    
    winf = np.zeros(filter_size)
    wsup = np.zeros(filter_size)
    
    winf[w_norm <= 0.0] = 1.0
    wsup[w_norm > 0.0] = 1.0

    profil += prof_add(prof['decroiss_1'], prof['alpha_1'], prof['B_1l'], prof['B_1r'])
    profil += prof_add(prof['decroiss_2'], prof['alpha_2'], prof['B_2l'], prof['B_2r'])
    profil += prof_add(prof['decroiss_3'], prof['alpha_3'], prof['B_3l'], prof['B_3r'])
    profil += prof_add(prof['decroiss_4'], prof['alpha_4'], prof['B_4l'], prof['B_4r'])
    
    return profil
