'''
Created on 17/01/2014

@author: morisset
'''
import numpy as np
import pickle
import pyssn
from pyssn.utils.misc import execution_path
from scipy.interpolate import interp1d

class CST(object):
    BOLTZMANN = 1.3806488e-16 # erg/K - NIST 2010
    CLIGHT = 2.99792458e10 # cm/s - NIST 2010
#    HPLANCK = 6.62606957e-27 # erg s - NIST 2010
    HPLANCK = 6.62607004e-27 # erg s - NIST 2014
    EMASS = 9.10938291e-28 # g - NIST 2010
    ECHARGE = 1.602176565e-19 # Electron charge in Coulomb - NIST 2010
    PI = 3.141592653589793238462643
    BOLTZMANN_ANGK = (BOLTZMANN) / (HPLANCK * CLIGHT * 1.e8) # Boltzmann constant in (Ang * Kelvin) ** -1
    RYD = 109737.31568508 # Rydberg constant in cm^-1 - NIST 2014     
    RYD_EV = HPLANCK * CLIGHT * RYD * 1.e-7 / ECHARGE # infinite mass Rydberg in eV
    RYD_ANG = 1.e8 / RYD # infinite mass Rydberg in A
    RYD_erg =  2.179872325e-11 # HPLANCK * CLIGHT * RYD
    KCOLLRATE = pow((2 * PI / BOLTZMANN), 0.5) * pow(HPLANCK / (2 * PI), 2) / pow(EMASS, 1.5) # constant of collisional rate equation
    BOLTZMANN_eVK = 8.617343e-5 # Boltzmann constant in eV/K
    ERG_S_EV = 1.6021765e-12
    
def Planck(w, T):
    
    """
    w in Angstrom
    T in Kelvin
    
    Gives the blackbody flux (i.e. PI*Intensity) ergs/cm2/s/A
    """
    w_cm = w / 1e8
    C1 = 3.7417749e-13 # 2 * np.pi * CST.HPLANCK * CST.CLIGHT**2 * 1e-8
    C2 = 1.4387687  # CST.HPLANCK *  CST.CLIGHT / CST.BOLTZMANN
    
    
    res = C1 / (w_cm**5.0 * (np.exp(C2 / w_cm / T) - 1.))
    res[~np.isfinite(res)] = 0.
    
    return res

def make_cont_Ercolano(T_in, case, lam):
    """
    Adapted from http://adsabs.harvard.edu/abs/2006MNRAS.372.1875E
    """
    n_lam = len(lam)
    hnu =  CST.CLIGHT * 1e8 / lam * CST.HPLANCK  #!phy.c_ang_sec/lam*!phy.h
    with open(execution_path('../data/coeff_ercolano.pickle'), 'rb') as handle:
        BE = pickle.load(handle)
    
    if case == 'H':
        tab_T = 10**BE['th']
        D = BE['dh']
    elif case == 'He1':
        tab_T = 10**BE['the1']
        D = BE['dhe1']
    elif case == 'He2':
        tab_T = 10**BE['the2']
        D = BE['dhe2']
    else:
        pyssn.log_.error('Invalid case {0}'.format(case), calling='make_cont_Ercolano')
        return None
    if (T_in < np.min(tab_T)) or (T_in > np.max(tab_T)):
        pyssn.log_.error('Invalid temperature {0}'.format(T_in), calling='make_cont_Ercolano')
        return None
    
    BE_E_Ry = D[:,1]
    BE_E_erg = BE_E_Ry * CST.RYD_erg
    BE_E_Thr = BE_E_erg[D[:,0] == 1]
    Delta_E = np.zeros(n_lam)
    for i in np.arange(n_lam):
        DE = hnu[i] - BE_E_Thr
        Delta_E[i] = np.min(DE[DE > 0])
        
    n_T_sup = np.min(np.where(tab_T >= T_in)[0])
    n_T_inf = n_T_sup - 1
    T_sup = tab_T[n_T_sup]
    T_inf = tab_T[n_T_inf]
    
    BE_coeff_sup = D[:, n_T_sup+2]
    BE_coeff_inf = D[:, n_T_inf+2]
    
    coeff_sup = interp1d(BE_E_erg, BE_coeff_sup)(hnu)
    coeff_inf = interp1d(BE_E_erg, BE_coeff_inf)(hnu)

    C_interp= (np.log10(T_in) - np.log10(T_inf)) / (np.log10(T_sup) - np.log10(T_inf))
    
    coeff = coeff_sup * C_interp + coeff_inf*(1. - C_interp)
    
    cont = coeff * 1e-34 * T_in**(-1.5) * np.exp(-Delta_E / T_in / CST.BOLTZMANN) / lam**2. * CST.CLIGHT * 1e8 # erg/s.cm3/A
    return cont


def gff(Z, T, lam):
    """
     Adaptated from http://adsabs.harvard.edu/abs/1991CoPhC..66..129S
    """
    D= np.array([8.986940175e+00, -4.009515855e+00,  8.808871266e-01,
        2.640245111e-02, -4.580645915e-02, -3.568055702e-03,   
        2.827798067e-03,  3.365860195e-04, -8.006936989e-01,
        9.466021705e-01,  9.043402532e-02, -9.608451450e-02,
        -1.885629865e-02,  1.050313890e-02,  2.800889961e-03, 
        -1.078209202e-03, -3.781305103e-01,  1.102726332e-01, 
        -1.543619180e-02,  8.310561114e-03,  2.179620525e-02, 
        4.259726289e-03, -4.181588794e-03, -1.770208330e-03,   
        1.877213132e-02, -1.004885705e-01, -5.483366378e-02,   
        -4.520154409e-03,  8.366530426e-03,  3.700273930e-03, 
        6.889320423e-04,  9.460313195e-05,  7.300158392e-02,   
        3.576785497e-03, -4.545307025e-03, -1.017965604e-02,   
        -9.530211924e-03, -3.450186162e-03,  1.040482914e-03, 
        1.407073544e-03, -1.744671550e-03,  2.864013856e-02,   
        1.903394837e-02,  7.091074494e-03, -9.668371391e-04,   
        -2.999107465e-03, -1.820642230e-03, -3.874082085e-04, 
        -1.707268366e-02, -4.694254776e-03,  1.311691517e-03, 
        5.316703136e-03,  5.178193095e-03,  2.451228935e-03,   
        -2.277321615e-05, -8.182359057e-04,  2.567331664e-04, 
        -9.155339970e-03, -6.997479192e-03, -3.571518641e-03, 
        -2.096101038e-04,  1.553822487e-03,  1.509584686e-03, 
        6.212627837e-04,  4.098322531e-03,  1.635218463e-03,   
        -5.918883504e-04, -2.333091048e-03, -2.484138313e-03, 
        -1.359996060e-03, -5.371426147e-05,  5.553549563e-04, 
        3.837562402e-05,  2.938325230e-03,  2.393747064e-03,   
        1.328839809e-03,  9.135013312e-05, -7.137252303e-04,   
        -7.656848158e-04, -3.504683798e-04, -8.491991820e-04, 
        -3.615327726e-04,  3.148015257e-04,  8.909207650e-04, 
        9.869737522e-04,  6.134671184e-04,  1.068883394e-04,   
        -2.046080100e-04 ])



    XLF = np.log10(CST.CLIGHT * 1e8 / lam)
    N_lam = len(lam)
    G = np.zeros(N_lam)
    D = D.reshape(11, 8)
    B = np.zeros(11)
    C = np.zeros(8)

    XLRKT = 5.1983649 - np.log10(T)
    TXG = 0.66666667 * (2.0 * np.log10(Z) + XLRKT)

    for j in np.arange(7):
        B[10] = D[10, j]
        B[9] = TXG * B[10] + D[9, j]
        for IR in np.arange(8)[::-1]:
            B[IR] = TXG * B[IR+1] - B[IR+2] + D[IR, j]
        C[j] = 0.25 * (B[0] - B[2])               


    CON=0.72727273 * XLRKT - 10.376127  

    for i in np.arange(N_lam):
        TXU = 0.72727273 * XLF[i] + CON 
        B[7] = C[7]
        B[6] = TXU * B[7] + C[6] 
        for IR in np.arange(5)[::-1]: 
            B[IR] = TXU * B[IR+1] - B[IR+2] + C[IR]
        G[i] = B[0] - B[2]

    return G

def make_ion_frac(N, co=None):

    import pandas as pd
    import pymysql
    from pyneb.utils.physics import IP, sym2name, Z
    import numpy as np
    
    sym2name['S'] = 'sulphur'
    if co is None:
        co = pymysql.connect(host='132.248.1.102', db='3MdB', user='OVN_user', passwd='oiii5007') 
        close_after=True
    else:
        close_after=False
    res = pd.read_sql("select * from abion where N = {}".format(N), con=co)
    if close_after:
        co.close()
    
    elems = np.unique([s.split('_')[1].lower() for s in res.keys()[2::]])
    f = open('{}_ionfrac.dat'.format(N), 'w')
    
    for sym in sorted(sym2name):
        if sym[0] not in '0123456789':
            if sym2name[sym] in elems:
                for i in range(30):
                    ion = '{}{}'.format(sym, i)
                    ion_frac_str = 'A_{}_vol_{}'.format(sym2name[sym].upper(), i)
                    if ion_frac_str in res.keys():
                        f.write('{} {}\n'.format(ion, res[ion_frac_str][0]))
    f.close()
            
            
            
            
            
            
            
    