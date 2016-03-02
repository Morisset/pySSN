"""
pySSN is available under the GNU licence providing you cite the developpers names:

    Ch. Morisset (Instituto de Astronomia, Universidad Nacional Autonoma de Mexico)

    D. Pequignot (Meudon Observatory, France)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from scipy import interpolate
import pyneb as pn

from pyssn import log_, config
from ..utils.physics import CST, Planck, make_cont_Ercolano, gff
from ..utils.misc import execution_path, change_size, convol, rebin, is_absorb, no_red_corr, gauss, carre, lorentz, convolgauss, vactoair, clean_label
from ..core.profiles import profil_instr
from ..utils.misc import get_parser

"""
ToDo:
1) Define the special lines in a table or dictionnary where ref, color, style are set up.
2) avoid reduced tick values for lambdas when zooming.
"""

def read_data(filename, NF=True):
    dtype = 'i8, a1, a9, f, f, f, f, a1, i8, i4, f, a25'
    if NF:
        delimiter = [14, 1, 9, 11, 6, 10, 7, 1, 14, 4, 7, 25]
    else:
        delimiter = [ 9, 1, 9, 11, 6, 10, 7, 1,  9, 4, 7, 25]
    names = ['num', 'foo', 'id', 'lambda','l_shift', 'i_rel', 'i_cor', 'foo2', 'ref', 'profile', 
             'vitesse', 'comment']
    usecols = (0, 2, 3, 4, 5, 6, 8, 9, 10, 11)
    dd = np.genfromtxt(filename, dtype=dtype, delimiter=delimiter, names = names, usecols = usecols)

    if np.isnan(dd['num']).sum() > 0:
        log_.error('Some line ID are not defined {}'.format(dd['id'][np.isnan(dd['num'])]))
    if np.isnan(dd['lambda']).sum() > 0:
        log_.error('Some wavelengths are not defined {}'.format(dd['num'][np.isnan(dd['lambda'])]))
    if np.isnan(dd['l_shift']).sum() > 0:
        log_.error('Some wavelengths shifts are not defined {}'.format(dd['num'][np.isnan(dd['l_shift'])]))
    if np.isnan(dd['i_cor']).sum() > 0:
        log_.error('Some intensity corrections are not defined {}'.format(dd['num'][np.isnan(dd['i_cor'])]))
    if np.isnan(dd['i_rel']).sum() > 0:
        log_.error('Some relative intensities are not defined {}'.format(dd['num'][np.isnan(dd['i_rel'])]))
    
    return dd.view(np.recarray)

class spectrum(object):
    
    def __init__(self, config_file=None, phyat_file=None, profil_instr=profil_instr, 
                 do_synth = True, do_read_liste = True, do_cosmetik = True, do_run = True, limit_sp = None,
                 spectr_obs=None, sp_norm=None, obj_velo=None):
        """
        Main pySSN object.
        It reads the configuration file given by the config_file parameter. 
        It reads the atomic data, model and cosmetik files. It reads the observation. It computes the reddening
        correction, the 
        """
        
        self.profil_instr = profil_instr
        self.do_cosmetik = do_cosmetik        
        self.calling = 'spectrum'
        self.full_config_file = config_file
        if '/' in self.full_config_file:
            file_name = self.full_config_file.split('/')[-1]
            dir = self.full_config_file.split(file_name)[0]
            if dir == '':
                dir = './'
            self.directory = dir
            self.config_file = file_name
        else:
            self.directory = './'
            self.config_file = self.full_config_file
        config.addDataFilePath(self.directory, inpySSN=False)
        
        self.init_vars()
        
        self.read_conf(self.config_file)
        self.init_obs(spectr_obs=spectr_obs, sp_norm=sp_norm, obj_velo=obj_velo, limit_sp=limit_sp)
        self.init_red_corr()
        self.make_continuum()
            
        if phyat_file is not None:
            self.phyat_file = phyat_file
        else:
            self.phyat_file = self.get_conf('phyat_file', 'liste_phyat.dat')
        if do_run:
            self.run(do_synth = do_synth, do_read_liste = do_read_liste)

    def init_vars(self):
        self.fig1 = None
        self.fig2 = None
        self.fig3 = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.cursor_width = 0.02
        self.aire_ref = 1.0
        self.zoom_fact = 0.1
        self._cid = None
        self.plot_magenta = None
        self.plot_cyan = None
        self.label_magenta = None
        self.label_cyan = None        

        self.hr = False
        self.split = True
        self.do_ax2 = True
        self.do_buttons = True
        self.do_ax3 = True
        self.cut_plot2 = 1.
        self.ax2_fontsize = 12
        self.legend_loc = 1

        self.x_plot_lims = None
        self.y1_plot_lims = None
        self.y2_plot_lims = None
        self.y3_plot_lims = None
 
    def init_obs(self, spectr_obs=None, sp_norm=None, obj_velo=None, limit_sp=None):
        
        if spectr_obs is not None:
            self.set_conf('spectr_obs', spectr_obs)
        if sp_norm is not None:
            self.set_conf('sp_norm', sp_norm)
        if obj_velo is not None:
            self.set_conf('obj_velo', obj_velo)
        if limit_sp is None:
            self.limit_sp = self.get_conf('limit_sp')
        else:
            self.limit_sp = limit_sp

        self.read_obs()
        
    def run(self, do_synth = True, do_read_liste = True, do_profiles=True):
        
        if do_profiles:
            self.do_profile_dict()
        
        if do_synth:
            if do_read_liste:
                self.phyat_arr, self.n_data = self.read_phyat(self.phyat_file)
                self.fic_model = self.get_conf('fic_modele', message='error')
                self.model_arr = self.read_model(self.fic_model)
                self.n_models = len(self.model_arr)
                
                self.fic_cosmetik = self.get_conf('fic_cosmetik', message='warn')
                if self.fic_cosmetik is None:
                    self.cosmetik_arr = []
                else:
                    self.cosmetik_arr = self.read_cosmetik(self.fic_cosmetik)
                self.n_cosmetik = len(self.cosmetik_arr)
        
                self.sp_theo, self.liste_totale, self.liste_raies = \
                    self.append_lists(self.phyat_arr, self.model_arr, self.cosmetik_arr)
        
            self.sp_theo, self.sp_synth = self.make_synth(self.liste_raies, self.sp_theo)
            self.n_sp_theo = len(self.sp_theo['spectr'])
        self.f *= self.aire_ref

        self.sp_abs = self.make_sp_abs(self.sp_theo)

        self.make_filter_instr()
        self.sp_synth_tot = self.convol_synth(self.cont, self.sp_synth)
        self.cont_lr, self.sp_synth_lr = self.rebin_on_obs()
                
    def do_profile_dict(self, return_res=False):
        
        self.fic_profs = self.get_conf('fic_profile', None)
        
        if self.fic_profs is None:
            self.fic_profs = execution_path('./')+'../data/default_profiles.dat'
        else:
            self.fic_profs = self.directory + self.fic_profs
        
        if not os.path.isfile(self.fic_profs):
            log_.error('File not found {}'.format(self.fic_profs), calling=self.calling)
        emis_profiles = {}
        emis_profiles['1'] = {'T4': 0.0, 'vel':0.0, 'params': [['G', '1.00', '0.0', '20.0']]}
        prof_params = None
        with open(self.fic_profs) as f:
            for l in f:
                if l[0] not in ('#', 'c', 'C', ';'):
                    if '#' in l:
                        l = l.split('#')[0]
                    if ';' in l:
                        l = l.split(';')[0]
                    if ':' in l:
                        if prof_params is not None:
                            emis_profiles[key] = {'T4': T4,
                                                  'vel': vel,
                                                  'params' : prof_params}
                        key, T4_vel = l.split(':')
                        T4_vel = T4_vel.split()
                        T4 = np.float(T4_vel[0].strip())
                        if len(T4_vel) == 2:
                            vel = np.float(T4_vel[1].strip())
                        else:
                            vel = 0.0
                        prof_params = []
                    else:
                        if l.split() != []:
                            params = l.split()
                            params[1::] = [np.float(p.strip()) for p in params[1::]]
                            prof_params.append(params)
        log_.message('profile read from {0}'.format(self.fic_profs), 
                                   calling = self.calling)
        emis_profiles[key] = {'T4': T4, 
                              'vel': vel,
                              'params' : prof_params}
        if return_res:
            return emis_profiles
        self.emis_profiles = emis_profiles
        
    def compare_profiles(self):
        
        ref_diff = []
        new_profile = self.do_profile_dict(return_res=True)
        for k in new_profile.keys():
            if k not in self.emis_profiles.keys():
                ref_diff.append(k)
            if new_profile[k]['T4'] != self.emis_profiles[k]['T4']:
                ref_diff.append(k)
            if new_profile[k]['vel'] != self.emis_profiles[k]['vel']:
                ref_diff.append(k)
            for lo, ln in zip(self.emis_profiles[k]['params'], new_profile[k]['params']):
                for llo,lln in zip(lo, ln):
                    if llo != lln:
                        ref_diff.append(k)
        return np.unique(ref_diff)

                
    def get_profile(self, raie):
        
        basic_profiles_dic = {'G': (3, gauss),
                              'C': (3, carre),
                              'L': (3, lorentz)}
        
        
        profile_key =  str(raie['profile'])
        if profile_key not in self.emis_profiles:
            profile_key = '1'
        T4 = self.emis_profiles[profile_key]['T4']
        vel = self.emis_profiles[profile_key]['vel']
        params_str = self.emis_profiles[profile_key]['params']

        lambda_0 = raie['lambda'] + raie['l_shift'] + self.get_conf('lambda_shift', 0.0) 
        w_norm = self.w - lambda_0 - vel * lambda_0 / CST.CLIGHT * 1e5
        
        profile = np.zeros_like(self.w)
        largeur = raie['vitesse'] * lambda_0 / CST.CLIGHT * 1e5
        masse =  2 * (raie['num'] - raie['num'] % 100000000000)/100000000000
        if (masse == 2 and (raie['num'] - raie['num'] % 101000000000)/100000000 == 1010) : 
            masse = 1
        
        for param in params_str:
            profile_type = param[0]
            if profile_type not in basic_profiles_dic:
                log_.error('Wrong number profile reference{}'.format(profile_type), calling=self.calling)
            params = param[1::]
            if len(params) != basic_profiles_dic[profile_type][0]:
                log_.error('Wrong number of parameters {} for profile {}'.format(len(params), profile_type), calling=self.calling)
            profile += basic_profiles_dic[profile_type][1](w_norm, params[0], params[1]*largeur, params[2]*largeur)
        if T4 > 0.0:
            fwhm_therm = 21.4721 * np.sqrt(T4 / masse) * lambda_0 / CST.CLIGHT * 1e5 #km/s
            profile = convolgauss(profile, self.w, lambda_0, fwhm_therm)
        profile[~np.isfinite(profile)] = 0.0
        return profile
                    
                    
    def read_conf(self, config_file=None):
        
        if config_file is None:
            config_file = self.config_file
        else:
            self.config_file = config_file
        self.conf = {}
        execfile(execution_path('./')+'init_defaults.py', self.conf)
        if self.config_file is not None:
            try:
                execfile(self.directory + self.config_file, self.conf)
                log_.message('configuration read from {0}'.format(self.config_file), 
                                   calling = self.calling)
            except:
                log_.warn('configuration NOT read from {0}'.format(self.config_file),
                                calling = self.calling)
            
    def get_conf(self, key=None, undefined=None, message=None):
        """
        Return the value of the key parameter in the configuration. 
        If key is not defined, return the value of the undefined keyword, default being None.
        """
        if key is None:
            for k in self.conf.keys():
                self.get_conf(k)
            return None
        
        if key not in self.conf:
            if message == 'warn':
                log_.warn('{0} not defined in configuration file'.format(key), calling=self.calling)
            elif message == 'message':
                log_.message('{0} not defined in configuration file'.format(key), calling=self.calling)
            elif message == 'error':
                log_.error('{0} not defined in configuration file'.format(key), calling=self.calling)
            else:
                pass
            return undefined
        else:
            return self.conf[key]
                    
        if 'fic_modele' not in self.conf:
            log_.warn('fic_model not defined in configuration file', calling=self.calling)
            return None        
                    
    def set_conf(self, key, value):
        """
        Set the value of the configuration "key" parameter to "value".
        """
        self.conf[key] = value
                    
    def read_phyat(self, phyat_file):
        
        self.phyat_file = phyat_file
        phyat_arr = []
        for dir_ in config.DataPaths:
            try:
                phyat_arr = read_data('{0}/{1}'.format(dir_, self.phyat_file))
                log_.message('phyat data read from {0}/{1}'.format(dir_, self.phyat_file),
                                calling = self.calling)
                break
            except:
                log_.debug('No phyat file found as {0}/{1}'.format(dir_, self.phyat_file),
                                calling = self.calling)
        if len(phyat_arr) == 0:
            log_.error('No phyat file read', calling = self.calling)
            return None
        return phyat_arr, len(phyat_arr)
        
    def read_model(self, model_file):
        
        model_arr = []
        if model_file == 'from phyat':
            mask = self.phyat_arr['ref'] == 999
            model_arr = self.phyat_arr.copy()[mask]
            model_arr['num'] -= 90000000000000
            model_arr['vitesse'] = 10
            log_.message('data initialized from phyat',
                                calling = self.calling) 
        else:
            try:
                model_arr = read_data('{0}'.format(self.directory + model_file))
                log_.message('data read from {0}'.format(self.directory + model_file),
                                    calling = self.calling)
            except:
                log_.warn('unable to read from {0}'.format(self.directory + model_file),
                                    calling = self.calling)
        model_arr['ref'] = 0
        return model_arr
        
    def read_cosmetik(self, cosmetik_file):
        
        cosmetik_arr = []
        try:
            cosmetik_arr = read_data('{0}'.format(self.directory + cosmetik_file))
            log_.message('cosmetik read from {0}'.format(self.directory + cosmetik_file),
                                calling = self.calling)
        except:
            log_.warn('unable to read from {0}'.format(self.directory + cosmetik_file),
                                calling = self.calling)
        return cosmetik_arr

    def read_obs(self, k_spline = 1):
        
        if self.get_conf('spectr_obs') is None:
            n_pix = (self.limit_sp[1] - self.limit_sp[0]) / self.conf['lambda_pix']
            self.w = np.linspace(self.limit_sp[0], self.limit_sp[1], n_pix)
            self.f = np.ones_like(self.w)
        else:
            obs_file = self.directory + self.conf['spectr_obs']+'.spr.gz'
            if not os.path.isfile(obs_file):
                obs_file = self.directory + self.conf['spectr_obs']+'.spr'
            try:                
                self.obs = np.loadtxt(obs_file)
                log_.message('Observations read from {0}'.format(obs_file),
                                    calling = self.calling)
            except:
                self.f = None
                log_.warn('Observations NOT read from {0}'.format(obs_file),
                                    calling = self.calling)
                self.n_lambda = 0.
                return None
        if bool(self.get_conf('data_incl_w', undefined = False)):
            self.w = self.obs[:,0]
            self.f = self.obs[:,1]
        else:
            self.f = self.obs 
        
        if bool(self.get_conf('reverse_spectra', undefined=False)):
                self.f = self.f[::-1]
        self.n_lambda = len(self.f)
        self.tab_pix = np.arange(self.n_lambda)
        
        self.f *= self.get_conf('sp_norm', undefined = 1.)
        
        if ("cal_lambda" in self.conf) and ("cal_pix" in self.conf):
            cal_lambda = np.array(self.conf["cal_lambda"])
            cal_pix = np.array(self.conf["cal_pix"])
            arg_sort = cal_lambda.argsort()
            cal_lambda = cal_lambda[arg_sort]
            cal_pix = cal_pix[arg_sort]
            interp_lam = interpolate.UnivariateSpline(cal_pix, cal_lambda, k=k_spline)
            self.w = interp_lam(self.tab_pix)
            log_.message('Wavelength table generated using spline of order {0}'.format(k_spline),
                                calling = self.calling)
                        
        self.obj_velo = self.get_conf("obj_velo", undefined=0.)
        self.w *= 1 - self.obj_velo/(CST.CLIGHT/1e5)
        log_.message('Wavelenghts shifted by Vel = {} km/s'.format(self.conf["obj_velo"]),
                                calling = self.calling)

        lims = ((self.w >= self.limit_sp[0]) & (self.w <= self.limit_sp[1]))
        log_.message('Observations resized from {0} to {1}'.format(len(self.w), lims.sum()), calling=self.calling)
        self.w = self.w[lims]
        self.f = self.f[lims]

        self.w_ori = self.w.copy()
        self.f_ori = self.f.copy()
        
        resol = self.get_conf('resol', undefined = 1, message=None)
        log_.message('Observations resized from {0} by a factor of {1}'.format(len(self.w), resol), 
                           calling=self.calling)
        self.w = change_size(self.w, resol)
        self.f = change_size(self.f, resol)
        self.n_lambda = len(self.f)
        self.tab_pix = change_size(self.tab_pix, resol)
        self.lambda_pix = (np.max(self.w) - np.min(self.w)) / self.n_lambda
        
    def renorm(self, new_norm):
        self.f /= self.get_conf('sp_norm', undefined = 1.)
        self.f_ori /= self.get_conf('sp_norm', undefined = 1.)
        self.set_conf('sp_norm', new_norm)
        self.f *= self.get_conf('sp_norm', undefined = 1.)
        self.f_ori *= self.get_conf('sp_norm', undefined = 1.)

    def init_red_corr(self):
        self.E_BV = self.get_conf('e_bv', 0.)
        if self.E_BV > 0:
            RC = pn.RedCorr(E_BV = self.E_BV, law=self.get_conf('red_corr_law', message='error'))
            self.red_corr = RC.getCorr(self.w, self.get_conf('lambda_ref_rougi', message='error'))
            log_.message('Reddening correction set to {0}'.format(self.E_BV), calling=self.calling)
        else:
            self.red_corr = np.ones_like(self.w)
        
    def make_continuum(self):
        
        self.conts = {}
        
        user_cont = np.zeros_like(self.w)
        if bool(self.get_conf("cont_in_lambda", False)):
            user_cont_int = interpolate.interp1d(self.conf["cont_lambda"], self.conf["cont_intens"])
            user_cont = user_cont_int(self.w)
        
        cont_pix = self.get_conf("cont_pix", 0.)
        if cont_pix != 0:
            arg_sort = np.array(cont_pix).argsort()
            user_cont_int = interpolate.interp1d(np.array(cont_pix)[arg_sort], 
                                                         np.array(self.get_conf("cont_intens", message='error'))[arg_sort])
            user_cont = user_cont_int(self.tab_pix)    
        self.conts['user'] = user_cont
        
        bb_cont = np.zeros_like(self.w)        
        if "cont_bb_t" in self.conf:
            if np.ndim(self.conf["cont_bb_t"]) == 0:
                tab_T = np.array([self.conf["cont_bb_t"]])
                tab_I = np.array([self.conf["cont_bb_i"]])
            else:
                tab_T = np.array(self.conf["cont_bb_t"])
                tab_I = np.array(self.conf["cont_bb_i"])
            for I, T in zip(tab_I, tab_T):
                bb_cont += I * Planck(self.w, T) / T**4
        self.conts['bb'] = bb_cont
           
        pl_cont = np.zeros_like(self.w)  # Power law
        if "cont_pl_alpha" in self.conf:
            if np.ndim(self.conf["cont_pl_alpha"]) == 0:
                tab_alpha = np.array([self.conf["cont_pl_alpha"]])
                tab_I = np.array([self.conf["cont_pl_i"]])
            else:
                tab_alpha = np.array(self.conf["cont_pl_alpha"])
                tab_I = np.array(self.conf["cont_pl_i"])
            for I, alpha in zip(tab_I, tab_alpha):
                pl_cont += I * (self.w / 5000.)**alpha
        self.conts['pl'] = pl_cont
        
        
        alfa = 1e-13 * 0.668 * (self.conf["cont_Thi"]/1e4)**(-0.507) / \
            (1. + 1.221*(self.conf["cont_Thi"]/1e4)**(0.653)) * 1.000 
        emis_Hi = alfa * CST.HPLANCK * CST.CLIGHT * 1e8 / 4861.3 # erg/s.cm3
        H_cont = self.conf["cont_ihi"] * make_cont_Ercolano(self.conf["cont_Thi"],'H',self.w) / emis_Hi 
        H_cont[~np.isfinite(H_cont)] = 0.
        self.conts['H'] = H_cont
        
        alfa = 1e-13 * 0.331 * (self.conf["cont_Thei"]/1e4)**(-0.615) / \
            (1. + 0.910*(self.conf["cont_Thei"]/1e4)**(0.780)) * 0.7986
        emis_Hei = alfa * CST.HPLANCK * CST.CLIGHT * 1e8 / 4471.5
        He1_cont = self.conf["cont_ihei"] * make_cont_Ercolano(self.conf["cont_Thei"],'He1',self.w) / emis_Hei 
        He1_cont[~np.isfinite(He1_cont)] = 0.
        self.conts['He1'] = He1_cont

        alfa = 2. * 1e-13 * 1.549 * (self.conf["cont_Theii"]/1e4/4.)**(-0.693) / \
            (1. + 2.884*(self.conf["cont_Theii"]/1e4/4.)**(0.609))*1.000
        emis_Heii = alfa * CST.HPLANCK * CST.CLIGHT * 1e8 / 4685.8
        He2_cont = self.conf["cont_iheii"] * make_cont_Ercolano(self.conf["cont_Theii"],'He2',self.w) / emis_Heii 
        He2_cont[~np.isfinite(He2_cont)] = 0.
        self.conts['He2'] = He2_cont
                
        gff_HI = gff(1., self.conf["cont_Thi"], self.w)
        gff_HeI = gff(1., self.conf["cont_Thei"], self.w)
        gff_HeII = gff(4., self.conf["cont_Theii"], self.w)
        
        #32.d0*!phy.e^4.*!phy.h/3./!phy.m_e^2./!phy.c^3.*sqrt(!dpi*13.6*!phy.erg_s_ev/3./!phy.k)= 6.8391014e-38

        FF_cont = (6.8391014e-38 * CST.CLIGHT * 1e8 / self.w**2. * (
                    self.conf["cont_ihi"] * 1.0**2. / np.sqrt(self.conf["cont_Thi"]) * np.exp(-CST.HPLANCK*CST.CLIGHT*1e8/self.w/CST.BOLTZMANN/self.conf["cont_Thi"]) * gff_HI/emis_Hi + 
                    self.conf["cont_ihei"] * 1.0**2./ np.sqrt(self.conf["cont_Thei"]) * np.exp(-CST.HPLANCK*CST.CLIGHT*1e8/self.w/CST.BOLTZMANN/self.conf["cont_Thei"]) * gff_HeI/emis_Hei  + 
                    self.conf["cont_iheii"] * 2.0**2. / np.sqrt(self.conf["cont_Theii"]) * np.exp(-CST.HPLANCK*CST.CLIGHT*1e8/self.w/CST.BOLTZMANN/self.conf["cont_Theii"]) * gff_HeII / emis_Heii))
        FF_cont[~np.isfinite(FF_cont)] = 0.
        self.conts['FF'] = FF_cont

        
        # 2-photons
        #http://adsabs.harvard.edu/abs/1984A%26A...138..495N
        y = 1215.7 / self.w
        A = 202.0 * (y * (1. - y) * (1. -(4. * y * (1 - y))**0.8) + 0.88 * ( y * (1 - y))**1.53 * (4. * y * (1 - y))**0.8)
        alfa_eff = 0.838e-13 * (self.conf["cont_Thi"] / 1e4)**(-0.728) # fit DP de Osterbrock
        q = 5.31e-4 * (self.conf["cont_Thi"] / 1e4)**(-0.17) # fit DP de Osterbrock
        n_crit = 8.226 / q
        twophot_cont = self.conf["cont_ihi"] * CST.HPLANCK * CST.CLIGHT * 1e8 / self.w**3. * 1215.7 * A / 8.226 * alfa_eff / (1. + self.conf["cont_edens"]/n_crit) / emis_Hi
        twophot_cont[~np.isfinite(twophot_cont)] = 0.
        self.conts['2photons'] = twophot_cont

        self.cont = np.zeros_like(self.w)  
        for key in self.conts:
            self.cont += self.conts[key]
        self.cont *= self.aire_ref
        self.cont /= self.red_corr
        
    def plot_conts(self):
        colors = {'bb': 'yellow', 'pl': 'green', '2photons': 'blue', 'FF': 'red',
                  'H': 'red', 'He1': 'yellow', 'He2': 'blue', 'user': 'black'}
        fig_conts = plt.figure()
        if self.ax1 is not None:
            self.ax_cont = fig_conts.add_subplot(111, sharex=self.ax1)
        else:
            self.ax_cont = fig_conts.add_subplot(111)
        
        for key in self.conts:
            if key[0] == 'H':
                style=':'
            else:
                style = '-'
            self.ax_cont.plot(self.w, self.conts[key], linestyle=style, label = key, color = colors[key])
        self.ax_cont.plot(self.w, self.cont, label = 'TOTAL', linestyle='--', linewidth = 2)
        
        self.ax_cont.legend()
        self.ax_cont.set_yscale('log')
            
        
    def append_lists(self, phyat_arr, model_arr, cosmetik_arr):
        
        n_models = len(model_arr)
        liste_totale = phyat_arr.copy()
        liste_totale.resize(len(phyat_arr) + len(model_arr))
        for i,j in enumerate(np.arange(len(phyat_arr), len(phyat_arr) + len(model_arr))):
            liste_totale[j] = model_arr[i]

        sp_theo = {}
        sp_theo['raie_ref'] = model_arr
        sp_theo['correc'] = np.zeros(n_models)
        sp_theo['spectr'] = np.zeros((n_models, len(self.w)))

        if "do_icor_outside_cosmetik" in self.conf:
            liste_totale.i_rel *= liste_totale.i_cor
            liste_totale.i_cor = 1.
            sp_theo['raie_ref'].i_rel *= sp_theo['raie_ref'].i_cor
            sp_theo['raie_ref'].i_cor = 1.
        
        if self.do_cosmetik:
            for line_cosmetik in cosmetik_arr:
                if (line_cosmetik['ref'] == 0) and not bool(self.conf['do_icor_on_ref']):
                    log_.warn('No cosmetik on {0}, reference line'.format(line_cosmetik['num']),
                                    calling = self.calling)
                else:
                    to_change = (liste_totale.num == line_cosmetik['num'])
                    if to_change.sum() == 1:
                        line_to_change = liste_totale[to_change][0]
                        if (line_to_change['lambda'] == line_cosmetik['lambda']) or (line_cosmetik['l_shift'] == 0.):
                            line_to_change['l_shift'] = line_cosmetik['l_shift']
                        if (line_to_change['i_rel'] == line_cosmetik['i_rel']) or (line_cosmetik['i_cor'] == 1.):
                            line_to_change['i_cor'] = line_cosmetik['i_cor']
                        liste_totale[to_change] = line_to_change
                        log_.debug('Cosmetik on {0}'.format([line_cosmetik]), calling=self.calling)
                    elif to_change.sum() == 0:
                        if self.get_conf('warn_on_no_cosmetik'):
                            log_.warn('No cosmetik on {0}, undefined line'.format(line_cosmetik['num']),
                                        calling = self.calling)
                    else:
                        log_.warn('No cosmetik on {0}, multiple defined line'.format(line_cosmetik['num']),
                                        calling = self.calling)

        liste_raies = self.restric_liste(liste_totale)
        log_.message('Size of the line list: {0}, size of the restricted line list: {1}'.format(len(liste_totale),
                                                                                                      len(liste_raies)), calling=self.calling)

        if self.do_cosmetik:
            for line_cosmetik in cosmetik_arr:
                if bool(self.conf['do_icor_on_ref']):
                    to_change = (liste_raies.num == line_cosmetik['num'])
                    if to_change.sum() == 1:
                        line_to_change = liste_raies[to_change][0]
                        line_to_change['vitesse'] *= line_cosmetik['vitesse']
                        if line_cosmetik['profile'] != -1:
                            line_to_change['profile'] = line_cosmetik['profile']
                        liste_raies[to_change] = line_to_change
        
        return sp_theo, liste_totale, liste_raies
        
    def restric_liste(self, liste_in):
        """
        This function changes liste_in
        """
        
        
        """
        We set ref=999 for all the lines depending on a 999 one.
        """
        while True:
            the_end = True
            non_affich = np.where(liste_in['ref'] == 999)[0]
            for i_999 in non_affich:
                this_num = liste_in['num'][i_999]
                dep_non_affich = ((liste_in['ref'] == this_num) & (liste_in['ref'] != 999))
                if dep_non_affich.sum() != 0:
                    before = liste_in['ref'][dep_non_affich]
                    liste_in['ref'][dep_non_affich] = 999
                    log_.message('Before = {0}, After = {1}'.format(before, liste_in['ref'][dep_non_affich]), calling=self.calling)
                    the_end = False
            if the_end:
                break

        where_restr = (((liste_in['lambda'] + liste_in['l_shift']) < np.max(self.w)) & 
                       ((liste_in['lambda'] + liste_in['l_shift']) > np.min(self.w)) &
                       (liste_in['ref'] != 999))
        liste_out = liste_in.copy()[where_restr]
        log_.message('Old size = {0}, new_size = {1}'.format(len(liste_in), len(liste_out)), calling=self.calling)
        
        last_loop = 0
        the_end = False
        
        while True:
            if last_loop == 1:
                the_end = True
            satellites = np.where(liste_out['ref'] != 0)[0]
            for i_satellite in satellites[::-1]:
                if liste_out['ref'][i_satellite] != -1:
                    raie_synth = liste_out[i_satellite].copy()
                    i_main_line = np.where(liste_in['num'] == raie_synth['ref'])[0]
                    if len(i_main_line) != 1:
                        if self.get_conf('warn_on_no_reference'):
                            log_.warn('Satellite sans raie de reference:{0} looking for {1}'.format(raie_synth['num'], raie_synth['ref']), 
                                        calling=self.calling)
                        raie_synth['i_rel'] = 0.0
                        raie_synth['comment'] = '!pas de ref' + raie_synth['comment']
                    else:
                        main_line = liste_in[i_main_line]
                        if main_line['ref'] != 0:
                            raie_synth['i_rel'] *= main_line['i_rel']
                            raie_synth['ref'] = main_line['ref']
                            if bool(self.conf['recursive_i_cor']):
                                raie_synth['i_cor'] *= main_line['i_cor']
                            raie_synth['profile'] = main_line['profile']
                            last_loop = 2
                            log_.debug('filling {0} with {1}'.format(liste_out[i_satellite]['num'],
                                                                             main_line['num']),
                                               calling = self.calling)
                        if last_loop == 1:
                            raie_synth['i_rel'] *= main_line['i_rel']
                            raie_synth['vitesse'] *= main_line['vitesse']
                            raie_synth['l_shift'] += main_line['l_shift']
                            if bool(self.conf['recursive_i_cor']):
                                raie_synth['i_cor'] *= main_line['i_cor']
                            raie_synth['profile'] = main_line['profile']
                            log_.debug('filling {0} with {1}, last loop'.format(liste_out[i_satellite]['num'],
                                                                             main_line['num']),
                                               calling = self.calling)
                    liste_out[i_satellite] = raie_synth
            if last_loop == 0:
                last_loop = 1
            else:
                last_loop = 0
                
            if the_end:
                break
        tt = (np.abs(liste_out['i_rel']) > 1e-50)
        log_.message('number of lines with i_rel > 1e-50: {0}'.format(tt.sum()), calling=self.calling)
        return liste_out[tt]
        
    def make_synth(self, liste_raies, sp_theo):
        
        if bool(self.conf['do_calcul_aire_ref']):
            self.aire_ref = 1.0
                   
        sp_synth = np.zeros_like(self.w)
        sp_theo['spectr'] *= 0.0 
        sp_theo['correc'] *= 0.0
        
        for raie in liste_raies:
            #sp_tmp = self.profil_emis(self.w, raie, self.conf['lambda_shift'])
            sp_tmp = self.get_profile(raie)
            aire = np.trapz(sp_tmp, self.w)
            if np.isfinite(aire) and (aire != 0.):
                max_sp = np.max(sp_tmp)
                if (np.abs(sp_tmp[0]/max_sp) > 1e-3) or (np.abs(sp_tmp[-1]/max_sp) > 1e-3):
                    log_.warn('Area of {0} {1} could be wrong'.format(raie['id'], raie['lambda']), 
                                      calling = self.calling)
                intens_pic = raie['i_rel'] * raie['i_cor'] * self.aire_ref / aire

                if raie['ref'] == 0:
                    tab_tmp = (sp_theo['raie_ref'].num == raie['num'])
                else:
                    tab_tmp = (sp_theo['raie_ref'].num == raie['ref'])
                this_line = intens_pic * sp_tmp
                if not no_red_corr(raie):
                    this_line /= self.red_corr
                if not is_absorb(raie):
                    sp_synth += this_line
                sp_theo['spectr'][tab_tmp] +=  this_line
                sp_theo['correc'][tab_tmp] = 1.0
        tt = (sp_theo['correc'] != 0.)
        for key in ('correc', 'raie_ref', 'spectr'):
            sp_theo[key] = sp_theo[key][tt]
        
        log_.message('Number of theoretical spectra: {0}'.format(len(self.sp_theo['correc'])), calling=self.calling)
        return sp_theo, sp_synth
        
    def make_sp_abs(self, sp_theo):
        
        sp_tau = np.zeros_like(self.w)
        """
        WARNING check also misc.is_absorb(raie)
        """
        index_abs = is_absorb(self.sp_theo['raie_ref'])        
        for i_abs in index_abs:
            sp_tau += self.sp_theo['spectr'][i_abs] * self.sp_theo['correc'][i_abs]
        
        sp_abs = np.exp(sp_tau)
        
        
        if self.get_conf('fic_atm') is not None:
            if type(self.get_conf('fic_atm')) not in (list, tuple):
                self.conf['fic_atm'] = (self.conf['fic_atm'],)
                self.conf['coeff_atm'] = (self.conf['coeff_atm'],)
            if len(self.get_conf('fic_atm')) != len(self.get_conf('coeff_atm')):
                log_.error('fic_atm number {} != coeff_atm number {}'.format(len(self.get_conf('fic_atm')), len(self.get_conf('coeff_atm'))), 
                                      calling = self.calling)
            for fic_atm, coeff_atm, shift_atm in zip(self.get_conf('fic_atm'), self.get_conf('coeff_atm'), self.get_conf('shift_atm')):
                try:
                    d = np.genfromtxt(fic_atm, dtype=None, names=('wl', 'abs'))
                    d['wl'] = vactoair(d['wl'])
                    if type(coeff_atm) not in (list, tuple):
                        coeff_atm = (coeff_atm, )
                    if type(shift_atm) not in (list, tuple):
                        shift_atm = (shift_atm, )
                    for c_atm, s_atm in zip(coeff_atm, shift_atm):
                        abs_interp = interpolate.interp1d(d['wl']*(1+s_atm/CST.CLIGHT*1e5), d['abs'])
                        sp_abs *= np.exp(np.log(abs_interp(self.w)) * c_atm)
                except:
                    log_.warn('Problem in using data from {}'.format(fic_atm), 
                                      calling = self.calling)
        
        # sp_abs /= self.red_corr
        return sp_abs

    def make_filter_instr(self):
        
        filter_size = 11
        increm = 1.1
        detect_limit = 1e-3 / np.max(self.sp_synth)
        
        while True:
            filter_size = int(filter_size * increm)
            if filter_size/2*2 == filter_size:
                filter_size += 1
            if filter_size > self.n_lambda:
                break
            self.filter_ = self.profil_instr(filter_size, self.conf['prof'], self.lambda_pix)
            if (abs(self.filter_[0]) < detect_limit) and (abs(self.filter_[-1]) < detect_limit):
                break
            
        self.filter_ /= self.filter_.sum()
    
    def convol_synth(self, cont, sp_synth):
        
        input_arr = (cont + sp_synth) * self.sp_abs 
        kernel = self.filter_
        sp_synth_tot = convol(input_arr, kernel)
        return sp_synth_tot
        
    def rebin_on_obs(self):
        resol = self.get_conf('resol', undefined = 1, message=None)
        cont_lr = rebin(self.cont, resol)
        sp_synth_lr = rebin(self.sp_synth_tot, resol)
        return cont_lr, sp_synth_lr 
                
    def adjust(self):

        new_model_arr = self.read_model(self.fic_model)        
        new_cosmetik_arr = self.read_cosmetik(self.fic_cosmetik)

        new_sp_theo, new_liste_totale, new_liste_raies = self.append_lists(self.phyat_arr, new_model_arr, new_cosmetik_arr)

        mask_diff = np.zeros(len(new_liste_raies), dtype=bool)
        for key in ('lambda', 'l_shift', 'i_rel', 'i_cor', 'vitesse', 'profile'):
            mask_diff = mask_diff | (new_liste_raies[key] != self.liste_raies[key])
        
        log_.debug('{} differences in lines from files'.format(mask_diff.sum()),
                           calling=self.calling + ' adjust')
        ref_diff = self.compare_profiles()
        log_.debug('{} differences in profile'.format(len(ref_diff)),
                           calling=self.calling + ' adjust')
        for im, l in enumerate(new_liste_raies.profile):
            if np.str(l) in ref_diff:
                mask_diff[im] = True
        
        if mask_diff.sum() > 0:
            old_sp_theo = self.sp_theo.copy()
            if len(new_sp_theo) != len(old_sp_theo):
                log_.error('The new list has different number of elements', 
                                 calling = self.calling+'.adjust')
            
            liste_old_diff = self.liste_raies[mask_diff]
            old_sp_theo, old_sp_synth = self.make_synth(liste_old_diff, old_sp_theo)
            if len(ref_diff) > 0:
                self.do_profile_dict()
                
            liste_new_diff = new_liste_raies[mask_diff]
            new_sp_theo, new_sp_synth = self.make_synth(liste_new_diff, new_sp_theo)
            if log_.level >= 3:
                print('Old values:')
                self.print_line(liste_old_diff)
                print('New values:')
                self.print_line(liste_new_diff)
            do_abs = False
            for i_change in np.arange(len(new_sp_theo['raie_ref'])):
                to_change = (self.sp_theo['raie_ref']['num'] == new_sp_theo['raie_ref'][i_change]['num'])
                #print to_change, i_change
                new_sp_theo['correc'][i_change] = self.sp_theo['correc'][to_change].copy()
                old_sp_theo['correc'][i_change] = self.sp_theo['correc'][to_change].copy()
                if (new_sp_theo['raie_ref'][i_change]['i_rel'] != old_sp_theo['raie_ref'][i_change]['i_rel']):
                    new_sp_theo['correc'][i_change] = 1.0
                self.sp_theo['spectr'][to_change] += new_sp_theo['spectr'][i_change] - old_sp_theo['spectr'][i_change]
                self.sp_theo['raie_ref'][to_change] = new_sp_theo['raie_ref'][i_change]
                self.sp_theo['correc'][to_change] = new_sp_theo['correc'][i_change]
                if is_absorb(new_sp_theo['raie_ref'][i_change]):
                    do_abs = True
                else:
                    self.sp_synth += (new_sp_theo['correc'][i_change] * new_sp_theo['spectr'][i_change] -
                                      old_sp_theo['correc'][i_change] * old_sp_theo['spectr'][i_change]) 
                log_.message('change line {0}'.format(new_sp_theo['raie_ref'][i_change]['num']),
                                   calling=self.calling + ' adjust')            
            if do_abs:
                self.sp_abs = self.make_sp_abs(self.sp_theo)
            self.liste_raies = new_liste_raies
            self.sp_synth_tot = self.convol_synth(self.cont, self.sp_synth)
            self.cont_lr, self.sp_synth_lr = self.rebin_on_obs()
        log_.message('{} differences'.format(mask_diff.sum()), calling=self.calling + ' adjust')

        return mask_diff.sum()
        
        #self.update_plot2()
            
#     def modif_intens(self, raie_num, fact):
#         
#         if fact <= 0.:
#             log_.error('fact must be >0. {0}'.format(fact))
#             return None
#         a_changer = (self.sp_theo['raie_ref']['num'] == raie_num)
#         if a_changer.sum() == 1:
#             old_correc = self.sp_theo['correc'][a_changer][0]
#             if is_absorb(self.sp_theo['raie_ref'][a_changer]):
#                 sp_abs_old = self.make_sp_abs(self.sp_theo[a_changer])
#                 self.sp_theo['correc'][a_changer] = fact
#                 sp_abs_new = self.make_sp_abs(self.sp_theo[a_changer])
#                 self.sp_abs =  self.sp_abs - sp_abs_old + sp_abs_new
#             else:
#                 self.sp_synth += (fact-old_correc) * self.sp_theo['spectr'][a_changer][0] 
#                 self.sp_theo['correc'][a_changer] = fact
#                 #self.sp_theo['spectr'][a_changer] *= fact/old_correc # adding this break the tool
#             self.sp_synth_tot = self.convol_synth(self.cont, self.sp_synth)
#             self.cont_lr, self.sp_synth_lr = self.rebin_on_obs()
#             self.update_plot2()
#             
    def print_line(self, line, sort='lambda'):
        if type(line) == np.core.records.recarray:
            sorts = np.argsort(line[sort])
            for isort in sorts:
                self.print_line(line[isort])
            return
        print('{0[num]:>14d} {0[id]:9s}{0[lambda]:11.3f}{0[l_shift]:6.3f}{0[i_rel]:10.3e}{0[i_cor]:7.3f}'\
              ' {0[ref]:>14d}{0[profile]:5d}{0[vitesse]:7.2f}{0[comment]:25s}'.format(line))
       
    def line_info(self, line_num, sat_info=True, print_header=True, sort='lambda'):
        
        if print_header:
            print('-------- INFO LINES --------')
        if type(line_num) == type(()) or type(line_num) == type([]):
            for line in line_num:
                self.line_info(line, sat_info=sat_info, print_header=False)
            return 
        to_print = (self.liste_raies['num'] == line_num)
        if to_print.sum() == 1:
            raie = self.liste_raies[to_print][0]
            self.print_line(raie, sort=sort)
            if raie['ref'] != 0 and sat_info:
                print('Satellite line of:')
                self.line_info(raie['ref'])
                
        to_print = (self.sp_theo['raie_ref']['num'] == line_num)
        if to_print.sum() > 0 and sat_info:
            raie = self.sp_theo['raie_ref'][to_print][0]
            self.print_line(raie, sort=sort)
            print('Reference line')
            satellites_tab = (self.liste_raies['ref'] == raie['num'])
            Nsat = satellites_tab.sum()
            if Nsat > 0:
                print('{0} satellites'.format(Nsat))
                for raie in self.liste_raies[satellites_tab]:
                    self.print_line(raie, sort=sort)
            if self.sp_theo['correc'][to_print][0] != 1.0:
                print('Intensity corrected by {0}'.format(self.sp_theo['correc'][to_print][0]))
            
    def plot1(self):
        f, ax = plt.subplots()
        ax.step(self.w_ori, self.f_ori, label='Obs')
        ax.step(self.w_ori, self.sp_synth_lr, label='Synth')
        ax.legend()
        return f, ax
        
    def plot2(self, hr=False, cut=None, split=False, do_ax2 = True, do_ax3 = True,
              do_buttons=True, xlims=None, fontsize=12, legend_loc=1, fig=None, 
              magenta_ref=None, magenta_lab=None,
              cyan_ref = None, cyan_lab=None, call_init_axes=True):
        
        
        log_.message('entering plots, ID(ax1)'.format(id(self.fig1)), calling=self.calling)
        self.hr = hr
        self.split = split
        self.do_ax2 = do_ax2
        self.do_buttons = do_buttons
        self.do_ax3 = do_ax3
        if cut is not None:
            self.cut_plot2 = cut
        self.ax2_fontsize = fontsize
        self.legend_loc = legend_loc

        if magenta_ref is not None:
            self.plot_magenta = magenta_ref
            self.label_magenta = magenta_lab
        else:
            self.plot_magenta = self.get_conf('plot_magenta')
            self.label_magenta = self.get_conf('label_magenta')      
        if cyan_ref is not None:
            self.plot_cyan   = cyan_ref
            self.label_cyan = cyan_lab
        else:
            self.plot_cyan = self.get_conf('plot_cyan')
            self.label_cyan = self.get_conf('label_cyan')      
        
        if fig is None:
            self.fig1 = plt.figure()
            log_.message('creating new figure ID {}'.format(id(self.fig1)), calling=self.calling)
        else:
            self.fig1 = fig
            self.fig1.clf()
            log_.message('using argument figure ID: {} {}'.format(id(self.fig1), id(fig)), calling=self.calling)
        if split:
            if do_ax2:
                self.fig2 = plt.figure()
            if do_ax3:
                self.fig3 = plt.figure()
            self.ax1 = self.fig1.add_subplot(111)
            if do_ax2:
                self.ax2 = self.fig2.add_subplot(111, sharex=self.ax1)
            if do_ax3:
                self.ax3 = self.fig3.add_subplot(111, sharex=self.ax1)
        else:
            n_subplots = 1
            i_ax2 = 2
            i_ax3 = 2
            if do_ax2:
                n_subplots += 1
                i_ax3 += 1
            if do_ax3:
                n_subplots += 1
            self.ax1 = self.fig1.add_subplot(n_subplots, 1, 1)
            if do_ax2:
                self.ax2 = self.fig1.add_subplot(n_subplots, 1, i_ax2, sharex=self.ax1)
            if do_ax3:
                self.ax3 = self.fig1.add_subplot(n_subplots, 1, i_ax3, sharex=self.ax1)
        
        self.plot_ax1(self.ax1, xlims=xlims)
        if do_ax2:
            self.plot_ax2(self.ax2)
        if do_ax3:
            self.plot_ax3(self.ax3)
        if do_buttons:
            self._make_buttons(split=split)
        if call_init_axes:
            self.init_axes()
        self.restore_axes()
        plt.subplots_adjust(hspace=0.0)
               
    def plot_ax1(self, ax, xlims=None):
        
        ax.step(self.w_ori, self.f_ori, label='Obs', c='red')
        self.ax1_line_synth = ax.step(self.w_ori, self.sp_synth_lr, label='Synth', c='blue', linewidth=2)[0]
        if self.hr:
            ax.step(self.w, self.sp_synth, c='green')
        if self.plot_magenta is not None:
            i_magenta = np.where(self.sp_theo['raie_ref']['num'] == self.plot_magenta)[0]
            if self.label_magenta is None:
                label_magenta = clean_label(self.sp_theo['raie_ref'][i_magenta]['id'][0])
            else:
                label_magenta = self.label_magenta
            if len(i_magenta) == 1:
                self.ax1_line_magenta = ax.step(self.w, self.cont+self.sp_theo['spectr'][i_magenta][0], c='magenta', 
                              label=label_magenta, linestyle='--')[0]
        if self.plot_cyan is not None:
            i_cyan = np.where(self.sp_theo['raie_ref']['num'] == self.plot_cyan)[0]
            if self.label_cyan is None:
                label_cyan = clean_label(self.sp_theo['raie_ref'][i_cyan]['id'][0])
            else:
                label_cyan = self.label_cyan
            if len(i_cyan) == 1:
                self.ax1_line_cyan = ax.step(self.w, self.cont+self.sp_theo['spectr'][i_cyan][0], c='cyan', 
                              label=label_cyan, linestyle='-')[0]
        ax.legend(loc=self.legend_loc)
        log_.debug('ax1 drawn on ax ID {}'.format(id(ax)), calling=self.calling)
        
    def plot_ax2(self, ax):        
        
        for line in self.liste_raies:
            wl = line['lambda'] + line['l_shift'] + self.conf['lambda_shift']
            i_rel = line['i_rel']
            if (abs(i_rel) > self.cut_plot2) & (wl > ax.get_xlim()[0]) & (wl < ax.get_xlim()[1]):
                ax.plot([wl, wl], [0, 1], color='blue')
                ax.text(wl, -0.2, '{0} {1:7.4f}'.format(line['id'], i_rel), 
                              rotation='vertical', fontsize=self.ax2_fontsize)
        """
        ax.set_xlim(self.get_conf('x_plot_lims'))
        ax.set_ylim(self.get_conf('y2_plot_lims'))
        """
        log_.debug('ax2 drawn on ax ID {}'.format(id(ax)), calling=self.calling)
        

    def plot_ax3(self, ax):     

        ax.step(self.w, self.f - self.cont, c = 'red', linestyle='--')
        ax.plot((0, 1e10), (0.0, 0.0), c='green')
        ax_line_diff = ax.step(self.w_ori, self.f_ori - self.sp_synth_lr, c='blue')[0]

        """    
        ax.set_xlim(self.get_conf('x_plot_lims'))
        if self.get_conf('y3_plot_lims') is None:
            mask = (self.w > ax.get_xlim()[0]) & (self.w < ax.get_xlim()[1]) 
            ax.set_ylim((np.min((self.f - self.cont)[mask]), np.max((self.f - self.cont)[mask])))            
        else:
            ax.set_ylim(self.get_conf('y3_plot_lims'))
        """
        log_.debug('ax3 drawn on ax ID {}'.format(id(ax)), calling=self.calling)
        
    def update_plot2(self):
        
        if self.ax1 is None:
            return
        self.ax1_line_synth.remove()
        self.ax1_line_synth = self.ax1.step(self.w_ori, self.sp_synth_lr, label='Synth', c='blue', linewidth=2)[0]
        self.ax1.legend(loc=self.legend_loc)
        if self.plot_magenta is not None:
            try:
                self.ax1_line_magenta.remove()
            except:
                pass
            i_magenta = np.where(self.sp_theo['raie_ref']['num'] == self.plot_magenta)[0]
            if self.label_magenta is None:
                self.label_magenta = self.sp_theo['raie_ref'][i_magenta]['id']
            if len(i_magenta) == 1:
                self.ax1_line_magenta = self.ax1.step(self.w, self.cont+self.sp_theo['spectr'][i_magenta][0], c='magenta', 
                                                      label=self.label_magenta, linestyle='--')[0]
        if self.plot_cyan is not None:
            try:
                self.ax1_line_cyan.remove()
            except:
                pass
            i_cyan = np.where(self.sp_theo['raie_ref']['num'] == self.plot_cyan)[0]
            if self.label_cyan is None:
                self.label_cyan = self.sp_theo['raie_ref'][i_cyan]['id']
            if len(i_cyan) == 1:
                self.ax1_line_cyan = self.ax1.step(self.w, self.cont+self.sp_theo['spectr'][i_cyan][0], c='cyan', 
                                                   label=self.label_cyan, linestyle='-')[0]
        
        for i in np.arange(len(self.ax2.texts)):
            self.ax2.texts.pop()
        for i in np.arange(len(self.ax2.lines)):
            self.ax2.lines.pop()
        
        i_max = np.max(self.liste_raies['i_rel'])
        for line in self.liste_raies:
            wl = line['lambda'] + line['l_shift'] + self.conf['lambda_shift']
            i_rel = line['i_rel']
            if (abs(i_rel) > self.cut_plot2) & (wl > self.ax2.get_xlim()[0]) & (wl < self.ax2.get_xlim()[1]):
                self.ax2.plot([wl, wl], [0, 1], color='blue')
                self.ax2.text(wl, -0.2, '{0} {1:7.4f}'.format(line['id'], i_rel), 
                                rotation='vertical', fontsize=self.ax2_fontsize)
                #self.ax2.set_ylim((-1.5, 1))
 
        if self.do_ax3:
            self.ax3_line_diff.remove()
            self.ax3_line_diff = self.ax3.step(self.w_ori, self.f_ori - self.sp_synth_lr, c='blue')[0]
        
        self.fig1.canvas.draw()
 
    def init_axes(self):
        self.x_plot_lims = self.get_conf('x_plot_lims')
        if self.x_plot_lims is None:
            self.x_plot_lims = (np.min(self.w), np.max(self.w))
            
        self.y1_plot_lims = self.get_conf('y1_plot_lims')
        if self.y1_plot_lims is None:
            mask = (self.w_ori > self.x_plot_lims[0]) & (self.w_ori < self.x_plot_lims[1]) 
            self.y1_plot_lims = (np.min(self.sp_synth_lr[mask]), np.max(self.sp_synth_lr[mask]))      
        
        self.y2_plot_lims = self.get_conf('y2_plot_lims')
        if self.y2_plot_lims is None:
            self.y2_plot_lims = (-1.5, 1)
        
        self.y3_plot_lims = self.get_conf('y3_plot_lims')
        if self.y3_plot_lims is None:
            mask = (self.w_ori > self.x_plot_lims[0]) & (self.w_ori < self.x_plot_lims[1])
            self.y3_plot_lims = (np.min((self.f - self.cont)[mask]), np.max((self.f - self.cont)[mask]))
        log_.message('Axes initialized', calling=self.calling)
        self.print_axes()
                         
    def save_axes(self):
        if self.ax1 is not None:
            self.x_plot_lims = self.ax1.get_xlim()
            self.y1_plot_lims = self.ax1.get_ylim()
        else:
            self.x_plot_lims = None
            self.y1_plot_lims = None
        if self.ax2 is not None:
            self.y2_plot_lims = self.ax2.get_ylim()
        else:
            self.y2_plot_lims = None
        if self.ax3 is not None:
            self.y3_plot_lims = self.ax3.get_ylim()
        else:
            self.y3_plot_lims = None
        log_.message('Axes saved', calling=self.calling)
        self.print_axes()
        
    def restore_axes(self):
        if self.x_plot_lims is not None:
            if self.ax1 is not None:
                self.ax1.set_xlim(self.x_plot_lims)
                log_.message('X-axes restored to {}'.format(self.ax1.get_xlim()), calling=self.calling)
            else:
                log_.message('ax1 is None', calling=self.calling)
        else:
            log_.message('x_plot_lims is None', calling=self.calling)
        if self.y1_plot_lims is not None:
            if self.ax1 is not None:
                self.ax1.set_ylim(self.y1_plot_lims)
        if self.y2_plot_lims is not None:
            if self.ax2 is not None:
                self.ax2.set_ylim(self.y2_plot_lims)
        if self.y3_plot_lims is not None:
            if self.ax3 is not None:
                self.ax3.set_ylim(self.y3_plot_lims)
        
        log_.message('Axes restored', calling=self.calling)
        self.print_axes()
        
    def print_axes(self):
        log_.debug('{} {} {} {}'.format(self.x_plot_lims, self.y1_plot_lims, self.y2_plot_lims, self.y3_plot_lims), calling=self.calling)
    
    def rerun(self):
        self.run(do_synth = True, do_read_liste = True, do_profiles=True)

    def replot2(self):    
        self.save_axes()
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.plot2(hr=self.hr, cut=self.cut_plot2, split=self.split, 
                   do_ax2=self.do_ax2, do_ax3=self.do_ax3, do_buttons=self.do_buttons, 
                   fontsize=self.ax2_fontsize, legend_loc=self.legend_loc, fig=self.fig1, call_init_axes=False)
        self.fig1.canvas.draw()
        
    def _make_buttons(self, split):
        if split:
            self.fig1.subplots_adjust(bottom=0.3)
        else:
            self.fig1.subplots_adjust(bottom=0.2)
            
        self.buttons = {}
        b_w = 0.06
        b_h = 0.06
        b_x0 = 0.05
        b_y0 = 0.02
        ax_zxm = self.fig1.add_axes([b_x0, b_y0, b_w, b_h]) 
        self.buttons['ZX-'] = Button(ax_zxm, 'ZX-')
        self.buttons['ZX-'].on_clicked(self._ZoomZXm)
        ax_zxp = self.fig1.add_axes([b_x0, b_y0 + b_h, b_w, b_h]) 
        self.buttons['ZX+'] = Button(ax_zxp, 'ZX+')
        self.buttons['ZX+'].on_clicked(self._ZoomZXp)
        ax_zym = self.fig1.add_axes([b_x0 + b_w, b_y0, b_w, b_h]) 
        self.buttons['Zy-'] = Button(ax_zym, 'ZY-')
        self.buttons['Zy-'].on_clicked(self._ZoomZYm)
        ax_zyp = self.fig1.add_axes([b_x0 + b_w, b_y0 + b_h, b_w, b_h]) 
        self.buttons['Zy+'] = Button(ax_zyp, 'ZY+')
        self.buttons['Zy+'].on_clicked(self._ZoomZYp)
        ax_sxm = self.fig1.add_axes([b_x0 + 2*b_w, b_y0, b_w, b_h]) 
        self.buttons['SX-'] = Button(ax_sxm, 'SX-')
        self.buttons['SX-'].on_clicked(self._ZoomSXm)
        ax_sxp = self.fig1.add_axes([b_x0 + 2*b_w, b_y0 + b_h, b_w, b_h]) 
        self.buttons['SX+'] = Button(ax_sxp, 'SX+')
        self.buttons['SX+'].on_clicked(self._ZoomSXp)
        ax_sym = self.fig1.add_axes([b_x0 + 3*b_w, b_y0, b_w, b_h]) 
        self.buttons['Sy-'] = Button(ax_sym, 'SY-')
        self.buttons['Sy-'].on_clicked(self._ZoomSYm)
        ax_syp = self.fig1.add_axes([b_x0 + 3*b_w, b_y0 + b_h, b_w, b_h]) 
        self.buttons['Sy+'] = Button(ax_syp, 'SY+')
        self.buttons['Sy+'].on_clicked(self._ZoomSYp)
        ax_curson = self.fig1.add_axes([b_x0 + 5*b_w, b_y0, 2*b_w, b_h]) 
        self.buttons['CursOn'] = Button(ax_curson, 'CursOn')
        self.buttons['CursOn'].on_clicked(self._cursOn)
        ax_curson = self.fig1.add_axes([b_x0 + 5*b_w, b_y0 + b_h, 2*b_w, b_h]) 
        self.buttons['CursOff'] = Button(ax_curson, 'CursOff')
        self.buttons['CursOff'].on_clicked(self._cursOff)
        ax_rerun = self.fig1.add_axes([b_x0 + 7*b_w, b_y0 + b_h, 2*b_w, b_h])
        self.buttons['Rerun'] = Button(ax_rerun, 'Rerun')
        self.buttons['Rerun'].on_clicked(self._call_rerun)
        ax_adjust = self.fig1.add_axes([b_x0 + 7*b_w, b_y0, 2*b_w, b_h])
        self.buttons['Adjust'] = Button(ax_adjust, 'Adjust')
        self.buttons['Adjust'].on_clicked(self._call_adjust)
        ax_readobs = self.fig1.add_axes([b_x0 + 9*b_w, b_y0 + b_h, 2*b_w, b_h])
        self.buttons['ReadObs'] = Button(ax_readobs, 'ReadObs')
        self.buttons['ReadObs'].on_clicked(self._call_readobs)
        ax_replot = self.fig1.add_axes([b_x0 + 9*b_w, b_y0, 2*b_w, b_h])
        self.buttons['RePlot'] = Button(ax_replot, 'RePlot')
        self.buttons['RePlot'].on_clicked(self._call_replot)
        
    def _ZoomZXm(self, event=None): 
        self._Zoom('ZX-')

    def _ZoomZXp(self, event=None): 
        self._Zoom('ZX+')
        
    def _ZoomZYm(self, event=None): 
        self._Zoom('ZY-')

    def _ZoomZYp(self, event=None): 
        self._Zoom('ZY+')
        
    def _ZoomSXm(self, event=None): 
        self._Zoom('SX-')

    def _ZoomSXp(self, event=None): 
        self._Zoom('SX+')
        
    def _ZoomSYm(self, event=None): 
        self._Zoom('SY-')

    def _ZoomSYp(self, event=None): 
        self._Zoom('SY+')
        
    def _Zoom(self, zoom_direction):
        
        """
        zoom_direction = 'ABC', with A in ['S', 'Z'], B in ['X', 'Y'], and C in ['+', '-']
        """
        xmin, xmax = self.ax1.get_xlim()
        dx = xmax - xmin
        
        ymin, ymax = self.ax1.get_ylim()
        dy = ymax - ymin
        
        if zoom_direction[0] == 'S':
            if zoom_direction[2] == '+':
                coeff = self.zoom_fact
            elif zoom_direction[2] == '-':
                coeff = -self.zoom_fact
            if zoom_direction[1] == 'X':
                xmin += coeff * dx
                xmax += coeff * dx
            elif zoom_direction[1] == 'Y':
                ymin += coeff * dy
                ymax += coeff * dy
        elif zoom_direction[0] == 'Z':
            if zoom_direction[2] == '+':
                coeff = self.zoom_fact
            elif zoom_direction[2] == '-':
                coeff = -self.zoom_fact
            if zoom_direction[1] == 'X':
                xmin += coeff * dx 
                xmax -= coeff * dx 
            elif zoom_direction[1] == 'Y':
                ymin += coeff * dy
                ymax -= coeff * dy 
        
        self.ax1.set_xlim((xmin, xmax))
        self.ax1.set_ylim((ymin, ymax))
        self.fig1.canvas.draw()
        
        
    def _cursOn(self, event=None):
        self._cid = self.fig1.canvas.mpl_connect('button_press_event', self._curs_onclick)
        log_.message('Cursor ON', calling=self.calling)

    def _cursOff(self, event=None):
        if self._cid is not None:
            self.fig1.canvas.mpl_disconnect(self._cid)
            log_.message('Cursor OFF', calling=self.calling)
                    
    def _curs_onclick(self, event):
        
        wl = event.xdata
        if (wl > self.ax1.get_xlim()[1]) or (wl < self.ax1.get_xlim()[0]) or (event.button == 2):
            self._cursOff()
            return None
        if event.button in (1,3):
            wl_lim = self.cursor_width * (self.ax1.get_xlim()[1] - self.ax1.get_xlim()[0])
            tt = (np.abs(self.liste_raies['lambda'] + self.liste_raies['l_shift'] + self.conf['lambda_shift'] - wl) < wl_lim)
            if tt.sum() > 0:
                print('')
                print('------------ CURSOR on {0:10.3f} -----------'.format(wl))
                self.print_line(self.liste_raies[tt])
                print('---------------------------------------------')

    def _call_adjust(self, event=None):
        self.adjust()
        self.update_plot2()

    def _call_rerun(self, event=None):
        self.rerun()
        self.replot2()
        
    def _call_readobs(self, event=None):
                
        self.init_obs(spectr_obs=None, sp_norm=None, obj_velo=None, limit_sp=self.limit_sp)
        self.init_red_corr()
        self.make_continuum()               
        self.sp_synth_tot = self.convol_synth(self.cont, self.sp_synth)
        self.cont_lr, self.sp_synth_lr = self.rebin_on_obs()
        self.replot2()

    def _call_replot(self, event=None):
        self.replot2()        

    def plot_indiv_sp(self, y_shift_coeff=None, legend_zoom=.115):
        """
        Seems buggy, loosing ax1.xlim when plotted.
        """
        if y_shift_coeff is None:
            y_shift_coeff = np.max(self.sp_theo['spectr'])/100.
        fig_indiv_spectra = plt.figure()
        if self.ax1 is not None:
            ax_is = fig_indiv_spectra.add_subplot(111, sharex=self.ax1)
        else:
            ax_is = fig_indiv_spectra.add_subplot(111)
        for i in np.arange(self.n_sp_theo):
            label = self.sp_theo['raie_ref']['id'][i]
            ax_is.plot(self.w, self.sp_theo['spectr'][i] + y_shift_coeff*(self.n_sp_theo - i), label=label)
        ax_is.set_ylim((0, y_shift_coeff*(i+2)))
        ax_is.legend(fontsize= i * legend_zoom )
        
        self.fig_indiv_spectra = fig_indiv_spectra
        self.ax_indiv_spectra = ax_is
        
    def plot_profile(self):
        
        self.fig_prof = plt.figure()
        self.ax_prof = plt.semilogy(self.filter_)
        
def main_loc(config_file):
    """
    In case of not having Qt4.
    Usage:
    from pyssn.core.spectrum import main_loc
    sp = main_loc('./s6302_n_c_init.py')
    """
    sp = spectrum(config_file=config_file)
    fig = plt.figure(figsize=(20, 7))
    sp.plot2(fig=fig)
    sp.save_axes()
    plt.show()
    return sp

def main():
    """
    
    """
    parser = get_parser()
    args = parser.parse_args()
    if args.file is None:
        log_.error('A file name is needed, use option -f')
    log_.level = args.verbosity
    sp = spectrum(config_file=args.file)
    fig = plt.figure(figsize=(20, 7))
    sp.plot2(fig=fig)
    plt.show()

