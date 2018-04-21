"""
pySSN is available under the GNU licence providing you cite the developpers names:

    Ch. Morisset (Instituto de Astronomia, Universidad Nacional Autonoma de Mexico)

    D. Pequignot (Meudon Observatory, France)
"""
import time
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from scipy import interpolate
from collections import OrderedDict

from pyssn import log_, config
if config.INSTALLED['PyNeb']:
    import pyneb as pn

from ..utils.physics import CST, Planck, make_cont_Ercolano, gff
from ..utils.misc import execution_path, change_size, convol, rebin, is_absorb, no_red_corr, gauss, carre, lorentz, convolgauss 
from ..utils.misc import vactoair, clean_label,  get_parser, read_data, my_execfile as execfile
from ..core.profiles import profil_instr

"""
ToDo:
1) Define the special lines in a table or dictionnary where ref, color, style are set up.
"""


def save_data(filename, array_to_save, NF=True):
#     #dtype = 'i8, a1, a9, f, f, f, f, a1, i8, i4, f, a25'
#     dtype = 'i8, a1, a9, f, f, f, f, a1, i8, i4, f, a100'
#     if NF:
#         #delimiter = [14, 1, 9, 11, 6, 10, 7, 1, 14, 4, 7, 25]
#         delimiter = [14, 1, 9, 11, 6, 10, 7, 1, 14, 4, 7, 100]
#     else:
#         delimiter = [ 9, 1, 9, 11, 6, 10, 7, 1,  9, 4, 7, 25]
#     names = ['num', 'foo', 'id', 'lambda','l_shift', 'i_rel', 'i_cor', 'foo2', 'ref', 'profile', 
#              'vitesse', 'comment']
#     usecols = (0, 2, 3, 4, 5, 6, 8, 9, 10, 11)
#     #np.savetxt(filename, array_to_save, dtype=dtype, delimiter=delimiter, names = names, usecols = usecols)
#     #np.savetxt(filename, array_to_save, delimiter=delimiter, fmt="%s")
# 
#     fmt="%s"
#     fmt = '{>14d} {9s}{11.3f}{6.3f}{10.3e}{7.3f}{>14d}{5d}{7.2f}{s}'
#             
#     #np.savetxt(filename, array_to_save, fmt=fmt)

    for line in array_to_save:
        print('{0[num]:>14d} {0[id]:9s}{0[lambda]:11.3f}{0[l_shift]:6.3f}{0[i_rel]:10.3e}{0[i_cor]:7.3f} {0[ref]:>14d}{0[profile]:5d}{0[vitesse]:7.2f}{1:1s}'.format(line, line['comment'].strip()))

def read_data(filename, NF=True):
    #dtype = 'i8, a1, a9, f, f, f, f, a1, i8, i4, f, a25'
    dtype = 'i8, a1, a9, f, f, f, f, a1, i8, i4, f, a100'
    if NF:
        #delimiter = [14, 1, 9, 11, 6, 10, 7, 1, 14, 4, 7, 25]
        delimiter = [14, 1, 9, 11, 6, 10, 7, 1, 14, 4, 7, 100]
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

    if dd.size == 1:
        dd = np.atleast_1d(dd)
    return dd.view(np.recarray)


class spectrum(object):
    
    def __init__(self, config_file=None, phyat_file=None, profil_instr=profil_instr, 
                 do_synth = None, do_read_liste = None, do_cosmetik = None, do_run = True, limit_sp = None,
                 spectr_obs=None, sp_norm=None, obj_velo=None, post_proc_file=None):
        """
        Main pySSN object.
        It reads the configuration file given by the config_file parameter. 
        It reads the atomic data, model and cosmetik files. It reads the observation. It computes the reddening
        correction, the 
        """
                
        self.cursor = None
        
        self.selected_ions_data = None
        
        self.fields = [ 'num', 'id', 'lambda', 'l_shift', 'i_rel', 'i_cor', 'ref', 'profile', 'vitesse', 'comment' ]

        self.field_width = { 'num'     : 14,
                             'id'      : 9, 
                             'lambda'  : 11, 
                             'l_shift' : 6, 
                             'l_tot'   : 11, 
                             'i_rel'   : 10, 
                             'i_cor'   : 7, 
                             'i_tot'   : 10,
                             'ref'     : 14,
                             'profile' : 4,
                             'vitesse' : 7, 
                             'comment' : 100 }

        self.field_align = { 'num'     : '>',
                             'id'      : '<', 
                             'lambda'  : '>', 
                             'l_shift' : '>', 
                             'l_tot'   : '>', 
                             'i_rel'   : '>', 
                             'i_cor'   : '>', 
                             'i_tot'   : '>',
                             'ref'     : '>',
                             'profile' : '>',
                             'vitesse' : '>', 
                             'comment' : '<' }

        self.field_pos = { 'num'     : 0,
                           'id'      : 15, 
                           'lambda'  : 24, 
                           'l_shift' : 35, 
                           'i_rel'   : 41, 
                           'i_cor'   : 51, 
                           'ref'     : 59,
                           'profile' : 73,
                           'vitesse' : 77, 
                           'comment' : 85 }

        self.field_format = { 'num'     : '{:>14d}',
                              'id'      : '{:9s}', 
                              'lambda'  : '{:11.3f}', 
                              'l_shift' : '{:6.3f}', 
                              'l_tot'   : '{:11.3f}', 
                              'i_rel'   : '{:10.3e}', 
                              'i_cor'   : '{:7.3f}', 
                              'i_tot'   : '{:10.3e}',
                              'ref'     : '{:>14d}',
                              'profile' : '{:>4d}',
                              'vitesse' : '{:7.2f}', 
                              'comment' : '{:s}' }

        self.field_tip = { 'num'     : 'line code number',
                           'id'      : 'ion', 
                           'lambda'  : 'wavelength in air', 
                           'l_shift' : 'wavelength additive correction', 
                           'l_tot'   : 'corrected wavelength', 
                           'i_rel'   : 'relative intensity',
                           'i_cor'   : 'intensity correction factor',
                           'i_tot'   : 'corrected intensity', 
                           'ref'     : 'reference line code number',
                           'profile' : 'line profile code number',
                           'vitesse' : 'natural line width', 
                           'comment' : 'comment' }

        self.field_abbr = { 'num'     : 'line number',
                            'id'      : 'ion', 
                            'lambda'  : 'wavelength', 
                            'l_shift' : 'w shift', 
                            'l_tot'   : 'corr wave', 
                            'i_rel'   : 'intensity',
                            'i_cor'   : 'i factor',
                            'i_tot'   : 'corr int', 
                            'ref'     : 'ref line',
                            'profile' : 'profile',
                            'vitesse' : 'v factor', 
                            'comment' : 'comment' }
        
        self.calling = 'spectrum'
        self.full_config_file = config_file
        if '/' in self.full_config_file:
            file_name = self.full_config_file.split('/')[-1]
            dir_ = self.full_config_file.split(file_name)[0]
            if dir_ == '':
                dir_ = './'
            self.directory = dir_
            self.config_file = file_name
        else:
            self.directory = './'
            self.config_file = self.full_config_file
        config.addDataFilePath(self.directory, inpySSN=False)
        
        self.init_vars()
         
        self.read_conf(self.config_file)
        log_.level = self.get_conf('log_level', 2)
        if not self.get_conf('do_synth'):
            self.set_conf('plot_ax2', False)
            self.set_conf('plot_ax3', False)
            self.set_conf('fic_cosmetik', 'NO_cosmetik.dat')
            self.set_conf('fic_modele', 'NO_modele.dat')
            self.set_conf('phyat_file', 'NO_phyat.dat')
        if self.get_conf('spectr_obs') is None:
            self.set_conf('plot_ax3', False)            
            
        if do_synth is None:
            self.do_synth = self.get_conf('do_synth')
        else:
            self.do_synth = do_synth
            
        if do_read_liste is None:
            do_read_liste = self.get_conf('do_read_liste')
        if do_cosmetik is None:
            do_cosmetik = self.get_conf('do_cosmetik')

        self.profil_instr = profil_instr
        self.do_cosmetik = do_cosmetik        
        self.post_proc_file = post_proc_file
        self.init_obs(spectr_obs=spectr_obs, sp_norm=sp_norm, obj_velo=obj_velo, limit_sp=limit_sp)
        self.init_red_corr()
        self.make_continuum()
            
        if phyat_file is not None:
            self.phyat_file = phyat_file
        else:
            self.phyat_file = self.get_conf('phyat_file', 'liste_phyat.dat')
        if do_run:
            self.run(do_synth = self.do_synth, do_read_liste = do_read_liste)
        
    def init_vars(self):
        self.fig1 = None
        self.fig2 = None
        self.fig3 = None
        self.ax1 = None
        self.ax2 = None
        self.ax3 = None
        self.cursor_width = 0.02
        self.cursor_w0 = None
        self.cursor_w1 = None
        self.cursor_w2 = None
        self.firstClick = True
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
                self.fic_model = self.get_conf('fic_modele', message='error')
                self.model_arr = self.read_model(self.fic_model)
                self.phyat_arr, self.n_data = self.read_phyat(self.phyat_file)
                self.n_models = len(self.model_arr)
                
                self.fic_cosmetik = self.get_conf('fic_cosmetik', message='warn')
                self.cosmetik_arr = self.read_cosmetik(self.fic_cosmetik, self.do_cosmetik)
                self.n_cosmetik = len(self.cosmetik_arr)

                self.sp_theo, self.liste_totale, self.liste_raies = \
                    self.append_lists(self.phyat_arr, self.model_arr, self.cosmetik_arr)
        
            self.sp_theo, self.sp_synth = self.make_synth(self.liste_raies, self.sp_theo)
            self.n_sp_theo = len(self.sp_theo['spectr'])
        else:
            self.sp_theo = None
            self.sp_synth = None
            self.n_sp_theo = 0
            
        self.f *= self.aire_ref

        self.sp_abs = self.make_sp_abs(self.sp_theo)

        self.make_filter_instr()
        self.sp_synth_tot = self.convol_synth(self.cont, self.sp_synth)
        self.cont_lr, self.sp_synth_lr = self.rebin_on_obs()
                
    def do_profile_dict(self, return_res=False):
        
        self.fic_profs = self.get_conf('fic_profile', None)
        self.instr_prof_file = self.get_conf('fic_instr_prof', None)
        
        if self.fic_profs is None:
            self.fic_profs = execution_path('./')+'../data/default_profiles.dat'
        else:
            self.fic_profs = self.directory + self.fic_profs

        if self.instr_prof_file is not None:
            try:
                user_module = {}
                execfile(self.instr_prof_file, user_module)
                prof = user_module['prof']
                self.conf['prof'] = prof
                log_.message('instrumental profile read from {}'.format(self.instr_prof_file))
            except:
                log_.warn('instrumental profile NOT read from {}. Check it contains only ASCII characters.'.format(self.instr_prof_file), 
                              calling = self.calling)
                self.instr_prof_file = None

        if not os.path.isfile(self.fic_profs):
            log_.error('File not found {}'.format(self.fic_profs), calling=self.calling)
        emis_profiles = {}
        emis_profiles['1'] = {'T4': 1e4, 'vel':0.0, 'params': [['G', '1.00', '0.0', '20.0']]}
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
        
        self.plot_magenta = self.get_conf('plot_magenta', None)
        self.label_magenta = self.get_conf('label_magenta', None)      
        self.plot_cyan = self.get_conf('plot_cyan', None)
        self.label_cyan = self.get_conf('label_cyan', None)
        
        # If you DON'T want an i_cor on a main line to affect the satellites, 
        # set the following variable to False
        self.set_conf('recursive_i_cor', True)
        
        # Is i_cor applied in the atomic physic database AND model database?
        # If this is the case, i_cor on phyat_database will be directly
        # applied on i_rel and will not appear as i_cor in the printed liste
        # of lines or with the cursor.
        self.set_conf('do_icor_outside_cosmetik', True)
        
        # If you want to perform cosmetik on reference lines (which have ref = 0):
        self.set_conf('do_icor_on_ref', True)

        # Here follow caracteristics of the reference line.
        # This line will be assumed to have a flux
        # at center of 1.00/A.  NO!!!
        
        self.set_conf('do_calcul_aire_ref', False)
        self.set_conf('raie_ref ', {"vitesse" : 25.0, "lambda" : 4861.0, "profile" : 1})     # depuis 25/10/01
            
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
            log_.error('No phyat file read from {0}/{1}'.format(dir_, self.phyat_file), calling = self.calling)
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
        
    def read_cosmetik(self, cosmetik_file, do_cosmetik):
        cosmetik_arr = []
        if cosmetik_file is None:
            return cosmetik_arr
        if not os.path.isabs(cosmetik_file):
            cosmetik_file = '{0}'.format(self.directory + cosmetik_file)
        if do_cosmetik:
            try:
                cosmetik_arr = read_data(cosmetik_file)
                log_.message('cosmetik read from {0}'.format(cosmetik_file), calling = self.calling)
            except:
                log_.warn('unable to read from {0}'.format(cosmetik_file), calling = self.calling)
        return cosmetik_arr

    def read_obs(self, k_spline = 1):
        
        if self.get_conf('spectr_obs') is None:
            n_pix = (self.limit_sp[1] - self.limit_sp[0]) / self.conf['lambda_pix']
            self.w = np.linspace(self.limit_sp[0], self.limit_sp[1], n_pix)
            self.f = np.ones_like(self.w)
        else:
            if self.conf['spectr_obs'].split('.')[-1] == 'fits':
                from astropy.io import fits
                try:
                    self.f, header = fits.getdata(self.conf['spectr_obs'], header=True)
                    dispersion_start = header['CRVAL1'] - (header['CRPIX1'] - 1) * header['CDELT1']
                    self.w = dispersion_start + np.arange(len(self.f)) * header['CDELT1']
                except:
                    log_.warn('Observations NOT read from {0}'.format(self.conf['spectr_obs']),
                                    calling = self.calling)
                    self.n_lambda = 0.
                    return None
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
                    self.w = None
                
                if bool(self.get_conf('reverse_spectra', undefined=False)) :
                        self.f = self.f[::-1]
        if self.get_conf('wave_unit') == 'mu':
            self.w *= 10000.
        self.n_lambda = len(self.f)
        self.tab_pix = np.arange(self.n_lambda)
        
        self.f *= self.get_conf('sp_norm', undefined = 1.)
        
        if ("cal_lambda" in self.conf) and ("cal_pix" in self.conf) and (self.w is None):
            cal_lambda = np.array(self.conf["cal_lambda"])
            cal_pix = np.array(self.conf["cal_pix"])
            arg_sort = cal_lambda.argsort()
            cal_lambda = cal_lambda[arg_sort]
            cal_pix = cal_pix[arg_sort]
            interp_lam = interpolate.UnivariateSpline(cal_pix, cal_lambda, k=k_spline)
            self.w = interp_lam(self.tab_pix)
            log_.message('Wavelength table generated using spline of order {0}'.format(k_spline),
                                calling = self.calling)
            if bool(self.get_conf('reverse_spectra', undefined=False)) :
                    self.f = self.f[::-1]
                        
        if self.limit_sp[0] < 0.01:
            self.limit_sp[0] = np.min(self.w) * (1. + self.get_conf("delta_limit_sp")/100.)
        if self.limit_sp[1] > 0.9e10:
            self.limit_sp[1] = np.max(self.w) * (1. - self.get_conf("delta_limit_sp")/100.)
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
        
        if self.conf["cont_hi_i"]  != 0.:
            alfa = 1e-13 * 0.668 * (self.conf["cont_hi_t"]/1e4)**(-0.507) / \
                (1. + 1.221*(self.conf["cont_hi_t"]/1e4)**(0.653)) * 1.000 
            emis_Hi = alfa * CST.HPLANCK * CST.CLIGHT * 1e8 / 4861.3 # erg/s.cm3
            H_cont = self.conf["cont_hi_i"] * make_cont_Ercolano(self.conf["cont_hi_t"],'H',self.w) / emis_Hi 
            H_cont[~np.isfinite(H_cont)] = 0.
            self.conts['H'] = H_cont
        else:
            self.conts['H'] = np.zeros_like(self.w)
        
        if self.conf["cont_hei_i"] != 0.0:
            alfa = 1e-13 * 0.331 * (self.conf["cont_hei_t"]/1e4)**(-0.615) / \
                (1. + 0.910*(self.conf["cont_hei_t"]/1e4)**(0.780)) * 0.7986
            emis_Hei = alfa * CST.HPLANCK * CST.CLIGHT * 1e8 / 4471.5
            He1_cont = self.conf["cont_hei_i"] * make_cont_Ercolano(self.conf["cont_hei_t"],'He1',self.w) / emis_Hei 
            He1_cont[~np.isfinite(He1_cont)] = 0.
            self.conts['He1'] = He1_cont
        else:
            self.conts['He1'] = np.zeros_like(self.w)

        if self.conf["cont_heii_i"] != 0.0:
            alfa = 2. * 1e-13 * 1.549 * (self.conf["cont_heii_t"]/1e4/4.)**(-0.693) / \
                (1. + 2.884*(self.conf["cont_heii_t"]/1e4/4.)**(0.609))*1.000
            emis_Heii = alfa * CST.HPLANCK * CST.CLIGHT * 1e8 / 4685.8
            He2_cont = self.conf["cont_heii_i"] * make_cont_Ercolano(self.conf["cont_heii_t"],'He2',self.w) / emis_Heii 
            He2_cont[~np.isfinite(He2_cont)] = 0.
            self.conts['He2'] = He2_cont
        else:
            self.conts['He2'] = np.zeros_like(self.w)
                
        gff_HI = gff(1., self.conf["cont_hi_t"], self.w)
        gff_HeI = gff(1., self.conf["cont_hei_t"], self.w)
        gff_HeII = gff(4., self.conf["cont_heii_t"], self.w)
        
        #32.d0*!phy.e^4.*!phy.h/3./!phy.m_e^2./!phy.c^3.*sqrt(!dpi*13.6*!phy.erg_s_ev/3./!phy.k)= 6.8391014e-38

        if self.conf["cont_hi_i"] != 0 and self.conf["cont_hei_i"] != 0 and self.conf["cont_heii_i"] != 0 :
            FF_cont = (6.8391014e-38 * CST.CLIGHT * 1e8 / self.w**2. * (
                        self.conf["cont_hi_i"] * 1.0**2. / np.sqrt(self.conf["cont_hi_t"]) * np.exp(-CST.HPLANCK*CST.CLIGHT*1e8/self.w/CST.BOLTZMANN/self.conf["cont_hi_t"]) * gff_HI/emis_Hi + 
                        self.conf["cont_hei_i"] * 1.0**2./ np.sqrt(self.conf["cont_hei_t"]) * np.exp(-CST.HPLANCK*CST.CLIGHT*1e8/self.w/CST.BOLTZMANN/self.conf["cont_hei_t"]) * gff_HeI/emis_Hei  + 
                        self.conf["cont_heii_i"] * 2.0**2. / np.sqrt(self.conf["cont_heii_t"]) * np.exp(-CST.HPLANCK*CST.CLIGHT*1e8/self.w/CST.BOLTZMANN/self.conf["cont_heii_t"]) * gff_HeII / emis_Heii))
            FF_cont[~np.isfinite(FF_cont)] = 0.
            self.conts['FF'] = FF_cont
        else:
            self.conts['FF'] = np.zeros_like(self.w)

        
        # 2-photons
        #http://adsabs.harvard.edu/abs/1984A%26A...138..495N
        if self.conf["cont_hi_i"] != 0:
            y = 1215.7 / self.w
            A = 202.0 * (y * (1. - y) * (1. -(4. * y * (1 - y))**0.8) + 0.88 * ( y * (1 - y))**1.53 * (4. * y * (1 - y))**0.8)
            alfa_eff = 0.838e-13 * (self.conf["cont_hi_t"] / 1e4)**(-0.728) # fit DP de Osterbrock
            q = 5.31e-4 * (self.conf["cont_hi_t"] / 1e4)**(-0.17) # fit DP de Osterbrock
            n_crit = 8.226 / q
            twophot_cont = self.conf["cont_hi_i"] * CST.HPLANCK * CST.CLIGHT * 1e8 / self.w**3. * 1215.7 * A / 8.226 * alfa_eff / (1. + self.conf["cont_edens"]/n_crit) / emis_Hi
            twophot_cont[~np.isfinite(twophot_cont)] = 0.
            self.conts['2photons'] = twophot_cont
        else:
            self.conts['2photons'] = np.zeros_like(self.w)
            
        self.cont = np.zeros_like(self.w)  
        for key in self.conts:
            self.cont += self.conts[key]
        self.cont *= self.aire_ref
        self.cont /= self.red_corr
        
    def plot_conts(self, ax):
        if self.sp_synth_lr is None:
            return
        colors = {'bb': 'cyan', 'pl': 'green', '2photons': 'blue', 'FF': 'red',
                  'H': 'red', 'He1': 'green', 'He2': 'blue', 'user': 'black'}        
        labels = {'bb': 'bb', 'pl': 'pl', '2photons': '2q', 'FF': 'ff',
                  'H': 'H I', 'He1': 'He I', 'He2': 'He II', 'user': 'interpol'}
        for key in self.conts:
            if key[0] == 'H':
                style=':'
            else:
                style = '-'
            ax.plot(self.w, self.conts[key], linestyle=style, label = labels[key], color = colors[key])
        ax.plot(self.w, self.cont, label = 'total cont', linestyle='--', linewidth = 2)
        ax.legend()
                    
        
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

    def make_synth_test(self, liste_raies):

        sp_theo = self.sp_theo.copy()
        
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
                    log_.message('Area of {0} {1} could be wrong'.format(raie['id'], raie['lambda']), 
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
        
        log_.message('Number of theoretical spectra: {0}'.format(len(sp_theo['correc'])), calling=self.calling)
        return sp_theo, sp_synth

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
                    log_.message('Area of {0} {1} could be wrong'.format(raie['id'], raie['lambda']), 
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
        
        log_.message('Number of theoretical spectra: {0}'.format(len(sp_theo['correc'])), calling=self.calling)
        return sp_theo, sp_synth
        
    def make_sp_abs(self, sp_theo):
        
        if sp_theo is None:
            return None
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
                self.conf['shift_atm'] = (self.conf['shift_atm'],)
            if len(self.get_conf('fic_atm')) != len(self.get_conf('coeff_atm')):
                log_.error('fic_atm number {} != coeff_atm number {}'.format(len(self.get_conf('fic_atm')), len(self.get_conf('coeff_atm'))), 
                                      calling = self.calling)
            for fic_atm, coeff_atm, shift_atm in zip(self.get_conf('fic_atm'), self.get_conf('coeff_atm'), self.get_conf('shift_atm')):
                try:
                    d = np.genfromtxt(fic_atm, dtype=None, names=('wl', 'abs'))
                    d['wl'] = vactoair(d['wl'], self.conf['vactoair_inf'], self.conf['vactoair_sup'])
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
        
        if self.sp_synth is None:
            self.filter_ = None
            return None
        
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
        
        if sp_synth is None:
            return None
        
        input_arr = (cont + sp_synth) * self.sp_abs 
        kernel = self.filter_
        sp_synth_tot = convol(input_arr, kernel)
        return sp_synth_tot
        
    def rebin_on_obs(self):
        if self.sp_synth_tot is None:
            return None, None
        
        resol = self.get_conf('resol', undefined = 1, message=None)
        cont_lr = rebin(self.cont, resol)
        sp_synth_lr = rebin(self.sp_synth_tot, resol)
        return cont_lr, sp_synth_lr 
                
    def adjust(self):
        spectr0 = self.sp_theo['spectr'].copy()
        new_model_arr = self.read_model(self.fic_model)        
        new_cosmetik_arr = self.read_cosmetik(self.fic_cosmetik, self.do_cosmetik)
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
        if mask_diff.sum() > 0 and len(mask_diff) == len(self.liste_raies):
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
            self.sp_theo['spectr'] = spectr0
            for i_change in np.arange(len(new_sp_theo['raie_ref'])):
                to_change = (self.sp_theo['raie_ref']['num'] == new_sp_theo['raie_ref'][i_change]['num'])
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
    def print_line(self, line, sort='lambda', reverse=False):
        if type(line) == np.core.records.recarray:
            sorts = np.argsort(line[sort])
            if reverse:
                sorts = sorts[::-1]
            for isort in sorts:
                self.print_line(line[isort])
            return
        print('{0[num]:>14d} {0[id]:9s}{0[lambda]:11.3f}{0[l_shift]:6.3f}{0[i_rel]:10.3e}{0[i_cor]:7.3f}'\
              ' {0[ref]:>14d}{0[profile]:5d}{0[vitesse]:7.2f}{1:1s}'.format(line, line['comment'].strip()))
       
    def get_line_info(self, line_num, sort='lambda', reverse=False):
        line = None
        refline = None
        satellites = None
        refline_num = -1
        to_select = (self.liste_raies['num'] == line_num)
        if to_select.sum() > 0:
            line = self.liste_raies[to_select][0]
            refline_num = line['ref']
        if line is None:
            refline_num = line_num
        to_select = (self.sp_theo['raie_ref']['num'] == refline_num)
        if to_select.sum() > 0:
            refline = self.sp_theo['raie_ref'][to_select][0]
            to_select = (self.liste_raies['ref'] == refline_num)
            satellites = self.liste_raies[to_select]
            order = np.argsort(satellites[sort])
            satellites = np.array(satellites)[order]
        return line, refline, satellites

    def read_satellites(self, filename, refline_num):
        with open(filename, 'r') as f:
            satellites = []
            for eachline in f:
                if int(self.fieldStrFromLine(eachline,'ref')) == refline_num:
                    satellites.append(eachline)
        return satellites

    def cosmetic_line_unchanged(self, line_c):
        if line_c == None:
            return None
        line_num = int(self.fieldStrFromLine(line_c,'num'))
        line = self.read_line(self.phyat_file, line_num)
        if line == None:
            log_.warn('Error in cosmetic file: line {0:} does not exist in the atomic database\n'.format(str(line_num)), calling=self.calling)
            return None
        else:    
            line = line.rstrip()
            keys = ['l_shift', 'i_cor', 'i_rel', 'profile', 'vitesse']
            v0 = {i: float(self.fieldStrFromLine(line, i)) for i in keys}
            v1 = {i: float(self.fieldStrFromLine(line_c, i)) for i in keys}
            if v0 == v1:
                return True
            else:
                return False

    def cosmetic_line_ok(self, line_c):
        if line_c == None:
            return None
        line_num = int(self.fieldStrFromLine(line_c,'num'))
        line = self.read_line(self.phyat_file, line_num)
        if line == None:
            log_.warn('Error in cosmetic file: line {0:} does not exist in the atomic database\n'.format(str(line_num)), calling=self.calling)
            return None
        else:    
            line = line.rstrip()
            keys = [ 'lambda', 'i_rel' ]
            v0 = {i: float(self.fieldStrFromLine(line, i)) for i in keys}
            v1 = {i: float(self.fieldStrFromLine(line_c, i)) for i in keys}
            if v0['i_rel'] != v1['i_rel'] or v0['lambda'] != v1['lambda']:
                log_.warn('Error in cosmetic file for line {}\n'.format(str(line_num)), calling=self.calling)
                log_.warn('(cosmetic)   ' + line_c, calling=self.calling)
                log_.warn('(database)   ' + line, calling=self.calling)
                return False
            else:
                return True

    def read_line(self, filename, line_num):
        line = None
        line_num_str = str(line_num)
        k = len(line_num_str)
        if not os.path.isfile(filename):
            return None
        else:
            with open(filename, 'r') as f:
                line = None
                for eachline in f:
                    s = self.fieldStrFromLine(eachline,'num')
                    #s = s.strip()
                    #if int(s) == line_num:
                    if s[-k:] == line_num_str:
                        line = eachline
                        break
        return line

    def replace_field(self, line, field, value):
        w = self.field_width[field]
        if len(value) > w:
            return None
        elif len(value) < w:
            a = self.field_align[field]
            value = '{:{a}{w}}'.format(value)
        j = self.field_pos[field]
        k = j + w
        line = line[:j] + value + line[k:]
        return line

    def remove_line(self, filename, line_num):
        line = self.read_line(filename, line_num)
        if line == None:
            return False
        if not os.path.isfile(filename):
            return False
        else:
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            i = lines.index(line)
            if i >= 0:
                del lines[i]
                with open(filename, 'w') as f:
                    f.writelines(lines)
                return True
            else:
                return False        

    def replace_line(self, filename, line):
        line_num = int(self.fieldStrFromLine(line,'num'))
        if os.path.isfile(filename):
            lineNotFound = True
            with open(filename, 'r') as f:
                lines = f.read().splitlines()
                for i in range(0, len(lines)):
                    curr_line = lines[i]
                    if int(self.fieldStrFromLine(curr_line,'num')) == line_num:
                        lines[i] = line + '\n'
                        lineNotFound = False
                    else:
                        lines[i] = lines[i] + '\n'
            if lineNotFound:    
                lines.append(line)
        else:
            lines = [line]                
        with open(filename, 'w') as f:
            f.writelines(lines)

    def fieldStrFromLine(self, lineOfFile, field):
        if lineOfFile == None:
            return None
        fieldStr = None
        if field in self.fields:
            i = self.field_pos[field]
            j = i+self.field_width[field]
            fieldStr = lineOfFile[i:j]
        return fieldStr

    def line_info(self, line_num, sat_info=True, print_header=True, sort='lambda', reverse=False):
        
        if print_header:
            print('\n{0:-^45}'.format(' INFO LINES '))           
        if type(line_num) == type(()) or type(line_num) == type([]):
            for line in line_num:
                self.line_info(line, sat_info=sat_info, print_header=False)
            return 
        to_print = (self.liste_raies['num'] == line_num)
        if to_print.sum() == 1:
            raie = self.liste_raies[to_print][0]
            self.print_line(raie)
            if raie['ref'] != 0 and sat_info:
                print('\nSatellite line of:')
                self.line_info(raie['ref'], print_header=False)
                
        to_print = (self.sp_theo['raie_ref']['num'] == line_num)
        if to_print.sum() > 0 and sat_info:
            raie = self.sp_theo['raie_ref'][to_print][0]
            self.print_line(raie)
            print('')
            satellites_tab = (self.liste_raies['ref'] == raie['num'])
            Nsat = satellites_tab.sum()
            if Nsat > 0:
                print('{0} satellites'.format(Nsat))
                self.print_line(self.liste_raies[satellites_tab], sort=sort, reverse=reverse)
            if self.sp_theo['correc'][to_print][0] != 1.0:
                print('Intensity corrected by {0}'.format(self.sp_theo['correc'][to_print][0]))
        if print_header:
            print('-'*45)

    def get_ref_list(self, ions):
        ref_list = []
        for ion in ions:
            i_ion = np.where(self.sp_theo['raie_ref']['id'] == ion.ljust(9))[0]
            if len(i_ion) == 0:
                ref_list.append(-1)
            for i in i_ion:
                ref_list.append(self.sp_theo['raie_ref'][i][0])
        return ref_list

    def get_element_and_int_ion(self, ion):
        
        def roman_to_int(s):
            x = {'C': 100, 'L': 50, 'X': 10, 'V': 5, 'I': 1}
            s = s.strip()
            if len(s) == 0:
                n = -1
            elif len(s) == 1:
                if s == 'I':
                    n = 1
                elif s == 'V':
                    n = 5
                elif s == 'X':
                    n = 10
            else:
                n = sum([x[i] if x[i] >= x[j] else -x[i] for i, j in zip(s, s[1:])]) + x[j]
            return n
        
        ion = self.true_ion(ion)
        k = ion.find('_')
        element = ion[:k]
        int_ion = roman_to_int(ion[k+1:])
        return element, int_ion

    def set_selected_ions_data(self):

        pos_label = 0
        pos_ion = 1
        pos_ref = 2
        pos_i_ion = 3
        pos_color = 4
        pos_linestyle =5
        
        selected_ions = self.get_conf('selected_ions')
        colors = self.get_conf('color_selected_ions')
        linestyles = [ 'solid', 'dashed', 'dashdot', 'dotted' ]
        label_list = []
        ref_list = []
        for ion in selected_ions:
            i_ion = np.where(self.sp_theo['raie_ref']['id'] == ion.ljust(9))[0]

            if len(i_ion) == 0:
                i_ion_raies = np.where(self.liste_raies['id'] == ion.ljust(9))[0]
                if len(i_ion_raies) > 0:
                    ref_list = self.liste_raies[i_ion_raies]['ref']     
                    ref_list = list(set(ref_list))
                    i_ion = []
                    for i in ref_list:
                        i_ion = i_ion + list(np.where(self.sp_theo['raie_ref']['num'] == i)[0])
                    i_ion = list(set(i_ion))

            if self.get_conf('diff_lines_by') == 1:
                ref_list = []
                for i in range(0, len(i_ion)):
                    ref_line = self.sp_theo['raie_ref'][i_ion[i]]['num']
                    ref_list.append(ref_line)
                k = len(label_list)
                color = colors[k%len(colors)]
                linestyle = linestyles[(k/len(colors))%len(linestyles)]
                label_list.append([ion, [ion], ref_list, list(i_ion), color, linestyle])
            else:
                for i in range(0, len(i_ion)):
                    k = len(label_list)
                    color = colors[k%len(colors)]
                    linestyle = linestyles[(k/len(colors))%len(linestyles)]
                    if len(i_ion) == 1:
                        label = ion
                    else:
                        label = ion + '-' + str(i)
                    ref_line = self.sp_theo['raie_ref'][i_ion[i]]['num']
                    refline_str = str(ref_line).strip('0')
                    #refline_str = str(ref_line)[:-6]
                    label = label + ' (' + refline_str + ')'
                    label_list.append([label, [ion], [ref_line], list(i_ion[i:i+1]), color, linestyle])
                
        # sorting
        for i in range(0,len(label_list)-1):
            label1 = label_list[i][pos_label]
            true_ion1 = self.true_ion(label1)
            element1, int_ion1 = self.get_element_and_int_ion(true_ion1)
            for j in range(i+1, len(label_list)):
                label2 = label_list[j][pos_label]
                true_ion2 = self.true_ion(label2)
                element2, int_ion2 = self.get_element_and_int_ion(true_ion2)
                if (element2 < element1) or ((element2 == element1) and (int_ion2 < int_ion1)) or ((true_ion2 == true_ion1) and (label2 < label1)): 
                    prov = label_list[i]
                    label_list[i] = label_list[j]
                    label_list[j] = prov
                    label1 = label2
                    true_ion1 = true_ion2
                    element1 = element2
                    int_ion1 = int_ion2

        if self.get_conf('diff_lines_by') == 2:
            i = 0
            while i < len(label_list)-1:
                ion = label_list[i][pos_ion]
                j = i+1
                while j < len(label_list):
                    ion2 = label_list[j][pos_ion]
                    if self.true_ion(ion2[0]) == self.true_ion(ion[0]):
                        ion = list(set(ion + ion2))
                        ref_list = list(set(label_list[i][pos_ref] + label_list[j][pos_ref]))
                        i_ion = list(set(label_list[i][pos_i_ion] + label_list[j][pos_i_ion]))
                        label_list.pop(j)
                        label_list[i][pos_ion] = ion
                        label_list[i][pos_ref] = ref_list
                        label_list[i][pos_i_ion] = i_ion
                    else:
                        j += 1
                i += 1

            for i in range(len(label_list)):
                ion_list = label_list[i][pos_ion]
                ion_list.sort()
                ion = ion_list[0]
                for j in range(1, len(ion_list)):
                    ion = ion + '+' + ion_list[j]
                label_list[i][pos_label] = ion
        self.selected_ions_data = label_list
        return

    def get_refline_lists(self, ions):
        ref_code_list = []
        ref_index_list = []
        ref_label_list = []
        for ion in ions:
            ion = self.true_ion(ion)
            i_ion = np.where(self.sp_theo['raie_ref']['id'] == ion.ljust(9))[0]
            if len(i_ion) == 0:
                ref_code_list.append(-1)
                ref_index_list.append(-1)
                ref_label_list.append(ion.replace('_',' ') + ' (no lines)')
            for i in i_ion:
                ref_code_list.append(self.sp_theo['raie_ref'][i][0])
                ref_index_list.append(i)
                if len(i_ion) == 1:
                    ref_label_list.append(ion.replace('_',' '))
                else:
                    ref_label_list.append(ion.replace('_',' ')+ ' - ' + str(np.where(i_ion==i)[0][0]))
#              ref_label_list.append(ion.replace('_',' ')+ ' - ' + str(self.sp_theo['raie_ref'][i][0])[:-8])
        return ref_code_list, ref_index_list, ref_label_list

    def get_ref_index_list(self, ions):
        ref_index_list = []
        for ion in ions:
            ion = self.true_ion(ion)
            i_ion = np.where(self.sp_theo['raie_ref']['id'] == ion.ljust(9))[0]
            for i in i_ion:
                ref_index_list.append(i)
        return ref_index_list

    def get_refline_from_code(self, code_str):
        s = ''
        for line in self.liste_raies:
            if code_str == str(line[0]):
                s = line['ref']
        return s
        
    def get_ion_from_code(self,code_str):
        s = ''
        for line in self.liste_raies:
            if code_str == str(line[0]):
                s = line['id']
        return s
    
    def true_ion(self, ion):
        while len(ion) > 0 and ion[-1] not in [ 'I', 'V', 'X' ]:
            ion = ion[:-1]  
        return ion
                
    def get_all_ions_from_ion(self, ion):
        ion = self.true_ion(ion)
        ion_list = [ion]
        k = len(ion)
        for line in self.liste_raies:
            s = str(line['id'])
            if len(s) > k and s[:k] == ion and s[k] not in [ 'I', 'V', 'X' ]:
                ion_list.append(s.strip())
        return list(set(ion_list))
        
    def get_ions_from_element(self, elem):
        ion_list = []
        for line in self.liste_raies:
            s = str(line['id'])
            k = s.find('_')
            if elem == s[:k]:
                ion_list.append(s.strip())
        return list(set(ion_list))

    def save_lines(self):
        if self.get_conf('show_selected_intensities_only'):
            cut = self.get_conf('cut_plot2')
        else:
            cut = 0.0
        ref_list = self.get_ref_list(self.get_conf('selected_ions'))        
        sort_list = [ 'lambda', 'i_rel', 'id' ]
        k = self.get_conf('line_saved_ordered_by')
        sort = sort_list[k/2]
        filename = self.get_conf('line_saved_filename')
        extension = os.path.splitext(filename)[1][1:].lower()
        sep = ' '
        end = '\n'
        if extension == 'tex':
            sep = ' & '
            end = ' {0}{0}{1}'.format('\\', '\n')
        elif extension == 'csv':
            sep = ' ; '
            end = '\n'
        sorts = np.argsort(self.liste_raies[sort])
        if k%2 == 1:
            sorts = sorts[::-1]
        with open(filename, 'w') as f:
            field_print = self.get_conf('line_field_print')
            n = len(field_print)
            
            if self.get_conf('line_saved_header'):
                s = ''
                for item in field_print:
                    f.write('{0:9s} : {1:>}\n'.format(item, self.field_tip[item]))
                  
                for item in field_print:
                    width = self.field_width[item]
                    thisformat = '{{:{a}{w}s}{}}'
                    if ( item == field_print[n-1] ):
                        add_s = end
                    else:
                        add_s = sep
                    if ( item == field_print[0] ):
                        align = ''
                    else:
                        align = '>'
                    s = s + str('{:{a}{w}s}{}'.format(item, add_s, a=align, w=width))
                f.write('\n'+s+'\n')
            
            for i_sort in sorts:
                line = self.liste_raies[i_sort]
                wl = line['lambda'] + line['l_shift'] + self.conf['lambda_shift']
                i_rel = line['i_rel']
                i_tot = line['i_rel'] * line['i_cor']
                #if (abs(i_rel) > cut) and ( not self.get_conf('show_selected_ions_only') or line['ref'] in ref_list):
                if (abs(i_tot) > cut) and ( not self.get_conf('show_selected_ions_only') or line['ref'] in ref_list):
                    s = ''
                    n = len(field_print)
                    for item in field_print:
                        thisformat = self.field_format[item]
                        if item == 'l_tot':
                            r = wl
                        elif item == 'i_tot':
                            r = i_tot
                        else:
                            r = line[item]
                        if item == 'l_shift': 
                            s = s + ' '
                        s = s + str(thisformat.format(r))
                        if ( item == field_print[n-1] ):
                            s = s + end
                        else:
                            s = s + sep
                    f.write(s)

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
            self.set_conf('cut_plot2', cut)   
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
        
        if self.sp_synth_lr is None:
            return
        self.ax1_line_synth = ax.step(self.w_ori, self.sp_synth_lr, label='Synth', c='blue', linewidth=2)[0]
        if self.hr:
            ax.step(self.w, self.sp_synth, c='green')

        selected_ions = self.get_conf('selected_ions')
        self.set_selected_ions_data()
        label_list = self.selected_ions_data
        if selected_ions != [] and self.get_conf('plot_lines_of_selected_ions'):
            j = self.get_conf('index_of_current_ion')
            if j in range(0, len(selected_ions)):
                ions = []
                ions.append( selected_ions[j] )
            else:
                ions = selected_ions

            if j in range(0, len(label_list)):
                label_list = label_list[j:j+1]

            for item in label_list:
                label = item[0].replace('_',' ')
                i_ion = item[3]
                color = item[4]
                linestyle = item[5]
                y = 0
                for i in range(0,len(i_ion)):
                    y = y + self.sp_theo['spectr'][i_ion][i]
                ax.step(self.w, self.cont+y, c=color, label=label, linestyle=linestyle )[0]

        ax.legend(loc=self.legend_loc)
        log_.debug('ax1 drawn on ax ID {}'.format(id(ax)), calling=self.calling)
    
    def plot_line_ticks(self, ax, y1, y2, wmin=0., wmax=20000.):

        if self.sp_synth_lr is None:
            return
        lcolor = self.get_conf('line_tick_color')
        label_list = self.selected_ions_data
        j = self.get_conf('index_of_current_ion')
        if j in range(0, len(label_list)):
            label_list = label_list[j:j+1]

        for line in self.liste_raies:
            wl = line['lambda'] + line['l_shift'] + self.conf['lambda_shift']
            i_rel = line['i_rel']
            i_tot = line['i_rel'] * line['i_cor']
            if (wmin < wl) and (wl < wmax) and ((abs(i_tot) > self.get_conf('cut_plot2')) or ( not self.get_conf('show_selected_intensities_only'))): 
                if not self.get_conf('show_selected_ions_only'):
                    ax.axvline( wl, ymin=y1, ymax=y2, color = lcolor, linestyle = 'solid' )
                refline = line['ref']
                ion = line['id'].strip()
                for item in label_list:
                    ion_list = item[1]
                    ref_list = item[2]
                    color = item[4]
                    linestyle = item[5]
                    if refline in ref_list and ion in ion_list:
                        ax.axvline( wl, ymin=y1, ymax=y2, color = color, linestyle = linestyle, linewidth = 1.5 )
        
        # To add ticks to the legend of the figure when the spectrum of the selected ions are not plotted
        if not self.get_conf('plot_lines_of_selected_ions'):
            for item in label_list:
                label = item[0].replace('_',' ')
                color = item[4]
                linestyle = item[5]
                ax.step( [0,0], [0,100], color = color, linestyle = linestyle, label = label )
            ax.legend(loc=self.legend_loc)
        log_.debug('Lines drawn on ax ID {}'.format(id(ax)), calling=self.calling)

    def plot_ax3(self, ax):     
        if self.sp_synth_lr is not None:
            ax.step(self.w, self.f - self.cont, c = 'red', linestyle='--')
            ax.plot((0, 1e10), (0.0, 0.0), c='green')
            ax_line_diff = ax.step(self.w_ori, self.f_ori - self.sp_synth_lr, c='blue')[0]
            log_.debug('ax3 drawn on ax ID {}'.format(id(ax)), calling=self.calling)
        
    def update_plot2(self):
        
        if self.ax1 is None:
            return
        if self.sp_synth_lr is None:
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
            if (abs(i_rel) > self.get_conf('cut_plot2')) & (wl > self.ax2.get_xlim()[0]) & (wl < self.ax2.get_xlim()[1]):
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
            if self.sp_synth_lr is None:
                self.y1_plot_lims = (np.min(self.f), np.max(self.f)) 
            else:
                mask = (self.w_ori > self.x_plot_lims[0]) & (self.w_ori < self.x_plot_lims[1])
                self.y1_plot_lims = (np.min(self.sp_synth_lr[mask]), np.max(self.sp_synth_lr[mask]))      
        
        self.y2_plot_lims = self.get_conf('y2_plot_lims')
        if self.y2_plot_lims is None:
            self.y2_plot_lims = (-0.5, 1,5)
        
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
        #log_.message('Axes saved', calling=self.calling)
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
    
    def apply_post_proc(self):
        
        if self.post_proc_file is not None and self.post_proc_file is not "":
            try:
                user_module = {}
                execfile(os.path.abspath(self.directory)+'/'+self.post_proc_file, user_module)
                self.post_proc = user_module['post_proc']
                log_.message('function post_proc read from {}'.format(self.post_proc_file), calling=self.calling)
            except:
                self.post_proc = None
                log_.warn('function post_proc NOT read from {}'.format(self.post_proc_file), calling=self.calling)
            
            if self.post_proc is not None:
                self.post_proc(self.fig1)
        
    def rerun(self):
        self.run(do_synth = True, do_read_liste = True, do_profiles=True)

    def replot2(self):    
        self.save_axes()
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.plot2(hr=self.hr, cut=self.get_conf('cut_plot2'), split=self.split, 
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

    def get_nearby_lines(self, w1, w2, do_print=True, sort='i_tot', reverse=True):
        if w1 == None or w2 == None:
            return  None
        w = (w1 + w2)/2
        w_lim = abs(w2 - w1)/2
        tt = (np.abs(self.liste_raies['lambda'] + self.liste_raies['l_shift'] + self.conf['lambda_shift'] - w) < w_lim)
        nearby_lines = self.liste_raies[tt]
        i_tot = nearby_lines['i_rel']*nearby_lines['i_cor']
        if sort == 'i_tot':
            sorts = np.argsort(i_tot)
        else:
            sorts = np.argsort(nearby_lines[sort])
        if reverse:
            sorts = sorts[::-1]
        nearby_lines = np.array(nearby_lines)[sorts]
        if tt.sum() > 0:
            if do_print:
                print('\n{0:-^45}'.format(' CURSOR on {0:.3f} '.format(w)))
                self.print_line(self.liste_raies[tt])
                print('-'*45)
        return nearby_lines

                    
    def nearby_lines(self, event, do_print=True, sort='i_tot', reverse=True):
        
        nearby_lines = None
        w = event.xdata
        try:
            if (w > self.ax1.get_xlim()[1]) or (w < self.ax1.get_xlim()[0]) or (event.button == 2):
                self._cursOff()
                return None
        except AttributeError:
            log_.warn('ax1 not defined', calling=self.calling)
            return None
        try:
            if event.button in (1,3):
                if self.firstClick:
                    self.cursor_w0 = w
                    self.firstClick = False
                else:
                    self.cursor_w1 = self.cursor_w0
                    self.cursor_w2 = w
                    self.firstClick = True
                    w = (self.cursor_w1+self.cursor_w2)/2
                    w_lim = self.cursor_width * (self.ax1.get_xlim()[1] - self.ax1.get_xlim()[0])
                    if abs(self.cursor_w2-w) > w_lim/10:
                        w_lim = abs(self.cursor_w2-w)
                    self.cursor_w1 = w - w_lim
                    self.cursor_w2 = w + w_lim
                    nearby_lines = self.get_nearby_lines(self.cursor_w1, self.cursor_w2, do_print=do_print)
                    return nearby_lines
                    
        except AttributeError:
            log_.warn('ax1 not defined', calling=self.calling)
            return None
                    
                    
    def _curs_onclick(self, event):
        
        wl = event.xdata
        try:
            if (wl > self.ax1.get_xlim()[1]) or (wl < self.ax1.get_xlim()[0]) or (event.button == 2):
                self._cursOff()
                return None
        except AttributeError:
            log_.warn('ax1 not defined', calling=self.calling)
            return None
        try:
            if event.button in (1,3):
                wl_lim = self.cursor_width * (self.ax1.get_xlim()[1] - self.ax1.get_xlim()[0])
                tt = (np.abs(self.liste_raies['lambda'] + self.liste_raies['l_shift'] + self.conf['lambda_shift'] - wl) < wl_lim)
                if tt.sum() > 0:
                    print('\n{0:-^45}'.format(' CURSOR on {0:.3f} '.format(wl)))
                    self.print_line(self.liste_raies[tt])
                    print('-'*45)
        except AttributeError:
            log_.warn('ax1 not defined', calling=self.calling)
            return None

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
    sp = spectrum(config_file=args.file, post_proc_file=args.post_proc)
    fig = plt.figure(figsize=(20, 7))
    sp.plot2(fig=fig)
    sp.apply_post_proc()
    plt.show()
