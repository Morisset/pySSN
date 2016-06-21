'''
Created on 17/01/2014

@author: morisset
'''
#-------------------------------------------------------------------
#   Initialisation
#-------------------------------------------------------------------

log_level = 2
#
# All the following variables (most of them boolean) can be
# set here, or in the init_file of the model you are doing,
# or sometime in the control panel.
#
   
# Observation filename
# If set to None, limit_sp must be given
spectr_obs = None 
   
# If you want to make a synthesis, let it to True.
# Set it to False if you don't want to make a synthesis 
# 
do_synth = True

# If you want to read the list of lines, let it to True.
# Set it to False if you just want to run a new synthesis with 
# the same lines (modification of the profiles, lambda-calibration...)
# Note that if do_synth = False, the following flag has no meaning. 
do_read_liste = True


# If you DON'T want an i_cor on a main line to affect the satellites, 
# set the following variable to False
recursive_i_cor = True

# Is i_cor applied in the atomic physic database AND model database?
# If this is the case, i_cor on phyat_database will be directly
# applied on i_rel and will not appear as i_cor in the printed liste
# of lines or with the cursor.
do_icor_outside_cosmetik = True

# If you want to perform cosmetik on reference lines (which have ref = 0):
do_icor_on_ref = True

# Want be sure that the listes of lines are correct? set to True 
do_debug_liste = False
warn_on_no_cosmetik = True
warn_on_no_reference = True

# Applied to the synthese
lambda_shift = 0.0 # w units

spectr_obs = None 
# Factor applied to the observations
sp_norm = 1.0 #  
lambda_pix = 0.1 #
limit_sp = [0, 1.0e+10]
# Velocity of the object, applied to the wavelengths in the
# calibration of the observations

obj_velo = 0. # km/s

# If observed spectra is to be reversed to let the wavelength growing
reverse_spectra = False

# If the data include the wavelengths (1st column), set the following to True
data_incl_w =  True

# The following is used to define the lambda calibration, only in case of data_incl_w = False
cal_lambda = (0, 512)
cal_pix = (5000, 6700)

# Here are the 2 filenames of the data. The filename 
# of the atomic database is hard-coded (in read_liste.pro)
fic_cosmetik = 'liste_cosmetik.dat'
fic_modele = 'liste_modele.dat'

# Atmospheric lines in the following files:
fic_atm = None
coeff_atm = 0.0
shift_atm = 0.0

# Here is the format of the listes

# Here follow caracteristics of the reference line.
# This line will be assumed to have a flux
# at center of 1.00/A. 

do_calcul_aire_ref = 0  
raie_ref =  {"vitesse" : 25.0, "lambda" : 4861.0, "profile" : 1}      # depuis 25/10/01

# Here is the name of the function used for the emission profiles:
profil_emis_name = 'profil_emis'

#
#
# If limit_sp is undefined, then the computing interval is
# the observed one increased by delta_limit_sp at
# each side.
#
delta_limit_sp =  10.

#
# Rebinning factor applied to the observations before the synthesis.
# Nouveau mode.
resol = 1  # Resol must be small if the obs resol is large

# Resol should be odd (otherwise, 1 will be automatically added).
# If you change this, verify that all the
# lines with the same profile have the same area
# Look at make_synth.pro, there is a line to uncomment to get 
# the list of the areas of the lines.
# Will be replaced by Resol (as soon as the bug is fixed...)
high_resolution = 1

# initialisation of the continuum flux 
cont_in_lambda = False #set to True if the continuum is in wavelength
cont_unred = True #set to True if redeening is to be applied to cont
cont_lambda = 0.
cont_pix = 0.
cont_intens = 0.
cont_bb_t = 11000.
cont_bb_i = 0.
cont_pl_alpha = 1.
cont_pl_i = 0.
cont_ihi = 1.e4
cont_ihei = 3.e2
cont_iheii = 0.e3
cont_Thi = 1e4
cont_Thei = 1e4
cont_Theii = 1e4
cont_edens = 1e3

# Wavelength of the reference line for the reddening.
lambda_ref_rougi = 4861.3
red_corr_law = 'S79 H83 CCM89'

# definition of the structure for the limits of the plots
x_plot_lims = None
y1_plot_lims = None
y2_plot_lims = (-1.5, 1.)
y3_plot_lims = None

# Color lines
plot_magenta = None
label_magenta = ''
plot_cyan = None
label_cyan = ''

# Plot the second and third panels
plot_ax2 = True
plot_ax3 = True

# Scale factor for the plot of the theoretical (individual) spectra:
fact_multi_synth =  1.

# norm_intens: value by which are divided the line intensities in
# the printed table
norm_intens = 1.00e4

line_select = 0

# instrument
prof = {'largeur':0.50,
    'B_1r':0.00,'B_1l':0.00,'decroiss_1':2.70,'alpha_1':0.50,
    'B_2r':0.00,'B_2l':0.00,'decroiss_2':1.50,'alpha_2':0.75,
    'B_3r':0.00,'B_3l':0.00,'decroiss_3':1.15,'alpha_3':1.00,
    'B_4r':0.00,'B_4l':0.00,'decroiss_4':0.45,'alpha_4':1.50,
    'comment':' Gauss'   }

ghost =  {"do_ghost":0, "delta_lambda" : 0. , "intens" : [ 0.00]}

#ghost =  {do_ghost:1, $
#          delta_lambda : [2.74,-2.75,5.2,-5.2] , $
#          largeur : [0.5,0.5,0.7,0.7],$
#          intens : [ 0.00310000, $
#                      0.000880000, $
#                      3.50000e-05, $
#                      8.00000e-06] }

