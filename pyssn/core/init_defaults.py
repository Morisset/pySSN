#-------------------------------------------------------------------
#   Initialisation
#-------------------------------------------------------------------

# All the following variables (most of them boolean) can be
# set here or sometime in the control panel.
# This is a python script that will be executed.

# Verbosity level:
# 1: only errors
# 2: errors and warnings
# 3: errors, warnings and comments
log_level = 2

#-------------------------------------------------------------------
#   Observations
#-------------------------------------------------------------------
# Observation filename, without the extension (which must be .spr[.gz])
# If set to None, limit_sp must be given.
spectr_obs = None 

# Factor applied to the observations
sp_norm = 1.0 #  
# Not sure it works... Better use obj_velo
lambda_shift = 0.0 # w units. 
# Velocity of the object, applied to the wavelengths in the
# calibration of the observations
obj_velo = 0. # km/s

# Only used if no observed spectrum given
lambda_pix = 0.1 

# If the data include the wavelengths (1st column), set the following to True
data_incl_w =  True

# The following is used to define the lambda calibration, only in case of data_incl_w = False
cal_lambda = (0, 512)
cal_pix = (5000, 6700)

# If observed spectra is to be reversed to let the wavelength growing. Used only with cal_lambda and cal_pix
reverse_spectra = False

# Rebinning factor applied to the observations before the synthesis.
# Nouveau mode.
resol = 1  # Resol must be small if the obs resol is large


#-------------------------------------------------------------------
#   Synthesis
#-------------------------------------------------------------------

# If you want to make a synthesis, let it to True.
# Set it to False if you don't want to make a synthesis. In this case an observation must be given and will be displayed
do_synth = True
do_cosmetik = True

# The limit of the spectrum to be computed
limit_sp = [0, 1.0e+10]

# The wavelengths are in "Angstrom" or "mu" 
wave_unit = 'Angstrom'

# Here are the 2 filenames of the data. The filename 
# of the atomic database is hard-coded 
fic_modele = 'liste_modele.dat'
fic_cosmetik = 'liste_cosmetik.dat'

# To take into account 'fic_cosmetik' data
do_cosmetik = True

# If you want to read the list of lines, let it to True.
# Set it to False if you just want to run a new synthesis with 
# the same lines (modification of the profiles, lambda-calibration...)
# Note that if do_synth = False, the following flag has no meaning. 
do_read_liste = True

warn_on_no_cosmetik = True
warn_on_no_reference = True

# Applied to the synthese

# Atmospheric lines in the following files:
fic_atm = None
coeff_atm = 0.0
shift_atm = 0.0

vactoair_inf = 2000.
vactoair_sup = 20000.

# Here is the name of the function used for the emission profiles:
profil_emis_name = 'profil_emis'

# If limit_sp is undefined, then the computing interval is
# the observed one increased by delta_limit_sp % at
# each side.
#
delta_limit_sp =  0.1 # %


# initialisation of the continuum flux 
cont_plot = False #Set True to plot continuum spectrum
cont_in_lambda = False #Set to True if the continuum is in wavelength
cont_unred = True #set to True if reddening is to be applied to cont
cont_lambda = 0.
cont_pix = 0.
cont_intens = 0.
cont_bb_t = 11000.
cont_bb_i = 0.
cont_pl_alpha = 1.
cont_pl_i = 0.
cont_hi_i = 1.e4
cont_hei_i = 3.e2
cont_heii_i = 0.e3
cont_hi_t = 1e4
cont_hei_t = 1e4
cont_heii_t = 1e4
cont_edens = 1e3

# Wavelength of the reference line for the reddening.
lambda_ref_rougi = 4861.3
red_corr_law = 'S79 H83 CCM89'
e_bv = 0.0
# definition of the structure for the limits of the plots
x_plot_lims = None
y1_plot_lims = None
y2_plot_lims = (-0.5, 1.5)
y3_plot_lims = None
cut_plot2 = 1.0 # used to select the line ID shown in plot2

# Color lines
plot_magenta = None
label_magenta = ''
plot_cyan = None
label_cyan = ''

# Plot the second and third panels
plot_ax2 = True
plot_ax3 = True

# Show line ticks 
show_line_ticks = True

# line_tick_ax: 0 - On plot of spectra, 1 - On plot of residuals, 2 - On separate plot
line_tick_ax = 0

# line_tick_pos: 0 - Top, 1 - Middle, 2 - Bottom
line_tick_pos = 0

# Show line ticks for selected ions only
show_selected_ions_only = True

# Show line ticks for selected ions only
show_selected_intensities_only = True

# Differentiate line sequences of a same ion according to the reference
differentiate_lines = True

# Allow editing line parameters in line info dialog
allow_editing_lines = True

# Automatically update synthesis after editing line parameters in line info dialog
update_after_editing_lines = True

# Color of line ticks 
line_tick_color = 'black'

# Plot lines of selected ions
plot_lines_of_selected_ions = True

# Color of selected line ticks 
color_selected_ions = [ 'DeepSkyBlue', 'magenta', 'BlueViolet', 'LightGreen' ]

# List of selected ions 
selected_ions = []

# Index to cycle the selected ions
index_of_selected_ions = -1

# Scale factor for the plot of the theoretical (individual) spectra. Still in devel
fact_multi_synth =  1.

# norm_intens: value by which are divided the line intensities in
# the printed table. Still in devel
norm_intens = 1.00e4

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



# Fields to be printed 
line_field_print   = [ 'id', 'lambda', 'l_shift', 'i_rel', 'i_cor' ]

# Printed header  
line_saved_header = True

# line_saved_ordered_by: 0 - Wavelength, 1 - Intensity, 2 - ion
line_saved_ordered_by = 0

line_saved_format = 'tex' # tex or csv
line_saved_filename = 'lines.dat'

plot_filename = 'plot.pdf'

show_dialogs = True
