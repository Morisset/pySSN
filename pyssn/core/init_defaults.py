#-------------------------------------------------------------------
#   Initialisation
#-------------------------------------------------------------------

# All the following variables (most of them boolean) can be
# set here or sometime in the control panel.
# This is a python script that will be executed.

#-------------------------------------------------------------------
#   Qt (Graphical user interface)
#-------------------------------------------------------------------

# Qt widget style (ex. Plastique, Windows, ...)
qt_style = None

# Enable tooltips of Qt widgets
qt_enable_tooltips = True

# Show line infomation on a modal dialog (True) or on the terminal (False)
qt_show_dialogs = True

# Show plot of residuals
qt_plot_residuals = True

# Allow editing line parameters in line info dialog
qt_allow_editing_lines = True

# Automatically update synthesis after editing line parameters in line info dialog
qt_update_after_editing_lines = True

# The fraction of height reserved for space between subplots
fig_hspace=0.05

# The bottom limit of the figure
fig_bottom=0.08 

# The left limit of the figure
fig_left=0.05

# The right limit of the figure
fig_right=0.995

# The top limit of the figure
fig_top=0.995 

# Automatically adjust limits of the figure
fig_adjust = False

#   Verbosity
#-------------------------------------------------------------------

# Verbosity level:
# 0: None
# 1: Only errors
# 2: Errors and warnings
# 3: Errors, warnings and comments
# 4: Debug messages
log_level = 2

warn_on_no_cosmetik = True
warn_on_no_reference = True

#-------------------------------------------------------------------
#   Observations
#-------------------------------------------------------------------

# Observed spectrum filename
# Supported formats: plain text files, zipped text files (.gz), fits files (.fits)
spectr_obs = None 

# Factor applied to the observations
sp_norm = 1.0 #  
# Not sure it works... Better use obj_velo
lambda_shift = 0.0 # w units. 

# Wavelength correction table
# ex. [(lambda1,delta1), (lambda2,delta2), ... ]
lambda_shift_table = None

# Velocity of the object, applied to the wavelengths in the calibration of the observations
obj_velo = 0. # km/s

# Object velocity interpolation table
# ex. [(lambda1,v1), (lambda2,v2), ... ]
obj_velo_table = None

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
fic_modele = 'liste_modele.dat'
fic_cosmetik = None
fic_profile = None
phyat_file = None

# To take into account 'fic_cosmetik' data
do_cosmetik = True

# If you want to read the list of lines, let it to True.
# Set it to False if you just want to run a new synthesis with 
# the same lines (modification of the profiles, lambda-calibration...)
# Note that if do_synth = False, the following flag has no meaning. 
do_read_liste = True

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
#cont_in_lambda = False #Set to True if the continuum is in wavelength
cont_unred = True #set to True if reddening is to be applied to cont
cont_user_table = None
cont_user_func = 'linear'
cont_pix = 0.
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
r_v = 3.1

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

# Show line ticks 
show_line_ticks = True

# line_tick_ax: 0 - On plot of spectra, 1 - On plot of residuals, 2 - On separate plot
line_tick_ax = 0

# line_tick_pos: 0 - Top, 1 - Middle, 2 - Bottom
line_tick_pos = 0

# line_tick_pos_selectedLine: 0 - Top, 1 - Middle, 2 - Bottom, 3 - Same position (=line_tick_pos), 4 - alternate position
line_tick_pos_selectedLine = 3

# Show line ticks for selected ions only
show_selected_ions_only = True

# Show line ticks for selected ions only
show_selected_intensities_only = True

# Differentiate lines by: 0 - ion and reference line, 1 - ion and process, 2 - ion, 3 - element 
diff_lines_by = 2

# Color of line ticks 
line_tick_color = 'black'

# Plot lines of selected ions
plot_lines_of_selected_ions = True

# Color of selected line ticks 
color_selected_ions = [ 'DeepSkyBlue', 'magenta', 'BlueViolet', 'LightGreen' ]

# List of selected ions 
selected_ions = []

# Sort selected ions 
selected_ions_sort = False

# Index to cycle the selected ions
index_of_selected_ions = -1

# Scale factor for the plot of the theoretical (individual) spectra. Still in devel
fact_multi_synth =  1.

# norm_intens: value by which are divided the line intensities in
# the printed table. Still in devel
norm_intens = 1.00e4

# instrument profile
instr_prof = {'largeur':0.50,
    'B_1r':0.00,'B_1l':0.00,'decroiss_1':2.70,'alpha_1':0.50,
    'B_2r':0.00,'B_2l':0.00,'decroiss_2':1.50,'alpha_2':0.75,
    'B_3r':0.00,'B_3l':0.00,'decroiss_3':1.15,'alpha_3':1.00,
    'B_4r':0.00,'B_4l':0.00,'decroiss_4':0.45,'alpha_4':1.50,
    'comment':'default profile'}

ghost = {"do_ghost":0, "delta_lambda" : 0. , "intens" : [ 0.00]}

#ghost =  {do_ghost:1, $
#          delta_lambda : [2.74,-2.75,5.2,-5.2] , $
#          largeur : [0.5,0.5,0.7,0.7],$
#          intens : [ 0.00310000, $
#                      0.000880000, $
#                      3.50000e-05, $
#                      8.00000e-06] }

# List of plot legend formats for different line processes = [ [ [list of process codes], format ], ... ] 
"""
code process
1    Rec
2    Dielectronic
3    Collisional 
4    Bowen (set to Fluorescence) 
5    FREE  (set to Rec)
6    FREE  (set to Rec)
7    Fluorescence
8    Charge Exchange
9    Rec. fixed ones.
"""
process_code_format = [ [[0,1,5,6,9], '{} (rec)' ], [[2], '{} (die)'], [[3], '[{}]'], [[7,4], '{} (fl)'], [[8], '{} (chEx)'] ]

#-------------------------------------------------------------------
#   Save parameters to file
#-------------------------------------------------------------------

# Automatically save all synthesis and plot parameters at exit
save_parameters_on_exit = False

#-------------------------------------------------------------------
#   Save lines to file
#-------------------------------------------------------------------

# Fields to be printed 
save_lines_fields = [ 'num', 'id', 'lambda', 'l_shift', 'i_rel', 'i_cor', 'ref', 'profile', 'vitesse', 'comment' ]

# Print header  
save_lines_header = True

# Lines sorted by: 
#   0 - wavelength, 
#   1 - decreasing wavelength 
#   2 - intensity
#   3 - decreasing intensity
#   4 - ion
#   5 - decreasing ion
save_lines_sort = 0

# Save line filename
save_lines_filename = 'lines.dat'

#-------------------------------------------------------------------
#   Save plot to file
#-------------------------------------------------------------------

# Plot filename
plot_filename = 'plot.pdf'

# Order the cosmetic_file by line number and remove duplicate lines.
order_cosmetic_file = False

# Remove from cosmetic file the unchanged lines
clean_cosmetic_file = False
