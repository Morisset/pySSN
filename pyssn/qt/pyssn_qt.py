"""
This is the window manager part of pySSN

pySSN is available under the GNU licence providing you cite the developpers names:

    Ch. Morisset (Instituto de Astronomia, Universidad Nacional Autonoma de Mexico)

    D. Pequignot (Meudon Observatory, France)

Inspired by a demo code by: 
Eli Bendersky (eliben@gmail.com)
"""
import sys, os

from PyQt4 import QtCore, QtGui

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import numpy as np
from pyssn import log_, __version__
from ..core.spectrum import spectrum
from ..utils.misc import get_parser

from collections import OrderedDict
from ..utils.physics import CST

log_.level = 4

#ToDo : 

class NavigationToolbar( NavigationToolbar2QT ):

    curs = QtCore.pyqtSignal(bool)

    def __init__(self, canvas, parent ):
        
        NavigationToolbar2QT.__init__(self,canvas,parent)
        self.clearButtons=[]
        # Search through existing buttons
        # next use for placement of custom button
        next=None
        for c in self.findChildren(QtGui.QToolButton):
            if next is None:
                next=c
            # Don't want to see subplots and customize
            """
            if str(c.text()) in ('Subplots', 'Customize'):
                c.defaultAction().setVisible(False)
                continue
            """
            # Need to keep track of pan and zoom buttons
            # Also grab toggled event to clear checked status of picker button
            if str(c.text()) in ('Pan','Zoom'):
                c.toggled.connect(self.clearCurs)
                self.clearButtons.append(c)
                next=None

        # create custom button
        pm=QtGui.QPixmap(32,32)
        pm.fill(QtGui.QApplication.palette().color(QtGui.QPalette.Normal,QtGui.QPalette.Button))
        painter=QtGui.QPainter(pm)
        painter.fillRect(6,6,20,20,QtCore.Qt.red)
        painter.fillRect(15,3,3,26,QtCore.Qt.blue)
        painter.fillRect(3,15,26,3,QtCore.Qt.blue)
        painter.end()
        icon=QtGui.QIcon(pm)
        
        ac = self.addAction(icon, "Toggle Curs") 
        ac.setCheckable(True)
        
        #Ver como inicializar
        #ac.setChecked(True)
        
        ac.toggled.connect(self.curs_toggle)        
        
        self.ac = ac
        
        #button=QtGui.QToolButton(self)
        #button.setDefaultAction(self.ac)

        # Add it to the toolbar, and connect up event
        #self.insertWidget(next.defaultAction(),button)

        # Grab the picked event from the canvas
        canvas.mpl_connect('pick_event',self.canvasPicked)

    def clearCurs(self, checked):
        if checked:
            self.ac.setChecked(False)

    def curs_toggle(self, checked):
        self.curs.emit(checked)
    
    def canvasPicked(self, event):
        if self.ac.isChecked():
            self.curs.emit(event.ind)


class AppForm(QtGui.QMainWindow):
    
    def __init__(self, parent=None, init_filename=None, post_proc_file=None, use_workspace=False):
        
        self.calling = 'pySSN GUI'
        self.use_workspace = use_workspace
        QtGui.QMainWindow.__init__(self, parent)
        self.setWindowTitle('pySSN')
        self.sp = None
        self.axes = None
        self.axes2 = None
        self.axes3 = None
        self.fig = None
        self.init_file_name = init_filename
        self.init_line_num = None
        self.init_ion = None
        self.init_xmin = None
        self.init_xmax = None
        self.init_y1min = None
        self.init_y1max = None
        self.init_y3min = None
        self.init_y3max = None
        self.init_legend_fontsize = None
        self.init_legend_loc = None

        self.init_nearby_line_num = None
        self.init_nearby_ion = None
        self.init_nearby_xmin = None
        self.init_nearby_xmax = None
        self.init_nearby_y1min = None
        self.init_nearby_y1max = None
        self.init_nearby_y3min = None
        self.init_nearby_y3max = None
        self.init_nearby_legend_fontsize = None
        self.init_nearby_legend_loc = None

        self.call_on_draw = True
        self.cursor_on = False
        self.line_info_ref = 0
        self.x_plot_lims = None
        self.y1_plot_lims = None
        self.y2_plot_lims = None
        self.y3_plot_lims = None
        self.xscale = None
        self.yscale = None        
        self.post_proc_file = post_proc_file
        self.tick_file = None
        self.save_parameters_file = None
        self.do_save = True
        self.cont_par_changed = False
        self.axes_fixed = False
        self.showErrorBox = True
        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        self.exec_init()
        self.cont_pars_dialog = None
        self.cursor_w1 = None
        self.cursor_w2 = None
        self.nearbyLines = None
        self.nearbyLines_sort_by = 'i_tot'
        self.nearbyLines_sort_reverse = True
        self.nearbyLines_dialog = None
        self.nearbyLines_selected_ions = None
        self.line_info_dialog = None
        self.instr_prof_dialog = None
        self.refine_wave_dialog = None
        self.refine_wave_as_table = False
        self.fig_prof = None 
        self.green_tick_shown = False
        self.magenta_tick_shown = False
        self.addGreenTickToLegend = True
        self.show_true_ions = False
        self.nearbyDialogFilterIsActive = False
        
    def closeEvent(self, evnt):
        if self.sp.get_conf('save_parameters_on_exit'):
            self.save_pars_as()
        if self.cont_pars_dialog is not None:
            self.cont_pars_dialog.close()
        if self.nearbyLines_dialog is not None:
            self.nearbyLines_dialog.close()
        if self.line_info_dialog is not None:
            self.line_info_dialog.close()
            self.line_info_table.close()
        if self.instr_prof_dialog is not None:
            self.instr_prof_dialog.close()
        if self.refine_wave_dialog is not None:
            self.refine_wave_dialog.close()

    def image_extension_list(self):
        filetypes = self.canvas.get_supported_filetypes()
        file_extensions = filetypes.keys()
        file_extensions.sort()
        return file_extensions

    def image_filter(self, fileExt=''):
        filetypes = self.canvas.get_supported_filetypes_grouped()
        imagetype_list = filetypes.keys()
        imagetype_list.sort()
        s = ''
        k = 0
        for imagetype in imagetype_list:
            extension_list = filetypes[ imagetype ]
            if fileExt in extension_list:
                k = imagetype_list.index(imagetype)
            s = s + str(imagetype)
            s1 = ' (*.' + str(extension_list[0])
            for extension in extension_list[1:]:
                s1 = s1 + ' *.' + str(extension)
            s1 = s1 + ')'
            s = s + s1 + s1 + ';;'
        filter_str = s[:-2]
        selectedFilter = s.split(';;')[k] 
        return filter_str, selectedFilter

    def save_plot(self):
        path = self.sp.get_conf('plot_filename')
        self.canvas.print_figure(path, dpi=self.dpi)
        self.statusBar().showMessage('Plot saved to file %s' % path, 2000)
      
    def save_plot_as(self):
        path = self.sp.get_conf('plot_filename')
        extension = os.path.splitext(path)[1][1:].lower()
        file_choices, selectedFilter = self.image_filter(extension)
        path = unicode(QtGui.QFileDialog.getSaveFileName(self, 'Save plot to file', path, file_choices, selectedFilter))
        if path:
            extension = os.path.splitext(path)[1][1:].lower()
            if extension in self.image_extension_list():
                self.sp.set_conf('plot_filename', path)
                self.canvas.print_figure(path, dpi=self.dpi)
                self.statusBar().showMessage('Plot saved to file %s' % path, 2000)
            else:
                title = 'Error saving plot'
                msg = 'Format "{0}" not supported.'.format(extension)
                msg = msg + '\nSupported formats: '
                extension_list = self.image_extension_list()
                n = len(extension_list)-1
                s = ''
                for i in range(0,n):
                    s = s + extension_list[i] + ', '
                s = s + extension_list[n] + '.'
                msg = msg + s
                QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
    
    def on_about(self):
        msg = """ pySSN (Spectral Synthesis for Nebulae):        
        """
        QtGui.QMessageBox.about(self, "About the demo", msg.strip())
    
    def set_cursor(self, checked):
        self.cursor_on = checked
        self.sp.firstClick = True
                
    def on_click(self, event):
        if self.cursor_on:
            do_print = not self.sp.get_conf('qt_show_dialogs', True)
            nearbyLines = self.sp.nearby_lines(event, do_print, sort='i_tot', reverse=True)
            if nearbyLines is None:
                return
            self.nearbyLines = nearbyLines
            if not do_print:
                self.show_nearbyLines_dialog()

    def sort_nearbyLines(self, sort, reverse=False):
        if self.nearbyLines is None:
            return
        if sort == 'proc':
            sorts = np.argsort([ self.sp.process[str(line_num)[-9]] for line_num in self.nearbyLines['num'] ])
        else:
            sorts = np.argsort(self.nearbyLines[sort])
        if reverse:
            sorts = sorts[::-1]
        self.nearbyLines = np.array(self.nearbyLines)[sorts]

    def create_main_frame(self):
        
        if self.use_workspace:
            self.main_frame = QtGui.QWorkspace()
        else:
            self.main_frame = QtGui.QWidget()
        # Create the mpl Figure and FigCanvas objects. 
        #
        self.dpi = 100
        #self.fig = plt.figure(figsize=(15,15))
        self.fig = plt.figure(figsize=(15,15))
        # self.fig = plt.figure(figsize=(20.0, 15.0), dpi=self.dpi)
        
        log_.debug('creating figure {}'.format(id(self.fig)), calling=self.calling)
                
        self.canvas = FigureCanvas(self.fig)
        if self.use_workspace:
            self.main_frame.addWindow(self.canvas)
            self.fig2 = Figure((20.0, 15.0), dpi=self.dpi)
            self.canvas2 = FigureCanvas(self.fig2)
            #self.main_frame.addWindow(self.canvas2)
        else:
            self.canvas.setParent(self.main_frame)
                
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('figure_leave_event', self.leave_fig)
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        self.mpl_toolbar.curs.connect(self.set_cursor)   

        # Other GUI controls
        # 
        
        self.fix_axes_cb = QtGui.QCheckBox("fix")
        self.fix_axes_cb.setChecked(False)
        self.connect(self.fix_axes_cb, QtCore.SIGNAL('stateChanged(int)'), self.fix_axes)

        self.xlim_min_box = QtGui.QLineEdit()
        self.xlim_min_box.setMinimumWidth(50)
        #self.connect(self.xlim_min_box, QtCore.SIGNAL('editingFinished()'), self.validate_xlim_min)
        self.connect(self.xlim_min_box, QtCore.SIGNAL('returnPressed()'), self.set_plot_limits_and_draw)

        self.xlim_max_box = QtGui.QLineEdit()
        self.xlim_max_box.setMinimumWidth(50)
        #self.connect(self.xlim_max_box, QtCore.SIGNAL('editingFinished()'), self.validate_xlim_max)
        #self.xlim_max_box.editingFinished.connect(self.validate_xlim_max)
        self.connect(self.xlim_max_box, QtCore.SIGNAL('returnPressed()'), self.set_plot_limits_and_draw)
        
        self.y1lim_min_box = QtGui.QLineEdit()
        self.y1lim_min_box.setMinimumWidth(50)
        #self.connect(self.y1lim_min_box, QtCore.SIGNAL('editingFinished()'), self.validate_y1lim_min)
        self.connect(self.y1lim_min_box, QtCore.SIGNAL('returnPressed()'), self.set_plot_limits_and_draw)
        
        self.y1lim_max_box = QtGui.QLineEdit()
        self.y1lim_max_box.setMinimumWidth(50)
        #self.connect(self.y1lim_max_box, QtCore.SIGNAL('editingFinished()'), self.validate_y1lim_max)
        self.connect(self.y1lim_max_box, QtCore.SIGNAL('returnPressed()'), self.set_plot_limits_and_draw)
        
        self.y3lim_min_box = QtGui.QLineEdit()
        self.y3lim_min_box.setMinimumWidth(50)
        #self.connect(self.y3lim_min_box, QtCore.SIGNAL('editingFinished()'), self.validate_y3lim_min)
        self.connect(self.y3lim_min_box, QtCore.SIGNAL('returnPressed()'), self.set_plot_limits_and_draw)
        
        self.y3lim_max_box = QtGui.QLineEdit()
        self.y3lim_max_box.setMinimumWidth(50)
        #self.connect(self.y3lim_max_box, QtCore.SIGNAL('editingFinished()'), self.validate_y3lim_max)
        self.connect(self.y3lim_max_box, QtCore.SIGNAL('returnPressed()'), self.set_plot_limits_and_draw)
        
        """
        self.select_init_button = QtGui.QPushButton("Init file")
        self.connect(self.select_init_button, QtCore.SIGNAL('clicked()'), self.select_init)
        """
        
        self.run_button = QtGui.QPushButton("Run")
        self.connect(self.run_button, QtCore.SIGNAL('clicked()'), self.rerun)
        
        self.draw_button = QtGui.QPushButton("Draw")
        self.connect(self.draw_button, QtCore.SIGNAL('clicked()'), self.on_draw)
        
        self.Command_GroupBox = QtGui.QGroupBox("Execute")
        self.Command_GroupBox.setCheckable(False)
        
        self.ObsSpec_GroupBox = QtGui.QGroupBox("Parameters of the synthetic spectrum")
        self.ObsSpec_GroupBox.setCheckable(False)

        self.SpecPlot_GroupBox = QtGui.QGroupBox("Plot of spectra")
        self.SpecPlot_GroupBox.setCheckable(False)

        self.lineIDs_GroupBox = QtGui.QGroupBox("Show lines")
        self.lineIDs_GroupBox.setCheckable(True)
        self.lineIDs_GroupBox.setChecked(True)
        
        self.connect(self.lineIDs_GroupBox, QtCore.SIGNAL('clicked()'), self.show_lines_clicked)
        self.lineIDs_GroupBox_ToolTip = 'Check to show ticks at the central positions of the spectral lines and plot the lines of selected ions'

        self.residual_GroupBox = QtGui.QGroupBox("Plot of residuals")
        self.residual_GroupBox.setCheckable(True)
        self.residual_GroupBox.setChecked(True)
        self.connect(self.residual_GroupBox, QtCore.SIGNAL('clicked()'), self.residual_box_clicked)
        self.residual_GroupBox_ToolTip = 'Check to display the residual plot'

        self.adjust_button = QtGui.QPushButton("Update")
        self.adjust_button.setChecked(False)
        self.connect(self.adjust_button, QtCore.SIGNAL('clicked()'), self.adjust)

        self.post_proc_button = QtGui.QPushButton("Post proc")
        self.post_proc_button.setChecked(False)
        self.connect(self.post_proc_button, QtCore.SIGNAL('clicked()'), self.apply_post_proc)

        self.update_profile_button = QtGui.QPushButton("Update profiles")
        self.update_profile_button.setChecked(False)
        self.connect(self.update_profile_button, QtCore.SIGNAL('clicked()'), self.update_profile)

        self.sp_min_box = QtGui.QLineEdit()
        self.sp_min_box.setMinimumWidth(50)
        #self.connect(self.sp_min_box, QtCore.SIGNAL('editingFinished()'), self.set_limit_sp)
        self.connect(self.sp_min_box, QtCore.SIGNAL('returnPressed()'), self.set_limit_sp_and_run)
        
        self.sp_max_box = QtGui.QLineEdit()
        self.sp_max_box.setMinimumWidth(50)
        #self.connect(self.sp_max_box, QtCore.SIGNAL('editingFinished()'), self.set_limit_sp)
        self.connect(self.sp_max_box, QtCore.SIGNAL('returnPressed()'), self.set_limit_sp_and_run)
        
        self.sp_norm_box = QtGui.QLineEdit()
        self.sp_norm_box.setMinimumWidth(50)
        self.connect(self.sp_norm_box, QtCore.SIGNAL('returnPressed()'), self.sp_norm)

        self.obj_velo_box = QtGui.QLineEdit()
        self.obj_velo_box.setMinimumWidth(50)
        self.connect(self.obj_velo_box, QtCore.SIGNAL('returnPressed()'), self.obj_velo)

        self.ebv_box = QtGui.QLineEdit()
        self.ebv_box.setMinimumWidth(50)
        self.connect(self.ebv_box, QtCore.SIGNAL('returnPressed()'), self.ebv)
        
        self.resol_box = QtGui.QLineEdit()
        self.resol_box.setMinimumWidth(50)
        self.connect(self.resol_box, QtCore.SIGNAL('returnPressed()'), self.resol)
        
        self.cut2_box = QtGui.QLineEdit()
        self.cut2_box.setMinimumWidth(50)
        self.connect(self.cut2_box, QtCore.SIGNAL('returnPressed()'), self.cut2)
        
        self.cut_cb = QtGui.QCheckBox('')
        self.cut_cb.setChecked(False)
        self.connect(self.cut_cb, QtCore.SIGNAL('clicked()'), self.cut_cb_changed)
        
        self.ion_box = QtGui.QLineEdit()
        self.ion_box.setMinimumWidth(70)
        self.connect(self.ion_box, QtCore.SIGNAL('returnPressed()'), self.draw_ion)
        
        self.ion_cb = QtGui.QCheckBox('')
        self.ion_cb.setChecked(False)
        self.connect(self.ion_cb, QtCore.SIGNAL('clicked()'), self.ion_cb_changed)
        
        self.line_info_box = QtGui.QLineEdit()
        self.line_info_box.setFixedWidth(130)
        self.connect(self.line_info_box, QtCore.SIGNAL('returnPressed()'), self.line_info)

        self.mpl_toolbar.addSeparator()
        self.mpl_toolbar.addWidget(QtGui.QLabel('   line number '))
        self.mpl_toolbar.addWidget(self.line_info_box)

        self.magenta_box = QtGui.QLineEdit()
        self.magenta_box.setMinimumWidth(50)
        self.connect(self.magenta_box, QtCore.SIGNAL('returnPressed()'), self.magenta_line)

        self.magenta_label_box = QtGui.QLineEdit()
        self.magenta_label_box.setMinimumWidth(50)
        self.connect(self.magenta_label_box, QtCore.SIGNAL('returnPressed()'), self.magenta_line)

        self.cyan_box = QtGui.QLineEdit()
        self.cyan_box.setMinimumWidth(50)
        self.connect(self.cyan_box, QtCore.SIGNAL('returnPressed()'), self.cyan_line)
        
        self.cyan_label_box = QtGui.QLineEdit()
        self.cyan_label_box.setMinimumWidth(50)
        self.connect(self.cyan_label_box, QtCore.SIGNAL('returnPressed()'), self.cyan_line)


        self.setStyleSheet("""QToolTip { 
                           background-color: black; 
                           color: lightgray; 
                           min-width: 20em;
                           font-size: 14px;
                           font-family: "sans-serif";
                           border: black solid 10px
                           }""")
        
        s = 'Click to execute the synthesis from the beginning.'
        self.run_button_ToolTip = s        
        
        s = 'Click to update synthesis with changes in line intensities, profiles, and continuum parameters.'
        self.adjust_button_ToolTip = s
        
        s = 'Enter line number to get information on\n' \
            'the reference line and on its satellites.'
        self.line_info_box_ToolTip = s        

        s =  'Color excess E(B-V)\n\n' \
             'Set with: \n' \
             '    e_bv = <float>\n\n' \
             'Comment: \n' \
            u'    E(B-V) \u2248 C(H\u03B2) / 1.5'
        self.ebv_box_ToolTip = s  
        
        s = 'Radial velocity in km/s\n\n' \
            'Set with: \n' \
            '    obj_velo = <float>'
        self.obj_velo_box_ToolTip = s 

        s = 'Minimum wavelength of the synthetic spectrum (in angstroms)\n\n' \
            'Set with:  \n' \
            '    limit_sp = (<xmin>, <xmax>)'
        self.sp_min_box_ToolTip = s        

        s = 'Maximum wavelength of the synthetic spectrum  (in angstroms)\n\n' \
            'Set with:  \n' \
            '    limit_sp = (<xmin>, <xmax>)'
        self.sp_max_box_ToolTip = s        

        s = 'Minimum wavelength in the plots of spectra and residuals  (in angstroms)\n\n' \
            'Set with:  \n' \
            '    x_plot_lims = (<xmin>, <xmax>)'
        self.xlim_min_box_ToolTip = s        

        s = 'Maximum wavelength in the plots of spectra and residuals  (in angstroms)\n\n' \
            'Set with:  \n' \
            '    x_plot_lims = (<xmin>, <xmax>)'
        self.xlim_max_box_ToolTip = s        

        s = 'Minimum ordinate in the plot of spectra, in units of relative intensity \n\n' \
            'Set with:  \n' \
            '    y1_plot_lims = (<ymin>, <ymax>)'
        self.y1lim_min_box_ToolTip = s        

        s = 'Maximum ordinate in the plot of spectra, in units of relative intensity\n\n' \
            'Set with:  \n' \
            '    y1_plot_lims = (<ymin>, <ymax>)'
        self.y1lim_max_box_ToolTip = s        

        s = 'Minimum ordinate in the plot of residuals, in units of relative intensity\n\n' \
            'Set with:  \n' \
            '    y3_plot_lims = (<ymin>, <ymax>)'
        self.y3lim_min_box_ToolTip = s        

        s = 'Maximum ordinate in the plot of residuals, in units of relative intensity\n\n' \
            'Set with:  \n' \
            '    y3_plot_lims = (<ymin>, <ymax>)'
        self.y3lim_max_box_ToolTip = s        
        
        s = 'Check to retain the current limits of the plots while zooming and panning.'
        self.fix_axes_cb_ToolTip = s        
        
        s = 'Check to show only lines with intensities above cut. \n\n' \
             'Set with: \n' \
             '    show_selected_intensities_only = <boolean>'
        self.cut_cb_ToolTip = s        

        s = 'Check to show only lines of selected ions. \n\n' \
             'Set with: \n' \
             '    show_selected_ions_only = <boolean>'
        self.ion_cb_ToolTip = s        
        
        s = 'Normalization factor, ratio between the intensity and the \n' \
            u'observed flux of the reference line, usually 10\u2074/F(H\u03B2)\n\n' \
             'Set with: \n' \
             '    sp_norm = <float>'
        self.sp_norm_box_ToolTip = s        

        s = 'Rebinning factor, the odd integer factor by which the number of points \n' \
            'of the original spectrum is multiplied in the rebinning process\n\n' \
            'Set with: \n' \
            '    resol = <integer>\n\n' \
            'Usage: \n' \
            '    Set to \'1\' if the resolution of the observed spectrum is large enough' 
        self.resol_box_ToolTip = s 
        
        s = 'Minimum relative intensity of lines to be shown. \n\n' \
             'Set with: \n' \
             '    cut_plot2 = <float>'
        self.cut2_box_ToolTip = s        
        
        s = 'Comma-separated list of selected ions, elements, or line numbers to be shown. \n\n' \
             'Set with: \n' \
             '    selected_ions = [<ion1>,<ion2>,...]\n\n' \
             'Examples: \n' \
             '    \'O III\' (or \'O_III\') to show the lines of O III\n' \
             '    \'O III*\' (or \'O_III*\') to show the lines of O III, O IIIfl, O III5g, etc\n' \
             '    \'O III, O IV\' to show the lines of O III and O IV\n' \
             '    \'O\' to show the lines of all O ions\n' \
             '    \'Fe, N\' to show the lines of all Fe and N ions\n' \
             '    <line number> to show the lines of that same ion'                
        self.ion_box_ToolTip = s        

        #
        # Layout with box sizers
        # 

        CommandLayout = QtGui.QGridLayout()

        wList = [self.run_button,self.adjust_button]
        Nrow = 2

        for w in wList:
            k = wList.index( w )
            i = k%Nrow
            j = 1+2*(k/Nrow)
            CommandLayout.addWidget(w,i,j)
            CommandLayout.setAlignment(w,QtCore.Qt.AlignCenter)

        self.Command_GroupBox.setLayout(CommandLayout)

        ObsSpecLayout = QtGui.QGridLayout()

        lList = ['xmin', 'xmax', u'10\u2074/F(H\u03B2)', 'radial vel.', 'E(B-V)', 'N']
        wList = [self.sp_min_box, self.sp_max_box, self.sp_norm_box, self.obj_velo_box, self.ebv_box, self.resol_box ]
        Nrow = 2

        for l in lList:
            w = QtGui.QLabel(l)
            k = lList.index( l )
            i = k%Nrow
            j = 2*(k/Nrow)
            ObsSpecLayout.addWidget(w,i,j)
            ObsSpecLayout.setAlignment(w,QtCore.Qt.AlignRight)

        for w in wList:
            k = wList.index( w )
            i = k%Nrow
            j = 1+2*(k/Nrow)
            ObsSpecLayout.addWidget(w,i,j)
            ObsSpecLayout.setAlignment(w,QtCore.Qt.AlignRight)

        self.ObsSpec_GroupBox.setLayout(ObsSpecLayout)

        SpecPlotLayout = QtGui.QGridLayout()
        SpecPlotLayout.addWidget(QtGui.QLabel('xmin'),0,0)
        SpecPlotLayout.addWidget(QtGui.QLabel('xmax'),1,0)
        SpecPlotLayout.addWidget(QtGui.QLabel('ymin'),0,2)
        SpecPlotLayout.addWidget(QtGui.QLabel('ymax'),1,2)
        SpecPlotLayout.addWidget(self.xlim_min_box,0,1)
        SpecPlotLayout.addWidget(self.xlim_max_box,1,1)
        SpecPlotLayout.addWidget(self.y1lim_min_box,0,3)
        SpecPlotLayout.addWidget(self.y1lim_max_box,1,3)
        SpecPlotLayout.addWidget(self.fix_axes_cb,0,4)

        self.SpecPlot_GroupBox.setLayout(SpecPlotLayout)

        LineIDLayout = QtGui.QGridLayout()
        LineIDLayout.addWidget(QtGui.QLabel('cut'),0,0)
        LineIDLayout.addWidget(self.cut2_box,0,1)
        LineIDLayout.addWidget(self.cut_cb,0,2)
        LineIDLayout.addWidget(QtGui.QLabel('ion'),1,0)
        LineIDLayout.addWidget(self.ion_box,1,1)
        LineIDLayout.addWidget(self.ion_cb,1,2)

        self.lineIDs_GroupBox.setLayout(LineIDLayout)

        ResidualLayout = QtGui.QGridLayout()            
        ResidualLayout.addWidget(QtGui.QLabel('ymin'),0,0)
        ResidualLayout.addWidget(QtGui.QLabel('ymax'),1,0)
        ResidualLayout.addWidget(self.y3lim_min_box,0,1)
        ResidualLayout.addWidget(self.y3lim_max_box,1,1)
 
        self.residual_GroupBox.setLayout(ResidualLayout)
        
        grid = QtGui.QGridLayout()

        grid.addWidget(self.Command_GroupBox, 0, 1 )
        grid.addWidget(self.ObsSpec_GroupBox, 0, 2 )
        grid.addWidget(self.SpecPlot_GroupBox, 0, 3 )
        grid.addWidget(self.residual_GroupBox, 0, 4 )
        grid.addWidget(self.lineIDs_GroupBox, 0, 5 )

        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(grid)        
        #vbox.setAlignment(QtCore.Qt.AlignBottom)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
        
    def create_status_bar(self):
        self.status_text = QtGui.QLabel("pySSN, v{}".format(__version__))
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("File")
        
        open_init_action = self.create_action("Open init file",
                                              shortcut="", 
                                              slot=self.select_init, 
                                              tip="Open the initialization file and run the synthesis")
        
        save_pars_action = self.create_action("Save parameters",
                                              shortcut="Ctrl+S", 
                                              slot=self.save_pars_as, 
                                              tip="Save synthesis and plot parameters to file")
        
        save_pars_as_action = self.create_action("Save parameters as",
                                              shortcut="Ctrl+Shift+S", 
                                              slot=self.save_pars_as, 
                                              tip="Select file name and save parameters of the synthesis")
        
        self.save_plot_action = self.create_action("Save plot",
                                              shortcut="Ctrl+P", 
                                              slot=self.save_plot_as,
                                              tip="Save plot to file")
        
        save_plot_as_action = self.create_action("Save plot as",
                                              shortcut="Ctrl+Shift+P", 
                                              slot=self.save_plot_as, 
                                              tip="Select file name and save plot")

        save_lines_action = self.create_action("Save lines",
                                              shortcut="Ctrl+L", 
                                              slot=self.save_lines_as, 
                                              tip="Save list of lines to file")

        save_lines_as_action = self.create_action("Save lines as",
                                              shortcut="Ctrl+Shift+L", 
                                              slot=self.save_lines_as, 
                                              tip="Select file name and save list of lines")

        self.add_actions(self.file_menu, 
            (open_init_action, save_pars_action, None, self.save_plot_action, None, save_lines_action))
       
        #(open_init_action, save_pars_action, save_pars_as_action, None, self.save_plot_action, save_plot_as_action, None, save_lines_action, save_lines_as_action))

        self.line_sort_list = ['wavelength', 'decreasing wavelength', 'intensity', 'decreasing intensity', 'ion' , 'decreasing ion' ]
        s = 'Sort lines by:\n'
        for i in range(len(self.line_sort_list)):
            s = s + '    ' + str(i) + ' - ' + self.line_sort_list[i] + '\n'
        s = s + '\nSet with:\n' + '    save_lines_sort = <integer>'
        self.line_sort_ag = QtGui.QActionGroup(self, exclusive=True)

        self.line_sort_menu = self.file_menu.addMenu("Sort lines by")
        self.line_sort_menu_ToolTip = ''       
        
        for i in range(len(self.line_sort_list)):
            a = self.line_sort_ag.addAction(QtGui.QAction(self.line_sort_list[i], self, checkable=True))
            self.line_sort_menu.addAction(a)

        self.line_sort_ag.triggered.connect(self.line_sort)
        
        self.line_print_dic = OrderedDict( [ 
                                ( 'num'     , 'line number' ),
                                ( 'id'      , 'ion' ), 
                                ( 'lambda'  , 'wavelength' ), 
                                ( 'l_shift' , 'wavelength shift' ), 
                                ( 'l_tot'   , 'corrected wavelength' ), 
                                ( 'i_rel'   , 'intensity' ), 
                                ( 'i_cor'   , 'intensity correction factor' ), 
                                ( 'i_tot'   , 'corrected intensity' ),
                                ( 'ref'     , 'reference line number' ),
                                ( 'profile' , 'line profile code number' ),
                                ( 'vitesse' , 'natural line width' ), 
                                ( 'comment' , 'comment' ) ])

        items = list(self.line_print_dic.values())
        
        s = 'Fields to be printed:\n'
        for i in range(len(items)):
            s = s + '    ' + str(i) + ' - ' + items[i] + '\n'
        s = s + '\nSet with:\n' + '    save_lines_fields = <list>'

        self.line_field_menu = self.file_menu.addMenu("Show fields")
        self.line_field_menu_ToolTip = ''      

        for i in range(len(items)):
            a = self.create_action(items[i], 
              shortcut='', slot=self.set_line_fields_to_print, checkable=True, 
              tip=None)
            self.line_field_menu.addAction(a)
        self.file_menu.addMenu(self.line_field_menu)

        self.show_header_action = self.create_action("Show header",
                                              slot=self.set_show_header, 
                                              shortcut="", 
                                              checkable=True, 
                                              tip="Show header in list of lines")

        self.file_menu.addAction(self.show_header_action)

        self.open_cosmetic_file_action = self.create_action("Open cosmetic file",
                                              slot=self.set_cosmetic_file, 
                                              shortcut="", 
                                              tip="Open the cosmetic file")

        self.clean_cosmetic_file_action = self.create_action("Clean cosmetic file",
                                              slot=self.clean_cosmetic_file, 
                                              shortcut="", 
                                              tip="Remove the unchanged lines from the cosmetic file")

        self.empty_cosmetic_file_action = self.create_action("Empty cosmetic file",
                                              slot=self.empty_cosmetic_file, 
                                              shortcut="", 
                                              tip="Remove all lines from the cosmetic file")

        self.order_cosmetic_file_action = self.create_action("Order cosmetic file",
                                              slot=self.order_cosmetic_file, 
                                              shortcut="", 
                                              tip="Order the cosmetic file by line number and remove duplicate lines")

        quit_action = self.create_action("&Quit", 
                                         slot=self.fileQuit, 
                                         shortcut="Ctrl+Q", 
                                         tip="Close the application")

        self.add_actions(self.file_menu, (None, self.open_cosmetic_file_action, self.clean_cosmetic_file_action, 
                                          self.order_cosmetic_file_action, self.empty_cosmetic_file_action, None, quit_action))

        self.run_menu = self.menuBar().addMenu("Execute")
        
        run_action = self.create_action("Run",
                                         shortcut="Ctrl+F9", 
                                         slot=self.rerun, 
                                         tip="Execute synthesis from the beginning")
        
        update_action = self.create_action("Update",
                                            shortcut="F9", 
                                            slot=self.adjust, 
                                            tip="Update synthesis with changes in line intensities, profiles, and continuum parameters")

        draw_action = self.create_action("Draw",
                                         shortcut="F8", 
                                         slot=self.set_plot_limits_and_draw, 
                                         tip="Redraw plots")
        
        post_proc_action = self.create_action("Post-process",
                                               shortcut="Ctrl+F8", 
                                               slot=self.apply_post_proc, 
                                               tip="Edit the plots with python commands defined in an external file")
        
        open_profile_action = self.create_action("Instrumental profile",
                                              shortcut="F7", 
                                              slot=self.apply_instr_prof, 
                                              tip="Open the instrumental profile file and run the synthesis")
                                              
        refine_wavelengths_action = self.create_action("Wavelength-refining", 
                                         slot=self.refine_wavelengths, 
                                         shortcut="F6", 
                                         tip="Refine the wavelength calibration")
                                         
        self.add_actions(self.run_menu, (update_action, run_action, draw_action, None, 
                                         post_proc_action, open_profile_action, refine_wavelengths_action))

        self.line_menu = self.menuBar().addMenu('Lines')
        
        self.show_line_ticks_action = self.create_action('Plot line ticks', 
            shortcut='Alt+L', slot=self.show_line_ticks_action_clicked, checkable=True, 
            tip='Check to show line ticks')
        
        self.plot_lines_action = self.create_action('Plot spectra of selected ions', 
            shortcut='Alt+P', slot=self.show_line_ticks_action_clicked, checkable=True, 
            tip='Check to plot spectra of selected ions')
        
        self.selected_intensities_action = self.create_action('Only above the cut', 
            shortcut='Alt+K', slot=self.selected_lines_clicked, checkable=True, 
            tip='Check to show the ticks for lines with intensities above cut only')
        
        self.selected_ions_action = self.create_action('Only for selected ions', 
            shortcut='Alt+I', slot=self.selected_lines_clicked, checkable=True, 
            tip='Check to show the line ticks for selected ions only')

        self.add_actions(self.line_menu, 
            (self.plot_lines_action, None, self.show_line_ticks_action, self.selected_intensities_action, self.selected_ions_action))

        self.diff_lines_list = ['ion and reference line', 'ion and process', 'ion', 'element' ]
        s = 'Differentiate lines by:\n'
        for i in range(len(self.diff_lines_list)):
            s = s + '    ' + str(i) + ' - ' + self.diff_lines_list[i] + '\n'
        s = s + '\nSet with:\n' + '    diff_lines_by = <integer>'
        self.diff_lines_ag = QtGui.QActionGroup(self, exclusive=True)

        self.diff_lines_menu = self.line_menu.addMenu("Differentiate lines by")
        self.diff_lines_menu_ToolTip = ''        

        for i in range(len(self.diff_lines_list)):
            a = self.diff_lines_ag.addAction(QtGui.QAction(self.diff_lines_list[i], self, checkable=True))
            a.setShortcut('Alt+' + str(i+1))
            self.diff_lines_menu.addAction(a)
        self.diff_lines_ag.triggered.connect(self.diff_lines)
        
        self.cycle_forwards_ions_action = self.create_action('Cycle forwards selected ions', 
            shortcut='Alt+0', slot=self.cycle_forwards_ions, checkable=False, 
            tip='Click to cycle forwards the selected ions')
        
        self.cycle_backwards_ions = self.create_action('Cycle backwards selected ions', 
            shortcut='Alt+9', slot=self.cycle_backwards_ions, checkable=False, 
            tip='Click to cycle backwards the selected ions')

        self.add_actions(self.line_menu, 
            (None, self.cycle_forwards_ions_action, self.cycle_backwards_ions, None))

        self.line_tick_ax_menu = self.line_menu.addMenu('Window of line ticks')
       
        self.line_tick_ax_list = ['Plot of spectra', 'Plot of residuals', 'Separate plot' ]
        s = 'Show line ticks on:\n'
        for i in range(len(self.line_tick_ax_list)):
            s = s + '    ' + str(i) + ' - ' + self.line_tick_ax_list[i] + '\n'
        s = s + '\nSet with:\n' + '    line_tick_ax = <integer>'
        self.line_tick_ax_ag = QtGui.QActionGroup(self, exclusive=True)
        self.line_tick_ax_menu_ToolTip = ''        

        for i in range(len(self.line_tick_ax_list)):
            a = self.line_tick_ax_ag.addAction(QtGui.QAction(self.line_tick_ax_list[i], self, checkable=True))
            self.line_tick_ax_menu.addAction(a)
        self.line_tick_ax_ag.triggered.connect(self.set_plot_ax2)

        self.line_tick_pos_menu = self.line_menu.addMenu('Position of line ticks')

        self.line_tick_pos_list = ['Top', 'Middle', 'Bottom' ]
        s = 'Position line ticks:\n'
        for i in range(len(self.line_tick_pos_list)):
            s = s + '    ' + str(i) + ' - ' + self.line_tick_pos_list[i] + '\n'
        s = s + '\nSet with:\n' + '    line_tick_pos = <integer>'
        self.line_tick_pos_ag = QtGui.QActionGroup(self, exclusive=True)
        self.line_tick_pos_menu_ToolTip = ''        

        for i in range(len(self.line_tick_pos_list)):
            a = self.line_tick_pos_ag.addAction(QtGui.QAction(self.line_tick_pos_list[i], self, checkable=True))
            self.line_tick_pos_menu.addAction(a)
        self.line_tick_pos_ag.triggered.connect(self.set_plot_ax2)

        self.line_tick_color_action = self.create_action('Color of line ticks', 
            shortcut=None, slot=self.line_tick_color_clicked, checkable=False, 
            tip='Set color of line ticks')

        self.toggle_legend_action = self.create_action('Toggle legend position and zoom', 
            shortcut='Alt+Shift+L', slot=self.toggle_legend_clicked, checkable=False, 
            tip='Toggle the legend position and zoom')
        self.line_menu.addAction(self.toggle_legend_action)
        
        self.editing_lines_action = self.create_action('Allow editing line parameters', 
            slot=self.editing_lines_clicked, checkable=True, 
            tip='Check to allow editing line parameters in line info dialog')
        
        self.update_lines_action = self.create_action('Update after editing line parameters', 
            shortcut='Alt+U', slot=self.update_lines_clicked, checkable=True, 
            tip='Check to update synthesis after editing line parameters in line info dialog')
       
        self.show_line_ticks_from_file_action = self.create_action('Plot line ticks from file', 
            shortcut='F4', slot=self.show_line_ticks_from_file,  
            tip='Check to show line ticks defined in an external file')
        
        self.ask_tickfile_action = self.create_action("Ask for file name",
            checkable=True, tip="Check to be always asked for the text file containing a list of wavelengths to be ticked")
             
        self.add_actions(self.line_menu, (None, self.show_line_ticks_from_file_action))

        self.cont_menu = self.menuBar().addMenu('Continuum')

        self.plot_cont_action = self.create_action('Plot continuum',
                                          shortcut="Alt+C", 
                                          slot=self.plot_cont_action_clicked, 
                                          checkable=True, 
                                          tip='Check to plot the different components of the continuum spectrum')
        
        self.cont_action = self.create_action('Parameters',
                                          shortcut="Shift+Alt+C", 
                                          slot=self.cont_dialog, 
                                          tip='Paremeters of the continuum spectrum')
        
        self.add_actions(self.cont_menu, 
             (self.plot_cont_action, self.cont_action,))


        self.settings_menu = self.menuBar().addMenu('Settings')

        self.verbosity_list = ['None', 'Errors', 'Errors and warnings', 'Errors, warnings, and comments', 'Debug messages' ]
        s = 'Verbosity level:\n'
        for i in range(len(self.verbosity_list)):
            s = s + '    ' + str(i) + ' - ' + self.verbosity_list[i] + '\n'
        s = s + '\nSet with:\n' + '    log_level = <integer>'
        self.verbosity_ag = QtGui.QActionGroup(self, exclusive=True)
        
        #self.verbosity_menu = self.menuBar().addMenu("Verbosity")
        self.verbosity_menu = self.settings_menu.addMenu("Verbosity")
        self.verbosity_menu_ToolTip = ''        

        for i in range(len(self.verbosity_list)):
            a = self.verbosity_ag.addAction(QtGui.QAction(self.verbosity_list[i], self, checkable=True))
            self.verbosity_menu.addAction(a)
        self.verbosity_ag.triggered.connect(self.verbosity)

        self.style_list = list(QtGui.QStyleFactory.keys())
        s = 'Widget styles:\n'
        for i in range(len(self.style_list)):
            s = s + '    ' + str(i) + ' - ' + self.style_list[i] + '\n'
        s = s + '\nSet with:\n' + '    qt_style = <integer>'
        self.style_ag = QtGui.QActionGroup(self, exclusive=True)

        self.style_menu = self.settings_menu.addMenu('Widget style')

        self.style_menu_ToolTip = ''        

        for i in range(len(self.style_list)):
            a = self.style_ag.addAction(QtGui.QAction(self.style_list[i], self, checkable=True))
            self.style_menu.addAction(a)
        self.style_ag.triggered.connect(self.style)
                
        self.enable_tooltips_action = self.create_action('Enable tooltips', 
            slot=self.enable_tooltips_action_clicked, checkable=True, 
            tip='Check to enable tooltips')
                
        self.adjust_fig_action = self.create_action('Adjust figure', 
            slot=self.adjust_fig_action_clicked, checkable=True, 
            tip='Automatically adjust figure to avoid overlaps and to minimize the empty borders.')
                
        self.show_uncor_obs_action = self.create_action('Show uncorrected spectrum', 
            slot=self.show_uncor_obs_action_clicked, checkable=True, 
            tip='Show observational spectrum without the wavelength refining.')
            
        self.add_actions(self.settings_menu, 
            (None, self.enable_tooltips_action, self.adjust_fig_action, None, self.editing_lines_action, self.update_lines_action, self.show_uncor_obs_action))         
        
        self.help_menu = self.menuBar().addMenu("&Help")
        
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About the demo')
       
        self.add_actions(self.help_menu, (about_action,))

    def fileQuit(self):
        self.close()
    
    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QtGui.QAction(text, self)
        if icon is not None:
            action.setIcon(QtGui.QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, QtCore.SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action          
    
    def isInteger(self, str_):
        try: 
            int(str_)
            return True
        except ValueError:
            return False
    
    def isPositiveInteger(self, str_):
        if self.isInteger(str_):
            n = int(str_)
            if n > 0:
                return True
            else:
                return False
        else:        
            return False
    
    def isPositiveOdd(self, str_):
        if self.isInteger(str_):
            n = int(str_)
            if n%2 == 1 and n > 0:
                return True
            else:
                return False
        else:        
            return False

    def isFloat(self, str_):
        try:
            np.float(str_)
            return True
        except ValueError:
            return False

    def floatFixFormat(self, r, fix_fmt, align='>'):
        """
        floatFixFormat(1.23456789, '{:7.3f}')     = '  1.234'
        floatFixFormat(-1.23456789, '{:7.3f}')    = ' -1.234'
        floatFixFormat(123.456789, '{:7.3f}')     = ' 1.23e2'
        floatFixFormat(-123.456789, '{:7.3f}')    = '-1.23e2'
        floatFixFormat(1.23456789e+04, '{:7.3f}') = ' 1.23e4'
        floatFixFormat(1.23456789e-04, '{:7.3f}') = ' 1.2e-4'
        floatFixFormat(1.23456789e+34, '{:7.3f}') = ' 1.2e34'
        floatFixFormat(99.999999, '{:7.3f}')      = ' 1.2e34'
        """
        if not ( 'f' in fix_fmt and self.isFloat(r) ):
            return None
        s = fix_fmt.strip('{')    
        s = s.strip('}')    
        s = s.strip(':')    
        s = s.strip('f')
        k = s.index('.')    
        w = int(s[:k])
        p = int(s[k+1:])
        s0 = '{:{align}{w}.{p}f}'.format(float(abs(r)), w=w-1, p=p, align=align)
        s = '{:0.{w}e}'.format(float(abs(r)), w=w)
        if r < 0:
            sgn = '-'
        else:
            sgn = ''
        k = s.index('e')
        mantissa = s[:k]
        mantissa = mantissa[:p+2]
        e = int(s[k+1:])
        if p+e+2>w-3-len(str(e)) and len(s0) < w:
            s = s0.strip()
        else:
            s = '{:0.{p}e}'.format(float(abs(r)), p=min(p,w-4-len(str(e))))
            k = s.index('e')
            mantissa = s[:k]
            exponent = str(int(s[k+1:]))
            s = mantissa + 'e' + exponent
        s = '{:{align}{w}}'.format(sgn+s, w=w, align=align)
        return s

    def rightFormat(self, s, field):
        if field == 'comment':
            output = s.strip()
            return output
        try:
            if field == 'profile':
                r = int(s)
            else:
                r = np.float(s)
            fmt = self.sp.field_format[field]
            if 'f' in fmt:
                s = self.floatFixFormat(r, fmt)
            else:
                s = fmt.format(r)
            if len(s) == self.sp.field_width[field] and not np.isinf(r):
                if field == 'vitesse' and (r < 0 or s.strip() == '0.00'):
                    output = None
                else:    
                    output = s
            else:
                output = None
        except:
            output = None
        return output
            
    def ConvStrToValidTypes(self, str_):
        str_ = str_.strip(' ')
        str_ = str_.replace('Error in ','')
        if str_ == '':
            result = None
        elif str_.isdigit():
            result = int(str_)
        elif self.isFloat(str_):
            result = np.float(str_)
        elif str_.capitalize() == 'True':
            result = True
        elif str_.capitalize() == 'False':
            result = False
        elif str_.find(',') >= 0:
            try:
                str_ = str_.replace('[','')
                str_ = str_.replace(']','')
                result = [float(i) for i in str_.split(',')]
            except:
                result = None
        else:      
            result = None
        return result

    def save_par_in_file(self, field, value, path, help_=None):
        if self.isValidFilename(path):
            if os.path.isfile(path):
                f = open(path, 'r')
                lines = f.readlines()[::-1]
                f.close()
            else:
                lines = []
            j = 0
            found = False
            while ( j < len(lines) ) and ( not found ):
                line = str(lines[j])
                if line.find(field) == 0:
                    if type(value) is str:
                        s0 = ' = \''
                        s1 = '\'\n'
                    else:
                        s0 = ' = '
                        s1 = '\n'
                    line = '# ' + line + field + s0 + value  + s1
                    lines[j] = line
                    found = True
                    break
                j += 1
            if not found:
                if help_ is not None:
                    lines.insert(0, '\n# ' + help_ + '\n')
                lines.insert(0, field + ' = ' + value + '\n')
            lines = lines[::-1]
            f = open(path, 'w')
            f.writelines(lines)
            f.close()

    def save_cont_pars(self):
        file_choices = "Python files (*.py) (*.py);;Text files (*.txt *.dat) (*.txt *.dat);;All Files (*) (*)"
        filename = self.sp.config_file.split('/')[-1]
        path = unicode(QtGui.QFileDialog.getSaveFileName(self, 'Save to file', filename, file_choices))
        if path:
            if os.path.isfile(path):
                f = open(path, 'r')
                lines = f.readlines()[::-1]
                f.close()
            else:
                lines = []
            for i in range(0, self.table.rowCount()):
                field = str(self.table.item(i,0).text())
                value = str(self.table.item(i,1).text())
                help_ = str(self.table.item(i,2).text().toUtf8())
                help_ = help_.replace('\xce\xb2', 'beta')
                help_ = help_.replace('\xe2\x81\xbb\xc2\xb3', '-3')
                help_ = help_.replace('\xce\xb1', 'alpha')
                help_ = help_.replace('\xce\xbb/5000 \xe2\x84\xab', 'lambda/5000 A')
                j = 0
                found = False
                while ( j < len(lines) ) and ( not found ):
                    line = str(lines[j])
                    if line.find(field) == 0:
                        k = line.find('#')
                        if k > 0:
                            comment = ' ' + line[k:]
                        else:
                            comment = '\n'
                        line = field + ' = ' + value + comment
                        lines[j] = line
                        found = True
                        break
                    j += 1
                if not found:
                    lines.insert(0, '\n# ' + help_ + '\n')
                    lines.insert(0, field + ' = ' + value + '\n')
                
            lines = lines[::-1]
            f = open(path, 'w')
            f.writelines(lines)
            f.close()
        
    def get_shifts_from_profile(self, profile_key):
        if profile_key not in self.sp.emis_profiles:
            profile_key = '1'
        vel = self.sp.emis_profiles[profile_key]['vel']
        par_list = self.sp.emis_profiles[profile_key]['params']
        shift_list = []
        for item in par_list:
            shift = np.float(item[2])
            intensity = np.float(item[1])
            if item[0]=='G' and ( intensity > 0.2 ):
                shift_list.append(shift)
        shift_list.sort()
        return shift_list, vel
    
    def plot_tick_at(self, wavelength, ion, line_num):
        if self.green_tick_shown:
            self.on_draw()
        color = 'green'
        ion = ion.replace('_',' ').strip()
        to_select = (self.sp.liste_raies['num'] == np.int(line_num))
        vitesse = self.sp.liste_raies[to_select]['vitesse']
        profile_key = str(self.sp.liste_raies[to_select]['profile'][0])
        shift_list, vel = self.get_shifts_from_profile(profile_key)
        line_num = line_num.strip().strip('0')
        # label = ion + ' (' + line_num.strip() + ')'
        label = ion + ' {:.2f}'.format(wavelength)
        posTick = self.getTickPosOfSelectedLine()
        y1, y2 = self.get_line_tick_lim(posTick)
        k = self.sp.get_conf('line_tick_ax')    
        if not (k == 1 and self.residual_GroupBox.isChecked()):
            k = 0
        if len(shift_list) > 0:
            if posTick == 0:
                ys1 = 2*y1-y2
                ys2 = y1
                ym = y1
            else:
                ys1 = y2
                ys2 = 2*y2-y1
                ym = y2
            if k == 0:
                yy1 = self.y1_plot_lims[0] + ym*(self.y1_plot_lims[1] - self.y1_plot_lims[0])
            else:
                yy1 = self.y3_plot_lims[0] + ym*(self.y3_plot_lims[1] - self.y3_plot_lims[0])
        current_legend_loc = self.sp.legend_loc
        f = 0.15
        r =  (self.x_plot_lims[1] - self.x_plot_lims[0])/2
        if wavelength - self.x_plot_lims[0] < 2*r*f:
            current_legend_loc = 1
        if self.x_plot_lims[1] - wavelength < 2*r*f:
            current_legend_loc = 2
        self.fig.axes[k].axvline( wavelength, y1, y2, color = color, linestyle = 'solid', linewidth = 2.5 ) 
        wave_shifts = -vitesse*wavelength*shift_list / CST.CLIGHT * 1e5 + wavelength*vel / CST.CLIGHT * 1e5
        if len(wave_shifts) > 0:
            max_wave_shift = max(abs(wave_shifts))
        else:
            max_wave_shift = 0

        # Ticks for the profiles components are not shown if they are within 1000*f percent of the x-axis width. 
        f = 0.001 
        if max_wave_shift > f*(self.x_plot_lims[1] - self.x_plot_lims[0]):
            x1 = (wavelength - self.x_plot_lims[0])/(self.x_plot_lims[1] - self.x_plot_lims[0])
            for shift in wave_shifts:            
                self.fig.axes[k].axvline( wavelength+shift, ys1, ys2, color = color, linestyle = '--', linewidth = 2.5 )                 
                x2 = (wavelength + shift - self.x_plot_lims[0])/(self.x_plot_lims[1] - self.x_plot_lims[0])
                self.fig.axes[k].axhline( yy1, x1, x2, color = color, linestyle = '-', linewidth = 1.0 ) 
            
        if self.addGreenTickToLegend:
            self.fig.axes[k].step( [0,0], [0,100], color = color, linestyle = 'solid', label = label, linewidth = 2.5 )
        self.fig.axes[k].legend(loc=current_legend_loc, fontsize=self.sp.legend_fontsize)
        self.fig.canvas.draw()
        self.green_tick_shown = True
        self.magenta_tick_shown = False

    def show_line_info_dialog(self):

        def get_window_size_and_position():
            if self.line_info_dialog is None:
                font = QtGui.QFont()
                width = QtGui.QFontMetrics(font).width('='*120)
                self.line_info_dialog_width = width
                self.line_info_dialog_height = 470
                sG = QtGui.QApplication.desktop().screenGeometry()
                self.line_info_dialog_x = sG.width()-self.line_info_dialog_width
                self.line_info_dialog_y = 0
            else:
                self.line_info_dialog_width = self.line_info_dialog.width()
                self.line_info_dialog_height = self.line_info_dialog.height()
                self.line_info_dialog_x = self.line_info_dialog.pos().x()
                self.line_info_dialog_y = self.line_info_dialog.pos().y()
    
        def save_initial_plot_pars():
            self.init_line_num = self.line_info_box.text()
            self.init_ion = self.ion_box.text()
            self.init_xmin = self.xlim_min_box.text()
            self.init_xmax = self.xlim_max_box.text()
            self.init_y1min = self.y1lim_min_box.text()
            self.init_y1max = self.y1lim_max_box.text()
            self.init_y3min = self.y3lim_min_box.text()
            self.init_y3max = self.y3lim_max_box.text()
            self.init_legend_fontsize = self.sp.legend_fontsize
            self.init_legend_loc = self.sp.legend_loc
      
        def toggle_statusbar():
            self.showStatusBar = not self.showStatusBar
            statusBar.setVisible(self.showStatusBar)
         
        def redo_initial_plot():
            self.line_info_box.setText(self.init_line_num)
            self.ion_box.setText(self.init_ion)
            self.xlim_min_box.setText(self.init_xmin)
            self.xlim_max_box.setText(self.init_xmax)
            self.y1lim_min_box.setText(self.init_y1min)
            self.y1lim_max_box.setText(self.init_y1max)
            self.y3lim_min_box.setText(self.init_y3min)
            self.y3lim_max_box.setText(self.init_y3max)
            self.sp.legend_fontsize = self.init_legend_fontsize
            self.sp.legend_loc = self.init_legend_loc
            self.set_plot_limits_and_draw()
            #self.save_from_lim_boxes()
            #self.draw_ion()

        def do_reset():
            self.curr_line_num = self.init_line_num
            get_info(self.curr_line_num)
            fill_line_info_table()
            redo_initial_plot()

        def toggle_show_satellites():
            self.show_satellites = (self.show_satellites + 1)%3
            fill_line_info_table()
         
        def on_click():
            item = self.line_info_table.currentItem()
            row = item.row()
            col = item.column()
            s = item.text()
            if col == col_ion:
                ion = self.line_info_table.item(row, col).text()
                self.ion_box.setText(ion)
                self.draw_ion()
            if not self.isFloat(s):
                return               
            if col in [col_num, col_ref] and int(s) != 0:
                self.curr_line_num = s
                get_info(self.curr_line_num)
                self.line_info_box.setText(self.curr_line_num)
                fill_line_info_table()
         
        def on_doubleClick():
            item = self.line_info_table.currentItem()
            row = item.row()
            col = item.column()
            s = item.text()
            if col == col_ion:
                ion = self.line_info_table.item(row, col).text()
                self.ion_box.setText(ion)
                self.draw_ion()
            if not self.isFloat(s):
                return               
            if col in [col_num, col_ref] and int(s) != 0:
                self.curr_line_num = s
                get_info(self.curr_line_num)
                self.line_info_box.setText(self.curr_line_num)
                fill_line_info_table()
         
        def on_itemClicked():
            # to avoid blinking with itemSelectionChanged 
            item = self.line_info_table.currentItem()
            if item == self.selected_item:
                on_itemSelectionChanged()

        def on_itemSelectionChanged():

            if self.green_tick_shown:
                self.on_draw()
            self.green_tick_shown = False
            
            item = self.line_info_table.currentItem()
            if item == None:
                self.draw_ion()
                return
            self.selected_item = item
            row = item.row()
            col = item.column()
            s = item.text()
            l_shift_refline = np.float(self.sp.fieldStrFromLine(self.refline,'l_shift'))
            if col == col_wave:
                wavelength = np.float(s)
                ion = str(self.line_info_table.item(row, col_ion).text())
                line_num = str(self.line_info_table.item(row, col_num).text())
                max_wave = np.float(self.sp_max_box.text())
                min_wave = np.float(self.sp_min_box.text())
                if wavelength > min_wave and wavelength < max_wave:
                    l_shift = np.float(self.line_info_table.item(row, col_lshift).text())
                    wavelength = wavelength + l_shift + l_shift_refline
                    r =  (self.x_plot_lims[1] - self.x_plot_lims[0])/2
                    f = 0.05
                    if (wavelength < self.x_plot_lims[0] + f*r) or (wavelength > self.x_plot_lims[1] - f*r):
                        if wavelength-r < min_wave:
                            self.x_plot_lims = (min_wave-r*f, min_wave-r*f+2*r)
                        elif wavelength+r > max_wave:
                            self.x_plot_lims = (max_wave+r*f-2*r , max_wave+r*f)
                        else:
                            self.x_plot_lims = (wavelength-r,wavelength+r)
                        if not self.axes_fixed:
                            self.update_lim_boxes()
                        self.restore_axes()
                    self.plot_tick_at(wavelength, ion, line_num)
                elif wavelength == 1:
                    if str(self.line_info_table.item(row, col_ref).text()) == '0000000000000':
                        satellites = self.satellites
                    else:
                        satellites = self.sp.read_satellites(self.sp.phyat_file, int(line_num))
                        satellites = add_satellites_of_subreferences(satellites)
                    SelectedSatellites = []
                    max_wave = np.float(self.sp_max_box.text())
                    min_wave = np.float(self.sp_min_box.text())
                    for i in range(0, len(satellites)):
                        wavelength = np.float(self.sp.fieldStrFromLine(satellites[i],'lambda'))
                        if (wavelength > min_wave) and (wavelength < max_wave):
                            SelectedSatellites.append(satellites[i])
                    satellites = SelectedSatellites
                    self.plot_line_ticks_for(satellites, ion, line_num, self.refline)

        def isRefLine(line):
            s = self.sp.fieldStrFromLine(line,'ref').strip()
            if s == '0000000000000':
                return True
            else:
                return False

        def isSubRefLine(line):
            wavelength = np.float(self.sp.fieldStrFromLine(line,'lambda'))
            if not isRefLine(line) and (wavelength < 2.0):
                return True
            else:
                return False

        def fill_data(i, line, cat=''):
            if line == None:
                return
            editableCols = []
            if self.sp.get_conf('qt_allow_editing_lines', False):
                if cat == 'sat':
                    if do_cosmetics:
                        editableCols = ['l_shift', 'i_cor', 'profile', 'vitesse', 'comment']
                    else:
                        editableCols = []
                elif cat == 'subref':
                    if do_cosmetics:
                        editableCols = ['i_cor', 'comment']
                    else:
                        editableCols = []
                elif cat == 'ref':
                    editableCols = ['l_shift', 'i_cor', 'i_rel', 'profile', 'vitesse', 'comment']
            for j in range(0,len(fieldItems)):
                s = self.sp.fieldStrFromLine(line, fieldItems[j])
                s = s.strip()
                if j == col_ion:
                    if self.show_true_ions: 
                        s = self.sp.true_ion(s).replace('_',' ').strip()
                    isPseudoIon = self.sp.isPseudoIon(s)
                if j == fieldItems.index('proc'):                    
                    if isRefLine(line):
                        s = ''
                    elif isPseudoIon:
                        s = ''
                    else:
                        s = self.sp.process[s]
                item = QtGui.QTableWidgetItem(s)
                if fieldItems[j] in editableCols:
                    item.setBackgroundColor(self.editableCells_bg_color)
                else:
                    item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)
                    item.setBackgroundColor(self.readOnlyCells_bg_color)
                self.line_info_table.setItem(i,j,item)

        def fill_text(i, text):
            item = QtGui.QTableWidgetItem(text)
            item.setFlags(item.flags() ^ (QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled))
            item.setBackgroundColor(self.readOnlyCells_bg_color)
            item.setTextAlignment(QtCore.Qt.AlignBottom)
            item.setTextColor(QtCore.Qt.blue)
            self.line_info_table.setItem(i,0,item)
            self.line_info_table.setSpan(i,0,2,len(fieldItems))

        def add_satellites_of_subreferences(satellites):
            subref_list = []
            all_satellites = satellites
            for sat_line in satellites:
                if isSubRefLine(sat_line):
                    subref_list.append(sat_line)
            i = 0                    
            while i < len(subref_list):
                sat_line_num = self.sp.fieldStrFromLine(subref_list[i],'num')
                new_satellites = self.sp.read_satellites(self.sp.phyat_file, int(sat_line_num))
                for line in new_satellites:
                    if isSubRefLine(line):
                        subref_list.append(line)
                i += 1
                for line in new_satellites:
                    if not line in all_satellites:
                        all_satellites.append(line) 
            return all_satellites        

        def get_info(line_num):  
            line = None    
            refline = None    
            subrefline = None
            LineList = []
            if int(line_num) == 0:
                return
            while refline == None:
                refline = self.sp.read_line(self.sp.fic_model, int(line_num))
                if refline is None:
                    if do_cosmetics:
                        curr_line = self.sp.read_line(self.sp.fic_cosmetik, int(line_num))
                    else:
                        curr_line = None
                    if self.sp.cosmetic_line_ok(curr_line) is not True:
                        curr_line = None
                    if curr_line == None:
                        curr_line = self.sp.read_line(self.sp.phyat_file, int(line_num))
                    LineList.append(curr_line)
                    line_num = self.sp.fieldStrFromLine(curr_line,'ref')
            if len(LineList) > 0: 
                if isSubRefLine(LineList[0]):
                    subrefline = LineList[:1]
                else:
                    line = LineList[0]
                    if len(LineList) > 1:
                        subrefline = LineList[1:]
            if subrefline is not None:
                n_subref = len(subrefline)
            else:
                n_subref = 0
            subsatellites = []
            for k in range(0, n_subref):
                subsat = []
                subrefline_num = self.sp.fieldStrFromLine(subrefline[k], 'num')             
                subsat = self.sp.read_satellites(self.sp.phyat_file, int(subrefline_num))
                n_subsat = len(subsat)
                if do_cosmetics:
                    for i in range(0,n_subsat):
                        sat_line = subsat[i]
                        sat_line_num = int(self.sp.fieldStrFromLine(sat_line,'num'))
                        cosmetic_line = self.sp.read_line(self.sp.fic_cosmetik, sat_line_num)
                        if cosmetic_line is not None:
                            subsat[i] = cosmetic_line
                subsatellites = subsatellites + subsat
            subsatellites = add_satellites_of_subreferences(subsatellites)
            n_subsat = len(subsatellites)
            if refline is not None:
                refline_num = self.sp.fieldStrFromLine(refline,'num')
                satellites = self.sp.read_satellites(self.sp.phyat_file, int(refline_num))
                satellites = add_satellites_of_subreferences(satellites)
                n_sat = len(satellites)
                if do_cosmetics:
                    for i in range(0,n_sat):
                        sat_line = satellites[i]
                        sat_line_num = int(self.sp.fieldStrFromLine(sat_line,'num'))
                        cosmetic_line = self.sp.read_line(self.sp.fic_cosmetik, sat_line_num)
                        if cosmetic_line is not None:
                            satellites[i] = cosmetic_line
            else:
                n_sat = 0
            if line is None and refline is None:
                title = 'Error in line info dialog'
                msg = 'Line number not found.'
                QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
            self.line = line
            self.subrefline = subrefline
            self.refline = refline
            self.subsatellites = subsatellites
            self.satellites = satellites
            self.n_sat = n_sat 
            self.n_subsat = n_subsat
            self.n_subref = n_subref

        def do_sort(lines):
            waves = []
            for i in range(0,len(lines)):
                waves.append(self.sp.fieldStrFromLine(lines[i], 'lambda'))
            lines = [x for _,x in sorted(zip(waves,lines))]
            return lines
               
        def fill_line_info_table():  
            self.line_info_table.blockSignals(True)
            line = self.line
            subrefline = self.subrefline
            refline = self.refline
            subsatellites = self.subsatellites
            satellites = self.satellites
            n_sat = self.n_sat 
            n_subsat = self.n_subsat
            n_subref = self.n_subref
            SelectedSatellites = []
            SelectedSubSatellites = []
            if self.show_satellites == 0:
                n_sat = 0
                n_subsat = 0
            else:
                max_wave = np.float(self.sp_max_box.text())
                min_wave = np.float(self.sp_min_box.text())
                for i in range(0, len(satellites)):
                    wavelength = np.float(self.sp.fieldStrFromLine(satellites[i],'lambda'))
                    if self.show_satellites == 2 or \
                       (self.show_satellites == 1 and (wavelength > min_wave) and (wavelength < max_wave)):
                        SelectedSatellites.append(satellites[i])
                for i in range(0, len(subsatellites)):
                    wavelength = np.float(self.sp.fieldStrFromLine(subsatellites[i],'lambda'))
                    if self.show_satellites == 2 or \
                      (self.show_satellites == 1 and (wavelength > min_wave) and (wavelength < max_wave)):
                        SelectedSubSatellites.append(subsatellites[i])
                n_sat = len(SelectedSatellites)
                n_subsat = len(SelectedSubSatellites)
            self.line_info_table.clearContents()
            self.line_info_table.setRowCount(n_sat+n_subsat+20)
            self.line_info_table.clearSpans()
            k = 0
            sat_list = []
            if line is not None:
                fill_text(k,'Line:')
                k += 2
                fill_data(k, line, 'sat')
                k += 1
            if subrefline is not None:
                fill_text(k,'Subreference line:')
                k += 2
                for i in range(0,n_subref):
                    fill_data(k, subrefline[i], 'subref')
                    k += 1
                if n_subsat > 0:
                    SelectedSubSatellites = do_sort(SelectedSubSatellites)
                    fill_text(k, str(n_subsat) + ' satellites:')
                    sat_list.append([k,n_subsat])
                    k += 2
                    for i in range(0,n_subsat):
                        if isSubRefLine(SelectedSubSatellites[i]):
                            fill_data(k+i, SelectedSubSatellites[i], 'subref')
                        else:
                            fill_data(k+i, SelectedSubSatellites[i], 'sat')
                    k += n_subsat
            fill_text(k,'Reference line:')
            k += 2
            fill_data(k, refline, 'ref')
            k += 1
            if n_sat > 0:
                SelectedSatellites = do_sort(SelectedSatellites)
                fill_text(k, str(n_sat) + ' satellites:')
                sat_list.append([k,n_sat])               
                k += 2
                for i in range(0,n_sat):
                    if isSubRefLine(SelectedSatellites[i]):
                        fill_data(k+i, SelectedSatellites[i], 'subref')
                    else:
                        fill_data(k+i, SelectedSatellites[i], 'sat')
                k += n_sat
            self.line_info_table.setRowCount(k)
            self.line_info_table.resizeColumnsToContents()
            self.line_info_table.resizeRowsToContents()
            self.line_info_table.blockSignals(False)
            
            self.line_info_table.blockSignals(True)
            if self.show_satellites == 1:
                s0 = ' (in the synthesis range)'
            elif self.show_satellites == 2:
                s0 = ' (in the entire database and including subreferences)'
            else:
                s0 = ''
            for i in sat_list:
                k = i[0]
                n = i[1]
                fill_text(k, str(n) + ' satellites:' + s0)
            self.line_info_table.blockSignals(False)
        
        def on_itemChanged():
            self.line_info_table.blockSignals(True)
            item = self.line_info_table.currentItem()
            if not (item.flags() & QtCore.Qt.ItemIsEditable):
                self.line_info_table.blockSignals(False)
                return
            row = item.row()
            col = item.column()
            s = str(item.text())
            value = self.rightFormat(s, fieldItems[col])
            if value != None:
                self.line_info_table.setItem(row, col, QtGui.QTableWidgetItem(value.strip()))
                self.line_info_table.item(row, col).setBackgroundColor(self.editableCells_bg_color)
                save_change(row,col)
            else:
                self.line_info_table.item(row, col).setBackgroundColor(QtGui.QColor('red'))
                title = 'Invalid format for the ' + self.sp.field_tip[fieldItems[col]]
                s0 = self.sp.field_format[fieldItems[col]]
                s0 = s0[2:-1]
                msg = "'" + s + "' can not be converted into the proper field format: " + s0
                if col == self.sp.fields.index('vitesse'):
                    msg = msg + '\nor it is not a positive number.'
                QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
            get_info(self.curr_line_num)
            fill_line_info_table()
            self.line_info_table.blockSignals(False)
       
        def get_line_from_table(row):
            line = ' '*85
            jList = range(0,len(fieldItems))
            jList.remove(col_proc)
            for j in jList:
                s = self.line_info_table.item(row,j).text()
                width = self.sp.field_width[fieldItems[j]]
                align = self.sp.field_align[fieldItems[j]]
                pos = self.sp.field_pos[fieldItems[j]]
                s = '{:{a}{w}s}'.format(s, a=align, w=width)
                line = line[:pos] + s + line[pos:]
            line = line.rstrip()
            return line
       
        def save_change(row, col):
            line = get_line_from_table(row)
            if isRefLine(line):
                filename = self.sp.fic_model
            else:
                filename = self.sp.fic_cosmetik
            self.sp.replace_line(filename, line)
            if col != self.sp.fields.index('comment') and \
               self.sp.get_conf('qt_update_after_editing_lines', False):
                self.adjust()
                self.nearbyLines = self.sp.get_nearby_lines(self.cursor_w1, self.cursor_w2, do_print=False)
                if self.nearbyLines is not None and self.nearbyLines_dialog.isVisible():
                    self.fill_nearbyLines_table()
                
        def init_lines():
            self.line = None
            self.subrefline = None
            self.refline = None
            self.subsatellites = []
            self.satellites = []
            self.n_sat = 0 
            self.n_subsat = 0
            self.n_subref = 0

        statusBar = QtGui.QStatusBar()
        s = 'Click on \"Satellites\" to cycle the tri-state display of satellite lines:\n' \
            '   1 - The satellite lines in the spectral range of the synthesis are shown; \n' \
            '   2 - All satellite lines (including subreference lines and lines outside the spectral range of the synthesis) are shown. \n' \
            '   3 - No satellite line is shown; \n' \
            'Double-click on a line number to show the data for that line. \n' \
            'Double-click on an ion to plot line ticks and spectrum for that single ion. \n' \
            'Select or click on a wavelength to draw a tick at that position and recenter the spectrum if necessary. \n' \
            'Click on \"Reset\" to return to the original line and plot settings. \n' \
            'The green fields are editable.'
        statusBar.addWidget(QtGui.QLabel(s),1)
        self.showStatusBar = False
        statusBar.setVisible(self.showStatusBar)
        self.show_satellites = 1
        get_window_size_and_position()

        if self.line_info_dialog is not None:
            self.line_info_dialog.close()
            self.line_info_table.close()

        self.line_info_dialog = QtGui.QDialog()
        self.line_info_dialog.resize(self.line_info_dialog_width,self.line_info_dialog_height)
        self.line_info_dialog.move(self.line_info_dialog_x,self.line_info_dialog_y)
        self.line_info_table = QtGui.QTableWidget()
        fieldItems = self.sp.fields
        fieldNames = [ self.sp.field_abbr[item] for item in fieldItems ]
        col_num = fieldItems.index('num')
        col_ion = fieldItems.index('id')
        col_wave = fieldItems.index('lambda')
        col_proc = fieldItems.index('proc')
        col_lshift = fieldItems.index('l_shift')
        col_irel = fieldItems.index('i_rel')
        col_icor = fieldItems.index('i_cor')
        col_ref = fieldItems.index('ref')
        col_prof =  fieldItems.index('profile')              
        col_vel =  fieldItems.index('vitesse')              
        col_comm = fieldItems.index('comment')                          
        self.line_info_table.setColumnCount(len(fieldItems))
        self.line_info_table.setHorizontalHeaderLabels(fieldNames)
        if self.enable_tooltips_action.isChecked():
            for j in range(0,len(fieldItems)):
                self.line_info_table.horizontalHeaderItem(j).setToolTip(self.sp.field_tip[fieldItems[j]])
        self.line_info_table.horizontalHeaderItem(col_vel).setText(u'\u0394v (factor)')
        if self.enable_tooltips_action.isChecked():
            s = 'For a reference line, it is the thermal broadening parameter, in km/s. \n' \
                'For satellite line, it is the dimensionless correction factor for the thermal broadening parameter with respect to the reference line.'
            self.line_info_table.horizontalHeaderItem(col_vel).setToolTip(s)
        self.line_info_table.horizontalHeaderItem(col_comm).setTextAlignment(QtCore.Qt.AlignLeft)
        self.line_info_table.horizontalHeaderItem(col_comm).setText('  comment')
        init_lines()

        do_cosmetics = self.sp.get_conf('do_cosmetik')
        save_initial_plot_pars()
        self.curr_line_num = self.line_info_box.text()
        get_info(self.curr_line_num)
        fill_line_info_table()
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Help|
                                                QtGui.QDialogButtonBox.Close|
                                                QtGui.QDialogButtonBox.Reset|
                                                QtGui.QDialogButtonBox.Apply)
        self.buttonBox.button(QtGui.QDialogButtonBox.Apply).setText("Satellites")
        if self.enable_tooltips_action.isChecked():
            self.buttonBox.button(QtGui.QDialogButtonBox.Apply).setToolTip("Click to toggle the satellite lines")
        self.buttonBox.button(QtGui.QDialogButtonBox.Apply).clicked.connect(toggle_show_satellites)
        s = "Click to return to the initial states of the line info dialog and figures"
        if self.enable_tooltips_action.isChecked():
            self.buttonBox.button(QtGui.QDialogButtonBox.Reset).setToolTip(s)
        self.buttonBox.button(QtGui.QDialogButtonBox.Reset).clicked.connect(do_reset)
        self.buttonBox.button(QtGui.QDialogButtonBox.Help).clicked.connect(toggle_statusbar)
        self.buttonBox.rejected.connect(self.line_info_dialog.close)
        self.line_info_table.doubleClicked.connect(on_doubleClick)
        self.line_info_table.itemChanged.connect(on_itemChanged)
        self.selected_item = None
        self.line_info_table.itemSelectionChanged.connect(on_itemSelectionChanged)
        self.line_info_table.itemClicked.connect(on_itemClicked)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.line_info_table)
        vbox.addWidget(self.buttonBox)
        vbox.addWidget(statusBar)
        self.line_info_dialog.setLayout(vbox)
        self.line_info_dialog.setWindowTitle('line info dialog')
        self.line_info_dialog.setWindowModality(QtCore.Qt.NonModal)
        self.line_info_dialog.show()
    
    def fill_nearbyLines_table(self):
        if self.nearbyLines is None or self.nearbyLines_table is None:
            return
        k = self.sp.get_conf('diff_lines_by')
        fieldItems = self.sp.fields
        jList = range(0,len(fieldItems))
        jProc = fieldItems.index('proc')
        jList.remove(jProc)
        if self.nearbyDialogFilterIsActive:
            #selected_ions = self.sp.get_conf('selected_ions')
            selected_ions = self.nearbyLines_selected_ions
            
            selected_true_ions = [self.sp.true_ion(ion) for ion in selected_ions]
            nearbyLines = []
            for line in self.nearbyLines:
                ion = str(line[fieldItems.index('id')]).strip()
                true_ion = self.sp.true_ion(ion)
                selectThisIon = (( ion in selected_ions or true_ion in selected_ions ) and k == 1) or (true_ion in selected_true_ions and k != 1)
                if selectThisIon:
                    nearbyLines.append(line)
        else:
            nearbyLines = self.nearbyLines
        self.nearbyLines_table.setRowCount(len(nearbyLines))                 
        for i in range(0,len(nearbyLines)):
            ion = self.sp.true_ion(nearbyLines[i][fieldItems.index('id')])
            for j in jList:
                if j > jProc:
                    k = j - 1
                else:
                    k = j
                fmt = self.sp.field_format[fieldItems[j]]
                s = fmt.format(nearbyLines[i][k])
                s = str(s).strip()
                if j == fieldItems.index('num'):
                    if self.sp.isPseudoIon(ion):
                        proc_str = ''
                    else:
                        proc_str = self.sp.process[s[-9]]
                if j == fieldItems.index('id'):
                    if self.show_true_ions: 
                        s = self.sp.true_ion(s).replace('_',' ').strip()    
                item = QtGui.QTableWidgetItem(s)
                item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)
                item.setBackgroundColor(self.readOnlyCells_bg_color)
                self.nearbyLines_table.setItem(i,j,item)
            item = QtGui.QTableWidgetItem(proc_str)
            item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)
            item.setBackgroundColor(self.readOnlyCells_bg_color)
            self.nearbyLines_table.setItem(i,jProc,item)
        self.nearbyLines_table.resizeColumnsToContents()
        self.nearbyLines_table.resizeRowsToContents()
        self.nearbyLines_table.clearSelection()
                
    def show_nearbyLines_dialog(self):

        def get_window_size_and_position():
            if self.nearbyLines_dialog is None:
                font = QtGui.QFont()
                width = QtGui.QFontMetrics(font).width('='*120)
                self.nearbyLines_dialog_width = width
                self.nearbyLines_dialog_height = 470
                sG = QtGui.QApplication.desktop().screenGeometry()
                self.nearbyLines_dialog_x = sG.width()-self.nearbyLines_dialog_width
                self.nearbyLines_dialog_y = sG.height()-self.nearbyLines_dialog_height
            else:
                self.nearbyLines_dialog_width = self.nearbyLines_dialog.width()
                self.nearbyLines_dialog_height = self.nearbyLines_dialog.height()
                self.nearbyLines_dialog_x = self.nearbyLines_dialog.pos().x()
                self.nearbyLines_dialog_y = self.nearbyLines_dialog.pos().y()

        def do_reset():
            self.curr_line_num = self.init_nearby_line_num
            #get_info(self.curr_line_num)
            #fill_line_info_table()
            self.nearbyDialogFilterIsActive = True
            #self.nearbyLines_selected_ions = []
            toggle_filter()
            redo_initial_plot()

        def toggle_filter():
            self.nearbyLines_selected_ions = []
            if not self.nearbyDialogFilterIsActive:
                get_selected_ions()
                if len(self.nearbyLines_selected_ions) > 0:
                    self.nearbyDialogFilterIsActive = True
                    self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.RestoreDefaults).setStyleSheet('background-color:red;')
                    self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.RestoreDefaults).setText('deactivate ion filter')        

                else:        
                    QtGui.QMessageBox.critical(self, 'nearby lines dialog: ion filter', 'No ion selected.', QtGui.QMessageBox.Ok )
            else:
                self.nearbyDialogFilterIsActive = False
                self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.RestoreDefaults).setStyleSheet('')
                self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.RestoreDefaults).setText('Filter selected ions')        
            self.fill_nearbyLines_table()

        def save_initial_plot_pars():
            self.init_nearby_line_num = self.line_info_box.text()
            self.init_nearby_ion = self.ion_box.text()
            self.init_nearby_xmin = self.xlim_min_box.text()
            self.init_nearby_xmax = self.xlim_max_box.text()
            self.init_nearby_y1min = self.y1lim_min_box.text()
            self.init_nearby_y1max = self.y1lim_max_box.text()
            self.init_nearby_y3min = self.y3lim_min_box.text()
            self.init_nearby_y3max = self.y3lim_max_box.text()
            self.init_nearby_legend_fontsize = self.sp.legend_fontsize
            self.init_nearby_legend_loc = self.sp.legend_loc

        def redo_initial_plot():
            #self.line_info_box.setText(self.init_line_num)
            self.ion_box.setText(self.init_nearby_ion)
            self.xlim_min_box.setText(self.init_nearby_xmin)
            self.xlim_max_box.setText(self.init_nearby_xmax)
            self.y1lim_min_box.setText(self.init_nearby_y1min)
            self.y1lim_max_box.setText(self.init_nearby_y1max)
            self.y3lim_min_box.setText(self.init_nearby_y3min)
            self.y3lim_max_box.setText(self.init_nearby_y3max)
            self.sp.legend_fontsize = self.init_nearby_legend_fontsize
            self.sp.legend_loc = self.init_nearby_legend_loc
            self.set_plot_limits_and_draw()

        def toggle_statusbar():
            self.showStatusBar = not self.showStatusBar
            statusBar.setVisible(self.showStatusBar)

        def on_doubleClick():
            item = self.nearbyLines_table.currentItem()
            row = item.row()
            col = item.column()
            if col in [col_num, col_ref]:
                self.line_info_box.setText(item.text())
                self.show_line_info_dialog()
            elif col == col_ion:
                self.ion_box.setText(item.text())
                self.draw_ion()

        def on_itemClicked():
            # to avoid blinking with itemSelectionChanged 
            item = self.nearbyLines_table.currentItem()
            if item == self.selected_item:
                on_itemSelectionChanged()

        def on_itemSelectionChanged():
            item = self.nearbyLines_table.currentItem()
            self.selected_item = item
            row = item.row()
            col = item.column()
            if col == col_wave:
                wavelength = np.float(item.text())
                l_shift = np.float(self.nearbyLines_table.item(row,col_lshift).text())
                wavelength = wavelength + l_shift
                line_num = str(self.nearbyLines_table.item(row,col_num).text())
                ion = str(self.nearbyLines_table.item(row,col_ion).text())
                max_wave = np.float(self.sp_max_box.text())
                min_wave = np.float(self.sp_min_box.text())
                r =  (self.x_plot_lims[1] - self.x_plot_lims[0])/2
                f = 0.05
                if (wavelength < self.x_plot_lims[0] + f*r) or (wavelength > self.x_plot_lims[1] - f*r):
                    if wavelength-r < min_wave:
                        self.x_plot_lims = (min_wave-r*f, min_wave-r*f+2*r)
                    elif wavelength+r > max_wave:
                        self.x_plot_lims = (max_wave+r*f-2*r , max_wave+r*f)
                    else:
                        self.x_plot_lims = (wavelength-r,wavelength+r)
                    if not self.axes_fixed:
                        self.update_lim_boxes()
                    self.restore_axes()

                self.plot_tick_at(wavelength, ion, line_num)
            else:
                if self.green_tick_shown:
                    self.on_draw()
                    self.green_tick_shown = False

        def do_header_clicked(col):
            if col == col_ion:
                self.toggle_show_true_ions()
                self.fill_nearbyLines_table()
        
        def do_header_doubleClicked(col):
            sort = fieldItems[col]
            if sort == self.nearbyLines_sort_by:
                self.nearbyLines_sort_reverse = not self.nearbyLines_sort_reverse
            else:
                self.nearbyLines_sort_reverse = False
                self.nearbyLines_sort_by = sort
            self.sort_nearbyLines(sort, self.nearbyLines_sort_reverse)
            self.fill_nearbyLines_table()

        def get_selected_ions():
            selectedItems = self.nearbyLines_table.selectedItems()
            selected_ions = []
            for item in selectedItems:
                col = item.column()
                if col == col_ion:
                    ion = str(item.text())
                    if not ion in selected_ions:
                        selected_ions.append(ion)
            if len(selected_ions) > 0:
                self.nearbyLines_selected_ions = selected_ions
            else:
                #self.nearbyLines_selected_ions = self.sp.get_conf('selected_ions')
                self.nearbyLines_selected_ions = []
                            
        def do_selection():
            selectedItems = self.nearbyLines_table.selectedItems()
            selected_ions = []
            selected_lines = []
            for item in selectedItems:
                col = item.column()
                if col == col_ion:
                    ion = str(item.text())
                    if not ion in selected_ions:
                        selected_ions.append(ion)
                if col in [col_num, col_ref]:
                    line = item.text()
                    selected_lines.append(line)
            if len(selected_ions) > 0:
                s = ''
                for ion in selected_ions:
                    s = s + ion + ', '
                s = s[:-2]
                self.ion_box.setText(s)
                self.draw_ion()
            if len(selected_lines) > 0:
                s = selected_lines[0]
                self.line_info_box.setText(s)
                self.line_info()

        get_window_size_and_position()
        self.nearbyLines_dialog = QtGui.QDialog()
        self.nearbyLines_dialog.resize(self.nearbyLines_dialog_width, self.nearbyLines_dialog_height)
        self.nearbyLines_dialog.move(self.nearbyLines_dialog_x,self.nearbyLines_dialog_y)
        statusBar = QtGui.QStatusBar()
        s = 'Double-click on a line number (or select the line number and press \"Apply\") to show line info dialog. \n' \
            'Double-click on an ion to plot line ticks and spectrum for that single ion. \n' \
            'Click or select a wavelength to draw a tick at that position. \n' \
            'Select multiple ions (using click, Shift+click, and Ctrl+click) and press \"Plot selected ions\" plot line ticks and spectra for a list of ions. \n' \
            'Click on the ion header to select all ions. \n' \
            'Double-click on a column header to sort the table; Double-click again to toggle between ascending and descending order. \n' \
            'Click on \"Reset\" to return to the original selected ions and plot settings. \n' \
            'Click on \"Filter selected ions\" to activate/deactivate ion selection.'
        statusBar.addWidget(QtGui.QLabel(s),1)
        self.showStatusBar = False
        statusBar.setVisible(self.showStatusBar)
        self.nearbyLines_table = QtGui.QTableWidget()   
        self.nearbyLines_table.setRowCount(len(self.nearbyLines))
        fieldItems = self.sp.fields
        fieldNames = [ self.sp.field_abbr[item] for item in fieldItems ]
        col_num = fieldItems.index('num')
        col_ion = fieldItems.index('id')
        col_wave = fieldItems.index('lambda')
        col_proc = fieldItems.index('proc')
        col_lshift = fieldItems.index('l_shift')
        col_irel = fieldItems.index('i_rel')
        col_icor = fieldItems.index('i_cor')
        col_ref = fieldItems.index('ref')
        col_prof =  fieldItems.index('profile')              
        col_vel =  fieldItems.index('vitesse')              
        col_comm = fieldItems.index('comment')                          
        self.nearbyLines_table.setColumnCount(len(fieldNames))
        self.nearbyLines_table.setHorizontalHeaderLabels(fieldNames)
        if self.enable_tooltips_action.isChecked():
            for j in range(0,len(fieldItems)):
                self.nearbyLines_table.horizontalHeaderItem(j).setToolTip(self.sp.field_tip[fieldItems[j]])
        self.nearbyLines_table.horizontalHeaderItem(col_comm).setTextAlignment(QtCore.Qt.AlignLeft)
        self.nearbyLines_table.horizontalHeaderItem(col_vel).setText(u'\u0394v')
        if self.enable_tooltips_action.isChecked():
            s = u'\u0394v is the thermal broadening parameter of the line, in km/s. \n' \
                 'For a single Gaussian profile, it is the half-width of the line at the level of 1/e of the peak, \n' \
                 'related to the full-width at half maximum and the Gaussian standard deviation by:\n\n' \
                u'     \u0394v = FWHM/(2(ln2)^\u00BD) = FWHM/1.665\n' \
                u'     \u0394v = \u221A2 \u03C3\n' 
            self.nearbyLines_table.horizontalHeaderItem(col_vel).setToolTip(s)
        self.nearbyLines_table.horizontalHeaderItem(col_comm).setText('  comment')
        #self.nearbyDialogFilterIsActive = False
        self.fill_nearbyLines_table()
        save_initial_plot_pars()
        self.buttonBox_nearbyLines = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Help|
                                                QtGui.QDialogButtonBox.Reset|
                                                QtGui.QDialogButtonBox.RestoreDefaults|
                                                QtGui.QDialogButtonBox.Apply|
                                                QtGui.QDialogButtonBox.Close)
        self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.RestoreDefaults).setText('Filter selected ions')        
        self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.Apply).setText('Plot selected ions')
        self.buttonBox_nearbyLines.rejected.connect(self.nearbyLines_dialog.close)
        self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.Apply).clicked.connect(do_selection)
        self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.Help).clicked.connect(toggle_statusbar)
        self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.Reset).clicked.connect(do_reset)
        self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.RestoreDefaults).clicked.connect(toggle_filter)
        self.nearbyLines_table.doubleClicked.connect(on_doubleClick)
        self.nearbyLines_table.itemSelectionChanged.connect(on_itemSelectionChanged)
        self.nearbyLines_table.itemClicked.connect(on_itemClicked)
        self.nearbyLines_table.verticalHeader().sectionDoubleClicked.connect(do_selection)
        #self.nearbyLines_table.horizontalHeader().sectionClicked.connect(do_header_clicked)
        self.nearbyLines_table.horizontalHeader().sectionDoubleClicked.connect(do_header_doubleClicked)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.nearbyLines_table)
        vbox.addWidget(self.buttonBox_nearbyLines)
        vbox.addWidget(statusBar)
        self.nearbyLines_dialog.setLayout(vbox)
        s = 'nearby line dialog: list of lines between {0:.2f} and {1:.2f} angstroms'.format(self.sp.cursor_w1, self.sp.cursor_w2)
        self.nearbyLines_dialog.setWindowTitle(s)
        self.nearbyLines_dialog.setWindowModality(QtCore.Qt.NonModal)
        self.cursor_w1 = self.sp.cursor_w1
        self.cursor_w2 = self.sp.cursor_w2
        if self.nearbyDialogFilterIsActive:
            self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.RestoreDefaults).setStyleSheet('background-color:red;')
        else:
            self.buttonBox_nearbyLines.button(QtGui.QDialogButtonBox.RestoreDefaults).setStyleSheet('')
        self.nearbyLines_dialog.show()

    def cont_dialog(self):

        Pars = [ ( 'cont_unred'     , 'Set to True if reddening is to be applied to the continuum' ),
                 ( 'cont_edens'     , u'Electron density, in cm\u207B\u00B3' ), 
                 ( 'cont_hi_t'      , 'Temperature for the H I continuum, in K' ),       
                 ( 'cont_hi_i'      , u'Intensity of the H I continuum (in theory, intensity of H\u03B2)' ),       
                 ( 'cont_hei_t'     , 'Temperature for the He I continuum, in K' ),      
                 ( 'cont_hei_i'     , 'Intensity of the He I continuum (in theory, intensity of He I 4471)' ),      
                 ( 'cont_heii_t'    , 'Temperature for the He II continuum, in K' ),     
                 ( 'cont_heii_i'    , 'Intensity of the He II continuum (in theory, intensity of He I 4686)' ),     
                 ( 'cont_bb_t'      , 'Temperature of the blackbody continuum, in K' ),      
                 ( 'cont_bb_i'      , 'Intensity of the blackbody continuum' ),      
                 ( 'cont_pl_alpha'  , u'Index \u03B1 of the power-law continuum F = I*(\u03BB/5000 \u212B)**\u03B1' ),   
                 ( 'cont_pl_i'      , 'Intensity I of the power-law continuum' ),      
                 ( 'cont_in_lambda' , 'True (False) to interpolate continuum using list of wavelengths (pixels)' ),  
                 ( 'cont_lambda'    , 'List of wavelenghs of the interpolated continuum' ),    
                 ( 'cont_pix'       , 'List of pixels of the interpolated continuum' ),       
                 ( 'cont_intens'    , 'List of intensities of the interpolated continuum' ) ]
      
        def toggle_statusbar():
            self.showStatusBar = not self.showStatusBar
            statusBar.setVisible(self.showStatusBar)

        def get_window_size_and_position():
            if self.cont_pars_dialog is None:
                self.cont_pars_dialog_width = 800
                self.cont_pars_dialog_height = 460
                sG = QtGui.QApplication.desktop().screenGeometry()
                self.cont_pars_dialog_x = sG.width()-self.cont_pars_dialog_width
                self.cont_pars_dialog_y = sG.height()-self.cont_pars_dialog_height
                self.cont_pars_dialog_x = 0
                self.cont_pars_dialog_y = 0
            else:
                self.cont_pars_dialog_width = self.cont_pars_dialog.width()
                self.cont_pars_dialog_height = self.cont_pars_dialog.height()
                self.cont_pars_dialog_x = self.cont_pars_dialog.pos().x()
                self.cont_pars_dialog_y = self.cont_pars_dialog.pos().y()

        def set_conf_from_table(row):
            s = str(self.table.item(row,1).text())
            value = self.ConvStrToValidTypes(s)
            if value != None:
                self.sp.set_conf(Pars[row][0], value)
                self.table.setItem(row, 1, QtGui.QTableWidgetItem(str(value)))
            else:
                self.table.setItem(row, 1, QtGui.QTableWidgetItem('Error in ' + s))

        def on_itemChanged():
            self.table.blockSignals(True)
            item = self.table.currentItem()
            row = item.row()
            s = str(item.text())
            value = self.ConvStrToValidTypes(s)
            if value != None:
                self.sp.set_conf(Pars[row][0], value)
                self.table.setItem(row, 1, QtGui.QTableWidgetItem(str(value)))
                self.table.item(row, 1).setBackgroundColor(self.editableCells_bg_color)
                self.cont_par_changed = True
            else:
                self.table.setItem(row, 1, QtGui.QTableWidgetItem('Error in ' + s))
                self.table.item(row, 1).setBackgroundColor(QtGui.QColor('red'))
            self.table.blockSignals(False)
 
        get_window_size_and_position()
        self.cont_pars_dialog = QtGui.QDialog()
        self.cont_pars_dialog.resize(self.cont_pars_dialog_width, self.cont_pars_dialog_height)
        self.cont_pars_dialog.move(self.cont_pars_dialog_x, self.cont_pars_dialog_y) 
        statusBar = QtGui.QStatusBar()
        s = 'Click on \"Save\" to write the continuum parameters to a file. \n' \
            'Click on \"Update\" to adjust the synthesis to the changes in the continuum parameters. \n' \
            'The green fields are editable.'
        statusBar.addWidget(QtGui.QLabel(s),1)
        self.showStatusBar = False
        statusBar.setVisible(self.showStatusBar)
        self.table = QtGui.QTableWidget()   
        self.table.setRowCount(len(Pars))
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels([ 'parameter', 'value', 'help' ])
        for j in range(0,len(Pars)):
            item = QtGui.QTableWidgetItem(Pars[j][0])
            item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)
            item.setBackgroundColor(self.readOnlyCells_bg_color)
            self.table.setItem(j,0,item)
            item = QtGui.QTableWidgetItem(str(self.sp.get_conf(Pars[j][0])))
            item.setBackgroundColor(self.editableCells_bg_color)
            self.table.setItem(j,1,item)
            item = QtGui.QTableWidgetItem(Pars[j][1])
            item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)
            item.setBackgroundColor(self.readOnlyCells_bg_color)
            self.table.setItem(j,2,item)
        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()
        if self.table.columnWidth(1) > 300:
            self.table.setColumnWidth(1,300) 
        self.table.itemChanged.connect(on_itemChanged)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Help|
                                                QtGui.QDialogButtonBox.Save|
                                                QtGui.QDialogButtonBox.Apply|
                                                QtGui.QDialogButtonBox.Close)
        self.buttonBox.button(QtGui.QDialogButtonBox.Help).setDefault(True)
        self.buttonBox.button(QtGui.QDialogButtonBox.Apply).setText('Update')
        if self.enable_tooltips_action.isChecked():
            self.buttonBox.button(QtGui.QDialogButtonBox.Apply).setToolTip('Click to update synthesis with changes in the continuum parameters.')
        self.buttonBox.button(QtGui.QDialogButtonBox.Apply).clicked.connect(self.adjust)
        self.buttonBox.rejected.connect(self.cont_pars_dialog.close)
        self.buttonBox.button(QtGui.QDialogButtonBox.Save).clicked.connect(self.save_cont_pars)
        self.buttonBox.button(QtGui.QDialogButtonBox.Help).clicked.connect(toggle_statusbar)
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(self.table)
        vbox.addWidget(self.buttonBox)
        vbox.addWidget(statusBar)
        self.cont_pars_dialog.setLayout(vbox)
        self.cont_pars_dialog.setWindowTitle('Continuum parameters')
        self.cont_pars_dialog.show()

    def get_line_tick_lim(self, line_tick_pos):
        if line_tick_pos == 1:
            y1 = 0.43
            y2 = 0.57
        else:
            if line_tick_pos == 2:
                y1 = 0.05
                y2 = 0.19
            else:
                y1 = 0.81
                y2 = 0.95
        return y1, y2

    def getTickPosOfSelectedLine(self):
        posTick = self.sp.get_conf('line_tick_pos_selectedLine',3)
        if posTick not in [0,1,2]:
            posOtherTicks = self.sp.get_conf('line_tick_pos')
            if posTick == 4:
                if posOtherTicks == 2:
                   posTick = 0
                else:
                    posTick = 2
            else:
                posTick = posOtherTicks
        return posTick    
    
    def plot_line_ticks_for(self, satellites, ion, line_num, refline):
        k = self.sp.get_conf('line_tick_ax')    
        if not (k == 1 and self.residual_GroupBox.isChecked()):
            k = 0
        posTick = self.getTickPosOfSelectedLine()
        y1, y2 = self.get_line_tick_lim(posTick)
        if len(satellites) > 0:
            if ( k == 0 ):
                self.sp.plot_line_ticks_for(satellites, ion, line_num, refline, self.axes, y1, y2, self.x_plot_lims[0], self.x_plot_lims[1], self.addGreenTickToLegend)
            elif ( k == 1 ):
                self.sp.plot_line_ticks_for(satellites, ion, line_num, refline, self.axes3, y1, y2, self.addGreenTickToLegend)
            elif ( k == 2 ):
                self.sp.plot_line_ticks_for(satellites, ion, line_num, refline, self.axes2, 0.2, 0.8, self.addGreenTickToLegend)     
            self.green_tick_shown = True
        self.canvas.draw()                
                             
    def on_draw(self, show_legend=True):
        log_.debug('Entering on_drawn', calling=self.calling)
        if self.sp is None:
            log_.debug('Np sp in on_drawn', calling=self.calling)
            return
        
        if self.axes is None:
            log_.debug('Calling make_axes from on_draw (self.axes is None)', calling=self.calling)
            self.call_on_draw=False
            self.make_axes()
            self.init_axes()
            log_.debug('back from make_axes from on_draw', calling=self.calling)
            self.call_on_draw=True
        
        if self.do_save:
            self.save_axes()
        self.axes.cla()
        self.sp.plot_ax1(self.axes, show_legend)

        k = self.sp.get_conf('line_tick_ax')
        if self.show_line_ticks_action.isChecked() and ( k == 0 ):
            y1, y2 = self.get_line_tick_lim(self.sp.get_conf('line_tick_pos'))
            self.sp.plot_line_ticks(self.axes, y1, y2, self.x_plot_lims[0], self.x_plot_lims[1], show_legend=show_legend)

        if self.sp.get_conf('cont_plot', False):
            self.sp.plot_conts(self.axes)
        
        if self.residual_GroupBox.isChecked():
            self.axes3.cla()
            self.sp.plot_ax3(self.axes3, show_legend)
            if self.show_line_ticks_action.isChecked() and ( k == 1 ):
                y1, y2 = self.get_line_tick_lim(self.sp.get_conf('line_tick_pos'))
                self.sp.plot_line_ticks(self.axes3, y1, y2)

        if self.show_line_ticks_action.isChecked() and ( k == 2 ):
            self.axes2.cla()
            # self.sp.plot_ax2(self.axes2)
            self.sp.plot_line_ticks(self.axes2, 0.2, 0.8)
        
        if self.residual_GroupBox.isChecked():
            self.axes3.set_xlabel(r'Wavelength ($\AA$)')
            self.axes3.set_ylabel(r'Residual')
        #elif self.show_line_ticks_action.isChecked() and self.sp.get_conf(') and self.axes2 is not None:
        elif self.show_line_ticks_action.isChecked() and ( k == 2 ):
            self.axes2.set_xlabel(r'Wavelength ($\AA$)')
        else:
            self.axes.set_xlabel(r'Wavelength ($\AA$)')
        self.axes.set_ylabel(r'F$_\lambda$')
        
        self.restore_axes()
        # self.update_lim_boxes()
        if self.adjust_fig_action.isChecked(): 
            plt.tight_layout(0.1)
        self.canvas.draw()
        self.statusBar().showMessage('Redraw is finished.', 4000) 
        log_.debug('Exit on_drawn', calling=self.calling)
        self.magenta_tick_shown = False
        
    def show_lines_clicked(self):
        if self.lineIDs_GroupBox.isChecked():
            self.show_line_ticks_action.setChecked(True)
            self.plot_lines_action.setChecked(True)
            self.sp.set_conf('plot_lines_of_selected_ions', True)
            self.set_ion()
        else:  
            self.show_line_ticks_action.setChecked(False)
            self.plot_lines_action.setChecked(False)
            self.sp.set_conf('plot_lines_of_selected_ions', False)
        self.make_axes()
        
    def line_tick_color_clicked(self):
        color = QtGui.QColorDialog.getColor()
        self.sp.set_conf('line_tick_color', str(color.name()))
        if self.show_line_ticks_action.isChecked():
            self.make_axes()

    def toggle_show_true_ions(self):
        self.show_true_ions = not self.show_true_ions
           
    def toggle_legend_clicked(self):
        fontsize_list = ['small', 'medium', 'large']
        i = fontsize_list.index(self.sp.legend_fontsize) + 1
        if i == len(fontsize_list):
            self.sp.legend_fontsize = fontsize_list[0] 
            self.sp.legend_loc = (self.sp.legend_loc)%2+1
        else:
             self.sp.legend_fontsize = fontsize_list[i]
        self.make_axes()

    def enable_tooltips_action_clicked(self):
        if self.enable_tooltips_action.isChecked():   
            self.enableToolTips()         
            self.sp.set_conf('qt_enable_tooltips', True)
            log_.debug('Tooltips enabled', calling=self.calling)
        else:
            self.disableToolTips()
            self.sp.set_conf('qt_enable_tooltips', False)
            log_.debug('Tooltips disabled', calling=self.calling)

    def adjust_fig_action_clicked(self):
        if self.adjust_fig_action.isChecked():   
            self.sp.set_conf('fig_adjust', True)
            log_.debug('Adjust figure enabled', calling=self.calling)
        else:
            self.fig.subplots_adjust(hspace=self.sp.get_conf('fig_hspace'), 
                                 bottom=self.sp.get_conf('fig_bottom'), 
                                 right=self.sp.get_conf('fig_right'), 
                                 top=self.sp.get_conf('fig_top'), 
                                 left=self.sp.get_conf('fig_left'))
 
            log_.debug('Adjust figure disabled', calling=self.calling)
        self.draw_ion()
        
    def show_uncor_obs_action_clicked(self):
        if self.show_uncor_obs_action.isChecked():
            self.sp.show_uncor_spec = True
        else:
            self.sp.show_uncor_spec = False
        self.set_plot_limits_and_draw()

    
    def disableToolTips(self):
        self.lineIDs_GroupBox.setToolTip('')        
        self.residual_GroupBox.setToolTip('')        
        self.run_button.setToolTip('')        
        self.adjust_button.setToolTip('')        
        self.line_info_box.setToolTip('')        
        self.ebv_box.setToolTip('')  
        self.obj_velo_box.setToolTip('') 
        self.sp_min_box.setToolTip('')        
        self.sp_max_box.setToolTip('')        
        self.xlim_min_box.setToolTip('')        
        self.xlim_max_box.setToolTip('')        
        self.y1lim_min_box.setToolTip('')        
        self.y1lim_max_box.setToolTip('')        
        self.y3lim_min_box.setToolTip('')        
        self.y3lim_max_box.setToolTip('')        
        self.fix_axes_cb.setToolTip('')        
        self.cut_cb.setToolTip('')        
        self.ion_cb.setToolTip('')        
        self.sp_norm_box.setToolTip('')        
        self.resol_box.setToolTip('') 
        self.cut2_box.setToolTip('')        
        self.ion_box.setToolTip('')        
        self.line_sort_menu.setToolTip('')        
        self.line_field_menu.setToolTip('')        
        self.line_tick_ax_menu.setToolTip('')        
        self.line_tick_pos_menu.setToolTip('')        
        self.diff_lines_menu.setToolTip('')        
        self.verbosity_menu.setToolTip('')        
        self.style_menu.setToolTip('')        
        
    def enableToolTips(self):
        self.lineIDs_GroupBox.setToolTip(self.lineIDs_GroupBox_ToolTip)        
        self.residual_GroupBox.setToolTip(self.residual_GroupBox_ToolTip)        
        self.run_button.setToolTip(self.run_button_ToolTip)        
        self.adjust_button.setToolTip(self.adjust_button_ToolTip)        
        self.line_info_box.setToolTip(self.line_info_box_ToolTip)        
        self.ebv_box.setToolTip(self.ebv_box_ToolTip)  
        self.obj_velo_box.setToolTip(self.obj_velo_box_ToolTip) 
        self.sp_min_box.setToolTip(self.sp_min_box_ToolTip)        
        self.sp_max_box.setToolTip(self.sp_max_box_ToolTip)        
        self.xlim_min_box.setToolTip(self.xlim_min_box_ToolTip)        
        self.xlim_max_box.setToolTip(self.xlim_max_box_ToolTip)        
        self.y1lim_min_box.setToolTip(self.y1lim_min_box_ToolTip)        
        self.y1lim_max_box.setToolTip(self.y1lim_max_box_ToolTip)        
        self.y3lim_min_box.setToolTip(self.y3lim_min_box_ToolTip)        
        self.y3lim_max_box.setToolTip(self.y3lim_max_box_ToolTip)        
        self.fix_axes_cb.setToolTip(self.fix_axes_cb_ToolTip)        
        self.cut_cb.setToolTip(self.cut_cb_ToolTip)        
        self.ion_cb.setToolTip(self.ion_cb_ToolTip)        
        self.sp_norm_box.setToolTip(self.sp_norm_box_ToolTip)        
        self.resol_box.setToolTip(self.resol_box_ToolTip) 
        self.cut2_box.setToolTip(self.cut2_box_ToolTip)        
        self.ion_box.setToolTip(self.ion_box_ToolTip)        
        self.line_sort_menu.setToolTip(self.line_sort_menu_ToolTip)        
        self.line_field_menu.setToolTip(self.line_field_menu_ToolTip)        
        self.line_tick_ax_menu.setToolTip(self.line_tick_ax_menu_ToolTip)        
        self.line_tick_pos_menu.setToolTip(self.line_tick_pos_menu_ToolTip)        
        self.diff_lines_menu.setToolTip(self.diff_lines_menu_ToolTip)        
        self.verbosity_menu.setToolTip(self.verbosity_menu_ToolTip)        
        self.style_menu.setToolTip(self.style_menu_ToolTip)        
    
    def show_line_ticks_action_clicked(self):
        self.set_ion()
        if self.plot_lines_action.isChecked():
            self.sp.set_conf('plot_lines_of_selected_ions', True)
        else:
            self.sp.set_conf('plot_lines_of_selected_ions', False)
        if self.show_line_ticks_action.isChecked() or self.plot_lines_action.isChecked():
            self.lineIDs_GroupBox.setChecked(True)
        else:  
            self.lineIDs_GroupBox.setChecked(False)
        self.make_axes()

    def plot_cont_action_clicked(self):
        if self.plot_cont_action.isChecked():
            self.sp.set_conf('cont_plot', True)
        else:
            self.sp.set_conf('cont_plot', False)
        self.on_draw()
        
    def ion_cb_changed(self):
        if self.ion_cb.isChecked():
            self.sp.set_conf('show_selected_ions_only', True) 
            self.selected_ions_action.setChecked(True)
        else:  
            self.sp.set_conf('show_selected_ions_only', False) 
            self.selected_ions_action.setChecked(False)
        self.make_axes()
        
    def cut_cb_changed(self):
        if self.cut_cb.isChecked():
            self.sp.set_conf('show_selected_intensities_only', True) 
            self.selected_intensities_action.setChecked(True)
        else:  
            self.sp.set_conf('show_selected_intensities_only', False) 
            self.selected_intensities_action.setChecked(False)
        self.make_axes()
          
    def selected_lines_clicked(self):
        if self.selected_ions_action.isChecked():
            self.sp.set_conf('show_selected_ions_only', True) 
            self.ion_cb.setChecked(True)
        else:  
            self.sp.set_conf('show_selected_ions_only', False) 
            self.ion_cb.setChecked(False)
        if self.selected_intensities_action.isChecked():
            self.sp.set_conf('show_selected_intensities_only', True) 
            self.cut_cb.setChecked(True)
        else:  
            self.sp.set_conf('show_selected_intensities_only', False) 
            self.cut_cb.setChecked(False)
        self.make_axes()

    def diff_lines_by_process_clicked(self):
        if self.diff_lines_by_process_action.isChecked():
            self.sp.set_conf('diff_lines_by_process', True) 
        else:  
            self.sp.set_conf('diff_lines_by_process', False) 
        self.make_axes()

    def editing_lines_clicked(self):
        if self.editing_lines_action.isChecked():
            self.sp.set_conf('qt_allow_editing_lines', True) 
        else:  
            self.sp.set_conf('qt_allow_editing_lines', False) 

    def update_lines_clicked(self):
        if self.update_lines_action.isChecked():
            self.sp.set_conf('qt_update_after_editing_lines', True) 
        else:  
            self.sp.set_conf('qt_update_after_editing_lines', False) 

    def cycle_forwards_ions(self):
        j = self.sp.get_conf('index_of_current_ion')
        selected_ions = self.sp.get_conf('selected_ions')   
        if j in range(-1, len(self.sp.selected_ions_data)-1):
            j += 1
        else:
            j = -1
        self.sp.set_conf('index_of_current_ion', j)
        self.set_refline_to_info_box(j)
        self.make_axes()

    def cycle_backwards_ions(self):
        j = self.sp.get_conf('index_of_current_ion')
        selected_ions = self.sp.get_conf('selected_ions')   
        if j in range(0, len(self.sp.selected_ions_data)):
            j -= 1
        else:
            j = len(self.sp.selected_ions_data)-1
        self.sp.set_conf('index_of_current_ion', j)
        self.set_refline_to_info_box(j)
        self.make_axes()
         
    def show_line_ticks_from_file(self):
        file_choices = "Text files (*.txt *.dat) (*.txt *.dat);;Tex files (*.tex) (*.tex);;CSV files (*.csv) (*.csv);;All Files (*) (*)"
        if self.tick_file is None:
            path = ''
        else:
            path = self.tick_file
        path = unicode(QtGui.QFileDialog.getOpenFileName(self, 'Open file', path, file_choices))
        if path:
            self.tick_file = path
        else:
            return
        f = open(self.tick_file, 'r')
        lines = f.readlines()
        f.close()
        color = 'darkmagenta'
        posTick = self.sp.get_conf('line_tick_pos')
        y1, y2 = self.get_line_tick_lim(posTick)
        k = self.sp.get_conf('line_tick_ax')    
        if k == 2: 
            k = 1
            y1 = 0.2
            y2 = 0.8
        elif k == 1 and self.residual_GroupBox.isChecked():
            k = 1
        else:
            k = 0
        dy = (y2-y1)*0.30
        if self.magenta_tick_shown == True:
            self.draw_ion()
        for line in lines:
            line = line.strip()
            line = line.split(' ')[0]
            if self.isFloat(line):
                wavelength = np.float(line)
                if wavelength > self.x_plot_lims[0] and wavelength < self.x_plot_lims[1]:
                    self.fig.axes[k].axvline( wavelength, y1+dy, y2-dy, color = color, linestyle = 'solid', linewidth = 1.5 ) 
        self.fig.axes[k].step( [0,0], [0,100], color = color, linestyle = 'solid', linewidth = 1.5, label = self.tick_file.split('/')[-1] )
        self.fig.axes[k].legend(loc=self.sp.legend_loc, fontsize=self.sp.legend_fontsize)
        self.fig.canvas.draw()         
        self.magenta_tick_shown = True
        
    def residual_box_clicked(self):
        if self.residual_GroupBox.isChecked():
            self.sp.set_conf('qt_plot_residuals', True)
        else:
            self.sp.set_conf('qt_plot_residuals', False)
        self.make_axes()
        
    def make_axes(self):
        log_.debug('Entering make_axes', calling=self.calling)
        if self.call_on_draw: 
            self.save_axes()
        self.fig.clf()

        i_ax1 = 0
        i_ax2 = 1
        i_ax3 = 2
        rspan_ax1 = 4
        rspan_ax2 = 1
        rspan_ax3 = 4
        n_subplots = rspan_ax1

        k = self.sp.get_conf('line_tick_ax')
        ShowAx2 = self.show_line_ticks_action.isChecked() and ( k == 2 )
        if ShowAx2:
            i_ax2 = n_subplots
            n_subplots += rspan_ax2
        if self.residual_GroupBox.isChecked():
            i_ax3 = n_subplots
            n_subplots += rspan_ax3
        if self.axes is not None:
            del(self.axes)
        self.axes = plt.subplot2grid((n_subplots,1), (i_ax1,0), rowspan=rspan_ax1)
        self.sp.ax1 = self.axes
        if ShowAx2:
            if self.axes2 is not None:
                del(self.axes2)
            self.axes2 = plt.subplot2grid((n_subplots,1), (i_ax2,0), rowspan=rspan_ax2, sharex=self.axes )
            self.axes2.tick_params( left='off',labelleft='off' ) 
            self.sp.ax2 = self.axes2
            self.axes.get_xaxis().set_visible(False)
        else:
            self.axes2 = None
            self.sp.ax2 = None
        if self.residual_GroupBox.isChecked():
            if self.axes3 is not None:
                del(self.axes3)
            self.axes3 = plt.subplot2grid((n_subplots,1), (i_ax3,0), rowspan=rspan_ax3, sharex=self.axes )
            self.sp.ax3 = self.axes3
            if ShowAx2:
                self.axes2.get_xaxis().set_visible(False)
            self.axes.get_xaxis().set_visible(False)
        else:
            self.axes3 = None
            self.sp.ax3 = self.axes3
        self.fig.subplots_adjust(hspace=self.sp.get_conf('fig_hspace'), 
                                 bottom=self.sp.get_conf('fig_bottom'), 
                                 right=self.sp.get_conf('fig_right'), 
                                 top=self.sp.get_conf('fig_top'), 
                                 left=self.sp.get_conf('fig_left'))
        if self.call_on_draw: 
            log_.debug('Calling on_draw from make_axes', calling=self.calling)
            self.do_save = False
            self.on_draw()
            self.do_save = True
        
        log_.debug('Exit make_axes', calling=self.calling)

    def init_axes(self):
        self.x_plot_lims = self.sp.get_conf('x_plot_lims')
        if self.x_plot_lims is None:
            self.x_plot_lims = (np.min(self.sp.w), np.max(self.sp.w))
            
        self.y1_plot_lims = self.sp.get_conf('y1_plot_lims')
        if self.y1_plot_lims is None:
            mask = (self.sp.w_ori > self.x_plot_lims[0]) & (self.sp.w_ori < self.x_plot_lims[1]) 
            r = 1.2
            if self.sp.sp_synth_lr is None:
                a = np.min(self.sp.f[mask])
                b = np.max(self.sp.f[mask])
            else:
                a = np.min(self.sp.sp_synth_lr[mask])
                b = np.max(self.sp.sp_synth_lr[mask])
            self.y1_plot_lims = ((a*(1+r)+b*(1-r))/2, (a*(1-r)+b*(1+r))/2)
        
        self.y2_plot_lims = self.sp.get_conf('y2_plot_lims')
        if self.y2_plot_lims is None:
            self.y2_plot_lims = (-0.5, 1.5)
        
        self.y3_plot_lims = self.sp.get_conf('y3_plot_lims')
        if self.y3_plot_lims is None:
            mask = (self.sp.w_ori > self.x_plot_lims[0]) & (self.sp.w_ori < self.x_plot_lims[1])
            r = 1.2
            if self.sp.sp_synth_lr is None:
                self.y3_plot_lims = (-1,1)
            else:
                a = np.min((self.sp.f_ori - self.sp.sp_synth_lr)[mask])
                b = np.max((self.sp.f_ori - self.sp.sp_synth_lr)[mask])
            self.y3_plot_lims = ((a*(1+r)+b*(1-r))/2, (a*(1-r)+b*(1+r))/2)

        log_.debug('Axes initialized. IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)
        self.print_axes()

    def save_axes(self):
        if self.axes is not None:
            self.x_plot_lims = self.axes.get_xlim()
            self.y1_plot_lims = self.axes.get_ylim()
            self.xscale = self.axes.get_xscale()
            self.yscale = self.axes.get_yscale()
        if self.axes2 is not None:
            self.y2_plot_lims = self.axes2.get_ylim()
        if self.axes3 is not None:
            self.y3_plot_lims = self.axes3.get_ylim()
        self.sp.save_axes()
        log_.debug('Axes saved. IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)
        self.print_axes()
        
    def restore_axes(self):
        if self.x_plot_lims is not None:
            if self.axes is not None:
                self.axes.set_xlim(self.x_plot_lims)
                log_.debug('X-axes restored to {}'.format(self.axes.get_xlim()), calling=self.calling)
            else:
                log_.debug('axes is None', calling=self.calling)
        else:
            log_.debug('x_plot_lims is None', calling=self.calling)
        if self.y1_plot_lims is not None:
            if self.axes is not None:
                self.axes.set_ylim(self.y1_plot_lims)
        if self.y2_plot_lims is not None:
            if self.axes2 is not None:
                self.axes2.set_ylim(self.y2_plot_lims)
        if self.y3_plot_lims is not None:
            if self.axes3 is not None:
                self.axes3.set_ylim(self.y3_plot_lims)
        if self.xscale is not None:
            self.axes.set_xscale(self.xscale)
            log_.debug('X scale set to {}'.format(self.xscale))
        if self.yscale is not None:
            self.axes.set_yscale(self.yscale)
            log_.debug('Y scale set to {}'.format(self.yscale))
        
        log_.debug('Axes restored. IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)
        self.print_axes()
        
    def print_axes(self):
        log_.debug('lims: {} {} {} {}'.format(self.x_plot_lims, self.y1_plot_lims, self.y2_plot_lims, self.y3_plot_lims), calling=self.calling)
        log_.debug('Axes IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)
        log_.debug(' IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)

    def exec_init(self):
        if self.init_file_name is None:
            self.get_init_filename()
        if self.init_file_name:
            self.statusBar().showMessage('Running synthesis ...') 
            QtGui.QApplication.processEvents() 
            self.start_spectrum()
            self.do_save = False
            self.on_draw()
            self.do_save = True
            self.restore_axes()
            self.update_lim_boxes()
            self.save_parameters_file = None
        else:
            log_.warn('A filename must be given', calling=self.calling)
            sys.exit('An initialization filename must be given')

    def get_init_filename(self):
        file_choices = "Python initialization files (*init.py) (*init.py);;Python files (*.py) (*.py);;All files (*) (*)"
        title = 'Open pySSN initialization file'
        init_file = str(QtGui.QFileDialog.getOpenFileName(self, title, self.init_file_name, file_choices))
        if init_file and os.path.isfile(init_file):
            self.init_file_name = init_file
        else:
            self.init_file_name = ''

    def select_init(self):
        old_name = self.init_file_name
        self.get_init_filename()
        if self.init_file_name:
            self.exec_init()
        else:
            self.init_file_name = old_name

    def save_pars(self):
        path = self.sp.get_conf('save_parameters_filename')
        keys = self.sp.default_keys
        if '__builtins__' in keys:
            keys.remove('__builtins__')
        keys.sort()
        with open(path, 'w') as f:
            for key in keys:
                value = self.sp.conf[key]
                if isinstance(value, basestring):
                    value = '\"{}\"'.format(value)
                f.write('{} = {}\n'.format(key, value))
        self.statusBar().showMessage('Parameters saved to file %s' % path, 4000)

    def save_pars_as(self):
        if self.save_parameters_file is None:
            path = self.init_file_name
        else:
            path = self.save_parameters_file        
        keys = self.sp.default_keys
        keys_to_be_removed = ['__builtins__', 'plot_magenta', 'label_magenta', 'plot_cyan', 'label_cyan']
        for key in keys_to_be_removed: 
            if key in keys:
                keys.remove(key)
        keys.sort()
        file_choices = "pySSN initialization files (*init.py) (*init.py);;Python files (*.py) (*.py);;All files (*) (*)"
        title = 'Save synthesis and plot parameters'
        selectedFilter = 'pySSN initialization files (*init.py) (*init.py)'
        path = unicode(QtGui.QFileDialog.getSaveFileName(self, title, path, file_choices, selectedFilter))
        if path:
            with open(path, 'w') as f:
                for key in keys:
                    if key == 'instr_prof':
                        value = self.sp.format_instr_prof()
                    else:
                        value = self.sp.conf[key]
                        if isinstance(value, basestring):
                            value = '\"{}\"'.format(value)
                    f.write('{} = {}\n'.format(key, value))
            self.save_parameters_file = path
            self.statusBar().showMessage('Parameters saved to file %s' % path, 4000)
                       
    def teste_instr_prof(self, prof):
        if prof is None:
            return 'not defined'
        keys = prof.keys()
        keys.remove('comment')
        if not 'width' in keys:
            return 'The parameter \'width\' was not found.'       
        if prof['width'] == 0.0:
            return 'The value of \'width\' can not be zero'
        if not (self.sp.get_key_indexes('Bb', prof)==self.sp.get_key_indexes('Br', prof)==
           self.sp.get_key_indexes('beta', prof)==self.sp.get_key_indexes('alpha', prof)):
           return 'Invalid indexes por the parameters \'Bb\', \'Br\', \'alpha\', or \'beta\''  
        if not all((type(prof[key])==float or type(prof[key])==int) for key in keys):
            return 'The values of parameters must be numbers.'
        return ''
                    
    def apply_instr_prof(self):

        def do_update():
            path = str(prof_box.toPlainText()).strip()
            try:
                user_module = {}
                exec(path) in user_module
                prof = user_module['instr_prof']
                self.sp.set_conf('instr_prof', prof)
                log_.message('new instrumental profile is ok', calling = self.calling)
            except:
                title = 'Error reading instrument profile'
                msg = 'Unable to read instrumental profile'
                path = None
                if self.showErrorBox:
                    QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok)
                else:
                    log_.warn(msg, calling = self.calling)
                return
            msg = self.teste_instr_prof(prof)
            if not msg:
                self.update_profile()
            else:
                title = 'Error in the instrument profile'
                if self.showErrorBox:
                    QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok)
                else:
                    log_.warn(msg, calling = self.calling)            
            
        def toggle_statusbar():
            self.showHelpBrowser = not self.showHelpBrowser
            helpBrowser.setVisible(self.showHelpBrowser)
            if self.showHelpBrowser:
                self.instr_prof_dialog.resize(self.instr_prof_dialog_width, 2.1*self.instr_prof_dialog_height)
            else:
                self.instr_prof_dialog.resize(self.instr_prof_dialog_width, self.instr_prof_dialog_height)
                
        def get_window_size_and_position():
            if self.instr_prof_dialog is None:
                font = QtGui.QFont("Courier")
                width = QtGui.QFontMetrics(font).width('='*80)
                height = 15*QtGui.QFontMetrics(font).height()
                self.instr_prof_dialog_width = width
                self.instr_prof_dialog_height = height
                sG = QtGui.QApplication.desktop().screenGeometry()
                self.instr_prof_dialog_x = sG.width()-self.instr_prof_dialog_width
                self.instr_prof_dialog_y = sG.height()
            else:
                if not self.showHelpBrowser:
                    self.instr_prof_dialog_width = self.instr_prof_dialog.width()
                    self.instr_prof_dialog_height = self.instr_prof_dialog.height()
                    self.instr_prof_dialog_x = self.instr_prof_dialog.pos().x()
                    self.instr_prof_dialog_y = self.instr_prof_dialog.pos().y()
            
        get_window_size_and_position()
        self.instr_prof_dialog = QtGui.QDialog()
        self.instr_prof_dialog.resize(self.instr_prof_dialog_width, self.instr_prof_dialog_height)
        self.instr_prof_dialog.move(self.instr_prof_dialog_x,self.instr_prof_dialog_y)
        self.instr_prof_dialog.setWindowTitle('instrument profile dialog')
        prof_box = QtGui.QTextEdit()
        prof_box.setFontFamily("Courier")
        prof_box.setText('instr_prof = ' + self.sp.format_instr_prof())
        linkLabel = QtGui.QLabel('<a href="https://github.com/Morisset/pySSN/wiki">More help online</a>')
        linkLabel.setOpenExternalLinks(True)
        helpBrowser = QtGui.QTextBrowser()
        
        # text=open('instr_prof.html').read()
        # This text should go to a file open with text=open('instr_prof.html').read()
        text = """<title> Instrumental profile help</title>
<p>The instrumental profile if defined by the <a href="https://en.wikibooks.org/wiki/Python_Programming/Dictionaries">python dictionary</a> <b>instr_prof</b>.

<p>The main component of the instrumental profile is set by the parameter <b>width</b>, which is the only indispensable parameters.</p>

<p>If <b>width</b> > 0, the main component profile follows a <a href="https://en.wikipedia.org/wiki/Normal_distribution">Gaussian distribution</a>, P &prop; exp(-(&lambda;/<b>width</b>)<sup>2</sup>).
In this case, <b>width</b> is related to the normal full-width at half maximum by <b>width</b> = FWHM/(2(ln2)<sup>1/2</sup>) = FWHM/1.665.</p>

<p>If <b>width</b> &lt; 0, the main component profile follows a <a href="https://en.wikipedia.org/wiki/rectangular_distribution">rectangular distribution</a>, P = 1 for -|<b>width</b>|/2  &lt; &lambda; &lt; |<b>width</b>|/2, and P = 0 for all other values of &lambda;.</p>

<p>A variable number of optional components can be included, each defined by four parameters, <b>Bb</b>, <b>Br</b>, <b>alpha</b>, and <b>beta</b>, and following P &prop;  <b>B</b>exp(-(&lambda;/<b>beta</b>)<sup><b>alpha</b></sup>).
<b>Bb</b> and <b>Br</b> are the intensity scale parameters for the bluish and reddish sides of the profile, respectively.</p>

<p>If more than one optional component is in use, the parameters must be indexed as <b>alpha_1</b> <b>alpha_2</b>, etc.</p>

Special cases for the optional components:
<ul>
  <li><b>alpha</b> = 2 produces a <a href="https://en.wikipedia.org/wiki/Normal_distribution">Gaussian distribution</a>.
  <li><b>alpha</b> = 2, <b>Bb</b> = 0 (or <b>Br</b> = 0) produces a <a href="https://en.wikipedia.org/wiki/Half_normal_distribution">half-Gaussian distribution</a>.
  <li><b>alpha</b> = 1 produces an <a href="https://en.wikipedia.org/wiki/Exponential_distribution">exponential distribution</a>.
</ul>

<p>A comment may be included in <b>instr_prof</b>.</p>
<p>Examples:</p>
<ol>
<li>instr_prof = {'width': 0.5}<br>
<li>instr_prof = {'width': 0.5, 'comment': 'Gaussian profle'}<br>
<li>Example: instr_prof = {'width': 0.5, 'Bb':0.00016, 'Br':9e-05, 'beta': 2.2, 'alpha': 0.45}<br>
<li>instr_prof = {'width': 0.5, 'Bb_1':0.00016, 'Br_1':9e-05, 'beta_1': 2.2, 'alpha_1': 0.45, 'Bb_2': 0.0014, 'Br_2':0.001, 'beta_2':  1.4, 'alpha_2':  0.75}<br>
</ol>"""
        helpBrowser.document().setHtml(text)
        helpBrowser.setOpenExternalLinks(True)
        self.showHelpBrowser = False
        helpBrowser.setVisible(self.showHelpBrowser)
        policy = helpBrowser.sizePolicy()
        policy.setVerticalStretch(20)
        helpBrowser.setSizePolicy(policy)
        vbox = QtGui.QVBoxLayout()
        buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Help|
                                           QtGui.QDialogButtonBox.Close|
                                           QtGui.QDialogButtonBox.Apply)
        buttonBox.button(QtGui.QDialogButtonBox.Apply).setText("Update")
        vbox.addWidget(prof_box,0)
        vbox.addWidget(buttonBox)
        vbox.addWidget(linkLabel)
        vbox.addWidget(helpBrowser)
        buttonBox.button(QtGui.QDialogButtonBox.Help).clicked.connect(toggle_statusbar)
        buttonBox.button(QtGui.QDialogButtonBox.Apply).clicked.connect(do_update)
        buttonBox.rejected.connect(self.instr_prof_dialog.close)
        self.instr_prof_dialog.setLayout(vbox)
        self.instr_prof_dialog.setWindowModality(QtCore.Qt.NonModal)
        self.instr_prof_dialog.show()
                   
    def refine_wavelengths(self):

        def table2list(text):
            text = str(text)
            text = text.splitlines()
            s = ''
            for i in range(len(text)):
                line = text[i].split()
                if len(line) == 2 and sum([self.isFloat(x) for x in line]) == 2: 
                    s += '({}, {}), '.format(line[0], line[1])
                else:
                    if len(line) > 0:
                        title = 'Error in table'
                        msg = 'Error in line \'{}\'.\nEach line must have two numbers separated by whitespaces.'.format(text[i])
                        if self.showErrorBox:
                            QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok)
                        else:
                            log_.warn(msg, calling = self.calling)
                        return ''

            s = s.strip(' ,')
            if s == '':
                return 'lambda_shift_table = None'
            else:
                return 'lambda_shift_table = [{}]'.format(s)
        
        def toggle_table():
            self.refine_wave_as_table = not self.refine_wave_as_table
            if self.refine_wave_as_table:
                text = str(edit_box.toPlainText()).strip()
                edit_box.clear()
                text = text.replace('lambda_shift_table','')
                text = text.strip(' =[]')
                text = text.split(')')
                for i in range(len(text)-1):
                    line = text[i].strip(' (,')
                    line = line.split(',')
                    line = '{:<7} {}'.format(line[0].strip(),line[1].strip())
                    edit_box.append(line)
                buttonBox.button(QtGui.QDialogButtonBox.RestoreDefaults).setText("Show as list")
            else:
                text = table2list(edit_box.toPlainText())
                if text == '':
                    self.refine_wave_as_table = True
                    return
                edit_box.clear()
                edit_box.setText(text)
                buttonBox.button(QtGui.QDialogButtonBox.RestoreDefaults).setText("Show as table")

        def do_update():
            old_value = self.sp.get_conf('lambda_shift_table')
            if self.refine_wave_as_table:
                path = table2list(edit_box.toPlainText())
                if path == 'error':
                    return
            else:
                path = str(edit_box.toPlainText()).strip()
            try:
                user_module = {}
                exec(path) in user_module
                value = user_module['lambda_shift_table']
                self.sp.set_conf('lambda_shift_table', value)
                log_.message('new \'lambda_shit_table\' is ok', calling = self.calling)
            except:
                title = 'Error'
                msg = 'Unable to read \'lambda_shit_table\''
                path = None
                if self.showErrorBox:
                    QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok)
                else:
                    log_.warn(msg, calling = self.calling)
                return
            self.sp.show_uncor_spec = True            
            self.sp.init_obs()
            if self.sp.read_obs_error:
                self.sp.set_conf('lambda_shift_table', old_value)
                if self.showErrorBox:
                    title = 'Error'
                    msg = self.sp.read_obs_error
                    QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok)
                else:
                    log_.warn(msg, calling = self.calling)
            else:
                """
                self.sp.init_red_corr()
                self.sp.make_continuum()
                self.sp.run()
                self.set_plot_limits_and_draw()
                """
                self.rerun()
            if not self.show_uncor_obs_action.isChecked():
                self.sp.show_uncor_spec = False
            
        def toggle_help():
            self.showHelpBrowser = not self.showHelpBrowser
            helpBrowser.setVisible(self.showHelpBrowser)
            if self.showHelpBrowser:
                self.refine_wave_dialog.resize(self.refine_wave_dialog_width, 2.5*self.refine_wave_dialog_height)
            else:
                self.refine_wave_dialog.resize(self.refine_wave_dialog_width, self.refine_wave_dialog_height)
                
        def get_window_size_and_position():
            if self.refine_wave_dialog is None:
                font = QtGui.QFont("Courier")
                width = QtGui.QFontMetrics(font).width('='*80)
                height = 15*QtGui.QFontMetrics(font).height()
                self.refine_wave_dialog_width = width
                self.refine_wave_dialog_height = height
                sG = QtGui.QApplication.desktop().screenGeometry()
                self.refine_wave_dialog_x = sG.width()-self.refine_wave_dialog_width
                self.refine_wave_dialog_y = sG.height()
            else:
                if not self.showHelpBrowser:
                    self.refine_wave_dialog_width = self.refine_wave_dialog.width()
                    self.refine_wave_dialog_height = self.refine_wave_dialog.height()
                    self.refine_wave_dialog_x = self.refine_wave_dialog.pos().x()
                    self.refine_wave_dialog_y = self.refine_wave_dialog.pos().y()
            
        get_window_size_and_position()
        self.refine_wave_dialog = QtGui.QDialog()
        self.refine_wave_dialog.resize(self.refine_wave_dialog_width, self.refine_wave_dialog_height)
        self.refine_wave_dialog.move(self.refine_wave_dialog_x,self.refine_wave_dialog_y)
        self.refine_wave_dialog.setWindowTitle('wavelength-refining dialog')
        edit_box = QtGui.QTextEdit()
        edit_box.setFontFamily("Courier")
        self.refine_wave_as_table = False
        edit_box.setText('lambda_shift_table = ' + str(self.sp.get_conf('lambda_shift_table')))
        linkLabel = QtGui.QLabel('<a href="https://github.com/Morisset/pySSN/wiki">More help online</a>')
        linkLabel.setOpenExternalLinks(True)
        helpBrowser = QtGui.QTextBrowser()
        
        # text=open('wave_refining.html').read()
        # This text should go to a file open with text=open('wave-refining').read()
        text = """<title> Wavelength-refining help</title>
<p>The wavelength calibration of the observational spectrum can be refined with the use of 
the <a href="https://en.wikibooks.org/wiki/Python_Programming/Lists">python list</a> <b>lambda_shift_table</b>. 
Each element of this list is an ordered pair of numbers (&lambda;, &Delta;&lambda;), where &Delta;&lambda; is the wavelength shift at the wavelength &lambda; needed to improve the calibration, after the Doppler correction.</p>

<p>The data in <b>lambda_shit_table</b> will be linearly interpolated to provide the corrected wavelengths. 
Outside the range of wavelenghts given in <b>lambda_shit_table</b>, the correction will be extrapolated to zero.</p>

<p>To set aside the wavelength-refining, set <b>lambda_shit_table</b> to None.</p>

<p>Examples:</p>
<ol>
<li><p>lambda_shift_table = [(4674, 0.05), (4690, 0.1), (9000, 1)]</p></li>
<li><p>lambda_shift_table = None (to set aside the wavelength-refining)</p></li>
</ol>

<p>Button functions:</p>
<ul>   
<li><p>Click on <b><span style="color:red">Show as table</span></b> to display and edit the data contained in <b>lambda_shit_table</b> as a two columns table.</p></li>

<li><p>Click on <b><span style="color:red">Show as list</span></b> to get back the <b>lambda_shit_table</b> list from the two columns table.</p></li>

<li><p>Click on <b><span style="color:red">Update</span></b> to refine the wavelength calibration and redo the synthesis.</p></li>
</ul>
"""
        helpBrowser.document().setHtml(text)
        helpBrowser.setOpenExternalLinks(True)
        self.showHelpBrowser = False
        helpBrowser.setVisible(self.showHelpBrowser)
        policy = helpBrowser.sizePolicy()
        policy.setVerticalStretch(20)
        helpBrowser.setSizePolicy(policy)
        vbox = QtGui.QVBoxLayout()
        buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Help|
                                           QtGui.QDialogButtonBox.RestoreDefaults|
                                           QtGui.QDialogButtonBox.Close|
                                           QtGui.QDialogButtonBox.Apply)
        buttonBox.button(QtGui.QDialogButtonBox.Apply).setText("Update")
        buttonBox.button(QtGui.QDialogButtonBox.RestoreDefaults).setText("Show as table")
        vbox.addWidget(edit_box,0)
        vbox.addWidget(buttonBox)
        vbox.addWidget(linkLabel)
        vbox.addWidget(helpBrowser)
        buttonBox.button(QtGui.QDialogButtonBox.Help).clicked.connect(toggle_help)
        buttonBox.button(QtGui.QDialogButtonBox.RestoreDefaults).clicked.connect(toggle_table)
        buttonBox.button(QtGui.QDialogButtonBox.Apply).clicked.connect(do_update)
        buttonBox.rejected.connect(self.refine_wave_dialog.close)
        self.refine_wave_dialog.setLayout(vbox)
        self.refine_wave_dialog.setWindowModality(QtCore.Qt.NonModal)
        self.refine_wave_dialog.show()

    def isValidFilename(self, filename):
        if filename is None:
            return False
        try:
            open(filename,'r')
            return True
        except IOError:
            try:
                open(filename, 'w')
                return True
            except IOError:
                return False        

    def set_cosmetic_file(self):
        file_choices = "Line cosmetic files (*cosm*.dat) (*cosm*.dat);;Data files (*.dat) (*.dat);;All files (*) (*)"
        title = 'Set the line cosmetic file'
        cosmetic_file = str(QtGui.QFileDialog.getSaveFileName(self, title, '', file_choices, options=QtGui.QFileDialog.DontConfirmOverwrite))
        msg = "Line cosmetic file '{}' not valid!".format(cosmetic_file)
        if cosmetic_file and not self.isValidFilename(cosmetic_file):
            QtGui.QMessageBox.critical(self, 'pySSN', msg, QtGui.QMessageBox.Ok )  
            cosmetic_file = None
        if cosmetic_file:
            self.sp.set_conf('do_cosmetik', True)
            dir_ = os.path.dirname(cosmetic_file)
            if dir_ == os.getcwd():
                cosmetic_file = cosmetic_file.split('/')[-1]
            self.sp.set_conf('fic_cosmetik', cosmetic_file)
            self.sp.fic_cosmetik = cosmetic_file
            
            if self.sp is not None:
                self.set_status_text()
            if self.axes is not None:
                self.adjust()

    def empty_cosmetic_file(self):
        if self.sp.fic_cosmetik is None or self.sp.phyat_file is None:
            return
        title = 'pySSN: cosmetic file'
        msg = 'All lines in the cosmetic file will be removed.\nConfirm?'
        ret = QtGui.QMessageBox.question(self, title, msg, QtGui.QMessageBox.Ok, QtGui.QMessageBox.Cancel )
        if ret == QtGui.QMessageBox.Ok:
            f = open(self.sp.fic_cosmetik, 'w')
            f.close()
    
    def order_lines(self, lines):
        if lines is None:
            return None
        numbers = []
        for line in lines:
            line_num = int(self.sp.fieldStrFromLine(line,'num'))
            numbers.append(line_num)
        lines = [x for _,x in sorted(zip(numbers, lines))]
        return lines
    
    def remove_duplicate_lines(self, lines):
        if lines is None:
            return None
        numbers = []
        output = []
        for line in lines:
            line_num = int(self.sp.fieldStrFromLine(line,'num'))
            if line_num not in numbers:
                numbers.append(line_num)
                output.append(line)
        return output
    
    def order_cosmetic_file(self):
        if self.sp.fic_cosmetik is None or not os.path.isfile(self.sp.fic_cosmetik):
            return
        f = open(self.sp.fic_cosmetik, 'r')
        cosmetic_lines = f.readlines()
        f.close()
        cosmetic_lines = self.order_lines(cosmetic_lines)
        n0 = len(cosmetic_lines)
        cosmetic_lines = self.remove_duplicate_lines(cosmetic_lines)
        n1 = len(cosmetic_lines)
        f = open(self.sp.fic_cosmetik, 'w')
        f.writelines(cosmetic_lines)
        f.close()
        if n0 > n1:
            s = ' and the duplicate lines removed'
        else:
            s = ''
        msg = 'The cosmetic \'{0:}\' file was ordered{1:}.'.format(self.sp.fic_cosmetik, s)
        self.statusBar().showMessage(msg, 4000) 
    
    def clean_cosmetic_file(self):

        def ShowCleanMessage(UnchangedLineList):
            nUL = len(UnchangedLineList)
            if nUL == 1:
                s1 = ''
                s2 = 'was'
                s3 = 'this line'
            elif nUL > 1:
                s1 = 's'
                s2 = 'were'
                s3 = 'these lines'
            msgBox = QtGui.QMessageBox()
            msgBox.setIcon(QtGui.QMessageBox.Question)
            msgBox.title = 'pySSN: cosmetic file'
            msg = '{0:} unchanged line{1:} in the cosmetic file {2:} found.'.format(nUL, s1, s2)
            msgBox.setText(msg)
            msgBox.setInformativeText('Do you want to delete {:}?\n'.format(s3))
            detailedText = 'Unchanged line{:}:\n\n'.format(s1)
            for i in UnchangedLineList:
                detailedText = detailedText + str(i) + '\n'
            msgBox.setDetailedText(detailedText)
            DelButton = msgBox.addButton(self.tr("Delete"), QtGui.QMessageBox.ActionRole)
            s = 'Delete from the cosmetic file all unchanged lines'
            if self.enable_tooltips_action.isChecked():
                DelButton.setToolTip(s)
            msgBox.addButton(QtGui.QMessageBox.Cancel)
            answer = msgBox.exec_()
            if msgBox.clickedButton() == DelButton:
                answer = True
            else:
                answer = False
            return answer

        if self.sp.fic_cosmetik is None or self.sp.phyat_file is None:
            return
        #if not self.sp.get_conf('clean_cosmetic_file'):
        #    return
        if not os.path.isfile(self.sp.fic_cosmetik):
            return
        f = open(self.sp.fic_cosmetik, 'r')
        cosmetic_lines = f.readlines()
        f.close()
        UnchangedLineList = []
        ChangedLines = []
        for i in range(len(cosmetic_lines)):
            line_c = cosmetic_lines[i].rstrip()
            line_num = int(self.sp.fieldStrFromLine(line_c,'num'))
            if self.sp.cosmetic_line_unchanged(line_c):
                UnchangedLineList.append(line_num)
            else:
                ChangedLines.append(line_c + '\n')
                
        if len(UnchangedLineList) > 0:
            ret = ShowCleanMessage(UnchangedLineList)
            if ret == True:
                f = open(self.sp.fic_cosmetik, 'w')
                f.writelines(ChangedLines)
                f.close()
        else:
            msg = 'No unchanged line in the cosmetic file {:}'.format(self.sp.fic_cosmetik)
            self.statusBar().showMessage(msg, 4000) 
                
            
    def match_cosmetic_phyat_files(self):
    
        def ShowErrorMessage():
            msg = 'The wavelength or intensity in the cosmetic file does not match that in the atomic database.\n\n' \
                  'Do you want to try to automatically correct the cosmetic file?'
            msgBox = QtGui.QMessageBox()
            msgBox.setText("Error in cosmetic file for line: " + str(line_num))
            msgBox.setInformativeText(msg)
            msgBox.addButton(QtGui.QMessageBox.Yes)
            msgBox.addButton(QtGui.QMessageBox.YesToAll)
            msgBox.addButton(QtGui.QMessageBox.No)
            msgBox.addButton(QtGui.QMessageBox.NoToAll)
            msgBox.setDefaultButton(QtGui.QMessageBox.Yes)
            answer = msgBox.exec_()
            return answer
    
        def ShowFinalMessage(nErr, nCor, nUnCor, nNfd, UnCorList, NotFound):
            msgBox = QtGui.QMessageBox()
            msgBox.setText('pySSN: error in cosmetic file')
            if nCor > 0:
                s0 = 'Rerun the synthesis to take into account the changes.\n\n'
            else:
                s0 = ''
            if nUnCor > 0:
                s1 = 'The cosmetic data for lines that still have problems will be ignored. ' \
                     'Do you want to delete them from the cosmetic file?'   
            else:
                s1 = ''
            msg = 'Number of lines with problems: {0:}\n' \
                  'Number of corrected lines: {1:}\n' \
                  'Number of uncorrected lines: {2:}\n' \
                  'Number of lines not found in the atomic database: {3:}\n\n' \
                  '{4:}{5:}'.format(nErr, nCor, nUnCor, nNfd, s0, s1)
            msgBox.setInformativeText(msg)
            if nNfd > 0:
                detailedText = 'Lines not found:\n\n'
                for i in NotFound:
                    detailedText = detailedText + i + '\n'
                detailedText =  detailedText + '\n'
            else:
                 detailedText = ''     
            if nUnCor > 0:
                detailedText = detailedText + 'Lines not corrected:\n\n'
                for i in UnCorList:
                    detailedText = detailedText + i + '\n'
            msgBox.setDetailedText(detailedText)
            DelAllButton = msgBox.addButton(self.tr("Delete all"), QtGui.QMessageBox.ActionRole)
            DelNotFndButton = msgBox.addButton(self.tr("delete not found"), QtGui.QMessageBox.ActionRole)
            DelUncorButton = msgBox.addButton(self.tr("delete uncorrected"), QtGui.QMessageBox.ActionRole)
            if self.enable_tooltips_action.isChecked():
                s = 'Delete from the cosmetic file all lines that still have problems'
                DelAllButton.setToolTip(s)
                s = 'Delete from the cosmetic file the lines not found in the atomic database'
                DelNotFndButton.setToolTip(s)
                s = 'Delete from the cosmetic file the uncorrected lines'
                DelUncorButton.setToolTip(s)
            msgBox.addButton(QtGui.QMessageBox.Cancel)
            msgBox.setMaximumHeight(16777215)
            msgBox.setMinimumHeight(800)
            # It does not expand! Why?
            msgBox.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
            msgBox.setSizeGripEnabled(True)
            if nUnCor == 0:
                DelUncorButton.setEnabled(False)
                DelAllButton.setEnabled(False)
            if nNfd == 0:
                DelNotFndButton.setEnabled(False)
                DelAllButton.setEnabled(False)
            answer = msgBox.exec_()
            if msgBox.clickedButton() == DelAllButton:
                answer = ['DelNotFnd', 'DelUncor']
            elif msgBox.clickedButton() == DelNotFndButton:
                answer = ['DelNotFnd']
            elif msgBox.clickedButton() == DelUncorButton:
                answer = ['DelUncor']
            else:
                answer = []    
            return answer
        
        if self.sp.fic_cosmetik is None:
            return
        if os.path.isfile(self.sp.fic_cosmetik):
            cosmetik_arr, errorMsg = self.sp.read_cosmetik()
            if len(errorMsg) > 0:
                self.sp.do_cosmetik = False
                self.sp.set_conf('do_cosmetik', False)        
                title = 'Error in cosmetic file: '
                msg = 'Unable to read cosmetic data from file \'{}\':{}\n\nLine cosmetics will be disabled!'.format(self.sp.get_conf('fic_cosmetik'), errorMsg)
                if self.showErrorBox:
                    QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
                else:
                    log_.warn('{}: {}'.format(title, msg), calling=self.calling)
                return
            ret = None
            f = open(self.sp.fic_cosmetik, 'r')
            cosmetic_lines = f.readlines()
            f.close()
            ErrorList = []
            CorrectedList = []
            UnCorList = []
            NotFound =[]
            k = self.sp.field_pos['id']
            keys = [ 'lambda', 'l_shift', 'i_rel', 'i_cor' ]
            for i in range(len(cosmetic_lines)):
                line_c = cosmetic_lines[i].rstrip()
                line_num = int(self.sp.fieldStrFromLine(line_c,'num'))
                cosmeticLineOk = self.sp.cosmetic_line_ok(line_c)
                if cosmeticLineOk == None:
                    NotFound.append(line_c[:k])
                    ErrorList.append(line_c[:k])
                elif cosmeticLineOk == False:
                    ErrorList.append(line_c[:k])
                    if ret != QtGui.QMessageBox.YesToAll and ret != QtGui.QMessageBox.NoToAll:
                            ret = ShowErrorMessage()
                    if ret == QtGui.QMessageBox.Yes or ret == QtGui.QMessageBox.YesToAll:
                        CorrectedList.append(line_c[:k])
                        line = self.sp.read_line(self.sp.phyat_file, line_num)
                        line = line.rstrip()
                        v0 = {i: np.float(self.sp.fieldStrFromLine(line, i)) for i in keys}
                        v1 = {i: np.float(self.sp.fieldStrFromLine(line_c, i)) for i in keys}
                        l_shift = v1['lambda'] + v1['l_shift'] - v0['lambda']
                        i_cor =  v1['i_cor'] * v1['i_rel'] / v0['i_rel']
                        l_shift_str = self.rightFormat(str(l_shift), 'l_shift')
                        i_cor_str = self.rightFormat(str(i_cor), 'i_cor')
                        line = self.sp.replace_field(line, 'l_shift', l_shift_str)
                        line = self.sp.replace_field(line, 'i_cor', i_cor_str)
                        log_.warn('(corrected) ' + line + '\n', calling=self.calling)
                        self.sp.replace_line(self.sp.fic_cosmetik, line)
                    else:
                        UnCorList.append(line_c[:k])
                        log_.warn('Not corrected.\n', calling=self.calling)
            nErr = len(ErrorList)
            nCor = len(CorrectedList)
            nUnCor = len(UnCorList)
            nNfd = len(NotFound)
            if nErr > 0:
                answer = ShowFinalMessage(nErr, nCor, nUnCor, nNfd, UnCorList, NotFound)

                if  'DelNotFnd' in answer:
                    for i in NotFound:
                        self.sp.remove_line(self.sp.fic_cosmetik, int(i))
                if  'DelUncor' in answer:
                    for i in UnCorList:
                        self.sp.remove_line(self.sp.fic_cosmetik, int(i))
    def set_status_text(self):
    
        if self.sp is None:
            return
            
        if self.sp.phyat_file == 'NO_phyat.dat':
            self.status_text.setText('pySSN, v {}. init file: {}, No synthesis'.format(__version__, 
                     self.sp.config_file.split('/')[-1]))
        elif self.sp.get_conf('do_cosmetik'):
            self.status_text.setText('pySSN, v {}. init file: {}, at. data: {}, model: {}, cosmetic: {}'.format(__version__, 
                     self.sp.config_file.split('/')[-1], 
                     self.sp.phyat_file.split('/')[-1],
                     self.sp.get_conf('fic_modele').split('/')[-1],
                     self.sp.get_conf('fic_cosmetik').split('/')[-1]))
        else:
            self.status_text.setText('pySSN, v {}. init file: {}, at. data: {}, model: {}, No cosmetic'.format(__version__, 
                     self.sp.config_file.split('/')[-1], 
                     self.sp.phyat_file.split('/')[-1],
                     self.sp.get_conf('fic_modele').split('/')[-1]))

    def test_init_file(self):   

        if self.sp == None:
            self.showErrorBox = False
        self.showErrorBox = True

        invalidCommands = []
        if os.path.isfile(self.init_file_name):
            f = open(self.init_file_name, 'r')
            lines = f.readlines()
            f.close()
        else:
            invalidCommands.append('\nFile not found')
            lines = []

        triple_quoted_string_found = False
        newlines = []
        rows = []
        for i in range(len(lines)):
            line = lines[i].split('#')[0].rstrip()
            k = line.find('=')
            if not (line.strip().startswith('#') or len(line.strip()) == 0):
                if '"""' in line:
                    triple_quoted_string_found = not triple_quoted_string_found
                    if triple_quoted_string_found:
                        newlines.append(line.split('#')[0].rstrip())
                        rows.append(i+1)
                    else:
                        s = line.split('#')[0].rstrip()
                        if len(s.strip()) > 0:
                            newlines[-1] += '\n' + s
                else:
                    if len(line) == len(line.lstrip()) and not triple_quoted_string_found:
                        newlines.append(line.split('#')[0].rstrip())
                        rows.append(i+1)
                    else:
                        s = line.split('#')[0].rstrip()
                        if len(s.strip()) > 0:
                            newlines[-1] += '\n' + s

        for i in range(len(newlines)):
            line = newlines[i]
            line_list = line.split('\n')
            if len(line_list) > 3:
                line_str = line_list[0] + '\n' + line_list[1] + '\n' + line_list[2] + '\n...'
            else:
                line_str = line
            
            try:
                exec(line)
            except IndentationError:
                invalidCommands.append('\nIndentation error, line {}:\n{}'.format(rows[i],line_str))
            except SyntaxError:
                if '"""' in line and triple_quoted_string_found:
                    invalidCommands.append('\nUnclosed triple-quotation mark, line {}:\n{}'.format(rows[i],line_str))
                else:
                    invalidCommands.append('\nInvalid syntax, line {}:\n{}'.format(rows[i],line_str))
            except(AttributeError, NameError): 
                invalidCommands.append('\nUndefined variable name or attribute, line {}:\n{}'.format(rows[i],line_str))
            except: 
                invalidCommands.append('\nUndefined error, line {}:\n{}'.format(rows[i],line_str))

        if len(invalidCommands) > 0:
            title = 'Fatal error'
            msg = 'Error in the initialization file \'{0}\': '.format(self.init_file_name)
            for line in invalidCommands:
                msg = msg + '\n' + line
            if self.showErrorBox:
                if self.sp == None:
                    buttom = QtGui.QMessageBox.Abort
                else:
                    buttom = QtGui.QMessageBox.Cancel
                QtGui.QMessageBox.critical(self, title, msg, buttom)
            else:
                log_.warn('{}: {}'.format(title, msg), calling=self.calling)
            return False

        return True

    def start_spectrum(self):
        init_file = self.init_file_name.split('/')[-1]
        dir_ = self.init_file_name.split(init_file)[0]
        if dir_ == '':
            dir_ = './'
        self.directory = dir_
        if not self.test_init_file():
            if self.sp == None:
                sys.exit() 
            else:
                return
        self.sp = spectrum(config_file=self.init_file_name)
        if self.sp.errorMsg:
            if self.showErrorBox:
                msg = 'Synthesis not possible. \n\n{}'.format(self.sp.errorMsg)
                msg = self.sp.errorMsg
                ret = QtGui.QMessageBox.critical(self, 'Critical Error', msg, QtGui.QMessageBox.Abort, QtGui.QMessageBox.Ignore)
                if ret == QtGui.QMessageBox.Abort:
                    sys.exit()
            self.sp.errorMsg = ''
        if len(self.sp.read_obs_error) > 0:
            title = 'Error reading observations'
            msg = self.sp.read_obs_error
            if self.showErrorBox:
                QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
            else:
                log_.warn('{}: {}'.format(title, msg), calling=self.calling)

        if ( self.sp.get_conf('fic_cosmetik') is None or
             self.sp.get_conf('fic_cosmetik') == '' ):
            self.sp.set_conf('do_cosmetik', False)
        
        if self.sp.get_conf('do_synth') and self.sp.get_conf('do_cosmetik'):
            self.match_cosmetic_phyat_files()
            if self.sp.get_conf('clean_cosmetic_file'):
                self.clean_cosmetic_file()
            if self.sp.get_conf('order_cosmetic_file'):
                self.order_cosmetic_file()
        
        self.set_status_text()
        self.axes = None
        self.sp.ax2_fontsize = 6
        self.sp_norm_box.setText('{}'.format(self.sp.get_conf('sp_norm')))   
        self.obj_velo_box.setText('{}'.format(self.sp.get_conf('obj_velo')))   
        self.ebv_box.setText('{}'.format(self.sp.get_conf('e_bv', 0)))   
        self.resol_box.setText('{}'.format(self.sp.get_conf('resol')))   
        self.cut2_box.setText('{}'.format(self.sp.get_conf('cut_plot2')))
        self.magenta_box.setText('{}'.format(self.sp.plot_magenta))
        self.magenta_label_box.setText('{}'.format(self.sp.label_magenta))
        self.cyan_box.setText('{}'.format(self.sp.plot_cyan))
        self.cyan_label_box.setText('{}'.format(self.sp.label_cyan))
        self.sp_min_box.setText('{}'.format(self.sp.get_conf('limit_sp')[0]))
        self.sp_max_box.setText('{}'.format(self.sp.get_conf('limit_sp')[1]))

        """
        if self.sp.get_conf('x_plot_lims') == None:
            self.sp.set_conf('x_plot_lims', self.sp.get_conf('limit_sp'))
        if self.sp.get_conf('y1_plot_lims') == None:
            self.sp.set_conf('y1_plot_lims', [0,1000])
        if self.sp.get_conf('y3_plot_lims') == None:
            self.sp.set_conf('y3_plot_lims', [-100,100])
        
        self.xlim_min_box.setText('{}'.format(self.sp.get_conf('x_plot_lims')[0]))
        self.xlim_max_box.setText('{}'.format(self.sp.get_conf('x_plot_lims')[1]))
        self.y1lim_min_box.setText('{}'.format(self.sp.get_conf('y1_plot_lims')[0]))
        self.y1lim_max_box.setText('{}'.format(self.sp.get_conf('y1_plot_lims')[1]))
        self.y3lim_min_box.setText('{}'.format(self.sp.get_conf('y3_plot_lims')[0]))
        self.y3lim_max_box.setText('{}'.format(self.sp.get_conf('y3_plot_lims')[1]))
        """

        self.init_axes()
        """
        log_.debug('x_plot_lims={}'.format(self.x_plot_lims), calling='start_spectrum')
        log_.debug('y1_plot_lims={}'.format(self.y1_plot_lims), calling='start_spectrum')
        log_.debug('y3_plot_lims={}'.format(self.y3_plot_lims), calling='start_spectrum')
        """
        self.xlim_min_box.setText('{}'.format(self.x_plot_lims[0]))
        self.xlim_max_box.setText('{}'.format(self.x_plot_lims[1]))
        self.y1lim_min_box.setText('{}'.format(self.y1_plot_lims[0]))
        self.y1lim_max_box.setText('{}'.format(self.y1_plot_lims[1]))
        self.y3lim_min_box.setText('{}'.format(self.y3_plot_lims[0]))
        self.y3lim_max_box.setText('{}'.format(self.y3_plot_lims[1]))

        self.verbosity_ag.actions()[self.sp.get_conf('log_level', 0)].setChecked(True)
        self.line_tick_ax_ag.actions()[self.sp.get_conf('line_tick_ax', 0)].setChecked(True)
        self.line_tick_pos_ag.actions()[self.sp.get_conf('line_tick_pos', 0)].setChecked(True)
        self.residual_GroupBox.setChecked(self.sp.get_conf('qt_plot_residuals', True))
        self.selected_ions_action.setChecked(self.sp.get_conf('show_selected_ions_only', False))
        self.ion_cb.setChecked(self.sp.get_conf('show_selected_ions_only', False))
        self.selected_intensities_action.setChecked(self.sp.get_conf('show_selected_intensities_only', False))
        self.cut_cb.setChecked(self.sp.get_conf('show_selected_intensities_only', False))
        self.diff_lines_ag.actions()[self.sp.get_conf('diff_lines_by', 0)].setChecked(True)
        self.line_tick_ax_ag.actions()[self.sp.get_conf('line_tick_ax', 0)].setChecked(True)
        self.editing_lines_action.setChecked(self.sp.get_conf('qt_allow_editing_lines', False))
        self.update_lines_action.setChecked(self.sp.get_conf('qt_update_after_editing_lines', False))
        self.plot_cont_action.setChecked(self.sp.get_conf('cont_plot', False))
        self.show_line_ticks_action.setChecked(self.sp.get_conf('show_line_ticks', False))
        self.plot_lines_action.setChecked(self.sp.get_conf('plot_lines_of_selected_ions', False))
        self.lineIDs_GroupBox.setChecked(self.sp.get_conf('show_line_ticks', False) or self.sp.get_conf('plot_lines_of_selected_ions', False))
        try:
            selected_ions = self.sp.get_conf('selected_ions')
            s = ''
            for ion in selected_ions:
                s = s + ion + ', '
            if not s == '':
                s = s[:-2]
            self.ion_box.setText(s)
            self.set_ion()
        except:
            self.ion_box.setText('')
        self.line_sort_ag.actions()[self.sp.get_conf('save_lines_sort', 0)].setChecked(True)
        self.show_header_action.setChecked(self.sp.get_conf('save_lines_header', False))
        self.get_line_fields_to_print()

        self.readOnlyCells_bg_color = QtGui.QColor('white')
        self.editableCells_bg_color = QtGui.QColor('lightgreen')

        if 'linux' in sys.platform and 'Plastique' in self.style_list:
            default_style = 'Plastique'
        elif 'darwin' in sys.platform and 'Macintosh (aqua)' in self.style_list:
            default_style = 'Macintosh (aqua)'
        else:
            default_style = self.style_list[0]
        
        if self.sp.get_conf('qt_style') not in self.style_list:
            
            if 'QT_STYLE' in os.environ:
                if os.environ['QT_STYLE'] in self.style_list:
                    self.sp.set_conf('qt_style', os.environ['QT_STYLE'])
                else:
                    log_.warn('Unknown Qt style {}, using {}'.format(os.environ['QT_STYLE'], default_style))
                    self.sp.set_conf('qt_style', default_style) 
            else:
                self.sp.set_conf('qt_style', default_style)  
        index_style = self.style_list.index(self.sp.get_conf('qt_style'))
        self.style_ag.actions()[index_style].setChecked(True)
        QtGui.qApp.setStyle(self.sp.get_conf('qt_style'))
        self.enable_tooltips_action.setChecked(self.sp.get_conf('qt_enable_tooltips', True))
        self.enable_tooltips_action_clicked()
        self.adjust_fig_action.setChecked(self.sp.get_conf('fig_adjust', True))

    def sp_norm(self):
        if self.sp is None:
            return
        if not self.validate_sp_norm():
            return
        old_sp_norm = self.sp.get_conf('sp_norm')
        new_sp_norm = np.float(self.sp_norm_box.text())
        if old_sp_norm == new_sp_norm:
            return
        log_.message('Changing sp_norm. Old: {}, New: {}'.format(old_sp_norm, new_sp_norm), calling=self.calling)
        self.statusBar().showMessage('Changing intensity scale of the observed spectrum ...') 
        QtGui.QApplication.processEvents() 
        self.sp.renorm(new_sp_norm)
        self.on_draw()

    def obj_velo(self):
        if self.sp is None:
            return
        if not self.validate_obj_velo():
            return
        old_obj_velo = self.sp.get_conf('obj_velo')
        new_obj_velo = np.float(self.obj_velo_box.text())
        if old_obj_velo == new_obj_velo:
            return
        self.sp.iterpolate_velocity = False
        self.sp.set_conf('obj_velo', new_obj_velo)
        log_.message('Changing obj_velo. Old: {}, New: {}'.format(old_obj_velo, new_obj_velo), calling=self.calling)
        self.statusBar().showMessage('Executing doppler correction of the observed spectrum ...') 
        QtGui.QApplication.processEvents() 
        self.sp.init_obs(obj_velo=new_obj_velo)
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run(do_synth = self.sp.do_synth, do_read_liste = True, do_profiles=False)
        self.on_draw()

    def ebv(self):
        if self.sp is None:
            return
        if not self.validate_ebv():
            return
        old_ebv = self.sp.get_conf('e_bv')
        new_ebv = np.float(self.ebv_box.text())
        if old_ebv == new_ebv and not self.cont_par_changed:
            return
        log_.message('Changing E B-V. Old: {}, New: {}'.format(old_ebv, new_ebv), calling=self.calling)
        self.statusBar().showMessage('Changing color excess E(B-V) ...', 4000) 
        self.statusBar().showMessage('Executing reddening correction of the synthetic spectrum ...') 
        QtGui.QApplication.processEvents() 
        self.sp.set_conf('e_bv', new_ebv)
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run(do_synth = self.sp.do_synth, do_read_liste = False, do_profiles=False)
        self.on_draw()
        self.cont_par_changed = False
            
    def rerun(self):
        if not self.validate_synthesis_parameters():
            return
        if ( self.x_plot_lims[0] < np.float(self.sp_min_box.text()) or
             self.x_plot_lims[1] > np.float(self.sp_max_box.text()) ):
            self.xlim_min_box.setText(self.sp_min_box.text())
            self.xlim_max_box.setText(self.sp_max_box.text())
        self.statusBar().showMessage('Rerunning synthesis ...') 
        QtGui.QApplication.processEvents() 
        self.sp.set_conf('limit_sp', (np.float(self.sp_min_box.text()), np.float(self.sp_max_box.text())))
        self.sp.set_conf('resol', np.int(self.resol_box.text()))
        self.sp.set_conf('obj_velo', np.float(self.obj_velo_box.text()))
        self.sp.set_conf('sp_norm', np.float(self.sp_norm_box.text()))
        self.sp.set_conf('e_bv', np.float(self.ebv_box.text()))
        self.sp.init_obs()
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run()
        self.set_plot_limits_and_draw()

    def adjust(self):
        if self.sp is None:
            return
        self.sp.errorMsg = ''
        self.statusBar().showMessage('Running update ...')
        QtGui.QApplication.processEvents() 
        self.sp_norm()
        self.obj_velo()
        self.ebv()
        if self.sp.errorMsg:
            if self.showErrorBox:
                msg = self.sp.errorMsg
                QtGui.QMessageBox.warning(self, 'Update error', msg, QtGui.QMessageBox.Ok)
                return 0
        ndiff, errorMsg = self.sp.adjust()
        if ndiff == -1:
            self.sp.do_cosmetik = False
            self.sp.set_conf('do_cosmetik', False)   
            self.sp.fic_cosmetik      
            self.set_status_text()
            title = 'Error in cosmetic file'
            msg = 'Unable to read from file \'{}\'\nChanging to \'no cosmetic\':\n{}'.format(self.sp.get_conf('fic_cosmetik'), errorMsg)
            if self.showErrorBox:
                QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
            else:
                log_.warn('{}: {}'.format(title, msg), calling=self.calling)
        if ndiff > 0:
            self.on_draw()

        """
        mvfc:  QObject::installEventFilter(): Cannot filter events for objects in a different thread.
               Segmentation fault

        self.line_info()
        """

        self.statusBar().showMessage('Update finished.', 4000)
        return ndiff 

    def apply_post_proc(self):
        path = str(self.post_proc_file or '')
        file_choices = "Python files (*.py) (*.py);;All files (*) (*)"
        title = 'Open post-process file'
        path = unicode(QtGui.QFileDialog.getOpenFileName(self, title, path, file_choices))
        path = path.split('/')[-1]
        if not path:
            return
        try:
            user_module = {}
            execfile(path, user_module)
            self.post_proc = user_module['post_proc']
            self.post_proc_file = path
            log_.message('function post_proc read from {}'.format(self.post_proc_file))
        except:
            self.post_proc = None
            title = 'Error reading post-process file'
            msg = 'Unable to read post-process file \'{}\''.format(path)
            if self.showErrorBox:
                QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok)
            else:
                log_.warn(msg, calling = self.calling)
            return
        try:
            self.post_proc(self.fig)
            self.canvas.draw()
        except:
            title = 'Error executing post-process'
            msg = 'Error in post-process file \'{}\''.format(self.post_proc_file)
            if self.showErrorBox:
                QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok)
            else:
                log_.warn(msg, calling = self.calling)

    def update_profile(self):
        if self.sp is None:
            return
        self.sp.run(do_synth = True, do_read_liste = False, do_profiles=True)
        self.on_draw()

    def cut2(self):
        if self.sp is None:
            return
        if not self.validate_cut():
            return    
        self.selected_intensities_action.setChecked(True)
        self.sp.set_conf('show_selected_intensities_only', True) 
        self.cut_cb.setChecked(True)
        self.draw_ion()

    def get_ion_str(self,s):
        s = s.strip()
        s = s.replace(' ', '_')
        if s.isdigit():
            line = self.sp.get_line_from_reduce_code(s)
            if line is None:
                s = ''
            else: 
                s = self.sp.fieldStrFromLine(line,'id').strip()
        return s

    def set_ion(self):
        if self.sp is None:
            return
        sList = []
        s = self.ion_box.text()
        k = s.indexOf(',')
        while k >= 0:
            s0 = self.get_ion_str(str(s[:k]))
            if s0 != '' and s0 != '*':
                sList.append(s0)
            s = s[k+1:]
            k = s.indexOf(',')
        s0 = self.get_ion_str(str(s))
        if s0 != '' and s0 != '*':
            sList.append(s0)
        s = ''
        for s0 in sList:
            s = s + s0 + ', '
        s = s[:-2]
        for item in sList[:]:
            sList.remove(item)
            if item[-1] == '*':
                item = item[:-1]
                this_ion_only = False
            else:
                this_ion_only = True
            self.sp.set_ion_list()
            if item.ljust(9) in self.sp.liste_raies['id']:
                if self.sp.true_ion(item) == item or this_ion_only:
                    sList = sList + [item]
                    if not this_ion_only:
                        sList = sList + self.sp.get_all_ions_from_ion(item)
            elif item.ljust(9) in self.sp.sp_theo['raie_ref']['id']:
                if self.sp.true_ion(item) == item or this_ion_only:
                    sList = sList + [item]
                    if not this_ion_only:
                        sList = sList + self.sp.get_all_ions_from_ion(item)
            else:
                ion_list = self.sp.get_ions_from_element(item)
                sList = sList + ion_list
        self.sp.set_conf('selected_ions', sList)
        self.ion_box.setText(s)
        
    def set_refline_to_info_box(self,j):
        if self.sp.get_conf('diff_lines_by') == 0 and len(self.sp.selected_ions_data) > 0: 
            if j == -1:
                j = 0
            s = str(self.sp.selected_ions_data[j][2][0])
            self.line_info_box.setText(s)
        
    def draw_ion(self):
        if self.cut_cb.isChecked():
            if self.validate_cut():
                self.sp.set_conf('cut_plot2', np.float(self.cut2_box.text()))
            else:
                return
        self.set_ion()
        self.sp.set_conf('index_of_current_ion', -1)
        self.sp.set_selected_ions_data()
        self.set_refline_to_info_box(-1)
        self.on_draw()
        
    def line_info(self):
        if self.sp is None:
            return
        msg = ''
        s = str(self.line_info_box.text())
        if s == '':
            return
        w = self.sp.field_width['num'] - 1
        s = s[-w:]
        if s[0] == '0':
            s = s[1:]
        self.line_info_box.setText(s)
        try:
            new_ref = int(s)
        except ValueError:
            msg = 'Invalid input.\n It is not an integer'
        if msg == '':
            line = self.sp.get_line_from_reduce_code(s)
            if line is None:
                msg = 'No line unambiguously associated with this number.'
        if msg == '':
            s = self.sp.fieldStrFromLine(line,'num').strip()
            self.line_info_box.setText(s) 
            self.line_info_ref = int(s)
            if self.sp.get_conf('qt_show_dialogs', True):
                self.show_line_info_dialog()
            else:
                self.sp.line_info(new_ref, sort='i_rel')
        else:
            title = 'Error in line number'
            QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
        
    def magenta_line(self):
        if self.sp is None:
            return
        ref_str = self.magenta_box.text()
        ref_txt = self.magenta_label_box.text()
        if ref_str == '':
            self.sp.plot_magenta = None
            self.sp.label_magenta = ''
            self.on_draw()
        else:
            new_ref = np.int(ref_str)
            self.sp.plot_magenta = new_ref
            self.sp.label_magenta = ref_txt
            self.on_draw()
        
    def cyan_line(self):
        if self.sp is None:
            return
        ref_str = self.cyan_box.text()
        ref_txt = self.cyan_label_box.text()
        if ref_str == '':
            self.sp.plot_cyan = None
            self.sp.label_cyan = ''
            self.on_draw()
        else:
            new_ref = np.int(ref_str)
            self.sp.plot_cyan = new_ref
            self.sp.label_cyan = ref_txt
            self.on_draw()
        
    def diff_lines(self):
        self.sp.set_conf('index_of_current_ion', -1)
        self.set_plot_ax2()
        if self.sp.get_conf('diff_lines_by') == 0 and len(self.sp.selected_ions_data) > 0:
            s = str(self.sp.selected_ions_data[0][2][0])
            self.line_info_box.setText(s)
        
    def set_plot_ax2(self):
        self.sp.set_selected_ions_data()
        k = self.line_tick_ax_list.index(self.line_tick_ax_ag.checkedAction().text())
        self.sp.set_conf('line_tick_ax',k)
        k = self.line_tick_pos_list.index(self.line_tick_pos_ag.checkedAction().text())
        self.sp.set_conf('line_tick_pos',k)
        k = self.diff_lines_list.index(self.diff_lines_ag.checkedAction().text())
        self.sp.set_conf('diff_lines_by',k)
        if self.show_line_ticks_action.isChecked():
            self.make_axes()
        
    def verbosity(self):
        verbosity = self.verbosity_list.index(self.verbosity_ag.checkedAction().text())
        if verbosity == log_.level:
            return
        log_.debug('Verbosity changed from {} to {}'.format(log_.level, verbosity), calling=self.calling)
        log_.level = verbosity
        self.sp.set_conf('log_level', verbosity)

    def style(self):
        new_style_str = str(self.style_ag.checkedAction().text())
        old_style_str = self.sp.get_conf('qt_style')
        if new_style_str == old_style_str:
            return
        self.sp.set_conf('qt_style', new_style_str)
        QtGui.qApp.setStyle(new_style_str)
        log_.debug('Widget style changed from {} to {}'.format(old_style_str, new_style_str), calling=self.calling)
            
    def update_lim_boxes(self):
        xformat = '{:.1f}'
        yformat = '{1:.{0}f}'
        min_diff = 2
        if abs(self.x_plot_lims[1] - self.x_plot_lims[0]) < min_diff:
            m = (self.x_plot_lims[0] + self.x_plot_lims[1])/2
            x_lims = (m - min_diff/2,m + min_diff/2)
        else:
            x_lims = self.x_plot_lims
        min_diff = 0.2
        if abs(self.y1_plot_lims[1] - self.y1_plot_lims[0]) < min_diff:
            m = (self.y1_plot_lims[0] + self.y1_plot_lims[1])/2
            y1_lims = (m - min_diff/2,m + min_diff/2)
        else:
            y1_lims = self.y1_plot_lims
        min_diff = 0.2
        if abs(self.y3_plot_lims[1] - self.y3_plot_lims[0]) < min_diff:
            m = (self.y3_plot_lims[0] + self.y3_plot_lims[1])/2
            y3_lims = (m - min_diff/2,m + min_diff/2)
        else:
            y3_lims = self.y3_plot_lims
        if self.x_plot_lims[0] != np.float(self.xlim_min_box.text()):
            self.xlim_min_box.setText(xformat.format(x_lims[0]))
        if self.x_plot_lims[1] != np.float(self.xlim_max_box.text()):
            self.xlim_max_box.setText(xformat.format(x_lims[1]))
        delta = abs(y1_lims[1]-y1_lims[0])
        if delta < 2:
            precision = 2
        else:
            precision = 1
        if self.y1_plot_lims[0] != np.float(self.y1lim_min_box.text()):
            self.y1lim_min_box.setText(yformat.format(precision, y1_lims[0]))
        if self.y1_plot_lims[1] != np.float(self.y1lim_max_box.text()):
            self.y1lim_max_box.setText(yformat.format(precision, y1_lims[1]))
        delta = abs(y3_lims[1]-y3_lims[0])
        if delta < 2:
            precision = 2
        else:
            precision = 1
        if self.y3_plot_lims[0] != np.float(self.y3lim_min_box.text()):
            self.y3lim_min_box.setText(yformat.format(precision, y3_lims[0]))
        if self.y3_plot_lims[1] != np.float(self.y3lim_max_box.text()):
            self.y3lim_max_box.setText(yformat.format(precision, y3_lims[1]))
        self.set_plot_limits_and_draw()

    def validate_input(self, editBox, field, title, varType = 'float', showError = True):
        value = editBox.text()
        if value == None:
            return False
        if ( ( varType == 'float' and not self.isFloat(value) ) or \
             ( varType == 'integer' and not self.isInteger(value) ) or \
             ( varType == 'positive integer' and not self.isPositiveInteger(value) ) or \
             ( varType == 'positive odd integer' and not self.isPositiveOdd(value) ) ):
            msg = '{} should be a {}'.format(field, varType)
            msg.replace('a integer', 'an integer')
            editBox.setFocus()
            if showError:
                if self.showErrorBox:
                    QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
                else:
                    log_.warn('{}: {}'.format(title, msg), calling=self.calling)
            return False
        else:
            return True

    def validate_sp_min(self):
        return self.validate_input(self.sp_min_box, 'xmin for the synthesis', 'Input error', 'float')

    def validate_sp_max(self):
        return self.validate_input(self.sp_max_box, 'xmax for the synthesis', 'Input error', 'float')

    def validate_sp_norm(self):
        return self.validate_input(self.sp_norm_box, 'normalization factor', 'Input error', 'float')

    def validate_ebv(self):
        return self.validate_input(self.ebv_box, 'color excess E(B-V)', 'Input error', 'float')

    def validate_obj_velo(self):
        return self.validate_input(self.obj_velo_box, 'radial velocity', 'Input error', 'float')

    def validate_resol(self):
        return self.validate_input(self.resol_box, 'rebinning factor', 'Input error', 'positive odd integer')

    def validate_xlim_min(self, showError = True):
        return self.validate_input(self.xlim_min_box, 'xmin', 'Invalid plot limit', 'float', showError)  

    def validate_xlim_max(self, showError = True):
        return self.validate_input(self.xlim_max_box, 'xmax', 'Invalid plot limit', 'float', showError)  

    def validate_y1lim_min(self):
        return self.validate_input(self.y1lim_min_box, 'ymin', 'Invalid plot limit', 'float')  

    def validate_y1lim_max(self):
        return self.validate_input(self.y1lim_max_box, 'ymax', 'Invalid plot limit', 'float')  

    def validate_y3lim_min(self):
        return self.validate_input(self.y3lim_min_box, 'residual ymin', 'Invalid plot limit', 'float')  

    def validate_y3lim_max(self):
        return self.validate_input(self.y3lim_max_box, 'residual ymax', 'Invalid plot limit', 'float')  

    def validate_cut(self):
        return self.validate_input(self.cut2_box, 'cut', 'Input error', 'float')  

    def sp_lim_in_range(self):
        xmin = np.float(self.sp_min_box.text())
        xmax = np.float(self.sp_max_box.text())
        if ( xmin < xmax - 9.999 ) and ( xmin > 0. ) and ( xmax < 200000000.):
            return True
        else:
            if self.showErrorBox:
                QtGui.QMessageBox.critical(self, 'Invalid synthesis limits', 'The acceptable values are:\n\n xmax - xmin > 10,\n xmin > 0,\n xmax < 200000000.', 
                                           QtGui.QMessageBox.Ok )
            else:
                log_.warn('Invalid synthesis limits', 'The acceptable values are:\n\n xmax - xmin > 10,\n xmin > 0,\n xmax < 200000000.', calling=self.calling)
            return False

    def validate_synthesis_parameters(self):
        return ( self.validate_sp_min() and
                 self.validate_sp_max() and
                 self.sp_lim_in_range() and
                 self.validate_sp_norm() and
                 self.validate_obj_velo() and
                 self.validate_ebv() and
                 self.validate_resol() )

    def validate_plot_parameters(self):
        return ( self.validate_xlim_min() and
                 self.validate_xlim_max() and
                 self.validate_y1lim_min() and
                 self.validate_y1lim_max() and
                 self.validate_y3lim_min() and
                 self.validate_y3lim_max() )

    def set_plot_limits_and_draw(self):
        if not self.validate_plot_parameters():
            return
        self.x_plot_lims = (np.float(self.xlim_min_box.text()), np.float(self.xlim_max_box.text()))
        self.y1_plot_lims = (np.float(self.y1lim_min_box.text()), np.float(self.y1lim_max_box.text()))
        self.y3_plot_lims = (np.float(self.y3lim_min_box.text()), np.float(self.y3lim_max_box.text()))
        self.sp.set_conf('x_plot_lims', self.x_plot_lims)
        self.sp.set_conf('y1_plot_lims', self.y1_plot_lims)
        self.sp.set_conf('y3_plot_lims', self.y3_plot_lims)
        self.restore_axes()
        self.draw_ion()

    def set_limit_sp(self):
        if not ( self.validate_sp_min() and 
                 self.validate_sp_max() and
                 self.sp_lim_in_range() ):
            return
        limit_sp = (np.float(self.sp_min_box.text()), np.float(self.sp_max_box.text()))
        self.sp.set_conf('limit_sp', limit_sp)
    
    def set_limit_sp_and_run(self):
        if str(self.sp_min_box.text()).strip() == '':
            self.sp_min_box.setText('{:.1f}'.format(self.sp.w_min))
        if str(self.sp_max_box.text()).strip() == '':
            self.sp_max_box.setText('{:.1f}'.format(self.sp.w_max))
        if not ( self.validate_sp_min() and 
                 self.validate_sp_max() and
                 self.sp_lim_in_range() ):
            return
        old_limit_sp = self.sp.get_conf('limit_sp')
        new_limit_sp = (np.float(self.sp_min_box.text()), np.float(self.sp_max_box.text()))
        if old_limit_sp == new_limit_sp:
            if not self.axes_fixed:
                self.xlim_min_box.setText(self.sp_min_box.text())
                self.xlim_max_box.setText(self.sp_max_box.text())
                self.set_plot_limits_and_draw()
            return
        if not self.validate_xlim_min(False):
            self.xlim_min_box.setText(self.sp_min_box.text())
        if not self.validate_xlim_max(False):
            self.xlim_max_box.setText(self.sp_max_box.text())            
        if ( np.float(self.xlim_min_box.text()) >= new_limit_sp[1] or
             np.float(self.xlim_max_box.text()) <= new_limit_sp[0] ):
            self.xlim_min_box.setText(self.sp_min_box.text())
            self.xlim_max_box.setText(self.sp_max_box.text())
        """
        else:
            if not self.axes_fixed:
                if np.float(self.xlim_min_box.text()) < new_limit_sp[0]:
                    self.xlim_min_box.setText(self.sp_min_box.text())
                if np.float(self.xlim_max_box.text()) > new_limit_sp[0]:
                    self.xlim_max_box.setText(self.sp_max_box.text())
        """
        self.sp.set_conf('limit_sp', new_limit_sp)
        log_.message('Changing limit_sp. Old: {}, New: {}'.format(old_limit_sp, new_limit_sp), calling=self.calling)
        self.statusBar().showMessage('Changing the synthesis wavelength limits ...') 
        QtGui.QApplication.processEvents() 
        self.sp.init_obs()
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run(do_synth = True, do_read_liste = True, do_profiles=False)
        self.set_plot_limits_and_draw()
        
    def resol(self):
        if self.sp is None:
            return
        if not self.validate_resol():
            return
        old_resol = self.sp.get_conf('resol')
        new_resol = np.int(self.resol_box.text())
        if old_resol == new_resol:
            return
        self.sp.set_conf('resol', new_resol)
        log_.message('Changing resol. Old: {}, New: {}'.format(old_resol, new_resol), calling=self.calling)
        self.statusBar().showMessage('Changing rebinning factor ...') 
        QtGui.QApplication.processEvents() 
        self.sp.set_conf('resol', new_resol)
        self.sp.init_obs()
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run(do_synth = True, do_read_liste = True, do_profiles=False)
        self.on_draw()
        
    def leave_fig(self, event):
        self.sp.firstClick = True
        if ( self.x_plot_lims != self.axes.get_xlim() or 
             self.y1_plot_lims != self.axes.get_ylim() or
             ( self.axes3 is not None and self.y3_plot_lims != self.axes3.get_ylim() ) ):
            limits_changed = True
        else:
            limits_changed = False
        if not self.axes_fixed and limits_changed:
            self.save_axes()
            self.update_lim_boxes()
        
    def fix_axes(self):
        if self.fix_axes_cb.isChecked():
            self.axes_fixed = True
        else:
            self.axes_fixed = False
            
    def get_line_fields_to_print(self):
        field_list = self.sp.get_conf('save_lines_fields')
        for i in range(0,len(self.line_field_menu.actions())):
            if self.line_print_dic.keys()[i] in field_list:
                self.line_field_menu.actions()[i].setChecked(True)
            else:
                self.line_field_menu.actions()[i].setChecked(False)
              
    def set_show_header(self):
        if self.show_header_action.isChecked():
            self.sp.set_conf('save_lines_header', True)
        else:
            self.sp.set_conf('save_lines_header', False)
      
    def set_line_fields_to_print(self):
        s = []
        for i in range(0,len(self.line_field_menu.actions())):
            if self.line_field_menu.actions()[i].isChecked():
                s.append( self.line_print_dic.keys()[i])
        self.sp.set_conf('save_lines_fields', s)
        
    def save_lines(self):
        self.sp.save_lines()
        path = self.sp.get_conf('save_lines_filename')
        self.statusBar().showMessage('Lines saved to file %s' % path, 4000)
    
    def save_lines_as(self):
        file_choices = "Text files (*.txt *.dat) (*.txt *.dat);;Tex files (*.tex) (*.tex);;CSV files (*.csv) (*.csv);;All Files (*) (*)"
        filename = self.sp.get_conf('save_lines_filename')
        extension = os.path.splitext(filename)[1][1:].lower()
        if extension in ['txt','dat']:
            selectedFilter = 'Text files (*.txt *.dat) (*.txt *.dat)'
        elif extension in ['tex']:
            selectedFilter = 'Tex files (*.tex) (*.tex)'
        elif extension in ['csv']:
            selectedFilter = 'CSV files (*.csv) (*.csv)'
        else:
            selectedFilter = 'All Files (*) (*)'
        path = unicode(QtGui.QFileDialog.getSaveFileName(self, 'Save lines to file', filename, file_choices, selectedFilter))
        if path:
            self.sp.set_conf('save_lines_filename', path)
            self.sp.save_lines()
            self.statusBar().showMessage('Lines saved to file %s' % path, 4000)

    def line_sort(self):
        k = self.line_sort_list.index(self.line_sort_ag.checkedAction().text())
        self.sp.set_conf('save_lines_sort',k)

def main_loc(init_filename=None, post_proc_file=None):
    app = QtGui.QApplication(sys.argv)
    form = AppForm(init_filename=init_filename, post_proc_file=post_proc_file)
    form.show()
    app.exec_()
    return form.fig
        
def main_loc_obj(init_filename=None, post_proc_file=None):
    app = QtGui.QApplication(sys.argv)
    form = AppForm(init_filename=init_filename, post_proc_file=post_proc_file)
    form.show()
    app.exec_()
    return form
        
def main():
    parser = get_parser()
    args = parser.parse_args()
    log_.level = args.verbosity    
    app = QtGui.QApplication(sys.argv)
    form = AppForm(init_filename=args.file, post_proc_file=args.post_proc)
    #import pdb
    #pdb.set_trace()
    form.show()
    app.exec_()
    
if __name__ == "__main__":
    main()
