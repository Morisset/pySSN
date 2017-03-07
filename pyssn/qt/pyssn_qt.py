"""
This is the window manager part of pySSN

pySSN is available under the GNU licence providing you cite the developpers names:

    Ch. Morisset (Instituto de Astronomia, Universidad Nacional Autonoma de Mexico)

    D. Pequignot (Meudon Observatory, France)

Inspired by a demo code by: 
Eli Bendersky (eliben@gmail.com)
"""
import sys, os
import argparse

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
        self.init_file_name = None
        self.call_on_draw = True
        self.cursor_on = False
        self.line_info_ref = 0
        self.x_plot_lims = None
        self.y1_plot_lims = None
        self.y2_plot_lims = None
        self.y3_plot_lims = None
        self.post_proc_file = post_proc_file
        self.do_save = True
        self.cont_par_changed = False
        self.axes_fixed = False

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        self.select_init(init_filename)
        self.cont_pars_dialog = None
        self.nearbyLines_dialog = None
        self.line_info_dialog = None

    def closeEvent(self, evnt):
        if self.cont_pars_dialog is not None:
          self.cont_pars_dialog.close()
        if self.nearbyLines_dialog is not None:
          self.nearbyLines_dialog.close()
        if self.line_info_dialog is not None:
          self.line_info_dialog.close()
        
    def image_extension_list(self):
        filetypes = self.canvas.get_supported_filetypes()
        file_extensions = filetypes.keys()
        file_extensions.sort()
        return file_extensions
        
    def image_filter(self):
        filetypes = self.canvas.get_supported_filetypes_grouped()
        imagetype_list = filetypes.keys()
        imagetype_list.sort()
        s = ''
        for imagetype in imagetype_list:
          extension_list = filetypes[ imagetype ]
          s = s + str(imagetype)
          s1 = ' (*.' + str(extension_list[0])
          for extension in extension_list[1:]:
            s1 = s1 + ' *.' + str(extension)
          s1 = s1 + ')'
          s = s + s1 + s1 + ';;'
        return s

    def set_save_plot_action_tip(self):
        plotFile = self.sp.get_conf('plot_filename')
        path, filename = os.path.split(plotFile)
        if path == os.getcwd():
          plotFile = filename
        s = "Save plot to file '" + plotFile + "' (initially set with 'plot_filename = <filename>'; " \
            "use option 'Save plot as' to change the file name and image format)" 
        self.save_plot_action.setStatusTip(s)

    def save_plot(self):
        path = self.sp.get_conf('plot_filename')
        self.canvas.print_figure(path, dpi=self.dpi)
        self.statusBar().showMessage('Saved to %s' % path, 2000)
      
    def save_plot_as(self):
        file_choices = self.image_filter()
        path = self.sp.get_conf('plot_filename')
        path = unicode(QtGui.QFileDialog.getSaveFileName(self, 'Save plot to file', path, file_choices))
        extension = os.path.splitext(path)[1][1:].lower()
        if path:
            if extension in self.image_extension_list():
              self.sp.set_conf('plot_filename', path)
              self.canvas.print_figure(path, dpi=self.dpi)
              self.statusBar().showMessage('Saved to %s' % path, 2000)
              self.set_save_plot_action_tip()
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
            do_print = not self.sp.get_conf('show_dialogs', True)
            self.nearbyLines = self.sp.nearby_lines(event, do_print)
            if not do_print:
              if self.nearbyLines is not None:
                self.show_nearbyLines_dialog()

    def create_main_frame(self):
        
        if self.use_workspace:
            self.main_frame = QtGui.QWorkspace()
        else:
            self.main_frame = QtGui.QWidget()
        # Create the mpl Figure and FigCanvas objects. 
        #
        self.dpi = 100
        self.fig = plt.figure(figsize=(15,15))
#        self.fig = plt.figure((figsize=(20.0, 15.0), dpi=self.dpi)
        
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
        self.connect(self.xlim_min_box, QtCore.SIGNAL('editingFinished()'), self.save_from_lim_boxes)
        self.connect(self.xlim_min_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes_and_draw)

        self.xlim_max_box = QtGui.QLineEdit()
        self.xlim_max_box.setMinimumWidth(50)
        self.connect(self.xlim_max_box, QtCore.SIGNAL('editingFinished()'), self.save_from_lim_boxes)
        self.connect(self.xlim_max_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes_and_draw)
        
        self.y1lim_min_box = QtGui.QLineEdit()
        self.y1lim_min_box.setMinimumWidth(50)
        self.connect(self.y1lim_min_box, QtCore.SIGNAL('editingFinished()'), self.save_from_lim_boxes)
        self.connect(self.y1lim_min_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes_and_draw)
        
        self.y1lim_max_box = QtGui.QLineEdit()
        self.y1lim_max_box.setMinimumWidth(50)
        self.connect(self.y1lim_max_box, QtCore.SIGNAL('editingFinished()'), self.save_from_lim_boxes)
        self.connect(self.y1lim_max_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes_and_draw)
        
        self.y3lim_min_box = QtGui.QLineEdit()
        self.y3lim_min_box.setMinimumWidth(50)
        self.connect(self.y3lim_min_box, QtCore.SIGNAL('editingFinished()'), self.save_from_lim_boxes)
        self.connect(self.y3lim_min_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes_and_draw)
        
        self.y3lim_max_box = QtGui.QLineEdit()
        self.y3lim_max_box.setMinimumWidth(50)
        self.connect(self.y3lim_max_box, QtCore.SIGNAL('editingFinished()'), self.save_from_lim_boxes)
        self.connect(self.y3lim_max_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes_and_draw)
        
        self.select_init_button = QtGui.QPushButton("Init file")
        self.connect(self.select_init_button, QtCore.SIGNAL('clicked()'), self.select_init)
        
        self.run_button = QtGui.QPushButton("Run")
        self.connect(self.run_button, QtCore.SIGNAL('clicked()'), self.rerun)
        
        self.draw_button = QtGui.QPushButton("Draw")
        self.connect(self.draw_button, QtCore.SIGNAL('clicked()'), self.on_draw)

#        self.savelines_button = QtGui.QPushButton("&Save")
#        self.savelines_button.setMinimumWidth(50)
#        self.connect(self.savelines_button, QtCore.SIGNAL('clicked()'), self.save_lines)
        
        self.Command_GroupBox = QtGui.QGroupBox("Execute")
        self.Command_GroupBox.setCheckable(False)
        
        self.ObsSpec_GroupBox = QtGui.QGroupBox("Parameters of the synthetic and observed spectra")
        self.ObsSpec_GroupBox.setCheckable(False)

        self.SpecPlot_GroupBox = QtGui.QGroupBox("Plot of spectra")
        self.SpecPlot_GroupBox.setCheckable(False)

        self.lineIDs_GroupBox = QtGui.QGroupBox("Show lines")
        self.lineIDs_GroupBox.setCheckable(True)
        self.lineIDs_GroupBox.setChecked(False)
        
#        self.lineIDs_GroupBox.setTristate(True)
#        self.lineIDs_GroupBox.setCheckState(QtCore.Qt.PartiallyChecked)
        
        self.connect(self.lineIDs_GroupBox, QtCore.SIGNAL('clicked()'), self.show_lines_clicked)
        self.lineIDs_GroupBox.setToolTip( 'Check to show ticks at the central positions of the spectral lines and plot the lines of selected ions' )        

        self.residual_GroupBox = QtGui.QGroupBox("Plot of residuals")
        self.residual_GroupBox.setCheckable(True)
        self.residual_GroupBox.setChecked(True)
        self.connect(self.residual_GroupBox, QtCore.SIGNAL('clicked()'), self.make_axes)
        self.residual_GroupBox.setToolTip( 'Check to display the residual plot' )        

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
        self.connect(self.sp_min_box, QtCore.SIGNAL('editingFinished()'), self.set_limit_sp)
        self.connect(self.sp_min_box, QtCore.SIGNAL('returnPressed()'), self.set_limit_sp_and_run)
        
        self.sp_max_box = QtGui.QLineEdit()
        self.sp_max_box.setMinimumWidth(50)
        self.connect(self.sp_max_box, QtCore.SIGNAL('editingFinished()'), self.set_limit_sp)
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
        self.line_info_box.setFixedWidth(100)
        self.connect(self.line_info_box, QtCore.SIGNAL('returnPressed()'), self.line_info)

        self.mpl_toolbar.addSeparator()
        self.mpl_toolbar.addWidget(QtGui.QLabel('   line code number '))
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
        self.run_button.setToolTip(s)        
        
        s = 'Click to update synthesis with changes in line intensities, profiles, and continuum parameters.'
        self.adjust_button.setToolTip(s)        
        
        s = 'Enter line code number to get information on\n' \
            'the reference line and on its satellites.'
        self.line_info_box.setToolTip(s)        

        s =  'Color excess E(B-V)\n\n' \
             'Set with: \n' \
             '    e_bv = <float>\n\n' \
             'Comment: \n' \
            u'    E(B-V) \u2248 C(H\u03B2) / 1.5'
        self.ebv_box.setToolTip(s)  
        
        s = 'Radial velocity in km/s\n\n' \
            'Set with: \n' \
            '    obj_velo = <float>'
        self.obj_velo_box.setToolTip(s) 

        s = 'Minimum wavelength of the synthetic spectrum (in angstroms)\n\n' \
            'Set with:  \n' \
            '    limit_sp = (<xmin>, <xmax>)'
        self.sp_min_box.setToolTip(s)        

        s = 'Maximum wavelength of the synthetic spectrum  (in angstroms)\n\n' \
            'Set with:  \n' \
            '    limit_sp = (<xmin>, <xmax>)'
        self.sp_max_box.setToolTip(s)        

        s = 'Minimum wavelength in the plots of spectra and residuals  (in angstroms)\n\n' \
            'Set with:  \n' \
            '    x_plot_lims = (<xmin>, <xmax>)'
        self.xlim_min_box.setToolTip(s)        

        s = 'Maximum wavelength in the plots of spectra and residuals  (in angstroms)\n\n' \
            'Set with:  \n' \
            '    x_plot_lims = (<xmin>, <xmax>)'
        self.xlim_max_box.setToolTip(s)        

        s = 'Minimum ordinate in the plot of spectra, in units of relative intensity \n\n' \
            'Set with:  \n' \
            '    y1_plot_lims = (<ymin>, <ymax>)'
        self.y1lim_min_box.setToolTip(s)        

        s = 'Maximum ordinate in the plot of spectra, in units of relative intensity\n\n' \
            'Set with:  \n' \
            '    y1_plot_lims = (<ymin>, <ymax>)'
        self.y1lim_max_box.setToolTip(s)        

        s = 'Minimum ordinate in the plot of residuals, in units of relative intensity\n\n' \
            'Set with:  \n' \
            '    y3_plot_lims = (<ymin>, <ymax>)'
        self.y3lim_min_box.setToolTip(s)        

        s = 'Maximum ordinate in the plot of residuals, in units of relative intensity\n\n' \
            'Set with:  \n' \
            '    y3_plot_lims = (<ymin>, <ymax>)'
        self.y3lim_max_box.setToolTip(s)        
        
        s = 'Check to retain the current limits of the plots while zooming and panning.'
        self.fix_axes_cb.setToolTip(s)        
        
        s = 'Check to show only lines with intensities above cut. \n\n' \
             'Set with: \n' \
             '    show_selected_intensities_only = <boolean>'
        self.cut_cb.setToolTip(s)        

        s = 'Check to show only lines of selected ions. \n\n' \
             'Set with: \n' \
             '    show_selected_ions_only = <boolean>'
        self.ion_cb.setToolTip(s)        
        
        s = 'Normalization factor, ratio between the intensity and the \n' \
            u'observed flux of the reference line, usually 10\u2074/F(H\u03B2)\n\n' \
             'Set with: \n' \
             '    sp_norm = <float>'
        self.sp_norm_box.setToolTip(s)        

        s = 'Rebinning factor, the integer factor by which the number of points \n' \
            'of the original spectrum is multiplied in the rebinning process\n\n' \
            'Set with: \n' \
            '    resol = <integer>\n\n' \
            'Usage: \n' \
            '    Set to \'1\' if the resolution of the observed spectrum is large enough' 
            
        self.resol_box.setToolTip(s) 
        
        s = 'Minimum relative intensity of lines to be shown. \n\n' \
             'Set with: \n' \
             '    cut_plot2 = <float>'
        self.cut2_box.setToolTip(s)        
        
        s = 'Comma-separated list of selected ions, elements, or line code numbers to be shown. \n\n' \
             'Set with: \n' \
             '    selected_ions = [<ion1>,<ion2>,...]\n\n' \
             'Examples: \n' \
             '    \'Fe III\' (or \'Fe_III\') to show the lines of Fe III\n' \
             '    \'Fe III, Fe IV\' to show the lines of Fe III and Fe IV\n' \
             '    \'Fe\' to show the lines of all Fe ions\n' \
             '    \'Fe, N\' to show the lines of all Fe and N ions\n' \
             '    <line code number> to show the lines of that same ion'                
        self.ion_box.setToolTip(s)        
              
        self.resol_box.setToolTip( s ) 

        self.verbosity_list = ['None', 'Errors', 'Errors and warnings', 'Errors, warnings, and comments', 'Debug messages' ]
        self.verbosity_button = QtGui.QPushButton('Verbosity')
        s = 'Verbosity level:\n'
        for i in range(len(self.verbosity_list)):
            s = s + '    ' + str(i) + ' - ' + self.verbosity_list[i] + '\n'
        s = s + '\nSet with:\n' + '    log_level = <integer>'
        self.verbosity_button.setToolTip( s )        
        self.verbosity_ag = QtGui.QActionGroup(self, exclusive=True)
        
        self.verbosity_menu = QtGui.QMenu()
        for i in range(len(self.verbosity_list)):
            a = self.verbosity_ag.addAction(QtGui.QAction(self.verbosity_list[i], self, checkable=True))
            self.verbosity_menu.addAction(a)
        self.verbosity_button.setMenu(self.verbosity_menu)
        self.verbosity_ag.triggered.connect(self.verbosity)
        
        #
        # Layout with box sizers
        # 
        hbox5 = QtGui.QHBoxLayout()

        for w in [self.line_info_box, self.magenta_box, self.magenta_label_box, self.cyan_box, self.cyan_label_box]:
            hbox5.addWidget(w)
            hbox5.setAlignment(w, QtCore.Qt.AlignVCenter)


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
#        vbox.addLayout(hbox5)
#        QtGui.QApplication.setStyle(QtGui.QStyleFactory.create(styleName))
#        self.changePalette()
        QtGui.qApp.setStyle('Windows')
        QtGui.qApp.setStyle('Cleanlooks')
        QtGui.qApp.setStyle('Plastique')
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
        
        self.save_plot_action = self.create_action("Save plot",
                                              shortcut="Ctrl+S", 
                                              slot=self.save_plot,
                                              tip="Save plot to default file")
        
        save_plot_as_action = self.create_action("Save plot as",
                                              shortcut="Ctrl+Shift+S", 
                                              slot=self.save_plot_as, 
                                              tip="Select file name and save plot")

        save_lines_action = self.create_action("Save lines",
                                              shortcut="Ctrl+L", 
                                              slot=self.save_lines, 
                                              tip="Save list of lines")

        save_lines_as_action = self.create_action("Save lines as",
                                              shortcut="Ctrl+Shift+L", 
                                              slot=self.save_lines_as, 
                                              tip="Select file name and save list of lines")

        self.add_actions(self.file_menu, 
            (open_init_action, None, self.save_plot_action, save_plot_as_action, None, save_lines_action, save_lines_as_action))

        self.line_sort_list = ['wavelength', 'decreasing wavelength', 'intensity', 'decreasing intensity', 'ion' , 'decreasing ion' ]
        s = 'Sort lines by:\n'
        for i in range(len(self.line_sort_list)):
            s = s + '    ' + str(i) + ' - ' + self.line_sort_list[i] + '\n'
        s = s + '\nSet with:\n' + '    line_saved_ordered_by = <integer>'
        self.line_sort_ag = QtGui.QActionGroup(self, exclusive=True)

        self.line_sort_menu = self.file_menu.addMenu("Sort lines by")
        self.line_sort_menu.setToolTip(s)        

        for i in range(len(self.line_sort_list)):
            a = self.line_sort_ag.addAction(QtGui.QAction(self.line_sort_list[i], self, checkable=True))
            self.line_sort_menu.addAction(a)

        self.line_sort_ag.triggered.connect(self.line_sort)
        
        self.line_print_dic = OrderedDict( [ ( 'id'      , 'ion' ), 
                                ( 'lambda'  , 'wavelength' ), 
                                ( 'l_shift' , 'wavelength shift' ), 
                                ( 'l_tot'   , 'corrected wavelength' ), 
                                ( 'i_rel'   , 'intensity' ), 
                                ( 'i_cor'   , 'intensity correction factor' ), 
                                ( 'i_tot'   , 'corrected intensity' ) ])
        
        items = self.line_print_dic.values()
        s = 'Fields to be printed:\n'
        for i in range(len(items)):
            s = s + '    ' + str(i) + ' - ' + items[i] + '\n'
        s = s + '\nSet with:\n' + '    line_field_print = <list>'

        self.line_field_menu = self.file_menu.addMenu("Show fields")
        self.line_field_menu.setToolTip(s)        

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
        
        quit_action = self.create_action("&Quit", 
                                         slot=self.fileQuit, 
                                         shortcut="Ctrl+Q", 
                                         tip="Close the application")

        self.add_actions(self.file_menu, (None, quit_action))

        self.run_menu = self.menuBar().addMenu("Execute")
        
        run_action = self.create_action("Run",
                                         shortcut="Ctrl+F9", 
                                         slot=self.rerun, 
                                         tip="Execute synthesis from the beginning")
        
        update_action = self.create_action("Update",
                                            shortcut="F9", 
                                            slot=self.adjust, 
                                            tip="Update synthesis with changes in line intensities, profiles, and continuum parameters")

#        self.draw_menu = self.menuBar().addMenu("Draw")
        
        draw_action = self.create_action("Draw",
                                         shortcut="F8", 
                                         slot=self.on_draw, 
                                         tip="Redraw plots")
        
        post_proc_action = self.create_action("Post-process",
                                               shortcut="Ctrl+F8", 
                                               slot=self.apply_post_proc, 
                                               tip="Edit the plots with python commands defined in an external file")
        
        self.ask_postprocfile_action = self.create_action("Ask for post-process file name",
                                               checkable=True, 
                                               tip="Check to be asked for the post-process file name before executing post-process")

        self.add_actions(self.run_menu, 
            (update_action, run_action, None, draw_action, post_proc_action, self.ask_postprocfile_action))

        self.line_menu = self.menuBar().addMenu('Lines')
        
        self.show_line_ticks_action = self.create_action('Show line ticks', 
            shortcut='Alt+L', slot=self.show_line_ticks_action_clicked, checkable=True, 
            tip='Check to show line ticks')
        
        self.plot_lines_action = self.create_action('Plot lines', 
            shortcut='Alt+P', slot=self.show_line_ticks_action_clicked, checkable=True, 
            tip='Check to plot spectra of selected ions')
        
        self.cycle_forwards_ions_action = self.create_action('Cycle forwards selected ions', 
            shortcut='Alt+0', slot=self.cycle_forwards_ions, checkable=False, 
            tip='Click to cycle forwards the selected ions')
        
        self.cycle_backwards_ions = self.create_action('Cycle backwards selected ions', 
            shortcut='Alt+9', slot=self.cycle_backwards_ions, checkable=False, 
            tip='Click to cycle backwards the selected ions')

        self.add_actions(self.line_menu, 
            (self.show_line_ticks_action, self.plot_lines_action, None, self.cycle_forwards_ions_action, self.cycle_backwards_ions, None))

        self.line_tick_ax_menu = self.line_menu.addMenu('Window of line ticks')
       
        self.line_tick_ax_list = ['Plot of spectra', 'Plot of residuals', 'Separate plot' ]
        s = 'Show line ticks on:\n'
        for i in range(len(self.line_tick_ax_list)):
            s = s + '    ' + str(i) + ' - ' + self.line_tick_ax_list[i] + '\n'
        s = s + '\nSet with:\n' + '    line_tick_ax = <integer>'
        self.line_tick_ax_ag = QtGui.QActionGroup(self, exclusive=True)
        self.line_tick_ax_menu.setToolTip(s)        

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
        self.line_tick_pos_menu.setToolTip(s)        

        for i in range(len(self.line_tick_pos_list)):
            a = self.line_tick_pos_ag.addAction(QtGui.QAction(self.line_tick_pos_list[i], self, checkable=True))
            self.line_tick_pos_menu.addAction(a)
        self.line_tick_pos_ag.triggered.connect(self.set_plot_ax2)

        self.line_tick_color_action = self.create_action('Color of line ticks', 
            shortcut=None, slot=self.line_tick_color_clicked, checkable=False, 
            tip='Set color of line ticks')
        self.line_menu.addAction(self.line_tick_color_action)
        
        self.selected_intensities_action = self.create_action('Only above the cut', 
            shortcut='Alt+K', slot=self.selected_lines_clicked, checkable=True, 
            tip='Check to show the ticks for lines with intensities above cut only')
        
        self.selected_ions_action = self.create_action('Only selected ions', 
            shortcut='Alt+I', slot=self.selected_lines_clicked, checkable=True, 
            tip='Check to show the line ticks for selected ions only')
        
        self.differentiate_lines_action = self.create_action('Differentiate line sequences', 
            shortcut='Alt+D', slot=self.differentiate_lines_clicked, checkable=True, 
            tip='Check to differentiate line sequences of a same ion according to the reference lines')
        
        self.editing_lines_action = self.create_action('Allow editing', 
            slot=self.editing_lines_clicked, checkable=True, 
            tip='Check to allow editing line parameters in line info dialog')
        
        self.update_lines_action = self.create_action('Update after editing', 
            shortcut='Alt+D', slot=self.update_lines_clicked, checkable=True, 
            tip='Check to update synthesis after editing line parameters in line info dialog')

        self.add_actions(self.line_menu, 
            (None, self.selected_intensities_action, self.selected_ions_action, self.differentiate_lines_action,
             None, self.editing_lines_action, self.update_lines_action))

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

        self.verbosity_list = ['None', 'Errors', 'Errors and warnings', 'Errors, warnings, and comments', 'Debug messages' ]
        s = 'Verbosity level:\n'
        for i in range(len(self.verbosity_list)):
            s = s + '    ' + str(i) + ' - ' + self.verbosity_list[i] + '\n'
        s = s + '\nSet with:\n' + '    log_level = <integer>'
        self.verbosity_ag = QtGui.QActionGroup(self, exclusive=True)
        
        self.verbosity_menu = self.menuBar().addMenu("Verbosity")
        self.verbosity_menu.setToolTip(s)        

        for i in range(len(self.verbosity_list)):
            a = self.verbosity_ag.addAction(QtGui.QAction(self.verbosity_list[i], self, checkable=True))
            self.verbosity_menu.addAction(a)
        self.verbosity_ag.triggered.connect(self.verbosity)

#        self.settings_menu = self.menuBar().addMenu("Settings")
        
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
      
    def ConvStrToValidTypes(self, str):
  
      def isFloat(str):
        try:
          float(str)
          return True
        except ValueError:
          return False
      
      def isBool(str):
        try:
          bool(str)
          return True
        except ValueError:
          return False

      str = str.strip(' ')
      str = str.replace('Error in ','')
      if str == '':
        result = None
      elif str.isdigit():
        result = int(str)
      elif isFloat(str):
        result = float(str)
      elif str.capitalize() == 'True':
        result = True
      elif str.capitalize() == 'False':
        result = False
      elif str.find(',') >= 0:
        try:
          str = str.replace('[','')
          str = str.replace(']','')
          result = [float(i) for i in str.split(',')]
        except:
          result = None
      else:      
        result = None
      return result
    
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
              help = self.table.item(i,2).text().toUtf8()
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
                lines.insert(0, '\n# ' + help + '\n')
                lines.insert(0, field + ' = ' + value + '\n')
                
            lines = lines[::-1]
            f = open(path, 'w')
            f.writelines(lines)
            f.close()

    # mvfc: I should have created classes for the dialogs to capture events like dragging, closing, resizing
    def close_line_info_dialog(self):
      self.line_info_dialog_width = self.line_info_dialog.width()
      self.line_info_dialog_height = self.line_info_dialog.height()
      self.line_info_dialog_x = self.line_info_dialog.pos().x()
      self.line_info_dialog_y = self.line_info_dialog.pos().y()
      self.line_info_dialog.close()
  
    def show_line_info_dialog(self):
         
       def on_doubleClick():
           item = self.line_info_table.currentItem()
           row = item.row()
           col = item.column()
           if col in [0,6]:
             self.line_info_box.setText(item.text())
             fill_line_info_table()

       def isRefLineNum(line_num_str):
           if line_num_str[-8:] == '00000000':
             return True
           else:
             return False

       def isRefLine(line):
           line_num = self.sp.fieldStrFromLine(line,'num')
           if line_num[-8:] == '00000000':
             return True
           else:
             return False

       def isSubRefLine(line):
           wavelength = float(self.sp.fieldStrFromLine(line,'lambda'))
           if not isRefLine(line) and (wavelength < 2.0):
             return True
           else:
             return False

       def isLine(line):
           if not (isRefLine(line) or isSubRefLine(line)):
             return True
           else:
             return False

       def fill_data(i, line, cat=''):
         if line == None:
           return
         editableCols = []
         if self.sp.get_conf('allow_editing_lines', False):
           if cat == 'sat':
             # mvfc: which are editable fields?
             editableCols = ['l_shift', 'i_cor']
           elif cat == 'ref':
             editableCols = ['i_rel']
         for j in range(0,len(fieldItems)):
           s = self.sp.fieldStrFromLine(line, fieldItems[j])
           s = s.strip()
           item = QtGui.QTableWidgetItem(s)
           if fieldItems[j] in editableCols:
             item.setBackgroundColor(self.editableCells_bg_color)
           else:
             item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)
             item.setBackgroundColor(self.readOnlyCells_bg_color)
           self.line_info_table.setItem(i,j,item)

       def fill_text(i, text):
           item = QtGui.QTableWidgetItem(text)
           item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)
           item.setBackgroundColor(self.readOnlyCells_bg_color)
           item.setTextAlignment(QtCore.Qt.AlignBottom)
           item.setTextColor(QtCore.Qt.blue)
           self.line_info_table.setItem(i,0,item)
           self.line_info_table.setSpan(i,0,2,len(fieldItems))

       def fill_line_info_table():  
           line_num = self.line_info_box.text()
           subrefline = None
           if isRefLineNum(line_num):
             line = None
             subrefline = None
             refline = self.sp.read_line(self.sp.fic_model, int(line_num))
             refline_num = line_num
           else:
             line = self.sp.read_line(self.sp.fic_cosmetik, int(line_num))
             if line == None:
               line = self.sp.read_line(self.sp.phyat_file, int(line_num))
             if line == None:
               subrefline = None
               refline = None
             else:
               if isSubRefLine(line):
                 subrefline = line
                 subrefline_num = line_num
                 line = None
                 line_num = None
                 refline_num = self.sp.fieldStrFromLine(subrefline, 'ref')
                 refline = self.sp.read_line(self.sp.fic_model, int(refline_num))
               else:
                 refline_num = self.sp.fieldStrFromLine(line, 'ref')
                 if isRefLineNum(refline_num):
                   refline = self.sp.read_line(self.sp.fic_model, int(refline_num))
                   subrefline = None
                 else:
                   subrefline_num = refline_num
                   subrefline = self.sp.read_line(self.sp.fic_cosmetik, int(subrefline_num))
                   if subrefline == None:
                     subrefline = self.sp.read_line(self.sp.phyat_file, int(subrefline_num))
                   refline_num = self.sp.fieldStrFromLine(subrefline, 'ref')
                   refline = self.sp.read_line(self.sp.fic_model, int(refline_num))
           if subrefline is not None:
             subsatellites = self.sp.read_satellites(self.sp.phyat_file, int(subrefline_num))
             n_subsat = len(subsatellites)
             for i in range(0,n_subsat):
               sat_line = subsatellites[i]
               sat_line_num = int(self.sp.fieldStrFromLine(sat_line,'num'))
               cosmetic_line = self.sp.read_line(self.sp.fic_cosmetik, sat_line_num)
               if cosmetic_line is not None:
                 subsatellites[i] = cosmetic_line
           else:
             n_subsat = 0
           if refline is not None:
             satellites = self.sp.read_satellites(self.sp.phyat_file, int(refline_num))
             n_sat = len(satellites)
             for i in range(0,n_sat):
               sat_line = satellites[i]
               sat_line_num = int(self.sp.fieldStrFromLine(sat_line,'num'))
               cosmetic_line = self.sp.read_line(self.sp.fic_cosmetik, sat_line_num)
               if cosmetic_line is not None:
                 satellites[i] = cosmetic_line
           else:
             n_sat = 0
           if line is None and refline is None:
             title = 'Error in line data display'
             msg = 'Line code number not found.'
             QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
             return
           self.line_info_table.setRowCount(n_sat+n_subsat+20)
           self.line_info_table.clearSpans()
           k = 0
           if line is not None:
             fill_text(k,'Line:')
             k += 2
             fill_data(k, line, 'sat')
             k += 1
           if subrefline is not None:
             fill_text(k,'Subreference line:')
             k += 2
             fill_data(k, subrefline)
             k += 1
             fill_text(k, str(n_subsat) + ' satellites:')
             k += 2
             for i in range(0,n_subsat):
               fill_data(k+i, subsatellites[i])
             k += n_subsat
           fill_text(k,'Reference line:')
           k += 2
           fill_data(k, refline, 'ref')
           k += 1
           fill_text(k, str(n_sat) + ' satellites:')
           k += 2
           for i in range(0,n_sat):
             fill_data(k+i, satellites[i])
           k += n_sat
           self.line_info_table.setRowCount(k)
           self.line_info_table.resizeColumnsToContents()
           self.line_info_table.resizeRowsToContents()
           return
 
       def rightFormat(s,col):
         try:
           r = float(s)
           s = self.sp.field_format[fieldItems[col]].format(r)
           if len(s) == self.sp.field_width[fieldItems[col]] and not np.isinf(r):
             output = s
           else:
             output = None
         except:
           output = None
         return output        
         
       def on_itemChanged():
         item = self.line_info_table.currentItem()
         if not (item.flags() & QtCore.Qt.ItemIsEditable):
           return
         self.line_info_table.blockSignals(True)
         row = item.row()
         col = item.column()
         s = item.text()
         value = rightFormat(s,col)
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
            QtGui.QMessageBox.critical(self, title, msg, QtGui.QMessageBox.Ok )
         self.line_info_table.blockSignals(False)
       
       def save_change(row, col):
         line = ' '*85
         for j in range(0,len(fieldItems)):
           s = self.line_info_table.item(row,j).text()
           width = self.sp.field_width[fieldItems[j]]
           align = self.sp.field_align[fieldItems[j]]
           pos = self.sp.field_pos[fieldItems[j]]
           s = '{:{a}{w}s}'.format(s, a=align, w=width)
           line = line[:pos] + s + line[pos:]
         line = line.rstrip()
         if col in [3,5]:
           filename = self.sp.fic_cosmetik
         else:
           filename = self.sp.fic_model
         self.sp.replace_line(filename, line)
         if self.sp.get_conf('update_after_editing_lines', False):
           self.adjust()
         
       self.line_info_dialog = QtGui.QDialog()
       self.line_info_dialog.resize(self.line_info_dialog_width,self.line_info_dialog_height)
       self.line_info_dialog.move(self.line_info_dialog_x,self.line_info_dialog_y)
       self.line_info_table = QtGui.QTableWidget()   
       self.line_info_table.setColumnCount(10)
       fieldItems = self.sp.fields
       fieldNames = [ self.sp.field_abbr[item] for item in fieldItems ]
       self.line_info_table.setHorizontalHeaderLabels(fieldNames)
       for j in range(0,len(fieldItems)):
         self.line_info_table.horizontalHeaderItem(j).setToolTip(self.sp.field_tip[fieldItems[j]])
       self.line_info_table.horizontalHeaderItem(9).setTextAlignment(QtCore.Qt.AlignLeft);
       self.line_info_table.horizontalHeaderItem(9).setText('  comment')
       fill_line_info_table()
       self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Close)
       #self.buttonBox.rejected.connect(self.line_info_dialog.close)
       self.buttonBox.rejected.connect(self.close_line_info_dialog)
       self.line_info_table.doubleClicked.connect(on_doubleClick)
       self.line_info_table.itemChanged.connect(on_itemChanged)
       vbox = QtGui.QVBoxLayout()
       vbox.addWidget(self.line_info_table)
       vbox.addWidget(self.buttonBox)
       self.line_info_dialog.setLayout(vbox)
       self.line_info_dialog.setWindowTitle('Line data')
       self.line_info_dialog.setWindowModality(QtCore.Qt.NonModal)
       self.line_info_dialog.show()
            
    def show_nearbyLines_dialog(self):

       def on_doubleClick():
         item = self.nearbyLines_table.currentItem()
         row = item.row()
         col = item.column()
         if col in [0,6]:
           self.line_info_box.setText(item.text())
           self.show_line_info_dialog()
         if col == 1:
           self.ion_box.setText(item.text())
           self.draw_ion()

       def do_selection():
         selectedItems = self.nearbyLines_table.selectedItems()
         selected_ions = []
         selected_lines = []
         for item in selectedItems:
           row = item.row()
           col = item.column()
           if col == 1:
             ion = item.text()
             if not ion in selected_ions:
               selected_ions.append(ion)
           if col in [0,6]:
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

       self.nearbyLines_dialog = QtGui.QDialog()
       self.nearbyLines_dialog.resize(800,300)
       sG = QtGui.QApplication.desktop().screenGeometry()
       x = (sG.width()-self.nearbyLines_dialog.width())
       y = (sG.height()-self.nearbyLines_dialog.height())
       self.nearbyLines_dialog.move(x,y)
       statusBar = QtGui.QStatusBar()
       statusBar.addWidget(QtGui.QLabel('Double-click on a line code number or on a reference line code number (or select and press ok) to edit parameters. '),1)
       self.nearbyLines_table = QtGui.QTableWidget()   
       self.nearbyLines_table.setRowCount(len(self.nearbyLines))
       self.nearbyLines_table.setColumnCount(10)
       fieldItems = self.sp.fields
       fieldNames = [ self.sp.field_abbr[item] for item in fieldItems ]
       self.nearbyLines_table.setHorizontalHeaderLabels(fieldNames)
       for j in range(0,len(fieldItems)):
         self.nearbyLines_table.horizontalHeaderItem(j).setToolTip(self.sp.field_tip[fieldItems[j]])
       self.nearbyLines_table.horizontalHeaderItem(9).setTextAlignment(QtCore.Qt.AlignLeft);
       self.nearbyLines_table.horizontalHeaderItem(9).setText('  comment')
       for i in range(0,len(self.nearbyLines)):
         for j in range(0,len(fieldItems)):
           s = self.sp.field_format[fieldItems[j]].format(self.nearbyLines[i][j])
           s = s.strip()
           item = QtGui.QTableWidgetItem(s)
           item.setFlags(item.flags() ^ QtCore.Qt.ItemIsEditable)
           item.setBackgroundColor(self.readOnlyCells_bg_color)
           self.nearbyLines_table.setItem(i,j,item)
       self.nearbyLines_table.resizeColumnsToContents()
       self.nearbyLines_table.resizeRowsToContents()
       self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Help|QtGui.QDialogButtonBox.Apply|QtGui.QDialogButtonBox.Close)
       self.buttonBox.rejected.connect(self.nearbyLines_dialog.close)
       self.buttonBox.clicked.connect(do_selection)
       self.nearbyLines_table.doubleClicked.connect(on_doubleClick)
       self.nearbyLines_table.verticalHeader().sectionDoubleClicked.connect(do_selection)
       vbox = QtGui.QVBoxLayout()
       vbox.addWidget(self.nearbyLines_table)
       vbox.addWidget(self.buttonBox)
       vbox.addWidget(statusBar)
       self.nearbyLines_dialog.setLayout(vbox)
       s = 'List of lines between {0:.2f} and {1:.2f} angstroms'.format(self.sp.cursor_w1, self.sp.cursor_w2)
       self.nearbyLines_dialog.setWindowTitle(s)
       self.nearbyLines_dialog.setWindowModality(QtCore.Qt.NonModal)
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
 
       self.cont_pars_dialog = QtGui.QDialog()
       self.cont_pars_dialog.resize(800,460)
       self.table = QtGui.QTableWidget()   
#       self.table.keyPressEvent = test
#       self.table.itemChanged(self.changeIcon)
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
       self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Save|QtGui.QDialogButtonBox.Close)
       self.buttonBox.rejected.connect(self.cont_pars_dialog.close)
       self.buttonBox.accepted.connect(self.save_cont_pars)
       vbox = QtGui.QVBoxLayout()
       vbox.addWidget(self.table)
       vbox.addWidget(self.buttonBox)
       self.cont_pars_dialog.setLayout(vbox)
       self.cont_pars_dialog.setWindowTitle('Continuum parameters')
#       self.cont_pars_dialog.setWindowModality(QtCore.Qt.ApplicationModal)
       self.cont_pars_dialog.setWindowModality(QtCore.Qt.NonModal)
#       self.cont_pars_dialog.exec_()
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
                          
    def on_draw(self):
        """ Redraws        for w in [self.sp_min_box, self.ebv_box, self.resol_box]:
            hbox3.addWidget(w)
            hbox3.setAlignment(w, QtCore.Qt.AlignVCenter)
 the figure
        """
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
        self.sp.plot_ax1(self.axes)
        
        k = self.sp.get_conf('line_tick_ax')
        if self.show_line_ticks_action.isChecked() and ( k == 0 ):
            y1, y2 = self.get_line_tick_lim(self.sp.get_conf('line_tick_pos'))
            self.sp.plot_line_ticks(self.axes, y1, y2)

        if self.sp.get_conf('cont_plot', False):
          self.sp.plot_conts(self.axes)
        
        if self.residual_GroupBox.isChecked():
            self.axes3.cla()
            self.sp.plot_ax3(self.axes3)
            if self.show_line_ticks_action.isChecked() and ( k == 1 ):
                y1, y2 = self.get_line_tick_lim(self.sp.get_conf('line_tick_pos'))
                self.sp.plot_line_ticks(self.axes3, y1, y2)

        if self.show_line_ticks_action.isChecked() and ( k == 2 ):
            self.axes2.cla()
#            self.sp.plot_ax2(self.axes2)
            self.sp.plot_line_ticks(self.axes2, 0.2, 0.8)
        
        if self.residual_GroupBox.isChecked():
            self.axes3.set_xlabel(r'Wavelength ($\AA$)')
        elif self.show_line_ticks_action.isChecked() and self.sp.get_conf('plot_ax2'):
            self.axes2.set_xlabel(r'Wavelength ($\AA$)')
        else:
            self.axes.set_xlabel(r'Wavelength ($\AA$)')
        
        self.restore_axes()
        self.update_lim_boxes()
        self.canvas.draw()
        self.statusBar().showMessage('Redraw is finished.', 4000) 
        log_.debug('Exit on_drawn', calling=self.calling)
        
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

    def differentiate_lines_clicked(self):
        if self.differentiate_lines_action.isChecked():
          self.sp.set_conf('differentiate_lines', True) 
        else:  
          self.sp.set_conf('differentiate_lines', False) 
        self.make_axes()

    def editing_lines_clicked(self):
        if self.editing_lines_action.isChecked():
          self.sp.set_conf('allow_editing_lines', True) 
        else:  
          self.sp.set_conf('allow_editing_lines', False) 

    def update_lines_clicked(self):
        if self.update_lines_action.isChecked():
          self.sp.set_conf('update_after_editing_lines', True) 
        else:  
          self.sp.set_conf('update_after_editing_lines', False) 

    def cycle_forwards_ions(self):
        j = self.sp.get_conf('index_of_current_ion')
        selected_ions = self.sp.get_conf('selected_ions')   
        if j in range(-1, len(selected_ions)-1):
          j += 1
        else:
          j = -1
        self.sp.set_conf('index_of_current_ion', j)
        self.make_axes()

    def cycle_backwards_ions(self):
        j = self.sp.get_conf('index_of_current_ion')
        selected_ions = self.sp.get_conf('selected_ions')   
        if j in range(0, len(selected_ions)):
          j -= 1
        else:
          j = len(selected_ions)-1
        self.sp.set_conf('index_of_current_ion', j)
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
        self.fig.subplots_adjust(hspace=0.0, bottom=0.11, right=0.98, top=0.99, left=0.1)
        if self.call_on_draw: 
            log_.debug('Calling on_draw from make_axes', calling=self.calling)
            self.do_save = False
            self.on_draw()
            self.do_save = True
        
        log_.debug('Exit make_axes', calling=self.calling)
        
    def make_axes_original(self):
        
        log_.debug('Entering make_axes', calling=self.calling)
        if self.call_on_draw: 
            self.save_axes()
        self.fig.clf()

        n_subplots = 1
        i_ax2 = 2
        i_ax3 = 2
        if self.show_line_ticks_action.isChecked():
            n_subplots += 1
            i_ax3 += 1
        if self.residual_GroupBox.isChecked():
            n_subplots += 1
        if self.axes is not None:
            del(self.axes)
        self.axes = self.fig.add_subplot(n_subplots, 1, 1)
        self.sp.ax1 = self.axes
        if self.show_line_ticks_action.isChecked():
            if self.axes2 is not None:
                del(self.axes2)
            self.axes2 = self.fig.add_subplot(n_subplots, 1, i_ax2, sharex=self.axes)
            self.sp.ax2 = self.axes2
            self.axes.get_xaxis().set_visible(False)
        else:
            self.axes2 = None
            self.sp.ax2 = None
        if self.residual_GroupBox.isChecked():
            if self.axes3 is not None:
                del(self.axes3)
            self.axes3 = self.fig.add_subplot(n_subplots, 1, i_ax3, sharex=self.axes)
            self.sp.ax3 = self.axes3
            if self.sp.get_conf('plot_ax2'):
                self.axes2.get_xaxis().set_visible(False)
            self.axes.get_xaxis().set_visible(False)
        else:
            self.axes3 = None
            self.sp.ax3 = self.axes3
        self.fig.subplots_adjust(hspace=0.0, bottom=0.11, right=0.97, top=0.97, left=0.1)
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
            if self.sp.sp_synth_lr is None:
                self.y1_plot_lims = (np.min(self.sp.f), np.max(self.sp.f))
            else:
                self.y1_plot_lims = (np.min(self.sp.sp_synth_lr[mask]), np.max(self.sp.sp_synth_lr[mask]))      
        
        self.y2_plot_lims = self.sp.get_conf('y2_plot_lims')
        if self.y2_plot_lims is None:
            self.y2_plot_lims = (-0.5, 1.5)
        
        self.y3_plot_lims = self.sp.get_conf('y3_plot_lims')
        if self.y3_plot_lims is None:
            mask = (self.sp.w_ori > self.x_plot_lims[0]) & (self.sp.w_ori < self.x_plot_lims[1])
            if self.sp.cont is None:
                self.y3_plot_lims = (-1,1)
            else:
                self.y3_plot_lims = (np.min((self.sp.f - self.sp.cont)[mask]), np.max((self.sp.f - self.sp.cont)[mask]))
        log_.debug('Axes initialized. IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)
        self.print_axes()

    def save_axes(self):
        if self.axes is not None:
            self.x_plot_lims = self.axes.get_xlim()
            self.y1_plot_lims = self.axes.get_ylim()
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
        
        log_.debug('Axes restored. IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)
        self.print_axes()
        
    def print_axes(self):
        log_.debug('{} {} {} {}'.format(self.x_plot_lims, self.y1_plot_lims, self.y2_plot_lims, self.y3_plot_lims), calling=self.calling)
        log_.debug('Axes IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)
        log_.debug(' IDs {} {} {}'.format(id(self.axes), id(self.axes2), id(self.axes3)), calling=self.calling)

    
    def select_init(self, init_file_name=None):
        file_choices = "Python initialization files (*init.py) (*init.py);;Python files (*.py) (*.py);;All files (*) (*)"
        title = 'Open pySSN initialization file'
        if init_file_name is None:
            self.init_file_name = str(QtGui.QFileDialog.getOpenFileName(self, title, '', file_choices))
        else:
            self.init_file_name = init_file_name
        while self.init_file_name and not os.path.isfile(self.init_file_name):
            msg = "Initialization file '{}' not found!".format(self.init_file_name)
            QtGui.QMessageBox.critical(self, 'pySSN', msg, QtGui.QMessageBox.Ok )  
            self.init_file_name = str(QtGui.QFileDialog.getOpenFileName(self, title, '', file_choices))
        
        if self.init_file_name:
            self.statusBar().showMessage('Running synthesis ...') 
            QtGui.QApplication.processEvents() 
            self.start_spectrum()
            self.do_save = False
            self.on_draw()
            self.do_save = True

            #self.lineIDs_GroupBox.setChecked(self.sp.get_conf('plot_ax2', True))
            self.lineIDs_GroupBox.setChecked(False)
            self.residual_GroupBox.setChecked(self.sp.get_conf('plot_ax3', True))

            self.restore_axes()
        else:
            if self.sp is None:
                raise ValueError('An initialization filename must be given')
            else:
                log_.warn('A filename must be given', calling=self.calling)
                return
        
    def start_spectrum(self):
        init_file = self.init_file_name.split('/')[-1]
        dir = self.init_file_name.split(init_file)[0]
        if dir == '':
            dir = './'
        self.directory = dir
        self.sp = spectrum(config_file=self.init_file_name)
        if self.sp.phyat_file == 'NO_phyat.dat':
            self.status_text.setText('pySSN, v {}. init file: {}, No synthesis'.format(__version__, 
                                                                                        self.sp.config_file.split('/')[-1]))
        else:
            self.status_text.setText('pySSN, v {}. init file: {}, at. data: {}, model: {}, cosmetic: {}'.format(__version__, 
                                                                                                      self.sp.config_file.split('/')[-1], 
                                                                                                      self.sp.phyat_file.split('/')[-1],
                                                                                                      self.sp.get_conf('fic_modele').split('/')[-1],
                                                                                                      self.sp.get_conf('fic_cosmetik').split('/')[-1]))
        self.axes = None
        self.sp.ax2_fontsize = 6
        self.sp_norm_box.setText('{}'.format(self.sp.get_conf('sp_norm')))   
        self.obj_velo_box.setText('{}'.format(self.sp.get_conf('obj_velo')))   
        self.ebv_box.setText('{}'.format(self.sp.get_conf('e_bv')))   
        self.resol_box.setText('{}'.format(self.sp.get_conf('resol')))   
        self.cut2_box.setText('{}'.format(self.sp.get_conf('cut_plot2')))
        self.magenta_box.setText('{}'.format(self.sp.plot_magenta))
        self.magenta_label_box.setText('{}'.format(self.sp.label_magenta))
        self.cyan_box.setText('{}'.format(self.sp.plot_cyan))
        self.cyan_label_box.setText('{}'.format(self.sp.label_cyan))
        self.sp_min_box.setText('{}'.format(self.sp.get_conf('limit_sp')[0]))
        self.sp_max_box.setText('{}'.format(self.sp.get_conf('limit_sp')[1]))
        self.verbosity_ag.actions()[self.sp.get_conf('log_level', 0)].setChecked(True)
        self.line_tick_ax_ag.actions()[self.sp.get_conf('line_tick_ax', 0)].setChecked(True)
        self.line_tick_pos_ag.actions()[self.sp.get_conf('line_tick_pos', 0)].setChecked(True)
        self.residual_GroupBox.setChecked(self.sp.get_conf('plot_ax3', True))
        self.selected_ions_action.setChecked(self.sp.get_conf('show_selected_ions_only', False))
        self.ion_cb.setChecked(self.sp.get_conf('show_selected_ions_only', False))
        self.selected_intensities_action.setChecked(self.sp.get_conf('show_selected_intensities_only', False))
        self.cut_cb.setChecked(self.sp.get_conf('show_selected_intensities_only', False))
        self.differentiate_lines_action.setChecked(self.sp.get_conf('differentiate_lines', False))
        self.editing_lines_action.setChecked(self.sp.get_conf('allow_editing_lines', False))
        self.update_lines_action.setChecked(self.sp.get_conf('update_after_editing_lines', False))
        self.plot_cont_action.setChecked(self.sp.get_conf('cont_plot', False))
        self.show_line_ticks_action.setChecked(self.sp.get_conf('show_line_ticks', False))
        self.plot_lines_action.setChecked(self.sp.get_conf('plot_lines_of_selected_ions', False))
        self.lineIDs_GroupBox.setChecked(self.sp.get_conf('show_line_ticks', False) or self.sp.get_conf('plot_lines_of_selected_ions', False))
        self.set_save_plot_action_tip()
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

        self.line_sort_ag.actions()[self.sp.get_conf('line_saved_ordered_by', 0)].setChecked(True)
        self.show_header_action.setChecked(self.sp.get_conf('line_saved_header', False))
        self.get_line_fields_to_print()

        self.readOnlyCells_bg_color = QtGui.QColor('white')
        self.editableCells_bg_color = QtGui.QColor('lightgreen')

        self.line_info_dialog_width = 800
        self.line_info_dialog_height = 470
        #???

        sG = QtGui.QApplication.desktop().screenGeometry()
        self.line_info_dialog_x = sG.width()-self.line_info_dialog_width
        self.line_info_dialog_y = 0
    
    def sp_norm(self):
        if self.sp is None:
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
        old_obj_velo = self.sp.get_conf('obj_velo')
        new_obj_velo = np.float(self.obj_velo_box.text())
        if old_obj_velo == new_obj_velo:
          return
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
        self.statusBar().showMessage('Rerunning synthesis ...') 
        QtGui.QApplication.processEvents() 
        self.set_limit_sp()
        self.sp.set_conf('resol', np.int(self.resol_box.text()))
        self.sp.set_conf('obj_velo', np.float(self.obj_velo_box.text()))
        self.sp.set_conf('sp_norm', np.float(self.sp_norm_box.text()))
        self.sp.set_conf('e_bv', np.float(self.ebv_box.text()))
        self.sp.init_obs()
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run()
        self.on_draw()

    def adjust(self):
        if self.sp is None:
            return
        self.statusBar().showMessage('Running update ...')
        QtGui.QApplication.processEvents() 
        self.sp_norm()
        self.obj_velo()
        self.ebv()
        N_diff = self.sp.adjust()
        if N_diff > 0:
            self.on_draw()
        self.statusBar().showMessage('Update finished.', 4000) 

    def apply_post_proc(self):
        if self.post_proc_file is None or self.ask_postprocfile_action.isChecked():
            file_choices = "Python files (*.py) (*.py);;All files (*) (*)"
            if self.post_proc_file is None:
              path = ''
            else:
              path = self.post_proc_file
            path = unicode(QtGui.QFileDialog.getOpenFileName(self, 'Open file', path, file_choices))
            if path:
              self.post_proc_file = path
            else:
              return
        try:
            user_module = {}
            execfile(self.post_proc_file, user_module)
            self.post_proc = user_module['post_proc']
            log_.message('function post_proc read from {}'.format(self.post_proc_file))
        except:
            self.post_proc = None
            log_.warn('function post_proc NOT read from {}'.format(self.post_proc_file), 
                          calling = self.calling)
        if self.post_proc is not None:
            try:
                self.post_proc(self.fig)
            except:
                log_.warn('Error in {}'.format(self.post_proc_file), 
                          calling = self.calling)

    def update_profile(self):
        if self.sp is None:
            return
        self.sp.run(do_synth = True, do_read_liste = False, do_profiles=True)
        self.on_draw()

    def cut2(self):
        if self.sp is None:
            return
        self.sp.set_conf('cut_plot2', np.float(self.cut2_box.text()))
        self.on_draw()

    def get_ion_str(self,s):
        s = s.strip()
        s = s.replace(' ', '_')
        if s.isdigit():
          self.sp.line_info(int(s), sort='i_rel')
          s = self.sp.get_ion_from_code(s)
          s = s.strip()
        return s

    def set_ion(self):
        if self.sp is None:
            return
        sList = []
        s = self.ion_box.text()
        k = s.indexOf(',')
        while k > 0:
          s0 = self.get_ion_str(str(s[:k]))
          sList.append(s0)
          s = s[k+1:]
          k = s.indexOf(',')
        s0 = self.get_ion_str(str(s))
        sList.append(s0)
        s = ''
        for s0 in sList:
          s = s + s0 + ', '
        s = s[:-2]
        for item in sList:
          if '_' not in item:
            ion_list = self.sp.get_ions_from_element(item)
            if len(ion_list) > 0:
              k = sList.index(item)
              for ion in ion_list[::-1]:
                sList.insert(k, ion)
            sList.remove(item)
        self.sp.set_conf('selected_ions', sList)
        self.ion_box.setText(s)
        
    def draw_ion(self):
        self.set_ion()
        self.sp.set_conf('index_of_current_ion', -1)
        self.on_draw()
        
    def line_info(self):
        if self.sp is None:
            return
        new_ref = np.int(self.line_info_box.text())
        self.line_info_ref = new_ref
        if self.sp.get_conf('show_dialogs', True):
          self.show_line_info_dialog()
        else:
          self.sp.line_info(new_ref, sort='i_rel')
        
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

    def set_plot_ax2(self):
        k = self.line_tick_ax_list.index(self.line_tick_ax_ag.checkedAction().text())
        self.sp.set_conf('line_tick_ax',k)
        k = self.line_tick_pos_list.index(self.line_tick_pos_ag.checkedAction().text())
        self.sp.set_conf('line_tick_pos',k)
        if self.show_line_ticks_action.isChecked():
          self.make_axes()
        
    def verbosity(self):
        verbosity = self.verbosity_list.index(self.verbosity_ag.checkedAction().text())
        log_.debug('Change verbosity from {} to {}'.format(log_.level, verbosity), calling=self.calling)
        log_.level = verbosity
        
    def update_lim_boxes(self):
        xformat = '{:.0f}'
        yformat = '{:.0f}'
        self.xlim_min_box.setText(xformat.format(self.x_plot_lims[0]))
        self.xlim_max_box.setText(xformat.format(self.x_plot_lims[1]))
        self.y1lim_min_box.setText(yformat.format(self.y1_plot_lims[0]))
        self.y1lim_max_box.setText(yformat.format(self.y1_plot_lims[1]))
        self.y3lim_min_box.setText(yformat.format(self.y3_plot_lims[0]))
        self.y3lim_max_box.setText(yformat.format(self.y3_plot_lims[1]))

    def save_from_lim_boxes(self):
        self.x_plot_lims = (float(self.xlim_min_box.text()), float(self.xlim_max_box.text()))
        self.y1_plot_lims = (float(self.y1lim_min_box.text()), float(self.y1lim_max_box.text()))
        self.y3_plot_lims = (float(self.y3lim_min_box.text()), float(self.y3lim_max_box.text()))
        self.restore_axes()

    def save_from_lim_boxes_and_draw(self):
        self.save_from_lim_boxes()
        self.on_draw()

    def set_limit_sp(self):
        limit_sp = (np.float(self.sp_min_box.text()), np.float(self.sp_max_box.text()))
        self.sp.set_conf('limit_sp', limit_sp)
    
    def set_limit_sp_and_run(self):
        old_limit_sp = self.sp.get_conf('limit_sp')
        new_limit_sp = limit_sp = (np.float(self.sp_min_box.text()), np.float(self.sp_max_box.text()))
        if old_limit_sp == new_limit_sp:
          return
        self.sp.set_conf('limit_sp', new_limit_sp)
        log_.message('Changing limit_sp. Old: {}, New: {}'.format(old_limit_sp, new_limit_sp), calling=self.calling)
        self.statusBar().showMessage('Changing the synthesis wavelength limits ...') 
        QtGui.QApplication.processEvents() 
        self.sp.init_obs()
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run(do_synth = True, do_read_liste = True, do_profiles=False)
        self.on_draw()

    def resol(self):
        if self.sp is None:
            return
        old_resol = self.sp.get_conf('resol')
        new_resol = np.int(self.resol_box.text())
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
        if not self.axes_fixed:
            self.save_axes()
            self.update_lim_boxes()
        
    def fix_axes(self):
        if self.fix_axes_cb.isChecked():
            self.axes_fixed = True
        else:
            self.axes_fixed = False
            
    def get_line_fields_to_print(self):
        field_list = self.sp.get_conf('line_field_print')
        for i in range(0,len(self.line_field_menu.actions())):
            if self.line_print_dic.keys()[i] in field_list:
              self.line_field_menu.actions()[i].setChecked(True)
            else:
              self.line_field_menu.actions()[i].setChecked(False)
    def set_show_header(self):
        if self.show_header_action.isChecked():
          self.sp.set_conf('line_saved_header', True)
        else:
          self.sp.set_conf('line_saved_header', False)
      
    def set_line_fields_to_print(self):
        s = []
        for i in range(0,len(self.line_field_menu.actions())):
            if self.line_field_menu.actions()[i].isChecked():
              s.append( self.line_print_dic.keys()[i])
        self.sp.set_conf('line_field_print', s)
        
    def save_lines(self):
        self.sp.save_lines()
        path = self.sp.get_conf('line_saved_filename')
        self.statusBar().showMessage('Saved to %s' % path, 4000)
    
    def save_lines_as(self):
        file_choices = "Text files (*.txt *.dat) (*.txt *.dat);;Tex files (*.tex) (*.tex);;CSV files (*.csv) (*.csv);;All Files (*) (*)"
        filename = self.sp.get_conf('line_saved_filename')
        path = unicode(QtGui.QFileDialog.getSaveFileName(self, 'Save lines to file', filename, file_choices))
        if path:
            self.sp.set_conf('line_saved_filename', path)
            self.sp.save_lines()
            self.statusBar().showMessage('Saved to %s' % path, 4000)

    def line_sort(self):
        k = self.line_sort_list.index(self.line_sort_ag.checkedAction().text())
        self.sp.set_conf('line_saved_ordered_by',k)

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
    form.show()
    app.exec_()
    
if __name__ == "__main__":
    main()
