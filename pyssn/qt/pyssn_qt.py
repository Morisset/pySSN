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

    def clearCurs( self, checked ):
        if checked:
            self.ac.setChecked(False)

    def curs_toggle(self, checked):
        self.curs.emit(checked)
    
    def canvasPicked(self, event ):
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
        self.axes_fixed = False

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        self.select_init(init_filename)
        
    def save_plot(self):
        file_choices = "PDF (*.pdf)|*.pdf"
        
        path = unicode(QtGui.QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def on_about(self):
        msg = """ pySSN (Spectral Synthesis for Nebulae):        
        """
        QtGui.QMessageBox.about(self, "About the demo", msg.strip())
    
    def set_cursor(self, checked):
        self.cursor_on = checked
                
    def on_click(self, event):
        if self.cursor_on:
            self.sp._curs_onclick(event)

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
        self.connect(self.xlim_min_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        
        self.xlim_max_box = QtGui.QLineEdit()
        self.xlim_max_box.setMinimumWidth(50)
        self.connect(self.xlim_max_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        
        self.y1lim_min_box = QtGui.QLineEdit()
        self.y1lim_min_box.setMinimumWidth(50)
        self.connect(self.y1lim_min_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        
        self.y1lim_max_box = QtGui.QLineEdit()
        self.y1lim_max_box.setMinimumWidth(50)
        self.connect(self.y1lim_max_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        
        self.y3lim_min_box = QtGui.QLineEdit()
        self.y3lim_min_box.setMinimumWidth(50)
        self.connect(self.y3lim_min_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        s = 'Minimum ordinate value of the residual plot\n\n' \
            'Set with:  \n' \
            '    y3_plot_lims = (<ymin>, <ymax>)'
        self.y3lim_min_box.setToolTip( s )        
        
        self.y3lim_max_box = QtGui.QLineEdit()
        self.y3lim_max_box.setMinimumWidth(50)
        self.connect(self.y3lim_max_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        s = 'Maximum ordinate value of the residual plot\n\n' \
            'Set with:  \n' \
            '    y3_plot_lims = (<ymin>, <ymax>)'
        self.y3lim_max_box.setToolTip( s )        
        
        self.select_init_button = QtGui.QPushButton("Init file")
        self.connect(self.select_init_button, QtCore.SIGNAL('clicked()'), self.select_init)
        
        self.rerun_button = QtGui.QPushButton("&Rerun")
        self.connect(self.rerun_button, QtCore.SIGNAL('clicked()'), self.rerun)
        
        self.draw_button = QtGui.QPushButton("&Draw")
        self.connect(self.draw_button, QtCore.SIGNAL('clicked()'), self.on_draw)

        self.savelines_button = QtGui.QPushButton("&Save")
        self.savelines_button.setMinimumWidth(50)
        self.connect(self.savelines_button, QtCore.SIGNAL('clicked()'), self.save_lines)
        
        self.Command_GroupBox = QtGui.QGroupBox("Main commands")
        self.Command_GroupBox.setCheckable(False)
        
        self.ObsSpec_GroupBox = QtGui.QGroupBox("Parameters related to the observed spectrum")
        self.ObsSpec_GroupBox.setCheckable(False)

        self.SpecPlot_GroupBox = QtGui.QGroupBox("Plot of the spectra")
        self.SpecPlot_GroupBox.setCheckable(False)

        self.lineIDs_GroupBox = QtGui.QGroupBox("Show lines")
        self.lineIDs_GroupBox.setCheckable(True)
        self.lineIDs_GroupBox.setChecked(False)
        
#        self.lineIDs_GroupBox.setTristate(True)
#        self.lineIDs_GroupBox.setCheckState(QtCore.Qt.PartiallyChecked)
        
        self.connect(self.lineIDs_GroupBox, QtCore.SIGNAL('clicked()'), self.make_axes)
        self.lineIDs_GroupBox.setToolTip( 'Check to display the central positions of the spectral lines' )        

        self.residual_GroupBox = QtGui.QGroupBox("Plot of residuals")
        self.residual_GroupBox.setCheckable(True)
        self.residual_GroupBox.setChecked(True)
        self.connect(self.residual_GroupBox, QtCore.SIGNAL('clicked()'), self.make_axes)
        self.residual_GroupBox.setToolTip( 'Check to display the residual plot' )        

        self.adjust_button = QtGui.QPushButton("&Update lines")
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
        self.connect(self.sp_min_box, QtCore.SIGNAL('returnPressed()'), self.change_sp)
        
        self.sp_max_box = QtGui.QLineEdit()
        self.sp_max_box.setMinimumWidth(50)
        self.connect(self.sp_max_box, QtCore.SIGNAL('returnPressed()'), self.change_sp)
        
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
        
        self.line_info_box = QtGui.QLineEdit()
        self.line_info_box.setMinimumWidth(50)
        self.connect(self.line_info_box, QtCore.SIGNAL('returnPressed()'), self.line_info)

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

        s =  'Color excess E(B-V)\n\n' \
             'Set with: \n' \
             '    e_bv = <float>\n\n' \
             'Comment: \n' \
            u'    E(B-V) \u2248 C(H\u03B2) / 1.5'
        self.ebv_box.setToolTip( s )  
        
        s = 'Radial velocity in km/s\n\n' \
            'Set with: \n' \
            '    obj_velo = <float>'
        self.obj_velo_box.setToolTip( s ) 
        
        s = 'Normalization factor, ratio between the intensity and the \n' \
            u'observed flux of the reference line, 10\u2074/F(H\u03B2)\n\n' \
             'Set with: \n' \
             '    sp_norm = <float>'
        self.sp_norm_box.setToolTip( s )        

        s = 'Rebinning factor, the integer factor by which the number of points \n' \
            'of the original spectrum is multiplied in the rebinning process\n\n' \
            'Set with: \n' \
            '    resol = <integer>\n\n' \
            'Usage: \n' \
            '    Set to \'1\' if the resolution of the observed spectrum is large enough' 
              
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

        wList = [self.select_init_button,self.rerun_button,self.draw_button,self.adjust_button,self.post_proc_button,self.verbosity_button]
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
        LineIDLayout.addWidget(self.savelines_button,1,1)

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
        self.file_menu = self.menuBar().addMenu("&File")
        
        save_plot_action = self.create_action("&Save Plot",
                                              shortcut="Ctrl+S", 
                                              slot=self.save_plot, 
                                              tip="Save the plot")
        quit_action = self.create_action("&Quit", 
                                         slot=self.fileQuit, 
                                         shortcut="Ctrl+Q", 
                                         tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (quit_action, save_plot_action))
        
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

        if self.lineIDs_GroupBox.isChecked() and self.sp.get_conf('plot_ax2'):
            self.axes2.cla()
            self.sp.plot_ax2(self.axes2)
#            self.sp.plot_lines(self.axes3)
        
        if self.residual_GroupBox.isChecked():
            self.axes3.cla()
            self.sp.plot_ax3(self.axes3)
            if self.lineIDs_GroupBox.isChecked() and not self.sp.get_conf('plot_ax2'):
                self.sp.plot_lines(self.axes3)
        
        if self.residual_GroupBox.isChecked():
            self.axes3.set_xlabel(r'Wavelength ($\AA$)')
        elif self.lineIDs_GroupBox.isChecked() and self.sp.get_conf('plot_ax2'):
            self.axes2.set_xlabel(r'Wavelength ($\AA$)')
        else:
            self.axes.set_xlabel(r'Wavelength ($\AA$)')
        
        self.restore_axes()
        self.update_lim_boxes()
        self.canvas.draw()
        log_.debug('Exit on_drawn', calling=self.calling)
        
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
        ShowAx2 = self.lineIDs_GroupBox.isChecked() and self.sp.get_conf('plot_ax2')
#        if self.lineIDs_GroupBox.isChecked():
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
#        if self.lineIDs_GroupBox.isChecked():
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
        if self.lineIDs_GroupBox.isChecked():
            n_subplots += 1
            i_ax3 += 1
        if self.residual_GroupBox.isChecked():
            n_subplots += 1
        if self.axes is not None:
            del(self.axes)
        self.axes = self.fig.add_subplot(n_subplots, 1, 1)
        self.sp.ax1 = self.axes
        if self.lineIDs_GroupBox.isChecked():
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
        if init_file_name is None:
            self.init_file_name = str(QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', '*init.py'))
        else:
            self.init_file_name = init_file_name
        if self.init_file_name:
            self.start_spectrum()
            self.do_save = False
            self.on_draw()
            self.do_save = True
            self.lineIDs_GroupBox.setChecked(self.sp.get_conf('plot_ax2', True))
            self.residual_GroupBox.setChecked(self.sp.get_conf('plot_ax3', True))
            self.restore_axes()
        else:
            if self.sp is None:
                raise ValueError('A filename must be given')
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
        self.verbosity_ag.actions()[self.sp.get_conf('log_level')].setChecked(True)
            
    def rerun(self):
        self.sp.read_obs()
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run()
        self.on_draw()
    
    def sp_norm(self):
        
        if self.sp is None:
            return
        old_sp_norm = self.sp.get_conf('sp_norm')
        new_sp_norm = np.float(self.sp_norm_box.text())
        log_.message('Changing sp_norm. Old: {}, New: {}'.format(old_sp_norm, new_sp_norm), calling=self.calling)
        
        self.sp.renorm(new_sp_norm)
        self.on_draw()

    def obj_velo(self):
        
        if self.sp is None:
            return
        old_obj_velo = self.sp.get_conf('obj_velo')
        new_obj_velo = np.float(self.obj_velo_box.text())
        log_.message('Changing obj_velo. Old: {}, New: {}'.format(old_obj_velo, new_obj_velo), calling=self.calling)
        
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
        log_.message('Changing E B-V. Old: {}, New: {}'.format(old_ebv, new_ebv), calling=self.calling)
        
        self.sp.set_conf('e_bv', new_ebv)
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run(do_synth = self.sp.do_synth, do_read_liste = False, do_profiles=False)
        self.on_draw()
            
    def adjust(self):
        if self.sp is None:
            return
        N_diff = self.sp.adjust()
        if N_diff > 0:
            self.on_draw()

    def apply_post_proc(self):
        if self.post_proc_file is None:
            file_ = str(QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', '*.py'))
            self.post_proc_file = file_.split('/')[-1]
        try:
            user_module = {}
            execfile(os.path.abspath(self.directory)+'/'+self.post_proc_file, user_module)
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
                
    def resol(self):
        if self.sp is None:
            return
        old_resol = self.sp.get_conf('resol')
        new_resol = np.int(self.resol_box.text())
        log_.message('Changing resol. Old: {}, New: {}'.format(old_resol, new_resol), calling=self.calling)
        self.sp.set_conf('resol', new_resol)
        self.sp.init_obs()
        self.sp.init_red_corr()
        self.sp.make_continuum()
        self.sp.run(do_synth = True, do_read_liste = True, do_profiles=False)
        self.on_draw()

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

    def line_info(self):
        if self.sp is None:
            return
        new_ref = np.int(self.line_info_box.text())
        self.line_info_ref = new_ref
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
    
    def verbosity(self):
        verbosity = self.verbosity_list.index(self.verbosity_ag.checkedAction().text())
        log_.debug('Change verbosity from {} to {}'.format(log_.level, verbosity), calling=self.calling)
        log_.level = verbosity
        
    def update_lim_boxes(self):
        self.xlim_min_box.setText('{}'.format(self.x_plot_lims[0]))
        self.xlim_max_box.setText('{}'.format(self.x_plot_lims[1]))
        self.y1lim_min_box.setText('{}'.format(self.y1_plot_lims[0]))
        self.y1lim_max_box.setText('{}'.format(self.y1_plot_lims[1]))
        self.y3lim_min_box.setText('{}'.format(self.y3_plot_lims[0]))
        self.y3lim_max_box.setText('{}'.format(self.y3_plot_lims[1]))

    def save_from_lim_boxes(self):
        self.x_plot_lims = (float(self.xlim_min_box.text()), float(self.xlim_max_box.text()))
        self.y1_plot_lims = (float(self.y1lim_min_box.text()), float(self.y1lim_max_box.text()))
        self.y3_plot_lims = (float(self.y3lim_min_box.text()), float(self.y3lim_max_box.text()))
        self.restore_axes()
        self.on_draw()
        
    def change_sp(self):
        limit_sp = (float(self.sp_min_box.text()), float(self.sp_max_box.text()))
        self.sp.set_conf('limit_sp', limit_sp)
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
    
    def save_lines(self):
        self.sp.save_lines()
    
    
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
