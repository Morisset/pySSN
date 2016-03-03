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
    
    def __init__(self, parent=None, init_filename=None, post_proc_file=None):
        self.calling = 'pySSN GUI'
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
        self.do_save = True
        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        self.select_init(init_filename)
        self.post_proc_file = post_proc_file
        
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
        self.main_frame = QtGui.QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        #
        self.dpi = 100
        self.fig = Figure((20.0, 15.0), dpi=self.dpi)
        log_.debug('creating figure {}'.format(id(self.fig)), calling=self.calling)
        
        
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
    
        self.canvas.mpl_connect('button_press_event', self.on_click)
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        self.mpl_toolbar.curs.connect(self.set_cursor)   
        # Other GUI controls
        # 
        
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
        self.connect(self.xlim_max_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        
        self.y3lim_min_box = QtGui.QLineEdit()
        self.y3lim_min_box.setMinimumWidth(50)
        self.connect(self.y3lim_min_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        
        self.y3lim_max_box = QtGui.QLineEdit()
        self.y3lim_max_box.setMinimumWidth(50)
        self.connect(self.y3lim_max_box, QtCore.SIGNAL('returnPressed()'), self.save_from_lim_boxes)
        
        self.select_init_button = QtGui.QPushButton("&Select init file")
        self.connect(self.select_init_button, QtCore.SIGNAL('clicked()'), self.select_init)
        
        self.draw_button = QtGui.QPushButton("&Draw")
        self.connect(self.draw_button, QtCore.SIGNAL('clicked()'), self.on_draw)
        
        self.pl_ax2_cb = QtGui.QCheckBox("Lines ID")
        self.pl_ax2_cb.setChecked(True)
        self.connect(self.pl_ax2_cb, QtCore.SIGNAL('stateChanged(int)'), self.make_axes)
        
        self.pl_ax3_cb = QtGui.QCheckBox("Residus")
        self.pl_ax3_cb.setChecked(True)
        self.connect(self.pl_ax3_cb, QtCore.SIGNAL('stateChanged(int)'), self.make_axes)

        self.adjust_button = QtGui.QPushButton("Update lines")
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

        self.cyan_box = QtGui.QLineEdit()
        self.cyan_box.setMinimumWidth(50)
        self.connect(self.cyan_box, QtCore.SIGNAL('returnPressed()'), self.cyan_line)
        
        self.verbosity_cb = QtGui.QComboBox()
        self.verbosity_cb.addItem('0: nothing')
        self.verbosity_cb.addItem('1: Errors')
        self.verbosity_cb.addItem('2: Errs + Warnings')
        self.verbosity_cb.addItem('3: Errs + Warns + Messages')
        self.verbosity_cb.addItem('4: All + Debug messages')
        self.verbosity_cb.setCurrentIndex(log_.level)
        self.connect(self.verbosity_cb, QtCore.SIGNAL('activated(QString)'), self.verbosity)
        
        #
        # Layout with box sizers
        # 
        hbox0t = QtGui.QHBoxLayout()
        hbox0 = QtGui.QHBoxLayout()
        hbox1 = QtGui.QHBoxLayout()
        hbox2 = QtGui.QHBoxLayout()
        hbox3 = QtGui.QHBoxLayout()
        hbox4 = QtGui.QHBoxLayout()
        hbox5 = QtGui.QHBoxLayout()
         
        for l in ['xmin', 'xmax', 'y1min', 'y1max', 'y3min', 'y3max']:
            w = QtGui.QLabel(l)
            hbox0t.addWidget(w)
            hbox0t.setAlignment(w, QtCore.Qt.AlignVCenter)

        for w in [self.xlim_min_box, self.xlim_max_box, self.y1lim_min_box, self.y1lim_max_box, self.y3lim_min_box, self.y3lim_max_box]:
            hbox0.addWidget(w)
            hbox0.setAlignment(w, QtCore.Qt.AlignVCenter)

        for w in [self.select_init_button, self.draw_button, self.pl_ax2_cb, self.pl_ax3_cb, 
                  self.adjust_button, self.post_proc_button]:
            hbox1.addWidget(w)
            hbox1.setAlignment(w, QtCore.Qt.AlignVCenter)

        for l in ['sp min', 'sp max', 'sp norm', 'obj velo', 'E_B-V', 'resol', 'cut2']:
            w = QtGui.QLabel(l)
            hbox2.addWidget(w)
            hbox2.setAlignment(w, QtCore.Qt.AlignVCenter)

        for w in [self.sp_min_box, self.sp_max_box, self.sp_norm_box, self.obj_velo_box, self.ebv_box, self.resol_box, self.cut2_box]:
            hbox3.addWidget(w)
            hbox3.setAlignment(w, QtCore.Qt.AlignVCenter)
          
        for l in ['line info', 'magenta', 'cyan', 'verbosity']:
            w = QtGui.QLabel(l)
            hbox4.addWidget(w)
            hbox4.setAlignment(w, QtCore.Qt.AlignVCenter)

        for w in [self.line_info_box, self.magenta_box, self.cyan_box, self.verbosity_cb]:
            hbox5.addWidget(w)
            hbox5.setAlignment(w, QtCore.Qt.AlignVCenter)

        vbox = QtGui.QVBoxLayout()
        
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(hbox0t)
        vbox.addLayout(hbox0)
        vbox.addLayout(hbox1)
        vbox.addLayout(hbox2)
        vbox.addLayout(hbox3)
        vbox.addLayout(hbox4)
        vbox.addLayout(hbox5)
        
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
        """ Redraws the figure
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
        
        if self.pl_ax2_cb.isChecked():
            self.axes2.cla()
            self.sp.plot_ax2(self.axes2)
        
        if self.pl_ax3_cb.isChecked():
            self.axes3.cla()
            self.sp.plot_ax3(self.axes3)
        
        if self.pl_ax3_cb.isChecked():
            self.axes3.set_xlabel(r'Wavelength ($\AA$)')
        elif self.pl_ax2_cb.isChecked():
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

        n_subplots = 1
        i_ax2 = 2
        i_ax3 = 2
        if self.pl_ax2_cb.isChecked():
            n_subplots += 1
            i_ax3 += 1
        if self.pl_ax3_cb.isChecked():
            n_subplots += 1
        if self.axes is not None:
            del(self.axes)
        self.axes = self.fig.add_subplot(n_subplots, 1, 1)
        self.sp.ax1 = self.axes
        if self.pl_ax2_cb.isChecked():
            if self.axes2 is not None:
                del(self.axes2)
            self.axes2 = self.fig.add_subplot(n_subplots, 1, i_ax2, sharex=self.axes)
            self.sp.ax2 = self.axes2
            self.axes.get_xaxis().set_visible(False)
        else:
            self.axes2 = None
            self.sp.ax2 = None
        if self.pl_ax3_cb.isChecked():
            if self.axes3 is not None:
                del(self.axes3)
            self.axes3 = self.fig.add_subplot(n_subplots, 1, i_ax3, sharex=self.axes)
            self.sp.ax3 = self.axes3
            if self.pl_ax2_cb.isChecked():
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
            self.y1_plot_lims = (np.min(self.sp.sp_synth_lr[mask]), np.max(self.sp.sp_synth_lr[mask]))      
        
        self.y2_plot_lims = self.sp.get_conf('y2_plot_lims')
        if self.y2_plot_lims is None:
            self.y2_plot_lims = (-1.5, 1)
        
        self.y3_plot_lims = self.sp.get_conf('y3_plot_lims')
        if self.y3_plot_lims is None:
            mask = (self.sp.w_ori > self.x_plot_lims[0]) & (self.sp.w_ori < self.x_plot_lims[1])
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
            self.pl_ax2_cb.setChecked(self.sp.get_conf('plot_ax2', True))
            self.pl_ax3_cb.setChecked(self.sp.get_conf('plot_ax3', True))
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
        self.cut2_box.setText('{}'.format(self.sp.cut_plot2))
        self.magenta_box.setText('{}'.format(self.sp.plot_magenta))
        self.cyan_box.setText('{}'.format(self.sp.plot_cyan))
        self.sp_min_box.setText('{}'.format(self.sp.get_conf('limit_sp')[0]))
        self.sp_max_box.setText('{}'.format(self.sp.get_conf('limit_sp')[1]))
            
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
        self.sp.run(do_synth = True, do_read_liste = True, do_profiles=False)
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
        self.sp.run(do_synth = True, do_read_liste = False, do_profiles=False)
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
        self.sp.cut_plot2 = np.float(self.cut2_box.text())
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
        new_ref = np.int(self.magenta_box.text())
        if new_ref != self.sp.plot_magenta:
            self.sp.plot_magenta = new_ref
            self.on_draw()
        
    def cyan_line(self):
        if self.sp is None:
            return
        new_ref = np.int(self.cyan_box.text())
        if new_ref != self.sp.plot_cyan:
            self.sp.plot_cyan = new_ref
            self.on_draw()
    
    def verbosity(self):
        verbosity = self.verbosity_cb.currentIndex()
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
        
        
def main_loc(init_filename=None, post_proc_file=None):
    app = QtGui.QApplication(sys.argv)
    form = AppForm(init_filename=init_filename, post_proc_file=post_proc_file)
    form.show()
    app.exec_()
    return form.fig
        
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


