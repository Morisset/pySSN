"""
This demo demonstrates how to embed a matplotlib (mpl) plot 
into a PyQt4 GUI application, including:

* Using the navigation toolbar
* Adding data to the plot
* Dynamically modifying the plot's properties
* Processing mpl events
* Saving the plot to a file from a menu

The main goal is to serve as a basis for developing rich PyQt GUI
applications featuring mpl plots (using the mpl OO API).

Eli Bendersky (eliben@gmail.com)
License: this code is in the public domain
Last modified: 19.01.2009
"""
import sys, os, random
import argparse

from PyQt4 import QtCore, QtGui

import matplotlib

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg
from matplotlib.figure import Figure
import numpy as np
import pyssn


#ToDo : 

pyssn.log_.level = 3

class NavigationToolbar( NavigationToolbar2QTAgg ):

    curs = QtCore.pyqtSignal(bool)

    def __init__(self, canvas, parent ):
        
        NavigationToolbar2QTAgg.__init__(self,canvas,parent)
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
    
    def __init__(self, parent=None):
        self.calling = 'AppForm'
        QtGui.QMainWindow.__init__(self, parent)
        self.setWindowTitle('pySSN')
        self.sp = None
        self.axes = None
        self.fig = None
        self.init_file_name = None
        self.call_on_draw = True
        self.cursor_on = False
        self.line_info_ref = 0
        
        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()

        self.select_init()
        
    def save_plot(self):
        file_choices = "PDF (*.pdf)|*.pdf"
        
        path = unicode(QtGui.QFileDialog.getSaveFileName(self, 
                        'Save file', '', 
                        file_choices))
        if path:
            self.canvas.print_figure(path, dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' % path, 2000)
    
    def on_about(self):
        msg = """ A demo of using PyQt with matplotlib:
        
         * Use the matplotlib navigation bar
         * Add values to the text box and press Enter (or click "Draw")
         * Show or hide the grid
         * Drag the slider to modify the width of the bars
         * Save the plot to a file using the File menu
         * Click on a bar to receive an informative message
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
        #self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)

        self.canvas.mpl_connect('button_press_event', self.on_click)
                
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        self.mpl_toolbar.curs.connect(self.set_cursor)   
        # Other GUI controls
        # 

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

        self.fix_axes_cb = QtGui.QCheckBox("Fix axes")
        self.fix_axes_cb.setChecked(False)
        self.connect(self.fix_axes_cb, QtCore.SIGNAL('stateChanged(int)'), self.fix_axes)

        self.adjust_button = QtGui.QPushButton("Update lines")
        self.adjust_button.setChecked(False)
        self.connect(self.adjust_button, QtCore.SIGNAL('clicked()'), self.adjust)

        self.update_profile_button = QtGui.QPushButton("Update profiles")
        self.update_profile_button.setChecked(False)
        self.connect(self.update_profile_button, QtCore.SIGNAL('clicked()'), self.update_profile)

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
        #
        # Layout with box sizers
        # 
        hbox = QtGui.QHBoxLayout()
        hbox2 = QtGui.QHBoxLayout()
        hbox3 = QtGui.QHBoxLayout()
         
        for w in [self.select_init_button, self.draw_button, self.pl_ax2_cb, self.pl_ax3_cb, self.fix_axes_cb, self.adjust_button]:
            hbox.addWidget(w)
            hbox.setAlignment(w, QtCore.Qt.AlignVCenter)

        for w in [self.sp_norm_box, self.obj_velo_box, self.ebv_box, self.resol_box, self.cut2_box]:
            hbox2.addWidget(w)
            hbox2.setAlignment(w, QtCore.Qt.AlignVCenter)
        
        for w in [self.line_info_box, self.magenta_box, self.cyan_box]:
            hbox3.addWidget(w)
            hbox3.setAlignment(w, QtCore.Qt.AlignVCenter)

        vbox = QtGui.QVBoxLayout()
        
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        vbox.addLayout(hbox)
        vbox.addLayout(hbox2)
        vbox.addLayout(hbox3)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QtGui.QLabel("pySSN, v{}".format(pyssn.__version__))
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
        pyssn.log_.message('Entering on_drawn', calling=self.calling)
        if self.sp is None:
            pyssn.log_.message('Np sp in on_drawn', calling=self.calling)
            return
        
        if self.axes is None:
            pyssn.log_.message('Calling make_axes from on_draw', calling=self.calling)
            self.call_on_draw=False
            self.make_axes()
            self.call_on_draw=True

        self.axes.cla()
        self.sp.plot_ax1(self.axes)
        
        if self.pl_ax2_cb.isChecked():
            self.axes2.cla()
            self.sp.plot_ax2(self.axes2)
        
        if self.pl_ax3_cb.isChecked():
            self.axes3.cla()
            self.sp.plot_ax3(self.axes3)
                 
        self.canvas.draw()
        pyssn.log_.message('Exit on_drawn', calling=self.calling)
        
    def make_axes(self):
        
        pyssn.log_.message('Entering make_axes', calling=self.calling)
        pyssn.log_.message('make_axes call_on_draw value:{}'.format(self.call_on_draw), calling=self.calling)
        self.fig.clf()

        n_subplots = 1
        i_ax2 = 2
        i_ax3 = 2
        if self.pl_ax2_cb.isChecked():
            n_subplots += 1
            i_ax3 += 1
        if self.pl_ax3_cb.isChecked():
            n_subplots += 1
        self.axes = self.fig.add_subplot(n_subplots, 1, 1)
        if self.pl_ax2_cb.isChecked():
            self.axes2 = self.fig.add_subplot(n_subplots, 1, i_ax2, sharex=self.axes)
            self.axes.get_xaxis().set_visible(False)
        if self.pl_ax3_cb.isChecked():
            self.axes3 = self.fig.add_subplot(n_subplots, 1, i_ax3, sharex=self.axes)
            if self.pl_ax2_cb.isChecked():
                self.axes2.get_xaxis().set_visible(False)
        if self.sp.fixed_axes:
            self.sp.restore_axes()
        self.fig.subplots_adjust(hspace=0.0)
        if self.call_on_draw: 
            pyssn.log_.message('Calling on_draw from make_axes', calling=self.calling)
            self.on_draw()
        pyssn.log_.message('Exit make_axes', calling=self.calling)
    
    def select_init(self, init_file_name=None):
        if init_file_name is None:
            self.init_file_name = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '', '*init.py')
        else:
            self.init_file_name = init_file_name
        self.start_spectrum()
        self.on_draw()
    
    def start_spectrum(self):
        self.sp = pyssn.spectrum(config_file=unicode(self.init_file_name))
        self.status_text.setText('pySSN, v {}. init file: {}, at. data: {}, model: {}, cosmetic: {}'.format(pyssn.__version__, 
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
    
    def fix_axes(self):
        if self.sp is None:
            return
        if self.fix_axes_cb.isChecked():
            self.sp.fix_axes()
        else:
            self.sp.unfix_axes()
            
    def sp_norm(self):
        
        if self.sp is None:
            return
        old_sp_norm = self.sp.get_conf('sp_norm')
        new_sp_norm = np.float(self.sp_norm_box.text())
        pyssn.log_.message('Entering sp_norm. Old: {}, New: {}'.format(old_sp_norm, new_sp_norm), calling=self.calling)
        if np.abs(new_sp_norm - old_sp_norm)/old_sp_norm > 1e-6:
            self.sp.renorm(new_sp_norm)
            self.on_draw()

    def obj_velo(self):
        
        if self.sp is None:
            return
        old_obj_velo = self.sp.get_conf('obj_velo')
        new_obj_velo = np.float(self.obj_velo_box.text())
        pyssn.log_.message('Entering obj_velo. Old: {}, New: {}'.format(old_obj_velo, new_obj_velo), calling=self.calling)
        if np.abs(new_obj_velo - old_obj_velo)/old_obj_velo > 1e-6:
            self.sp.init_obs(obj_velo=new_obj_velo)
            self.sp.make_continuum()
            self.sp.run()
            self.on_draw()

    def ebv(self):
        if self.sp is None:
            return
        old_ebv = self.sp.get_conf('e_bv')
        new_ebv = np.float(self.ebv_box.text())
        pyssn.log_.message('Entering ebv. Old: {}, New: {}'.format(old_ebv, new_ebv), calling=self.calling)
        if np.abs(new_ebv - old_ebv)/old_ebv > 1e-6:
            self.sp.set_conf('e_bv', new_ebv)
            self.sp.init_red_corr()
            self.sp.make_continuum()
            self.sp.run(do_synth = True, do_read_liste = False, do_profiles=False)
            self.on_draw()
            
    def adjust(self):
        if self.sp is None:
            return
        self.sp.adjust()
        self.on_draw()

    def resol(self):
        if self.sp is None:
            return
        old_resol = self.sp.get_conf('resol')
        new_resol = np.int(self.resol_box.text())
        pyssn.log_.message('Entering resol. Old: {}, New: {}'.format(old_resol, new_resol), calling=self.calling)
        if np.abs(new_resol - old_resol) > 0.5:
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
        if new_ref != self.line_info_ref:
            self.line_info_ref = new_ref
            self.sp.line_info(new_ref)
        
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
        
        
        
def main():
    app = QtGui.QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()
    
if __name__ == "__main__":
    main()


