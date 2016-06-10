'''
Created on 16/01/2014

@author: morisset
'''
import os
from .misc import execution_path
from .logging import my_logging

class _Config(object):
    """
    This is the place where to put stuff that any module may need, a kind of COMMON.
    An instantiation is done in the main __init__ file using the "config" name.

    """


    def __init__(self):
        """
        Define variables that will be known from everywhere:
        - INSTALLED: a dictionary whose keys are the libraries that PyNeb may need
            e.g. 'plt' for matplotlib.pyplot, 'scipy' for scipy.interpolate. and whose values
            are Boolean
        
        
        Parameters: none
        
        """
        self.log_ = my_logging()
        self.calling = '_Config'
        
        self.INSTALLED = {}
        try:
            import matplotlib.pyplot as plt
            self.INSTALLED['plt'] = True
        except:
            self.INSTALLED['plt'] = False
            self.log_.message('matplotlib not available', calling=self.calling)
        try:
            from scipy import interpolate
            self.INSTALLED['scipy'] = True
        except:
            self.INSTALLED['scipy'] = False
            self.log_.message('scipy not available', calling=self.calling)
        try:
            import multiprocessing as mp
            self.INSTALLED['mp'] = True
            self.Nprocs = mp.cpu_count()
        except:
            self.INSTALLED['mp'] = False
            self.log_.message('multiprocessing not available', calling=self.calling)
            self.Nprocs = 1
        try:
            from PyQt4 import QtCore, QtGui
            self.INSTALLED['Qt4'] = True
        except:
            self.INSTALLED['Qt4'] = False
            self.log_.message('pySSN Qt not available', calling=self.calling)
        try:
            import pyneb
            self.INSTALLED['PyNeb'] = True
        except:
            self.INSTALLED['PyNeb'] = False
        try:
            import yaml
            self.INSTALLED['YAML'] = True
        except:
            self.INSTALLED['YAML'] = False
            self.log_.message('YAML config files not available', calling=self.calling)
        
        self.DataPaths = []
        self.addDataFilePath('../data/', inpySSN=True)
        self.addDataFilePath('./', inpySSN=False)
        #self.addDataFilePath('/', inpySSN=False)
        self.unuse_multiprocs()
                    
    def use_multiprocs(self):
        self._use_mp = True
    
    def unuse_multiprocs(self):
        self._use_mp = False    
        
    def addDataFilePath(self, dir_=None, inpySSN=False):
        """
        Add a directory to the list of directories where atomic data files are searched for.
        
        Parameters:
           - dir_    directory
           - inpySSN     Boolean.
    
        """
        if dir_ is None:
            self.DataPaths = []
        else:
            try:
                if inpySSN:
                    self.DataPaths.insert(0, execution_path(dir_))
                else:
                    self.DataPaths.insert(0, os.path.abspath(dir_))
            except:
                self.log_.warn('{0} could not be added to the path list', calling='addDataFilePath')
                
