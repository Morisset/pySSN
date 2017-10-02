'''
Created on 1 oct. 2017

@author: christophemorisset
'''
from shutil import copyfile
from ..phyat_lists.entries import execution_path

def run_XSSN(outputcond_file):
    copyfile(outputcond_file, execution_path('outputcond.dat'))