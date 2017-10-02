'''
Created on 1 oct. 2017

@author: christophemorisset
'''
from shutil import copyfile
from ..phyat_lists.entries import execution_path
import subprocess
import sys, os

def run_XSSN(outputcond_file):
    copyfile(outputcond_file, execution_path('outputcond.dat'))
    to_run = 'cd {} ; XSSN_Phyat.exe'.format(os.path.dirname(sys._getframe(1).f_code.co_filename))
    stdout = open('XSSN.out', 'w')   
    subprocess.Popen(to_run, shell=True, stdout=stdout)
    