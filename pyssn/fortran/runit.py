'''
Created on 1 oct. 2017

@author: christophemorisset
'''
import shutil
import subprocess
import sys, os
from ..phyat_lists.entries import execution_path

def run_XSSN(outputcond_file):
    shutil.copyfile(outputcond_file, execution_path('outputcond.dat'))
    to_run = 'cd {} ; ./XSSN_Phyat.exe'.format(os.path.join(os.path.dirname(sys._getframe(1).f_code.co_filename), '../fortran/'))
    stdout = open('XSSN.out', 'w')   
    process = subprocess.Popen(to_run, shell=True, stdout=stdout, close_fds=True)
    stdout.close()
    