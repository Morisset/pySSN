'''
Created on 2 oct. 2017

@author: christophemorisset
'''
import os
import sys
import argparse
from shutil import copyfile

def execution_path(filename, extra=''):
    return os.path.join(os.path.dirname(sys._getframe(1).f_code.co_filename), extra, filename)

def print_phy_cond():
    with open('phy_cond.dat', 'w') as f:
        f.write("""IP 1000000   CR 1e4    1e2
IP 2000000 CR 1e4 1e2        
""")    

def print_asplund(file_out='asplund_2009.dat'):
    file_in = execution_path('asplund_2009.dat')
    copyfile(file_in, file_out)

def print_outputcond(file_out='outputcond.dat'):
    
    file_in = execution_path('outputcond.dat', extra='../fortran/')
    copyfile(file_in, file_out)

def print_ionfrac():
    parser = argparse.ArgumentParser()
    parser.add_argument("-T", "--Teff", help="Effective temperature")
    parser.add_argument("-logU", "--logU", help="Ionization parameter")
    parser.add_argument("-B", "--B", help="Matter/Radiation Bounded [R or M60]", default=True)
    
    args = parser.parse_args()

    if args.Teff not in ['25', '50', '75', '100', '150', '300']:
        raise Exception("Teff must be in [25, 50, 75, 100, 150, 300]")
    if args.logU not in ['-1', '-2', '-3']:
        raise Exception("Log U must be in [-1, -2, -3]")
    if args.B not in ['R', 'M60']:
        raise Exception('B must by R (Radiation bounded) or M60 (60% matter Bounded)')
    filename = '{}_{}_{}_ionfrac.dat'.format(args.Teff, args.logU, args.B)
    file_in = execution_path(filename, extra='ionfracs/')
    copyfile(file_in, filename)