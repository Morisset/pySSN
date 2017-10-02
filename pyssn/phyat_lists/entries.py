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

def print_asplund(file_out='abunds.dat'):
    file_in = execution_path('asplund_2009.dat')
    copyfile(file_in, file_out)

def print_outputcond(file_out='outputcond.dat'):
    
    file_in = execution_path('outputcond_ex.dat', extra='../fortran/')
    copyfile(file_in, file_out)

def print_ionfrac():
    parser = argparse.ArgumentParser()
    parser.add_argument("-T", "--Teff", help="Effective temperature")
    parser.add_argument("-logU", "--logU", help="Ionization parameter")
    parser.add_argument("-B", "--B", help="Matter/Radiation Bounded [R or M60]", default=True)
    args = parser.parse_args()
    print_ionfrac(args)

def print_ionfrac_args(args):

    if args.Teff not in ['25', '50', '75', '100', '150', '300']:
        raise Exception("Teff must be in [25, 50, 75, 100, 150, 300]")
    if args.logU not in ['-1', '-2', '-3']:
        raise Exception("Log U must be in [-1, -2, -3]")
    if args.B not in ['R', 'M60']:
        raise Exception('B must by R (Radiation bounded) or M60 (60% matter Bounded)')
    filename = '{}_{}_{}_ionfrac.dat'.format(args.Teff, args.logU, args.B)
    file_in = execution_path(filename, extra='ionfracs/')
    copyfile(file_in, filename)
    
def print_files():
    parser = argparse.ArgumentParser()
    parser.add_argument("-P", "--print_phy_cond", action="store_true", help="Write on disk the phy_cond.dat file")
    parser.add_argument("-I", "--print_ionfrac", action="store_true", help="Write on disk the ionfrac file")
    parser.add_argument("-A", "--print_asplund", action="store_true", help="Write on disk the asplund_2009.dat file")
    parser.add_argument("-O", "--print_outputcond", action="store_true", help="Write on disk the outputcond.dat file")
    
    parser.add_argument("-T", "--Teff", help="Effective temperature")
    parser.add_argument("-logU", "--logU", help="Ionization parameter")
    parser.add_argument("-B", "--B", help="Matter/Radiation Bounded [R or M60]", default=True)
    
    args = parser.parse_args()
    if args.print_phy_cond:
        print_phy_cond()
    if args.print_ionfrac:
        print_ionfrac_args(args)
    if args.print_asplund:
        print_asplund()
    if args.print_outputcond:
        print_outputcond()
