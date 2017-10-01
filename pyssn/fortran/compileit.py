'''
Created on 2 oct. 2017

@author: christophemorisset
'''
import os
import subprocess
import argparse

def compileit(fcompiler= 'gfortran', fname='XSSN_Phyat.f', oname='XSSN_Phyat.exe', foptions=''):
    fortran_path = os.path.dirname(os.path.abspath(__file__))
    to_run = 'cd {0} ; {1} {2} -o {3} {4}'.format(fortran_path, fcompiler,fname, oname, foptions)
    try:
        proc = subprocess.call(to_run, shell=True)
    except:
        raise NameError('Fortran compilation error')

def compile_XSSN():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--compiler", help="fortran compiler", default="gfortran")
    args = parser.parse_args()
    compileit(fcompiler=args.compiler)
