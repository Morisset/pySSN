'''
Created on 30 sept. 2017

@author: christophemorisset
'''
import os
import sys
import subprocess

def locate_fortran():
    path = os.path.dirname(os.path.abspath(__file__))
    return path+'/../fortran/'
 
def compile_XSSN(fcompiler= 'gfortran', fname='XSSN_Phyat.f', oname='XSSN_Phyat.exe', foptions=''):
    to_run = 'cd {0} ; {1} {2} -o {3} {4}'.format(locate_fortran(), fcompiler,fname, oname, foptions)
    try:
        proc = subprocess.call(to_run, shell=True)
    except:
        raise NameError('Fortran compilation error')

def make_link2data(datadir):
    link_name = os.path.dirname(os.path.abspath(__file__))
    link_name = '/'.join(link_name.split('/')[:-1]) + '/fortran/data'
    try:
        os.symlink(datadir, link_name)
    except OSError, e:
        if e.errno != os.errno.EEXIST:
            raise 
        pass
    if not os.path.exists(link_name+'/d1'):
        print('!!! The link does not seems to contain d1, d2, d3 and d4 subdirectories.')
    else:
        print('Link OK')

def run_XSSN(fphycond='phy_cond.dat', fionfrac='1789409_ionfrac.dat', fabund='asplund_2009.dat'):
    pass