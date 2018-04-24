'''
Created on 2 oct. 2017

@author: christophemorisset
'''
import os
import argparse

def make_link2data(datadir):
    link_name = os.path.dirname(os.path.abspath(__file__)) + '/data'
    if os.path.exists(link_name):
        os.unlink(link_name)
    try:
        os.symlink(datadir, link_name)
    except OSError, e:
        if e.errno != os.errno.EEXIST:
            raise 
        pass
    if not os.path.exists(link_name+'/d1'):
        print('!!! The link does not seem to contain d1 subdirectory.')
    else:
        print('Link OK')

def link2data():
    parser = argparse.ArgumentParser()
    parser.add_argument("datadir", help="data directory")
    args = parser.parse_args()
    
    make_link2data(args.datadir)
