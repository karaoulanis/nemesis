#!/usr/bin/env python

from __future__ import print_function
import subprocess as sp
import argparse
import os

def svn_checkout(url, path):
    #  
    command = 'svn checkout {0} {1}'.format(url, path)
    # output  = sp.Popen([command], shell=True, stdout=sp.PIPE).communicate()[0]
    proc = sp.Popen([command], shell=True)
    proc.wait()

def configure():
    command = './configure.py'
    proc = sp.Popen([command], shell=True)
    proc.wait()

def make():
    command = 'make -k'
    proc = sp.Popen([command], shell=True)
    proc.wait()

if __name__ == '__main__':
    # command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--url', metavar='url',
        dest = 'url',          help = 'svn url',
        default = 'http://nemesis-code.googlecode.com/svn/trunk/')
    parser.add_argument('-r', '--root', metavar='dir',
        dest = 'path',         help = 'define root directory',
        default = 'bld/auto')
    args = parser.parse_args()
    path = args.path
    url  = args.url
    
    saved_path = os.getcwd()
    svn_checkout(url, path)
    os.chdir(path)
    configure()
    make()
    os.chdir(saved_path)
   
