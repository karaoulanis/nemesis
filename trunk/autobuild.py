#!/usr/bin/env python

from __future__ import print_function
import subprocess as sp
import argparse
import os

def svn_checkout(url, path):
    #  
    command = 'svn checkout {0} {1}'.format(url, path)
    output  = sp.Popen([command], shell=True, stdout=sp.PIPE).communicate()[0]

if __name__ == '__main__':
    # command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--url', metavar='url',
        dest = 'url',          help = 'svn url',
        default = 'http://nemesis-code.googlecode.com/svn/trunk/')
    parser.add_argument('-r', '--root', metavar='dir',
        dest = 'path',         help = 'define root directory',
        default = 'path')
    args = parser.parse_args()
    path = args.path
    url  = args.url
    
    svn_checkout(url, path) 
