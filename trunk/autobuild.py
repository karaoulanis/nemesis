#!/usr/bin/env python

from __future__ import print_function
import subprocess as sp
import argparse
import os
import shutil 
import time
import re

def print_title(title):
     print('# '+'-'*77+'\n'+'# '+title+'\n'+'# '+'-'*77)

def get_revision(url):
    command = 'svn info {0}'.format(url)
    proc = sp.Popen([command], shell=True, stdout=sp.PIPE)
    proc.wait()
    m = re.search(r'Revision\:\s*([0-9]+)',proc.stdout.read())
    revision = '00000'
    if m != None:
        revision = '{0:0>5s}'.format(m.group(1))
    return revision

def exec_command(command, logfile, title):
    print_title(title)
    proc = sp.Popen([command], shell=True, stdout=sp.PIPE)
    while True:
        line = proc.stdout.readline()
        if line == '' and proc.poll() != None:
            break
        else:
            print(line.rstrip())
            logfile.write(line)
    print()
    logfile.write('\n')
    
if __name__ == '__main__':
    # command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--url', metavar='url',
        dest = 'url',          help = 'svn url',
        default = 'http://nemesis-code.googlecode.com/svn/trunk/')
    parser.add_argument('-r', '--root', metavar='dir',
        dest = 'path',         help = 'define root directory',
        default = 'bld')

    args        = parser.parse_args()
    build_path  = args.path.rstrip('/')
    temp_path   = build_path + '/temp'
    url         = args.url
    
    revision    = get_revision(url)
    timestamp   = time.strftime("%Y%m%dT%H%M%S")
    filename    = '{0}/build_{1}_{2}.log'.format(build_path, revision, timestamp)
    logfile     = open(filename, 'w')

    # save current path
    saved_path = os.getcwd()
    # svn checkout
    command = 'svn checkout {0} {1}'.format(url, temp_path)
    exec_command(command, logfile, 'checkout')
    # move to project path
    os.chdir(temp_path)
    # configure
    command = './configure.py'
    exec_command(command, logfile, 'configure')
    # make 
    command = 'make -k'
    exec_command(command, logfile, 'make')
    # return to initial path
    os.chdir(saved_path)
    # delete path
    shutil.rmtree(temp_path)
    

    logfile.close()
    

