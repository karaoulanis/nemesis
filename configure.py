#!/usr/bin/env python

from __future__ import print_function
import subprocess as sp
import os
import re
import argparse

def get_target(project, path, item, num_items, verbal, color, src_ext = '.cc', obj_ext = '.o'):
    # run gcc and get dependeny data 
    objname = path.replace(src_ext, obj_ext)
    command = 'g++ {0} -Isrc -MM -MT \'{1}\''.format(path, objname)
    output  = sp.Popen([command], shell=True, stdout=sp.PIPE).communicate()[0]
    num_deps= len(output.replace('\\', '').split())-1
    size    = os.path.getsize(path)
    sign    = '' if verbal else '@'
    color0  = '\e[00m'    if color else ''
    color1  = '\e[01;32m' if color else ''
    color2  = '\e[01;34m' if color else ''
    # code analysis
    print('{0:<60}: {1:3d} dependencies'.format(objname, num_deps))
    # create target string and return
    s  = '\n# file : {0}\n'.format(path)
    s += output
    s += '\t@echo -e \'{0}compiling '.format(color1)
    s += '[{0:3d}/{1:3d}| '.format(item+1, num_items)
    s += '{0:5d} bytes|{1:3d} deps]{2}: '.format(size, num_deps, color0)
    s += '{0}{1}{2}\'\n'.format(color2, path, color0)
    s += '\t{0}$(CXX) $(CXXFLAGS) $({1}_CPPFLAGS) \\\n'.format(sign, project)
    s += '\t-c {0} \\\n\t-o {1}\n'.format(path, objname)
    return s

def get_objects(project, paths, verbal, color, src_ext = '.cc', obj_ext = '.o'):
    sign    = '' if verbal else '@'
    color0  = '\e[00m'    if color else ''
    color1  = '\e[01;32m' if color else ''
    color2  = '\e[01;34m' if color else ''
    s  = '\n{0}_OBJS = '.format(project)
    for path in paths:
        objname = os.path.basename(path).replace(src_ext, obj_ext)
        dirname = os.path.dirname(path) 
        s += ' \\\n\t{0}/{1}'.format(dirname, objname)
    s += '\n\n{0}: $({0}_OBJS)\n'.format(project)
    s += '\t@echo -e \'{0}linking   '.format(color1)
    s += '[{0:16d} object files]{1}: '.format(len(paths), color0)
    s += '{0}{1}{2}\'\n'.format(color2, project, color0)
    s += '\t{0}$(CXX) -o {1} $(LDFLAGS) $(LDADD) $({1}_LDFLAGS) $({1}_LDADD) $({1}_OBJS)\n'.format(sign, project)
    return s

def get_paths(source_dir, include_pat, exclude_pat = None):
    paths = []
    for root, dirs, files in os.walk(source_dir):
        for name in files:
            if ((re.search(include_pat, name) != None) and 
                (exclude_pat == None or re.search(exclude_pat, name) == None)):
                paths.append(os.path.join(root, name))
    paths.sort()
    return paths

if __name__ == '__main__':
    # project specific 
    project        = 'nemesis'
   
    # command line arguments 
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--def', metavar='file',
        dest = 'makefile_def', help = 'makefile definitions file',
        default = 'Makefile.def')
    parser.add_argument('-r', '--root', metavar='dir',
        dest = 'src',          help = 'define source root directory',
        default = 'src')
    parser.add_argument('-m', '--monochrome', action = 'store_false',
        dest = 'color',        help = 'disable colored output')
    parser.add_argument('-s', '--silent', action = 'store_false',
        dest = 'verbal',       help = 'disable additional compilation info')
    args = parser.parse_args()
    src           = args.src
    verbal        = args.verbal 
    color         = args.color
    makefile_def  = args.makefile_def
    
    makefile      = ''
    # Makefile.def
    with open(makefile_def, 'r') as f:
        makefile += f.read()

    # project: nemesis
    paths = get_paths(src, r'\.cc$', r'_test\.cc$')
    makefile += get_objects(project, paths, verbal, color)
    for i, path in enumerate(paths):
        makefile += get_target(project, path, i, len(paths), verbal, color)
    
    # project: nemesis-tests
    paths = get_paths(src, r'\.cc$', None)
    makefile += get_objects('{0}-tests'.format(project), paths, verbal, color)
    paths = get_paths(src, r'_test\.cc$', None)
    for i, path in enumerate(paths):
         makefile += get_target('{0}-tests'.format(project), path, i, len(paths), verbal, color)
    
    ofile=open('Makefile','w')
    ofile.write(makefile)
    ofile.close()

