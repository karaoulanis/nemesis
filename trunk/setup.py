#!/usr/bin/env python2.7

from __future__ import print_function
import subprocess as sp
import os
import re
import argparse

makefile_def_default ="""
CXX             = g++
CXXFLAGS        = -g -O2 -Wno-deprecated -Wall -Wextra    \\
                  -Wc++0x-compat -pedantic -ansi -Weffc++
CPPFLAGS        = -Isrc -I/usr/include                    \\
                  -I/usr/include/python2.7                \\
                  -I/opt/acml4.4.0/gfortran32_mp/include
LIBFLAGS        = -L/usr/lib                              \\
                  -lpython2.7                             \\
                  -L/opt/acml4.4.0/gfortran32_mp/lib      \\
                  -lacml_mp -lgfortran -lgomp -lrt -lm    
TEST_CPPFLAGS   = 
TEST_LIBFLAGS   = -lgtest
"""

def get_target(project, path, item, num_items,
               cxx, cxxflags, cppflags,
               verbal, color,
               src_ext = '.cc', obj_ext = '.o'):
    # get object name 
    obj_name = path.replace(src_ext, obj_ext)
    # run gcc and get dependency data 
    command  = 'g++ {0} -Isrc -MM -MT \'{1}\''.format(path, obj_name)
    output   = sp.Popen([command], shell=True, stdout=sp.PIPE).communicate()[0]
    num_deps = len(output.replace('\\', '').split())-1
    size     = os.path.getsize(path)
    echo     = '' if verbal else '@'
    color0   = '\e[00m'    if color else ''
    color1   = '\e[01;32m' if color else ''
    color2   = '\e[01;34m' if color else ''
    print('{0:<60}: {1:3d} dependencies'.format(obj_name, num_deps))
    # create target string and return
    s  = '\n# file : {0}\n'.format(path)
    s += output
    s += '\t@echo -e \'{0}compiling '.format(color1)
    s += '[{0:3d}/{1:3d}| '.format(item+1, num_items)
    s += '{0:5d} bytes|{1:3d} deps]{2}: '.format(size, num_deps, color0)
    s += '{0}{1}{2}\'\n'.format(color2, path, color0)
    s += '\t{0}{1} {2} {3} '.format(echo, cxx, cxxflags, cppflags)
    s += '-c {0} -o {1}\n'.format(path, obj_name)
    return s

def get_objects(project, paths,
                libflags,
                verbal, color,
                src_ext = '.cc', obj_ext = '.o'):
    echo    = '' if verbal else '@'
    color0  = '\e[00m'    if color else ''
    color1  = '\e[01;32m' if color else ''
    color2  = '\e[01;34m' if color else ''
    s  = '\n{0}_OBJS = '.format(project)
    for path in paths:
        obj_name = os.path.basename(path).replace(src_ext, obj_ext)
        dir_name = os.path.dirname(path) 
        s += ' \\\n\t{0}/{1}'.format(dir_name, obj_name)
    s += '\n\n{0}: $({0}_OBJS)\n'.format(project)
    s += '\t@echo -e \'{0}linking   '.format(color1)
    s += '[{0:16d} object files]{1}: '.format(len(paths), color0)
    s += '{0}{1}{2}\'\n'.format(color2, project, color0)
    s += '\t{0}$(CXX) -o {1} {2} $({1}_OBJS)\n'.format(echo, project, libflags)
    return s

def get_paths(src, endswith):
    paths = []
    for root, dirs, files in os.walk(src):
        for name in files:
            if name.endswith(endswith):
                paths.append(os.path.join(root, name))
    paths.sort()
    return paths

def get_project(project_name, verbal, color):
    all_paths = get_paths('src', '.cc')
    s = ''
    # project objects
    object_paths = [i for i in all_paths if not i.endswith('_test.cc')]
    s += get_objects(project_name, object_paths,
                        '$(LIBFLAGS)',
                        verbal, color)
    # project targets
    target_paths = [i for i in all_paths if not i.endswith('_test.cc')]
    for i, path in enumerate(target_paths):
        s += get_target(project_name, path, i, len(target_paths),
                        '$(CXX)', '$(CXXFLAGS)', '$(CPPFLAGS)',
                        verbal, color)
    # project-test name
    project_test_name = '{0}-tests'.format(project_name)
    # project-test objects
    object_paths = [i for i in all_paths if os.path.basename(i)!='main.cc']
    s += get_objects(project_test_name, object_paths,
                        '$(TEST_LIBFLAGS)',
                        verbal, color)
    # project targets
    target_paths = [i for i in all_paths if i.endswith('_test.cc')]
    for i, path in enumerate(target_paths):
        s += get_target(project_test_name, path, i, len(target_paths),
                        '$(CXX)', '$(TEST_CXXFLAGS)', '$(TEST_CPPFLAGS)',
                        verbal, color)
    # project: clean
    s += '\n# clean\n'
    s += 'clean:\n'
    s += '\trm src/*/*.o\n'
    s += '\trm nemesis\n'
    # return 
    return s

if __name__ == '__main__':
    
    # command line arguments 
    parser = argparse.ArgumentParser()
    
    parser.add_argument('project',
                        nargs   = '?',
                        default = os.path.basename(os.getcwd()),
                        help    = '')
                        
    parser.add_argument('--def',
                        metavar='file',
                        dest = 'makefile_def',
                        default = 'Makefile.def',
                        help = 'makefile def filename (default: Makefile.def)')

    parser.add_argument('--nocolor',
                        action = 'store_false',
                        dest   = 'color',
                        help = 'suppress colored output')

    parser.add_argument('--silent',
                        action = 'store_false',
                        dest = 'verbal',
                        help = 'suppress additional compilation info')
    
    
    args = parser.parse_args()
    project       = args.project
    makefile_def  = args.makefile_def
    verbal        = args.verbal 
    color         = args.color
    
    
    
    # makefile
    with open('Makefile','w') as makefile:
        makefile.write('include {0}'.format(makefile_def))
        makefile.write(get_project(project, verbal, color))
