#!/usr/bin/env python

from __future__ import print_function
import subprocess as sp
import argparse
import os
import re

def get_paths(root, extension):
    # scan directories
    paths = []
    for root, dirs, files in os.walk(root):
        for name in files:
            if name.endswith(extension): 
                paths.append(os.path.join(root, name))
    paths.sort()
    return paths

def print_title(title):
     print('# '+'-'*77+'\n'+'# '+title+'\n'+'# '+'-'*77)

def dependencies(root, src_ext = '.cc'):
    print_title('find dependencies')
    # get paths
    paths = get_paths(root, src_ext)
    total_deps = 0
    total_size = 0
    for path in paths:
        # run gcc and get dependeny data 
        command = 'g++ {0} -I{1} -MM'.format(path, root)
        output  = sp.Popen([command], shell=True, stdout=sp.PIPE).communicate()[0]
        deps    = len(output.replace('\\', '').split())-1
        size    = os.path.getsize(path)
        # code analysis
        print('{0:<60}: {1:6d}b {2:4d} deps'.format(path, size, deps))
        total_deps += deps
        total_size += size
    print('Total {0:48d} files: {1:6d}b {2:4d} deps'.format(len(paths), total_size, total_deps))

def trailing_whitespaces(root, extensions):
    print_title('check for trailing whitespaces')
    paths = []
    for ext in extensions:
        paths += get_paths(root, ext)
    paths.sort()
    num = 0
    for path in paths:
        ifile = open(path, 'r')
        lines = ifile.readlines()
        ifile.close()
        for linenum, line in enumerate(lines):
            line = line[:-1] # strip newline 
            if line and line[-1].isspace():
                print('{0:<60}\tline {1:3d}'.format(path, linenum+1)) 
                num += 1
    print('Total trailing whitespaces: {0:5d}'.format(num)) 

def function_names(root, extensions):
    print_title('check if function names stast with lowercase')
    paths = []
    for ext in extensions:
        paths += get_paths(root, ext)
    paths.sort()
    num   = 0
    total = 0
    for path in paths:
        with open(path, 'r') as f:
            for line in f.readlines():
                m = re.search(r'(\w+)\(', line)
                if m != None and m.group(1) != 'Author':
                    total += 1
                    g = m.group(1)
                    if g[0].islower():
                        num += 1
                        print('{0:4d} {1:<24s} {2:<48s}'.format(num, g, line[0:m.end(0)]))
    perc = float(num)/float(total)*100.
    print('Functions  : {0:6d}\nLowercase  : {1:6d}\nPercentage : {2:6.2f}%'.format(num, total, perc))    

def grep(root, extensions, pattern):
    print_title('grep pattern')
    paths = []
    for ext in extensions:
        paths += get_paths(root, ext)
    paths.sort()
    total = 0
    for path in paths:
        with open(path, 'r') as f:
            num = 0
            print_path = False
            for line in f.readlines():
                m = re.search(pattern, line)
                if m != None:
                    total += 1
                    if not print_path:
                        print('{0}'.format(path))
                        print_path = True
                    print('  {0}'.format(line.strip()))
    print('Found {0:5d} occurencies of {1}'.format(total, pattern))    

def self_contained_headers(root, extensions, flags):
    print_title('check if header files are self contained')
    paths = []
    for ext in extensions:
        paths += get_paths(root, ext)
    paths.sort()
    num = 0
    for path in paths:
        num += 1
        f = open('tmp.cc', 'w')
        f.write('#include \"{0}\"\n'.format(path))
        f.close()
        print('checking: {0}'.format(path))
        command = 'g++ {0} -c tmp.cc'.format(flags)
        output  = sp.Popen([command], shell=True, stdout=sp.PIPE).communicate()[0]
        os.remove('tmp.cc')
        os.remove('tmp.o')

if __name__ == '__main__':
    flags  = '-Isrc -I/usr/include -I/usr/include/python2.6 -I/opt/acml4.4.0/gfortran32_mp/include -Wno-deprecated'
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--deps',        action = 'store_true', 
                             dest = 'deps',    help = 'grep pattern')
    parser.add_argument('-n', '--functions',   action = 'store_true', 
                             dest = 'fn',      help = 'functions lowercase')
    parser.add_argument('-w', '--whitespaces', action = 'store_true', 
                             dest = 'w',       help = 'trailing whitespaces')
    parser.add_argument('-g', '--grep',        default = None, 
        metavar = 'pattern', dest = 'greppat', help = 'grep pattern')
    parser.add_argument('-s', '--self-cont',   action = 'store_true',
                             dest = 'sc',      help = 'self contained headers')
    parser.add_argument('-c', '--cxx-ext',     default = '.cc', 
        metavar = 'src_ext', dest = 'src_ext', help = 'source file extension')
    parser.add_argument('-i', '--inc-ext',     default = '.h',
        metavar = 'inc_ext', dest = 'inc_ext', help = 'header file extension')
    parser.add_argument('-r', '--root',        default = 'src',
        metavar = 'src',     dest = 'src',     help = 'root directory')
    parser.add_argument('-f', '--flags',       default = flags, 
        metavar = 'flags',   dest = 'flags',   help = 'compilation flags')
    args         = parser.parse_args()
    src          = args.src
    src_ext      = args.src_ext
    inc_ext      = args.inc_ext
    
    if args.deps:
        dependencies(args.src)
    if args.w:
        trailing_whitespaces(args.src, [src_ext, inc_ext])
    if args.fn:
        function_names(args.src, [inc_ext])
    if args.greppat != None:
        grep(args.src, [src_ext, inc_ext], args.greppat)
    if args.sc:
        self_contained_headers(args.src, [inc_ext], args.flags)

