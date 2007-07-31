
import os,sys,getopt
from time import time,ctime

#==============================================================================
# Get project files
#==============================================================================
def getProjectFiles():
	""" Get files/paths from a given directory.
		1.	First walk down the directory tree and gather all paths/files 
			matching a pattern.
		2.	Then apply some platform specific changes (eg. replace \\ by \ in 
			win32.)
		3. 	Rearrange so as 'main.cpp' to be the first in the lists.
		4. 	Gather dependencies for each file found. This involves opening the
			file and reading includes.
		5. 	If any includes are not in the given directory, (eg. std, boost 
			etc.) remove them from the lists.
		6. 	Run iteratively and include to each dependency with its 
			dependencies. This should stop when no other inclusions take place.
		7.	Now replace each dependency with path+dependency.
		8. 	Create an array as [c_path,c_filename,[depedencies...],...].
		9. 	Create an array as [h_path,...].
		10. Return the above arrays.
	"""
	#--------------------------------------------------------------------------
	# Find paths/names at a directory according to some pattern
	#--------------------------------------------------------------------------
	print "Searching in the src directory    :",
	paths=[]
	names=[]
	for path, folders, files in os.walk('../..'+'/src'):
		for name in files:
			if os.path.splitext(name)[1] in ['.cpp','.h']:
				paths.append(path)
				names.append(name)
	nFiles=len(paths)
	print nFiles,"files found."
	#--------------------------------------------------------------------------
	# Change '\\' in paths (win32)
	#--------------------------------------------------------------------------
	for i in range(nFiles):
		paths[i]=paths[i].replace('\\','/')
	#--------------------------------------------------------------------------
	# Bring main.cpp in the first row
	#--------------------------------------------------------------------------
	index=names.index("main.cpp")
	tmp=paths[index];	paths[index]=paths[0];	paths[0]=tmp
	tmp=names[index];	names[index]=names[0];	names[0]=tmp
	#--------------------------------------------------------------------------
	# Find first level dependencies (open files)
	#--------------------------------------------------------------------------
	depds=[]
	for i in range(nFiles):
		ifile=open(paths[i]+'/'+names[i],'r')
		tmp=[]
		for line in ifile.readlines():
			line=line.replace('<',' ')
			line=line.replace('>',' ')
			line=line.replace('"',' ')
			words=line.split()
			if len(words)>0 and words[0]=='#include':
				tmp.append(words[1])
		ifile.close()
		depds.append(tmp)
	#--------------------------------------------------------------------------
	# Remove dependencies not in given directories
	#--------------------------------------------------------------------------
	for i in range(nFiles):
		tmp=[]
		for dependency in depds[i]:
			if dependency in names:
				tmp.append(dependency)
		depds[i]=tmp
	#--------------------------------------------------------------------------
	# Find dependencies iteratively
	#--------------------------------------------------------------------------
	print "Finding dependencies (passes)     :",
	k=1
	while k!=0:
		k=0
		print "#",
		for i in range(nFiles):
			tmp=[]
			for dependency in depds[i]:
				tmp=tmp+depds[names.index(dependency)]
			for dependency in tmp:
				if dependency not in depds[i]:
					depds[i].append(dependency)
					k=1
	print
	#--------------------------------------------------------------------------
	# Find dependencies paths
	#--------------------------------------------------------------------------
	for i in range(nFiles):
		tmp=[]
		for dependency in depds[i]:
			tmp.append(paths[names.index(dependency)]+'/'+dependency)
		depds[i]=tmp
	#--------------------------------------------------------------------------
	# Find cpp files
	#--------------------------------------------------------------------------
	cFiles=[]
	for i in range(nFiles):
		if os.path.splitext(names[i])[1] in ['.cpp']:
			cFiles.append([paths[i],names[i],depds[i]])
	#--------------------------------------------------------------------------
	# Find included files
	#--------------------------------------------------------------------------
	hPaths=[]
	for i in range(nFiles):
		if os.path.splitext(names[i])[1] in ['.h']:
			if paths[i] not in hPaths:
				hPaths.append(paths[i])
	return cFiles,hPaths
#==============================================================================
# Get configure options 
#==============================================================================
def getConfOptions(optlist):
	""" Get options from the command line.
		Four arrays are returned:
		hPaths	: Holds paths to external included directories.	[path,...]
		lPaths	: Holds paths to library directories.			[path,...]
		lFiles	: Holds library names.							[lib,... ]
		compiler: Holds compiler and compiler options.	 		[compiler,opts]
		Default values are needed for the compiler list, in case nothing is
		passed through the options list.
	"""
	#--------------------------------------------------------------------------
	# List and defaults
	#--------------------------------------------------------------------------
	hPaths=[]
	lPaths=[]
	lFiles=[]
	compiler=[]
	compiler.append('g++')
	compiler.append('-O2 -fWall')
	#--------------------------------------------------------------------------
	# Read options
	#--------------------------------------------------------------------------
	for i in optlist:
		if   i[0]=='--with-inc-path':
			print "Configuring with include path     :",i[1]
			hPaths.append('\''+i[1]+'\'')	# may have spaces
		elif i[0]=='--with-lib-path':
			print "Configuring with library path     :",i[1]
			lPaths.append('\''+i[1]+'\'')	# may have spaces
		elif i[0]=='--with-lib-name':
			print "Configuring with library name     :",i[1]
			lFiles.append(i[1])
		elif i[0]=='--with-cxx-comp':
			print "Configuring with compiler         :",i[1]
			compiler[0]=i[1]
		elif i[0]=='--with-cxx-opts':
			print "Configuring with compiler options :",i[1]
			compiler[1]=i[1]
		else:
			continue
	#--------------------------------------------------------------------------
	# Return
	#--------------------------------------------------------------------------
	return hPaths,lPaths,lFiles,compiler

#==============================================================================
# Create Makefile
#==============================================================================
def createMakefile(cFiles,hPaths,lPaths,lFiles,compiler):
	""" Create the Makefile.
		Parameters:
		cFiles	: [path,file,[depedencies],...]
		hPaths	: [path,...]
		lPaths	: [path,...]
		lFiles	: [lib,...]
		compiler: [compiler,options]
	"""
	mak = open("Makefile", 'w+')
	mak.write("#===============================================================\n")
	mak.write("# Makefile for nemesis                                          \n")
	mak.write("# Maintained by : fk._	               	                       \n")
	mak.write("# Created by    : slvMake.py [Automatically]                    \n")
	mak.write('# Date          : %s\n'%ctime(time()))
	mak.write("# Version       : 1.0d                                          \n")
	mak.write("#===============================================================\n")

	mak.write("#===============================================================\n")
	mak.write("# options\n")
	mak.write("#===============================================================\n")
	mak.write("EDIR    = ../../bin\n")
	mak.write("EXEC    = nemesis.exe\n")
	mak.write("ODIR    = obj\n")
	mak.write("CC      = %s\n"%(compiler[0]))
	mak.write("CCFLAGS = %s\n"%(compiler[1]))
	mak.write("\n")

	mak.write("#===============================================================\n")
	mak.write("# library directories\n")
	mak.write("#===============================================================\n")
	mak.write("LDIR = ")
	dirs=[]
	for i in lPaths:
		mak.write(' \\\n\t-I%s'%i)
	mak.write("\n\n")

	mak.write("#===============================================================\n")
	mak.write("# libraries\n")
	mak.write("#===============================================================\n")
	mak.write("LIBS = ")
	for i in lFiles:
		mak.write(' \\\n\t%s'%i)
	mak.write("\n\n")

	mak.write("#===============================================================\n")
	mak.write("# include directories\n")
	mak.write("#===============================================================\n")
	mak.write("INCS = ")
	dirs=[]
	for i in hPaths:
		mak.write(' \\\n\t-I%s'%i)
	mak.write("\n\n")

	mak.write("#===============================================================\n")
	mak.write("# object files\n")
	mak.write("#===============================================================\n")
	mak.write("OBJS = ")
	for i in cFiles:
		mak.write('\\\n\t$(ODIR)/%s.o'%i[1].replace('.cpp',''))
	mak.write("\n\n")

	mak.write("#===============================================================\n")
	mak.write("# link\n")
	mak.write("#===============================================================\n")
	mak.write("link: $(OBJS) \n")
	mak.write("\t@echo Linking...\n")
	mak.write("\t@$(CC) $(CCFLAGS) -o $(EDIR)/$(EXEC) $(OBJS) $(LDIR) $(LIBS)\n")
	mak.write("\n")

	mak.write("#===============================================================\n")
	mak.write("# compile\n")
	mak.write("#===============================================================\n")
	k=0
	for i in cFiles:
		k=k+1
		cpath=i[0]
		cfile=i[1].replace('.cpp','')
		cdeps=cpath+'/'+cfile+'.cpp'
		for i in i[2]:
			cdeps=cdeps+' '+i
		mak.write('$(ODIR)/%s.o: %s/%s.cpp %s\n'%(cfile,cpath,cfile,cdeps))
		mak.write('\t@echo Compiling [% 4i/% 4i]: %s.cpp\n'%(k,len(cFiles),cfile))
		mak.write('\t@$(CC) $(CCFLAGS) $(INCS) -c %s/%s.cpp -o $(ODIR)/%s.o\n'%(cpath,cfile,cfile))
		mak.write("\n")
	mak.write("\n")
	
	mak.write("#===============================================================\n")
	mak.write("# clean\n")
	mak.write("#===============================================================\n")
	mak.write("clean: \n")
	mak.write("\t@del \"$(ODIR)\*.o\"\n")
	mak.write("\t@del \"$(EDIR)\$(EXEC)\"\n")
	mak.write("\n")
	mak.close()
	
#==============================================================================
# Main function
#==============================================================================
if __name__=='__main__':
	argc = len(sys.argv)
	args=sys.argv
	optlist, args = getopt.getopt(sys.argv[1:],'',
		[
		'with-inc-path=',
		'with-lib-path=',
		'with-lib-name=',
		'with-cxx-comp=',
		'with-cxx-opts='
		])
	cFiles,hPaths                =getProjectFiles()
	hExtra,lPaths,lFiles,compiler=getConfOptions(optlist)
	hPaths                       =hPaths+hExtra
	createMakefile(cFiles,hPaths,lPaths,lFiles,compiler)

