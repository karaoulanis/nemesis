
import os
from time import time,ctime

def setProjectFiles(headers,sources):
	for path, folders, files in os.walk('../..'+'/src'):
		newpath=path.replace("\\","/")
		for name in files:
			if os.path.splitext(name)[1]=='.cpp':
				sources.append([path.replace("\\","/"),name])
			if os.path.splitext(name)[1]=='.h':
				headers.append([path.replace("\\","/"),name])
	# main.cpp should be the first to be compiled
	for index in range(len(sources)):
		if sources[index][1]=='main.cpp':
			break
	if index==len(sources):
		raise StandardError,"No file main.cpp found!"
	tmp=sources[0]
	sources[0]=sources[index]
	sources[index]=tmp

def setExternalFiles(headers,libraries,optlist):
	for i in optlist:
		if   i[0]=='--with-inc-dir':
			print "Configuring with include dir      :",i[1]
			headers.append(["\""+i[1]+"\"",'NULL'])
		elif i[0]=='--with-inc-nam':
			print "Configuring with include file     :",i[1]
			headers.append(['NULL',i[1]])
		elif i[0]=='--with-lib-dir':
			print "Configuring with library dir      :",i[1]
			libraries.append(["\""+i[1]+"\"",'NULL'])
		elif i[0]=='--with-lib-nam':
			print "Configuring with library file     :",i[1]
			libraries.append(['NULL',i[1]])
		else:
			continue

def setCompiler(compiler,optlist):
	for i in optlist:
		if   i[0]=='--with-cxx-cmp':
			print "Configuring with compiler         :",i[1]
			compiler[0]=i[1]
		elif i[0]=='--with-cxx-opt':
			print "Configuring with compiler options :",i[1]
			compiler[1]=i[1]
		else:
			continue
	if compiler[0]=='NULL':
		compiler[0]='g++'
	if compiler[1]=='NULL':
		compiler[1]='-O2 -fWall'

def createMakefile(headers,sources,libraries,compiler):
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
	for i in libraries:
		if i[0]!='NULL' and i[0] not in dirs:
			dirs.append(i[0])
	for i in dirs:
			mak.write(' \\\n\t-I%s'%i)
	mak.write("\n\n")

	mak.write("#===============================================================\n")
	mak.write("# libraries\n")
	mak.write("#===============================================================\n")
	mak.write("LIBS = ")
	for i in libraries:
		if i[1]!='NULL':
			lib=i[1].replace('lib','-l').replace('.a','')
			mak.write(' \\\n\t%s'%(lib))
	mak.write("\n\n")

	mak.write("#===============================================================\n")
	mak.write("# include directories\n")
	mak.write("#===============================================================\n")
	mak.write("INCS = ")
	dirs=[]
	for i in headers:
		if i[0]!='NULL' and i[0] not in dirs:
			dirs.append(i[0])
	for i in dirs:
		mak.write(' \\\n\t-I%s'%i)
	mak.write("\n\n")

	mak.write("#===============================================================\n")
	mak.write("# header files\n")
	mak.write("#===============================================================\n")
	mak.write("HDRS = ")
	for i in headers:
		if i[1]!='NULL':
			mak.write(' \\\n\t%s/%s'%(i[0],i[1]))
	mak.write("\n\n")

	mak.write("#===============================================================\n")
	mak.write("# object files\n")
	mak.write("#===============================================================\n")
	mak.write("OBJS = ")
	for i in sources:
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
	for i in sources:
		k=k+1
		cfile=i[1].replace('.cpp','')
		cdir =i[0]
		mak.write('$(ODIR)/%s.o: %s/%s.cpp $(HDRS)\n'%(cfile,cdir,cfile))
		mak.write('\t@echo Compiling [% 4i/% 4i]: %s.cpp\n'%(k,len(sources),cfile))
		mak.write('\t@$(CC) $(CCFLAGS) $(INCS) -c %s/%s.cpp -o $(ODIR)/%s.o\n'%(cdir,cfile,cfile))
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
	


import sys,getopt,os



if __name__=='__main__':
	# lists	
	headers=[]
	sources=[]
	libraries=[]
	compiler=['NULL','NULL']
	#
	argc = len(sys.argv)
	args=sys.argv
	optlist, args = getopt.getopt(sys.argv[1:],'',
		[
		'with-inc-dir=',
		'with-inc-nam=',
		'with-lib-dir=',
		'with-lib-nam=',
		'with-cxx-cmp=',
		'with-cxx-opt='
		])
	setProjectFiles(headers,sources)
	setExternalFiles(headers,libraries,optlist)
	setCompiler(compiler,optlist)
	createMakefile(headers,sources,libraries,compiler)
