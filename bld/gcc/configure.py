import os
import re
from time import time,ctime

def getFiles(root_path):
	""" Get files/paths from a given directory.
	"""
	print "Searching in the src directory..."
	paths=[]
	names=[]
	for path, folders, files in os.walk(root_path):
		for name in files:
			filename,extension=os.path.splitext(name)
			if extension in ['.cpp','.h']:
				head,tail=os.path.split(path)
				paths.append(head+'/'+tail)
				names.append(name)
	return dict(zip(names, paths))

def dfs(graph, start):
	"""iterative depth first search from start
	http://code.activestate.com/recipes/576723/"""
	#~ print "Finding dependencies for file     :",start
	path=[]
	q=[start]
	while q:
		v=q.pop(0)
		if v not in path:
			path=path+[v]
			q=graph[v]+q
	return path

def getIncludes(filename,names):
	#~ print "Finding included files for file   :",filename
	ifile=open(filename,'r')
	s=ifile.read()
	ifile.close()
	pat=re.compile('\s*\#include\s+["<]([^">]+[.]h)*[">]')
	tmp=pat.findall(s)
	inc=[]
	for i in tmp:
		if i in names:
			inc.append(i)
	return inc

def getGraph(paths):
	print "Creating graph..."
	graph={}
	for name in paths:
		filename=paths[name]+'/'+name
		graph[name]=getIncludes(filename,paths)
	return graph

def plotGraph(graph):
	print "Writting graph..."
	ofile=open("graph.dot",'w')
	ofile.write('digraph G {\n')
	for i in graph:
		for j in graph[i]:
			ofile.write('\t\"%s\" -> \"%s\";\n'%(i,j))
	ofile.write('}\n');
	ofile.close()
	#~ os.system("C:\\PROGRA~1\\Graphviz2.16\\bin\\dot.exe -Tpng graph.dot -o graph.png")


print "Gathering data..."
paths=getFiles('../..'+'/src')
graph=getGraph(paths)
plotGraph(graph)
names=paths.keys()
names.sort()

print "Opening file for writting..."
ofile=open("Makefile",'w')

print "Writing info..."
ofile.write("#==================================================================\n")
ofile.write("# Makefile for nemesis                                          	\n")
ofile.write("# Maintained by : Fotios E. Karaoulanis                       		\n")
ofile.write("# Created by    : configure.py [Automatically]                  	\n")
ofile.write('# Date          : %s\n'%ctime(time()))
ofile.write("# Version       : 2.0                                          	\n")
ofile.write("#==================================================================\n")
ofile.write("\n")

print "Including platfrom options..."
ofile.write("include platform.in												\n")
ofile.write("\n")

incs=[]
for name in names:
	if paths[name] not in incs:
		incs.append(paths[name])
	incs.sort()
print "Writing include files..."
ofile.write("#==================================================================\n")
ofile.write("# include directories												\n")
ofile.write("#==================================================================\n")
ofile.write("INCS = ")
for inc in incs:
	ofile.write('\\\n\t-I%s'%inc)
ofile.write("\n\n")

print "Writing object files..."
ofile.write("#==================================================================\n")
ofile.write("# object files														\n")
ofile.write("#==================================================================\n")
ofile.write("OBJS = ")
for name in names:
	filename,extension=os.path.splitext(name)
	if extension=='.cpp':
		ofile.write('\\\n\t$(ODIR)/%s.o'%filename)
ofile.write("\n\n")

print "Writing link target..."
ofile.write("#==================================================================\n")
ofile.write("# link																\n")
ofile.write("#==================================================================\n")
ofile.write("link: $(OBJS) 														\n")
ofile.write("	@echo Linking...												\n")
ofile.write("	@$(CC) $(CCFLAGS) -o $(EDIR)/$(EXEC) $(OBJS) $(LDIR) $(LIBS)	\n")
ofile.write("\n")

print "Writing compile target..."
ofile.write("#==================================================================\n")
ofile.write("# compile															\n")
ofile.write("#==================================================================\n")
for name in names:
	filename,extension=os.path.splitext(name)
	if extension=='.cpp':
		depends=dfs(graph, name)
		line='$(ODIR)/'+filename+'.o'+':'
		for d in depends:
			line+=' %s/%s'%(paths[d],d)
		line+='\n'
		ofile.write(line)
		line='@echo Compiling : %s\n'%(name)
		ofile.write(line)
		line='@$(CC) $(CCFLAGS) $(INCS) $(INCS_EXT) -c %s/%s -o $(ODIR)/%s.o\n\n'%(paths[name],name,filename)
		ofile.write(line)

print "Makefile created successfully..."
raw_input("Press enter to continue...")