
import os
import string
import node
import element

path="D:\\solver\\nemesis\\bin"

def run(options,file):
	file=os.getcwd().replace("\\","\\\\")+"\\\\"+file
	command=path+"\\tetgen.exe"+" "+options+" "+file+">tetgen.log"
	os.system(command)
	
def getNodes(file,mark=-1):
	ifile=open(file,'r')
	line=ifile.readline()
	tokens=line.split()
	nNodes=int(line.split()[0])
	nAttrs=int(line.split()[2])
	nMarks=int(line.split()[3])
	nodeList=[]
	for i in range(nNodes):
		line=ifile.readline()
		tokens=line.split()
		id=tokens[0]
		x =tokens[1]
		y =tokens[2]
		z =tokens[3]
		m=0
		if nMarks==1 and mark<0:
			m=int(tokens[4+nAttrs])
		elif nMarks==1 and mark>=0:
			m=int(tokens[4+nAttrs])
			if m!=mark: continue
		elif nMarks==0 and mark<0:
			m=0
		else: continue
		nodeList.append([id,x,y,z,m])	
	return nodeList
	
def createNodes(file):
	print "Creating nodes... Start"
	ifile=open(file,'r')
	line=ifile.readline()
	tokens=line.split()
	nNodes=int(line.split()[0])
	nAttrs=int(line.split()[2])
	nMarks=int(line.split()[3])
	for i in range(nNodes):
		line=ifile.readline()
		tokens=line.split()
		id=int(tokens[0])
		x =float(tokens[1])
		y =float(tokens[2])
		z =float(tokens[3])
		#if nMarks==1:
		#	m=int(tokens[4+nAttrs])
		node.add(id,x,y,z)
	print "Creating nodes... Finish"

def createElems(file):
	print "Creating elems... Start"
	ifile=open(file,'r')
	line=ifile.readline()
	tokens=line.split()
	nElems=int(line.split()[0])
	for i in range(nElems):
		line=ifile.readline()
		tokens=line.split()
		id=int(tokens[0])
		n1=int(tokens[1])
		n2=int(tokens[2])
		n3=int(tokens[3])
		n4=int(tokens[4])
		mat=int(tokens[5])
		element.tetra4d(id,n1,n2,n3,n4,mat)
	print "Creating elems... Finish"
def show(file):
	file=os.getcwd().replace("\\","\\\\")+"\\\\"+file
	command=path+"\\tetview.exe"+" "+file
	os.system(command)


