
import os
import string
import node
import element

path="D:\\solver\\nemesis\\bin"

def run(options,file):
	file=os.getcwd().replace("\\","\\\\")+"\\\\"+file
	command=path+"\\trigen.exe"+" "+options+" "+file+">trigen.log"
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
		m=0
		if nMarks==1 and mark<0:
			m=int(tokens[3+nAttrs])
		elif nMarks==1 and mark>=0:
			m=int(tokens[3+nAttrs])
			if m!=mark: continue
		elif nMarks==0 and mark<0:
			m=0
		else: continue
		nodeList.append([id,x,y,m])	
		ifile.close()
	return nodeList

def getEdges(file,flag):
	ifile=open(file,'r')
	line=ifile.readline()
	tokens=line.split()
	nEdges=int(line.split()[0])
	edgeList=[]
	for i in range(nEdges):
		line=ifile.readline()
		tokens=line.split()
		id=tokens[0]
		n1=tokens[1]
		n2=tokens[2]
		m=int(tokens[3])
		if m==flag:
			edgeList.append([int(n1),int(n2)])	
	ifile.close()
	return edgeList

def createNodes(file):
	ifile=open(file,'r')
	line=ifile.readline()
	tokens=line.split()
	nNodes=int(line.split()[0])
	for i in range(nNodes):
		line=ifile.readline()
		tokens=line.split()
		id=int(tokens[0])
		x =float(tokens[1])
		y =float(tokens[2])
		node.add(id,x,y)
	ifile.close()
	return nNodes

def createElems(file):
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
		mat=int(tokens[4])
		if mat==0: mat=100
		element.tria3d(id,n1,n2,n3,mat)
	ifile.close()
	return nElems
	
def createFixes(file):
	ifile=open(file,'r')
	line=ifile.readline()
	tokens=line.split()
	nNodes=int(line.split()[0])
	nFixes=0
	for i in range(nNodes):
		line=ifile.readline()
		tokens=line.split()
		id=int(tokens[0])
		mark=int(tokens[3])
		if mark==1 or mark==2:
			node.fix(id,mark)
			nFixes=nFixes+1
		if mark==3:
			node.fix(id,1)
			node.fix(id,2)
			nFixes=nFixes+2
	ifile.close()
	return nFixes
