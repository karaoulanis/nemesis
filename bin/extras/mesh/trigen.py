
import os
import node
import element
import math
import load

#==============================================================================
# Variables 
#==============================================================================
path  = "D:\\solver\\nemesis\\bin"
file  =''
info  = 1
tol   = 0.0001
bmId  =  100000
brId  = 1000000
lines =[]
nodes =[]
attrs =[]

#==============================================================================
# Report 
#==============================================================================
def report(str,flag):
	if flag<=info:
		print str
		
#==============================================================================
# Report 
#==============================================================================
def getNodesBetween(nodeA,nodeB):
	# Get all nodes from file
	allNodes=[]
	ifile=open(file+'.1.node','r')
	nNodes=int(ifile.readline().split()[0])
	for i in range(nNodes):
		w=ifile.readline().split()
		allNodes.append([int(w[0]),float(w[1]),float(w[2])])
	ifile.close()
	# Get nA and nB
	nodeA=nodeA-1	# zero index arrays
	nodeB=nodeB-1	# zero index arrays
	xA=allNodes[nodeA][1]
	yA=allNodes[nodeA][2]
	xB=allNodes[nodeB][1]
	yB=allNodes[nodeB][2]
	xL=min(xA,xB)	
	xU=max(xA,xB)
	yL=min(yA,yB)
	yU=max(yA,yB)
	# Find nodes
	outNodes=[]
	if xB!=xA:
		a=(yB-yA)/(xB-xA)
		for i in range(nNodes):
			id = allNodes[i][0]
			x  = allNodes[i][1]
			y  = allNodes[i][2]
			if abs(y-yA-a*(x-xA))<tol and x>=xL and x<=xU:
				outNodes.append([id,x,y])
	else:
		for i in range(nNodes):
			id = allNodes[i][0]
			x  = allNodes[i][1]
			y  = allNodes[i][2]
			if x==xA and y>=yL and y<=yU:
				outNodes.append([id,x,y])
	# Sort
	if xB!=xA: 	sortBy=1
	else:		sortBy=2
	# Depending on the sort by x or y, bring that to 0 position
	for i in range(len(outNodes)):
		tmp=outNodes[i][0]
		outNodes[i][0]=outNodes[i][sortBy]
		outNodes[i][sortBy]=tmp
	# Now sort
	outNodes.sort()
	# Bring id in the 0 position again
	for i in range(len(outNodes)):
		tmp=outNodes[i][0]
		outNodes[i][0]=outNodes[i][sortBy]
		outNodes[i][sortBy]=tmp
	# If nA less than nB reverse order
	if outNodes[0][0]!=allNodes[nodeA][0]:
		outNodes.reverse()
	# Check for errors
	if outNodes[0][0]!=allNodes[nodeA][0] or outNodes[len(outNodes)-1][0]!=allNodes[nodeB][0]:
		print outNodes
		raise IOError,"Undefined error in handling nodes."
	# Return.
	return outNodes

#==============================================================================
# Find if point on line 
#==============================================================================
def isPointOnLine(x,y,x1,y1,x2,y2):
	ret=0
	if x2!=x1:
		xL=min(x1,x2)
		xU=max(x1,x2)
		a=(y2-y1)/(x2-x1)
		if abs(y-y1-((y2-y1)/(x2-x1))*(x-x1))<tol and x>=xL and x<=xU:
			ret=1
	else:
		yL=min(y1,y2)
		yU=max(y1,y2)
		if x==x1 and y>=yL and y<=yU:
			ret=1
	return ret

#==============================================================================
# Find if point on line 
#==============================================================================
def findMiddleNode(n1,n2,nodeCrds):
	x=0.5*(nodeCrds[n1-1][0]+nodeCrds[n2-1][0])
	y=0.5*(nodeCrds[n1-1][1]+nodeCrds[n2-1][1])
	found=0
	for i in range(len(nodeCrds)):
		if abs(x-nodeCrds[i][0])<tol and abs(y-nodeCrds[i][1])<tol:
			n3=i+1
			found=1
			break
	if found==0:
		raise StandardError, "Cannot find internal node!"
	return n3

#==============================================================================
# Find section for bar/beam
#==============================================================================
def findSection(n1,n2):
	sec=0
	x1,y1,x2,y2=nodes[n1-1][0],nodes[n1-1][1],nodes[n2-1][0],nodes[n2-1][1]
	for s in attrs:
		if s[0]=='section':
			x0,y0=s[1],s[2]
			if isPointOnLine(x0,y0,x1,y1,x2,y2)==1:
				sec=int(s[3])
				break
	if sec==0:
		raise StandardError, "Cannot find beam/bar section!"
	return sec

#==============================================================================
# Get spectific objects from lines/attrs
#==============================================================================
def getEntities(myList,myTypes):
	entities=[]
	for i in myList:
		if i[0] in myTypes:
			entities.append(i)
	return entities

#==============================================================================
# Read dxf block 
#==============================================================================
def readBlock(pairs):
	i=0
	block=pairs.pop(0)[1]
	data={}
	while len(pairs)>0 and pairs[i][0]!=0:
		p=pairs.pop(0)
		data[p[0]]=p[1]
	return block,data
	
#==============================================================================
# Parse dxf file 
#==============================================================================
def iDxf():
	pairs=[]
	# Open file and start reading
	ifile=open(file+'.dxf','r')
	lineA=ifile.readline().strip()
	lineB=ifile.readline().strip()
	# Skip the file until entities are found 
	while lineB!="ENTITIES":
		lineA=ifile.readline()
		lineB=ifile.readline().strip()
	# Read until ENDSEC is found and keep data in pairs
	while lineB!="ENDSEC":
		lineA=int(ifile.readline().strip())
		lineB=      ifile.readline().strip()
		pairs.append([int(lineA),lineB])
	# Read blocks of interest
	while len(pairs)>0:
		b,d=readBlock(pairs)
		# Read LINE and store [layer,x1,y1,x2,y2,marker] in lines
		if b=='LINE' and d[8] in ['geom','bar','beam','fixX','fixY','fixXY']:
			lines.append([d[8],		\
			[float(d[10]),float(d[20])],	\
			[float(d[11]),float(d[21])],-1])
		# Read INSERT and store [layer,id,x,y] in attrs
		elif b=='INSERT' and d[2] in ['node','region','hole','section']:
			id=readBlock(pairs)[1][1]
			attrs.append([d[2],float(d[10]),float(d[20]),id])
	# Close file
	ifile.close()
	# Add nodes to nodes list
	geomNodes=getEntities(attrs,['node'])
	for i in range(len(geomNodes)):
		tmp=[int(geomNodes[i][3]),geomNodes[i][1],geomNodes[i][2]]
		geomNodes[i]=tmp
	geomNodes.sort()
	for n in geomNodes:
		nodes.append([n[1],n[2]])
	# Add nodal coordinates to nodes list
	for l in lines:
		if l[1] not in nodes:
			nodes.append(l[1])
		if l[2] not in nodes:
			nodes.append(l[2])
	# Replace coordinates by nodal ids and find beam/bar sections
	for l in lines:
		l[1]=nodes.index(l[1])+1
		l[2]=nodes.index(l[2])+1
		if l[0] in ['bar','beam']:
			l[3]=findSection(l[1],l[2])

#==============================================================================
# Write poly file 
#==============================================================================
def writePoly():
	ofile=open(file+".poly",'w')
	# Nodes
	k=0
	ofile.write('%i 2 0 0\n'%len(nodes))
	for i in nodes:
		k=k+1
		ofile.write('%4i %12.5f %12.5f \n'%(k,i[0],i[1]))
	report('msh: Added %4i nodes'%len(nodes),1)
	# Lines
	k=0
	array=getEntities(lines,['geom','beam','fixX','fixY','fixXY'])
	ofile.write('\n')
	ofile.write('%i 1\n'%len(array))
	for i in array:
		k=k+1
		ofile.write('%4i %4i %4i %3i\n'%(k,i[1],i[2],i[3]))
	report('msh: Added %4i lines'%len(array),1)
	# Holes
	k=0
	array=getEntities(attrs,['hole'])
	ofile.write('\n')
	ofile.write('%i\n'%len(array))
	for i in array:
		k=k+1
		ofile.write('%4i %12.5f %12.5f \n'%(k,i[1],i[2]))
	report('msh: Added %4i holes'%len(array),1)
	# Regions
	k=0
	array=getEntities(attrs,['region'])
	ofile.write('\n')
	ofile.write('%i\n'%len(array))
	for i in array:
		k=k+1
		ofile.write('%4i %12.5f %12.5f %12.5f 1.0\n'%(k,i[1],i[2],float(i[3])))
	report('msh: Added %4i regions'%len(array),1)
	# Close file
	ofile.close()

#==============================================================================
# Create nodes
#==============================================================================
def createNodes():
	ifile=open(file+'.1.node','r')
	nNodes=int(ifile.readline().split()[0])
	for i in range(nNodes):
		w=ifile.readline().split()
		node.add(int(w[0]),float(w[1]),float(w[2]))
	ifile.close()
	report('msh: Created %i nodes.'%nNodes,1)
	return nNodes

#==============================================================================
# Create triangles
#==============================================================================
def createTrigs():
	ifile=open(file+'.1.ele','r')
	nElems,elType=map(int,ifile.readline().split()[0:2])
	for i in range(nElems):
		w=map(int,ifile.readline().split())
		if(elType==3):
			element.tria3d(w[0],w[1],w[2],w[3],w[4])
		else:
			element.tria6d(w[0],w[1],w[2],w[3],w[6],w[4],w[5],w[7])
	ifile.close()
	report('msh: Created %i triangles.'%nElems,1)
	return nElems

#==============================================================================
# Create Beams
#==============================================================================
def createBeams():
	# Get triangle type (3-node or 6-node)
	ifile=open(file+'.1.ele','r')
	elType=int(ifile.readline().split()[1])
	ifile.close()
	# if triangle type=6 then the middle node needs to be found
	nodeCrds=[]
	if(elType==6):
		ifile=open(file+'.1.node','r')
		for i in range(int(ifile.readline().split()[0])):
			nodeCrds.append(map(float,ifile.readline().split())[1:3])
		ifile.close()
	# Open edge file
	ifile=open(file+'.1.edge','r')
	nElems=0
	for i in range(int(ifile.readline().split()[0])):
		id,n1,n2,sec=map(int,ifile.readline().split())
		id+=bmId
		if sec>0:
			if elType==3:
				element.beam2t(id,n1,n2,sec,sec,1)
			else:
				n3=findMiddleNode(n1,n2,nodeCrds)
				element.beam3t(id,n1,n2,n3,sec,sec,2)
			nElems=nElems+1
	ifile.close()
	report('msh: Created %i beams.'%nElems,1)
	return nElems

#==============================================================================
# Create Beams
#==============================================================================
def createBars():
	k=0
	array=getEntities(lines,['bar'])
	for i in array:
		k+=1
		element.bar2s(k+brId,i[1],i[2],i[3],i[3])
	report('msh: Created %i bars.'%len(array),1)
	return k

#==============================================================================
# Create Fixes
#==============================================================================
def createFixes():
	k=0
	array=getEntities(lines,['fixX','fixY','fixXY'])
	for l in array:
		nodes=getNodesBetween(l[1],l[2])
		for n in nodes:
			if l[0]=='fixX':
				node.fix(n[0],1)
				k+=1
			elif l[0]=='fixY':
				node.fix(n[0],2)
				k+=1
			elif l[0]=='fixXY':
				node.fix(n[0],1)
				node.fix(n[0],2)
				k+=1
	report('msh: Fixed %i nodes.'%k,1)
	return k
#==============================================================================
# Define line load
#==============================================================================
def lineLoad(nodeA,nodeB,dof,pA,pB):
	d13=0.333333333333333333333333333333
	d16=0.166666666666666666666666666667
	# Determine type of elements
	ifile=open(file+'.1.ele','r')
	line=ifile.readline()
	tokens=line.split()
	nElems=int(tokens[0])
	elType=int(tokens[1])
	ifile.close()
	# Get nodes
	nodes=getNodesBetween(nodeA,nodeB)
	# Linear triangles
	if(elType==3):
		for k in range(0,len(nodes)-1):
			n1=nodes[k  ][0]
			x1=nodes[k  ][1]
			y1=nodes[k  ][2]
			n2=nodes[k+1][0]
			x2=nodes[k+1][1]
			y2=nodes[k+1][2]
			L=math.sqrt((y2-y1)**2+(x2-x1)**2)
			p1= d16*(2*pA+   pB)*L
			p2= d16*(   pA+2*pB)*L
			load.node(n1,dof,p1)
			load.node(n2,dof,p2)
	# Quadratic triangles
	if(elType==6):
		if (len(nodes)-3)%2!=0 or len(nodes)<3:
			print len(nodes)
			raise RuntimeError,"Input error"
		for k in range(0,len(nodes)-1,2):
			n1=nodes[k  ][0]
			x1=nodes[k  ][1]
			y1=nodes[k  ][2]
			n2=nodes[k+1][0]
			x2=nodes[k+1][1]
			y2=nodes[k+1][2]
			n3=nodes[k+2][0]
			x3=nodes[k+2][1]
			y3=nodes[k+2][2]
			L=math.sqrt((y3-y1)**2+(x3-x1)**2)
			p1= d16*pA*L
			p2= d13*(pA+pB)*L
			p3= d16*pB*L
			load.node(n1,dof,p1)
			load.node(n2,dof,p2)
			load.node(n3,dof,p3)

#==============================================================================
# Mesh
#==============================================================================
def mesh(type,size):
	writePoly()
	if type==1:
		options='-peAVqa%f'%size
	elif type==2:
		options='-peAVqa%fo2'%size
	command=path+"\\trigen.exe"+" "+options+" "+file+".poly >trigen.log"
	os.system(command)
	if info>1:
		ifile=open('trigen.log')
		for line in ifile.readlines():
			print line.strip()
		ifile.close()
	createNodes()
	createTrigs()
	createBeams()
	createBars()
	createFixes()

