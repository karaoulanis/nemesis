
import os
import string
import node
import element
import math
import load

trigenPath="D:\\solver\\nemesis\\bin"
verbose='on'
tol=0.00000000000001
file=''
loads=[]

def run(options):
	command=trigenPath+"\\trigen.exe"+" "+options+" "+file+".poly >trigen.log"
	os.system(command)
def getNodes(nodeA,nodeB):
	###########################################
	# Get all nodes from file
	###########################################
	ifile=open(file+'.1.node','r')
	allNodes=[]
	line=ifile.readline()
	tokens=line.split()
	nNodes=int(line.split()[0])
	for i in range(nNodes):
		line=ifile.readline()
		tokens=line.split()
		id=int(tokens[0])
		x=float(tokens[1])
		y=float(tokens[2])
		allNodes.append([id,x,y])
	ifile.close()
	###########################################
	# Get nA and nB
	###########################################
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
	###########################################
	# Find nodes
	###########################################
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
	###########################################
	# Sort
	###########################################
	if xB!=xA: 	sortBy=1
	else:		sortBy=2
	# Depending on the sort by x or y, bring that to 0 position.
	for i in range(len(outNodes)):
		tmp=outNodes[i][0]
		outNodes[i][0]=outNodes[i][sortBy]
		outNodes[i][sortBy]=tmp
	# Now sort.
	outNodes.sort()
	# Bring id in the 0 position again.
	for i in range(len(outNodes)):
		tmp=outNodes[i][0]
		outNodes[i][0]=outNodes[i][sortBy]
		outNodes[i][sortBy]=tmp
	# If nA less than nB reverse order.
	if outNodes[0][0]!=allNodes[nodeA][0]:
		outNodes.reverse()
	# Check for errors
	if outNodes[0][0]!=allNodes[nodeA][0] or outNodes[len(outNodes)-1][0]!=allNodes[nodeB][0]:
		print outNodes
		raise IOError,"Undefined error in handling nodes."
	# Return.
	return outNodes
def reportNodes(nodes):
	print "     id          x          y"
	print " ----------------------------"
	for i in range(len(nodes)):
		print ' %6i % 10.4f % 10.4f'%(nodes[i][0],nodes[i][1],nodes[i][2])
def report(str):
	if verbose=='on':
		print str
#
#
#
def createNodes():
	ifile=open(file+'.1.node','r')
	line=ifile.readline()
	tokens=line.split()
	nNodes=int(line.split()[0])
	for i in range(nNodes):
		line=ifile.readline()
		tokens=line.split()
		id=int(tokens[0])
		x=float(tokens[1])
		y=float(tokens[2])
		node.add(id,x,y)
		if verbose=='on':
			print 'msh: Created node : %i %8.4f %8.4f'%(id,x,y)
	ifile.close()
	print 'msh: Created %i nodes.'%nNodes
	return nNodes
#
#
#
def createElems():
	ifile=open(file+'.1.ele','r')
	line=ifile.readline()
	tokens=line.split()
	nElems=int(tokens[0])
	elType=int(tokens[1])
	for i in range(nElems):
		line=ifile.readline()
		tokens=line.split()
		id=int(tokens[0])
		if(elType==3):
			n1=int(tokens[1])
			n2=int(tokens[2])
			n3=int(tokens[3])
			mat=int(tokens[4])
			element.tria3d(id,n1,n2,n3,mat)
		else:
			n1=int(tokens[1])
			n2=int(tokens[2])
			n3=int(tokens[3])
			n4=int(tokens[6])
			n5=int(tokens[4])
			n6=int(tokens[5])
			mat=int(tokens[7])
			element.tria6d(id,n1,n2,n3,n4,n5,n6,mat)
	ifile.close()
	print 'msh: Created %i elements.'%nElems
	return nElems
#
#
#
def mesh(options):
	run(options)
	createNodes()
	createElems()
#
#
#
def fix(nodeA,nodeB,fix):
	nodes=getNodes(nodeA,nodeB)
	nNodes=len(nodes)
	for i in range(nNodes):
		node.fix(nodes[i][0],fix)
		if verbose=='on':
			print 'msh: Fixed node %i dof %i.'%(nodes[i][0],fix)
	print 'msh: Fixed %i nodes.'%nNodes
	return nNodes
#
#
#
def fixX(nodeA,nodeB):
	fix(nodeA,nodeB,1)
#
#
#
def fixY(nodeA,nodeB):
	fix(nodeA,nodeB,2)
#
#
#
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
	nodes=getNodes(nodeA,nodeB)
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
#
#
#
def preLoad(id,dof,pA,pB):
	lineLoad(loads[id-1][1],loads[id-1][2],dof,pA,pB)
#
#
#
def importDxf(filename,options):
	nodes=[]
	lines=[]
	holes=[]
	regions=[]
	fixesX=[]
	fixesY=[]
	#======================================================================
	# Open dxf
	#======================================================================
	ifile=open(filename+'.dxf','r')
	while 1:
		words=ifile.readline().split()
		if len(words)==0:
			continue
		if words[0]=="EOF":
			break
		#==============================================================
		# Read nodes and lines
		#==============================================================
		elif words[0]=="LINE":
			layer=""
			x1=0; y1=0; z1=0
			x2=0; y2=0; z2=0
			while 1:
				words=ifile.readline().split()
				if words[0]=="8":
					layer=ifile.readline().split()[0]
				elif words[0]=="10":
					x1=float('%12.5f'%float(ifile.readline()))
				elif words[0]=="20":
					y1=float('%12.5f'%float(ifile.readline()))
				elif words[0]=="30":
					z1=float('%12.5f'%float(ifile.readline()))
				elif words[0]=="11":
					x2=float('%12.5f'%float(ifile.readline()))
				elif words[0]=="21":
					y2=float('%12.5f'%float(ifile.readline()))
				elif words[0]=="31":
					z2=float('%12.5f'%float(ifile.readline()))
					break
			if [x1,y1,z1] not in nodes:
				nodes.append([x1,y1,z1])
			if [x2,y2,z2] not in nodes:
				nodes.append([x2,y2,z2])
			n1=nodes.index([x1,y1,z1])+1
			n2=nodes.index([x2,y2,z2])+1
			if layer=="geometry":
				lines.append([n1,n2])
			elif layer=="fixX":
				fixesX.append([n1,n2])
			elif layer=="fixY":
				fixesY.append([n1,n2])
			elif layer=="load_1":
				loads.append([1,n1,n2])
			elif layer=="load_2":
				loads.append([2,n1,n2])
			elif layer=="load_3":
				loads.append([3,n1,n2])
			elif layer=="load_4":
				loads.append([4,n1,n2])
			elif layer=="load_5":
				loads.append([5,n1,n2])
			elif layer=="load_6":
				loads.append([6,n1,n2])
			elif layer=="load_7":
				loads.append([7,n1,n2])
			elif layer=="load_8":
				loads.append([8,n1,n2])
			elif layer=="load_9":
				loads.append([9,n1,n2])
		#==============================================================
		# Read holes and regions
		#==============================================================
		elif words[0]=="TEXT":
			layer=""
			x=0; y=0; z=0
			n=""
			while 1:
				words=ifile.readline().split()
				if words[0]=="8":
					layer=ifile.readline().split()[0]
				elif words[0]=="10":
					x=float('%12.5f'%float(ifile.readline()))
				elif words[0]=="20":
					y=float('%12.5f'%float(ifile.readline()))
				elif words[0]=="30":
					z=float('%12.5f'%float(ifile.readline()))
				elif words[0]=="1":
					n=ifile.readline().split()[0]
					break
			if layer=="holes":
				holes.append([x,y])
			elif layer=="regions":
				regions.append([x,y,int(n)])
	ifile.close()
	#======================================================================
	# Save to poly
	#======================================================================
	ofile=open(filename+".poly",'w')
	ofile.write('%i 2 0 0\n'%len(nodes))
	# Nodes
	for i in range(len(nodes)):
		id=i+1
		x=nodes[i][0]
		y=nodes[i][1]
		ofile.write('%4i %12.5f %12.5f \n'%(id,x,y))
		report('msh: Added node  : %4i %12.5f %12.5f '%(id,x,y))
	# Lines
	ofile.write('\n')
	ofile.write('%i 0\n'%len(lines))
	for i in range(len(lines)):
		id=i+1
		n1=lines[i][0]
		n2=lines[i][1]
		ofile.write('%4i %4i %4i \n'%(id,n1,n2))
		report('msh: Added line  : %4i %4i %4i '%(id,n1,n2))
	# Holes
	ofile.write('\n')
	ofile.write('%i\n'%len(holes))
	for i in range(len(holes)):
		id=i+1
		x=holes[i][0]
		y=holes[i][1]
		ofile.write('%4i %12.5f %12.5f \n'%(id,x,y))
		report('msh: Added hole  : %4i %12.5f %12.5f'%(id,x,y))
	# Regions
	ofile.write('\n')
	ofile.write('%i\n'%len(regions))
	for i in range(len(regions)):
		id=i+1
		x=regions[i][0]
		y=regions[i][1]
		n=regions[i][2]
		ofile.write('%4i %12.5f %12.5f %12.5f 1.0\n'%(id,x,y,n))
		report('msh: Added region: %4i %12.5f %12.5f %12.5f 1.0'%(id,x,y,n))
	# Close file
	ofile.close()
	#======================================================================
	# Run mesh generator
	#======================================================================
	run(options)
	createNodes()
	createElems()
	#======================================================================
	# Create fixes
	#======================================================================
	for i in range(len(fixesX)):
		fixX(fixesX[i][0],fixesX[i][1])
	for i in range(len(fixesY)):
		fixY(fixesY[i][0],fixesY[i][1])
#
#
#
def exportToVtk():
	try:
		nodeFile = open(file+".1.node",'r')
		elemFile = open(file+".1.ele",'r')
		vtkFile   = open(file+".1.vtk",'w+')
	except IOError:
		print "Error while reading file!"
		sys.exit(0)

	#################################################
	# Read nodes
	#################################################
	line=nodeFile.readline()
	words=string.split(line)
	nNodes=int(words[0])
	nodes=[]
	for i in range(nNodes): 
		line=nodeFile.readline()
		words=string.split(line)
		nodes.append([int(words[0]),float(words[1]),float(words[2])])
	nodeFile.close()

	#################################################
	# Read elems
	#################################################
	line=elemFile.readline()
	words=string.split(line)
	nElems=int(words[0])
	elType=int(words[1])
	elems=[]
	for i in range(nElems): 
		line=elemFile.readline()
		words=string.split(line)
		if(elType==3):
			elems.append([int(words[0]),int(words[1]),int(words[2]),int(words[3])])
		else:
			elems.append([int(words[0]),int(words[1]),int(words[2]),int(words[3]),int(words[4]),int(words[5]),int(words[6])])
	elemFile.close()

	#################################################
	# Write vtk
	#################################################
	vtkFile.write('# vtk DataFile Version 3.0		\n')
	vtkFile.write('nemesis vtk output				\n')
	vtkFile.write('ASCII						\n')
	vtkFile.write('DATASET UNSTRUCTURED_GRID		\n')
	vtkFile.write('\n')
	vtkFile.write('POINTS %i FLOAT\n'%(nNodes))
	for i in range(nNodes):
		vtkFile.write('% 14.6f % 14.6f % 14.6f\n'%(nodes[i][1],nodes[i][2],0.))
	vtkFile.write('\n')
	if(elType==3):
		vtkFile.write('CELLS %i %i\n'%(nElems,4*nElems))
		for i in range(nElems):
			vtkFile.write('3 %i %i %i\n'%(elems[i][1]-1,elems[i][2]-1,elems[i][3]-1))
	else:
		vtkFile.write('CELLS %i %i\n'%(nElems,7*nElems))
		for i in range(nElems):
			vtkFile.write('6 %i %i %i %i %i %i\n'%(elems[i][1]-1,elems[i][2]-1,elems[i][3]-1,elems[i][6]-1,elems[i][4]-1,elems[i][5]-1))
	vtkFile.write('\n')
	vtkFile.write('CELL_TYPES %i \n'%(nElems))
	for i in range(nElems):
		if(elType==3):
			vtkFile.write('5\n')
		else:
			vtkFile.write('22\n')
	vtkFile.write('\n')
	vtkFile.write('CELL_DATA %i \n'%(nElems))
	vtkFile.write('SCALARS materials int\n')
	vtkFile.write('LOOKUP_TABLE DEFAULT\n')
	for i in range(nElems):
		vtkFile.write('%d\n'%1.0)
	vtkFile.close()