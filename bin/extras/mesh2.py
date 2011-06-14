import node
import element
import group

def struct(Lx,Ly,nElemsX,nElemsY):
	nNodesX=nElemsX+1
	nNodesY=nElemsY+1
	dx=float(Lx)/float(nElemsX)
	dy=float(Ly)/float(nElemsY)

	# Nodes
	id=0
	for i in range(0,nElemsY+1):
		for j in range(0,nElemsX+1):
			id=id+1
			x=j*dx
			y=Ly-i*dy
			#print 'node.add(%3i, %4f, %4f)'%(id,x,y)
			node.add(id,x,y)
	print '%12s %6i %s'%("mesh.struct :",(nElemsY+1)*(nElemsX+1),"nodes created.")
	# Elements
	for j in range(0,nElemsY):
		for i in range(0,nElemsX):
			id=j*nNodesX+i+1
			n1=j*nNodesX+i+1
			n2=i+(j+1)*nNodesX+1
			n3=i+(j+1)*nNodesX+2
			n4=j*nNodesX+i+2
			#print 'element.quad4StdDisp(%3i,%3i,%3i,%3i,%3i, 1, 1.)'%(id,n1,n2,n3,n4)
			element.quad4StdDisp(id,n1,n2,n3,n4,1,1.)
	print '%12s %6i %s'%("mesh.struct :",nElemsY*nElemsX,"elements created.")

def fixX(nStart,nEnd,dn):
	k=0
	for i in range(nStart,nEnd+1,dn):
		#print 'node.fix(%3i, 1, 0.)'%i;
		node.fix(i,1,0.)
		k+=1
	print '%12s %6i %s'%("mesh.fixX   :",k,"constraints created.")

def fixY(nStart,nEnd,dn):
	k=0
	for i in range(nStart,nEnd+1,dn):
		#print 'node.fix(%3i, 2, 0.)'%i;
		node.fix(i,2,0.)
		k+=1
	print '%12s %6i %s'%("mesh.fixY   :",k,"constraints created.")

def standard(Lx,Ly,nElemsX,nElemsY):
	struct(Lx,Ly,nElemsX,nElemsY)
	nNodesX=nElemsX+1
	nNodesY=nElemsY+1
	fixX(                    1, nNodesX*(nNodesY-2)+1, nNodesX) # left   - x
	fixX(nNodesX,               nNodesX*(nNodesY-1),   nNodesX) # right  - x
	fixX(nNodesX*(nNodesY-1)+1, nNodesX*nNodesY,       1) # lower - x
	fixY(nNodesX*(nNodesY-1)+1, nNodesX*nNodesY,       1) # lower - y	