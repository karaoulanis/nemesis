import math
import node
import load
import domain

def line(nA,nB,dn,dof,pA,pB):
	d13=0.333333333333333333333333333333
	d16=0.166666666666666666666666666667
	# Find grad
	
	x1=node.data(nA)['crds'][0]
	y1=node.data(nA)['crds'][1]
	x2=node.data(nB)['crds'][0]
	y2=node.data(nB)['crds'][1]
	L=math.sqrt((y2-y1)**2+(x2-x1)**2)
	grad=(pB-pA)/L
	p1=pA
	for i in range(nA,nB,dn):
		n1=i
		n2=i+dn
		x1=node.data(n1)['crds'][0]
		y1=node.data(n1)['crds'][1]
		x2=node.data(n2)['crds'][0]
		y2=node.data(n2)['crds'][1]
		dL=math.sqrt((y2-y1)**2+(x2-x1)**2)
		p2=p1+grad*dL
		if domain.type()==21 or domain.type()==22:
			load.node(n1,dof,-d13*p1*x1+d13*p1*x2-d16*p2*x1+d16*p2*x2)
			load.node(n2,dof,-d16*p1*x1+d16*p1*x2-d13*p2*x1+d13*p2*x2)
		elif domain.type()==23:
			load.node(n1,dof,(3*x1*x1*p1+x1*x1*p2-x2*x2*p1-  x2*x2*p2-2*x1*x2*p1)/(-12.))
			load.node(n2,dof,(  x1*x1*p1+x1*x1*p2-x2*x2*p1-3*x2*x2*p2+2*x1*x2*p2)/(-12.))
		else:
			exit(0)
		p1=p2
#~ def line(nA,nB,dn,dof,pA,pB):
	#~ dL=[]
	#~ DL=0
	#~ pI=[]
	#~ # Find lengths
	#~ for i in range(nA,nB,dn):
		#~ x1=node.data(i)[2]
		#~ y1=node.data(i)[3]
		#~ x2=node.data(i+1)[2]
		#~ y2=node.data(i+1)[3]
		#~ dL.append(math.sqrt((y2-y1)**2+(x2-x1)**2))
		#~ DL+=dL[len(dL)-1]
	#~ # Find load increments
	#~ grad=(pB-pA)/DL
	#~ pI.append(pA)
	#~ k=0
	#~ for i in range(nA+1,nB,dn):
		#~ pI.append(pI[k]+grad*dL[k])
		#~ k+=1
	#~ pI.append(pB)
	#~ # Find nodal loads
	#~ k=0
	#~ for i in range(nA,nB,dn):
		#~ load.node(i,dof,dL[k]*(0.5*pI[k]+0.1666666666666667*(pI[k+1]-pI[k])))
		#~ load.node(i+1,dof,dL[k]*(0.5*pI[k]+0.3333333333333333*(pI[k+1]-pI[k])))
		#~ k+=1