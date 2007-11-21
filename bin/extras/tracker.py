
import node
import element
import pylab
from numpy import array

def getNodePath(nodeID,tag,fac):
	# tags
	tags ={ 'lambda' :[0],
			'time'   :[1],
			'u_x'    :[2,'disp',0],
			'u_y'    :[2,'disp',1],
			'u_z'    :[2,'disp',2],
			'v_x'    :[2,'velc',0],
			'v_y'    :[2,'velc',1],
			'v_z'    :[2,'velc',2],
			'a_x'    :[2,'accl',0],
			'a_y'    :[2,'accl',1],
			'a_z'    :[2,'accl',2]}
	# data
	x=[]
	if tags[tag][0]==0 or tags[tag][0]==1:
		index=tags[tag][0]
		for i in node.path(nodeID):
			x.append(fac*i[index])
	elif tags[tag][0]==2:
		type=tags[tag][1]
		dof =tags[tag][2]
		for i in node.path(nodeID):
			x.append(fac*i[2][type][dof])
	else:
		raise InputError, 'Unknown tag :%s'%tag
	# return
	return x

def plotNodePath(nodeID,xTag,yTag,xFac=1.0,yFac=1.0,filename='NULL'):
	# data
	x=getNodePath(nodeID,xTag,xFac)
	y=getNodePath(nodeID,yTag,yFac)
	# plot
	pylab.plot(x,y,'k-o')
	pylab.xlabel(xTag)
	pylab.ylabel(yTag)
	pylab.grid(True)
	pylab.show()
	# file
	if filename!='NULL':
		ofile=open(filename,'w')
		for i in range(len(x)):
			ofile.write('%13.6f %13.6f\n'%(x[i],y[i]))


def getElemPath(elemID,matIndex,tag,fac):
	# tags
	tags ={ 'lambda'    :[0],
			'time'      :[1],
			'sigma'     :[2,'sigm',-1],
			'eps'       :[2,'epst',-1],
			'p'     	:[2,'epsv',-1],
			'q'     	:[2,'epsv',-1],
			'eps_v'     :[2,'epsv',-1],
			'sigma_x'   :[2,'sigm', 0],
			'sigma_y'   :[2,'sigm', 1],
			'sigma_z'   :[2,'sigm', 2],
			'tau_xy'    :[2,'sigm', 3],
			'tau_xz'    :[2,'sigm', 4],
			'tau_yz'    :[2,'sigm', 5],
			'eps_x'     :[2,'epst', 0],
			'eps_y'     :[2,'epst', 1],
			'eps_z'     :[2,'epst', 2],
			'eps_xy'    :[2,'epst', 3],
			'eps_xz'    :[2,'epst', 4],
			'eps_yz'    :[2,'epst', 5],
			'epsP_x'    :[2,'epsp', 0],
			'epsP_y'    :[2,'epsp', 1],
			'epsP_z'    :[2,'epsp', 2],
			'epsP_xy'   :[2,'epsp', 3],
			'epsP_xz'   :[2,'epsp', 4],
			'epsP_yz'   :[2,'epsp', 5]}
	# data
	x=[]
	if tags[tag][0]==0 or tags[tag][0]==1:
		index=tags[tag][0]
		for i in element.path(elemID,matIndex):
			x.append(fac*i[index])
	elif tags[tag][0]==2:
		type=tags[tag][1]
		index =tags[tag][2]
		if index>=0:
			for i in element.path(elemID,matIndex):
				x.append(fac*i[2][type][index])
		else:
			for i in element.path(elemID,matIndex):
				x.append(fac*i[2][type])
	else:
		raise InputError, 'Unknown tag :%s'%tag
	# return
	return x

def plotElemPath(elemID,matIndex,xTag,yTag,xFac=1.0,yFac=1.0,filename='NULL'):
	# data
	x=getElemPath(elemID,matIndex,xTag,xFac)
	y=getElemPath(elemID,matIndex,yTag,yFac)
	# plot
	pylab.plot(x,y,'k-')
	pylab.xlabel(xTag)
	pylab.ylabel(yTag)
	pylab.grid(True)
	pylab.show()
	# file
	if filename!='NULL':
		ofile=open(filename,'w')
		for i in range(len(x)):
			ofile.write('%13.6f %13.6f\n'%(x[i],y[i]))
