/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 3, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#include <PyParser.h>
#include <SException.h>
#include <sstream>

/******************************************************************************
* Static variables
******************************************************************************/
static Domain*   pD=0;
static Analysis* pA=0;
static int currentLC=0;

PyObject* buildTuple(const Packet& p)
{
	static PyObject* pTuple;
	int pos=0;
	pTuple=PyTuple_New(p.size);
	PyTuple_SetItem(pTuple,pos++,PyInt_FromLong(0));
	PyTuple_SetItem(pTuple,pos++,PyInt_FromLong(0));
	for(int i=0;i<p.dblSize;i++)	PyTuple_SetItem(pTuple,pos++,PyFloat_FromDouble(p.dblArray[i]));
	for(int i=0;i<p.intSize;i++)	PyTuple_SetItem(pTuple,pos++,PyInt_FromLong(p.intArray[i]));
	return pTuple;
}
PyObject* buildList(const Vector& v)
{
	static PyObject* pyList;
	pyList=PyList_New(v.size());
	for(int i=0;i<v.size();i++)
		PyList_SetItem(pyList,i,PyFloat_FromDouble(v[i]));
	return pyList;
}
PyObject* buildList(const Matrix& m)
{
	static PyObject* pyList;
	pyList=PyList_New(m.rows());
	for(int i=0;i<m.rows();i++)
	{
		PyObject* pyRow=PyList_New(m.cols());
		for(int j=0;j<m.cols();j++)
			PyList_SetItem(pyRow,j,PyFloat_FromDouble(m(i,j)));
		PyList_SetItem(pyList,i,pyRow);
	}
	return pyList;
}
static PyObject* buildDict(istream& s)
{
	static PyObject* pyDict;
	pyDict=PyDict_New();
	char name[128];
	int tag;
	s>>name;	// Read label, i.e. NODE
	s>>name;	// Read first tag name
	while(strcmp(name,"END"))
	{
		s>>tag;
		PyObject* pyKey=PyString_FromString(name);
		if(tag==1000)
		{
			int n;
			s>>n;
			PyDict_SetItem(pyDict,pyKey,PyInt_FromLong(n));
		}
		else if(tag==1020)
		{
			double d;
			s>>d;
			PyDict_SetItem(pyDict,pyKey,PyFloat_FromDouble(d));
		}
		else if(tag==1100)
		{
			int n;
			s>>n;
			Vector v(n);
			for(int i=0;i<n;i++) s>>v[i];
			PyDict_SetItem(pyDict,pyKey,buildList(v));
		}
		else if(tag==1200)
		{
			int rows,cols;
			s>>rows;
			s>>cols;
			Matrix m(rows,cols,0.);
			for(int i=0;i<rows*cols;i++) s>>m.data()[i];
			PyDict_SetItem(pyDict,pyKey,buildList(m));
		}
		else 
			throw SException("[nemesis:%d] %s",9999,"Internal error: Unknown tag.");
		s>>name;
	}
	return pyDict;
}
PyObject* buildList(istream& s)
{
	static PyObject* pyList;
	char name[128];
	int tag;
	s>>name;	// TRACKER
	s>>name;	// steps
	s>>tag;		// 1000
	int steps;
	s>>steps;
	pyList=PyList_New(steps);
	for(int i=0;i<steps;i++)
	{
		double lambda,time;
		s>>name;	// lambda
		s>>tag;		// 1010
		s>>lambda;
		s>>name;	// time
		s>>tag;		// 1010
		s>>time;
		s>>name;	// data
		s>>tag;		// 1020
		PyObject* pyRow=PyList_New(3);
		PyList_SetItem(pyRow,0,PyFloat_FromDouble(lambda));
		PyList_SetItem(pyRow,1,PyFloat_FromDouble(time));
		PyList_SetItem(pyRow,2,buildDict(s));
		PyList_SetItem(pyList,i,pyRow);
	}
	return pyList;
}
/******************************************************************************
* Database Commands
******************************************************************************/
static PyObject* pyDatabase_SQLite(PyObject *self, PyObject *args)
{
	const char* s;
	if(!PyArg_ParseTuple(args,"s",&s))	return NULL;
	try
	{
		Database* pDB=new SQLiteDatabase(s);
		pD->setDatabase(pDB);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDatabase_Store(PyObject *self, PyObject *args)
{
 	const char* s;
	if(!PyArg_ParseTuple(args,"s",&s))	return NULL;
	try
	{
		if(pD->getDatabase()==0) throw SException("[nemesis:%d] %s",1110,"No database set.");
		pD->storeState(s);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDatabase_Close(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,""))	return NULL;
	try
	{
		if(pD->getDatabase()==0) throw SException("[nemesis:%d] %s",1110,"No database set.");
		pD->getDatabase()->closeDB();
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDatabase_ExportToVtk(PyObject *self, PyObject *args)
{
 	const char* s;
	if(!PyArg_ParseTuple(args,"s",&s))	return NULL;
	try
	{
		if(pD->getDatabase()==0) throw SException("[nemesis:%d] %s",1110,"No database set yet.");
		pD->getDatabase()->exportToVtk(s);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyMethodDef DatabaseMethods[] = 
{
	{"SQLite",	pyDatabase_SQLite,	
		METH_VARARGS,	"Create and utilize an SQLite database."},
	{"close",	pyDatabase_Close,	
		METH_VARARGS,	"Close current database."},
	{"store",	pyDatabase_Store,	
		METH_VARARGS,	"Store domain's state to the database."},
	{"exportToVtk",	pyDatabase_ExportToVtk,	
		METH_VARARGS,	"Export a stored domain's state to vtk file format."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Domain Commands
******************************************************************************/
static PyObject* pyDomain_Dim(PyObject *self, PyObject *args)
{
	int n;
    if(!PyArg_ParseTuple(args,"i",&n))	return NULL;
	try
	{
		pD->setnDim(n);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDomain_PlaneStress(PyObject *self, PyObject *args)
{
	double t=1.0;
    if(!PyArg_ParseTuple(args,"|d",&t))	return NULL;
	try
	{
		pD->setTag(TAG_DOMAIN_PLANE_STRESS);
		pD->setFac(t);
		pD->setnDim(2);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDomain_PlaneStrain(PyObject *self, PyObject *args)
{
	double t=1.0;
    if(!PyArg_ParseTuple(args,"|d",&t))	return NULL;
	try
	{
		pD->setTag(TAG_DOMAIN_PLANE_STRAIN);
		pD->setFac(t);
		pD->setnDim(2);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDomain_Axisymmetric(PyObject *self, PyObject *args)
{
	double t=1.0;
    if(!PyArg_ParseTuple(args,"|d",&t))	return NULL;
	try
	{
		pD->setTag(TAG_DOMAIN_AXISYMMETRIC);
		pD->setFac(t);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDomain_RayleighDamping(PyObject *self, PyObject *args)
{
	double aK=0;
	double bM=0;
    if(!PyArg_ParseTuple(args,"dd",&aK,&bM))	return NULL;
	Vector facs(2);
	facs[0]=aK;
	facs[1]=bM;
	pD->setRayleighFactors(facs);
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDomain_Clear(PyObject *self, PyObject *args)
{
    if(!PyArg_ParseTuple(args,""))	return NULL;
	pD->clear();
	pA->clear();
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDomain_State(PyObject *self, PyObject *args)
{
	double facD;
    if(!PyArg_ParseTuple(args,"d",&facD))	return NULL;
	pD->state(facD);
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyDomain_EigenValues(PyObject *self, PyObject *args)
{
    if(!PyArg_ParseTuple(args,""))	return NULL;
    return buildList(pD->getEigenValues());
}
static PyObject* pyDomain_Type(PyObject *self, PyObject *args)
{
    if(!PyArg_ParseTuple(args,""))	return NULL;
	return Py_BuildValue("i",pD->getTag());
}
static PyObject* pyDomain_Gravity(PyObject *self, PyObject *args)
{
	double g,xG,yG,zG;
    if(!PyArg_ParseTuple(args,"d(ddd)",&g,&xG,&yG,&zG))	return NULL;
	try
	{
		pD->setGravity(g,xG,yG,zG);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
    return Py_None;
}
static PyMethodDef DomainMethods[] = 
{
	{"dim",		pyDomain_Dim,	
		METH_VARARGS,	"Set domain dimensions."},
	{"planeStress",		pyDomain_PlaneStress,	
		METH_VARARGS,	"Define a plain stress domain."},
	{"planeStrain",		pyDomain_PlaneStrain,	
		METH_VARARGS,	"Define a plain strain domain."},
	{"axisymmetric",		pyDomain_Axisymmetric,	
		METH_VARARGS,	"Define an axisymmetric domain."},
	{"clear",	pyDomain_Clear,	
		METH_VARARGS,	"Clear the domain."},
	{"state",	pyDomain_State,	
		METH_VARARGS,	"Set domain state."},
	{"RayleighDamping",		pyDomain_RayleighDamping,
		METH_VARARGS,	"Set Rayleigh damping factors."},
	{"eigenvalues",		pyDomain_EigenValues,
		METH_VARARGS,	"Return the eigenvalues of the domain."},
	{"type",				pyDomain_Type,	
		METH_VARARGS,	"Domain type."},
	{"gravity",				pyDomain_Gravity,	
		METH_VARARGS,	"Set gravity direction."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Node Commands
******************************************************************************/
static PyObject* pyNode_Add(PyObject *self, PyObject *args)
{
	int id;
	double x1,x2=0,x3=0;
	if(!PyArg_ParseTuple(args,"id|dd",&id,&x1,&x2,&x3))	return NULL;
	try
	{
		Node* pNode=new Node(id,x1,x2,x3);
		pD->add(pD->getNodes(),pNode);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
    return Py_None;
}
static PyObject* pyNode_Fix(PyObject *self, PyObject *args)
{
	int nodeID,dofID;
	double c=0.;
    if(!PyArg_ParseTuple(args,"ii|d",&nodeID,&dofID,&c))	return NULL;
	Constraint* pConstraint=new Constraint;
	try
	{
		pConstraint->setcDof(nodeID,dofID,1.0);
	}
	catch(SException e)
	{
		delete pConstraint;
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	pConstraint->setcVal(c);
	pD->add(pD->getConstraints(),pConstraint);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyNode_Data(PyObject *self, PyObject *args)
{
	int id;
    if(!PyArg_ParseTuple(args,"i",&id))	return NULL;
 	ostringstream os;
	try
	{
		pD->get<Node>(pD->getNodes(),id)->save(os);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	static istringstream is;
	is.str(os.str());
	return buildDict(is);
}
static PyObject* pyNode_Track(PyObject *self, PyObject *args)
{
	int id;
    if(!PyArg_ParseTuple(args,"i",&id))	return NULL;
	try
	{
		pD->get<Node>(pD->getNodes(),id)->addTracker();
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyNode_Path(PyObject *self, PyObject *args)
{
	int id;
    if(!PyArg_ParseTuple(args,"i",&id))	return NULL;
 	ostringstream os;
	try
	{
		pD->get<Node>(pD->getNodes(),id)->getTracker()->save(os);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	istringstream is(os.str());
	return buildList(is);
}
static PyMethodDef NodeMethods[] = 
{
	{"add",		pyNode_Add,		
		METH_VARARGS,	"Adds a node to the domain."},
	{"fix",		pyNode_Fix,		
		METH_VARARGS,	"Fix a nodal degree of freedom."},
	{"data",	pyNode_Data,	
		METH_VARARGS,	"Returns a tuple containing nodal data."},
	{"track",	pyNode_Track,	
		METH_VARARGS,	"Add a Tracker to the Node."},
	{"path",	pyNode_Path,	
		METH_VARARGS,	"Return the data recorded to the tracker."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* CrossSection Commands
******************************************************************************/
static PyObject* pySection_Rect(PyObject *self, PyObject *args)
{
	int id;
	double w;
	double h;
    if(!PyArg_ParseTuple(args,"idd",&id,&w,&h))
		return NULL;
	CrossSection* pSection=
		new RectangularCrossSection(id,w,h);
	pD->add(pD->getCrossSections(),pSection);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pySection_UserDefined(PyObject *self, PyObject *args)
{
	int id;
	double A=0,As2=0,As3=0,J1=0,J2=0,J3=0,h2=0,h3=0;
    if(!PyArg_ParseTuple(args,"id|ddddddd",&id,&A,&As2,&As3,&J1,&J2,&J3,&h2,&h3))
		return NULL;
	CrossSection* pSection=
		new UserDefinedCrossSection(id,A,As2,As3,J1,J2,J3,h2,h3);
	pD->add(pD->getCrossSections(),pSection);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pySection_Rem(PyObject *self, PyObject *args)
{
	return Py_None;
}
static PyObject* pySection_Info(PyObject *self, PyObject *args)
{
	int id;
    if(!PyArg_ParseTuple(args,"i",&id))	return NULL;
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef SectionMethods[] = 
{
	{"user",	pySection_UserDefined,	
		METH_VARARGS,"Adds a user defined section to the domain."},
	{"rect",	pySection_Rect,	
		METH_VARARGS,"Adds a rectangular section to the domain."},
	{"rem",		pySection_Rem, 
		METH_VARARGS,"Removes a section from the domain."},
	{"info",	pySection_Info,
		METH_VARARGS,"Prints section info."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Material Commands
******************************************************************************/
void createGroupByMaterial(int groupId)
{
	try
	{
		Group* pGroup=new Group(groupId);
		pD->add(pD->getGroups(),pGroup);
		pD->setCurrentGroup(groupId);
	}
	catch(SException e)
	{
	}
}
static PyObject* pyMaterial_SDof(PyObject *self, PyObject *args)
{
	int id;
	double E,rho=0;
    if(!PyArg_ParseTuple(args,"id|d",&id,&E,&rho))return NULL;
	Material* pMaterial=new SDofMaterial(id,E,rho);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_SpringElastic(PyObject *self, PyObject *args)
{
	int id;
	double Kn,Ks2=0.,Ks3=0.;
    if(!PyArg_ParseTuple(args,"id|dd",&id,&Kn,&Ks2,&Ks3))return NULL;
	Material* pMaterial=new SpringElastic(id,Kn,Ks2,Ks3);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* pyMaterial_SpringContact(PyObject *self, PyObject *args)
{
	int id;
	double Kn,Ks,mu,gap;
    if(!PyArg_ParseTuple(args,"idddd",&id,&Kn,&Ks,&mu,&gap))return NULL;
	Material* pMaterial=new SpringContact(id,Kn,Ks,mu,gap);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_UniaxialElastic(PyObject *self, PyObject *args)
{
	int id;
	double E=0,nu=0,rho=0,aT=0;
    if(!PyArg_ParseTuple(args,"id|ddd",&id,&E,&nu,&rho,&aT))return NULL;
	Material* pMaterial=new UniaxialElastic(id,E,nu,rho,aT);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_UniaxialElastoPlastic(PyObject *self, PyObject *args)
{
	int id;
	double E=0,nu=0,rho=0,aT=0,sy=0,Hiso=0,Hkin=0,eta=0;
    if(!PyArg_ParseTuple(args,"iddddddd|d",&id,&E,&nu,&rho,&aT,&sy,&Hiso,&Hkin,&eta)) 
		return NULL;
	Material* pMaterial=new UniaxialElastoPlastic(id,E,nu,rho,aT,sy,Hiso,Hkin,eta);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_UniaxialGap(PyObject *self, PyObject *args)
{
	int id;
	double E=0,nu=0,rho=0,aT=0,sy=0,gap=0;
    if(!PyArg_ParseTuple(args,"idddddd",&id,&E,&nu,&rho,&aT,&sy,&gap)) 
		return NULL;
	Material* pMaterial=new UniaxialGap(id,E,nu,rho,aT,sy,gap);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_Elastic(PyObject *self,PyObject *args)
{
	int id;
	double E=0,nu=0,rho=0,aT=0;
    if(!PyArg_ParseTuple(args,"idd|dd",&id,&E,&nu,&rho,&aT)) 
		return NULL;
	Material* pMaterial;
	if(pD->getTag()==TAG_DOMAIN_PLANE_STRESS)
		pMaterial=new PlaneStress(id,E,nu,rho,aT);
	else
		pMaterial=new MultiaxialElastic(id,E,nu,rho,aT);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_DuncanChang(PyObject *self,PyObject *args)
{
	int id;
	double E,nu,c,phi,m,Rf,pa,rho=0,aT=0;
    if(!PyArg_ParseTuple(args,"iddddddd|dd",&id,&E,&nu,&c,&phi,&m,&Rf,&pa,&rho,&aT)) 
		return NULL;
	Material* pMaterial=new DuncanChang(id,E,nu,c,phi,m,Rf,pa,rho,aT);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_UniaxialCyclic(PyObject *self,PyObject *args)
{
	int id;
	double E,nu,rho,aT,tmax,gmax;
    if(!PyArg_ParseTuple(args,"idddddd",&id,&E,&nu,&rho,&aT,&tmax,&gmax)) return NULL;
	Material* pMaterial=new UniaxialCyclic(id,E,nu,rho,aT,tmax,gmax);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_VonMises(PyObject *self,PyObject *args)
{
	int id,elasticId;
	double s0,K;
    if(!PyArg_ParseTuple(args,"iidd",&id,&elasticId,&s0,&K)) 
		return NULL;
	Material* pMaterial=new VonMises(id,elasticId,s0,K);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_MohrCoulomb(PyObject *self,PyObject *args)
{
	int id,elasticId;
	double c,phi,T;
    if(!PyArg_ParseTuple(args,"iiddd",&id,&elasticId,&c,&phi,&T)) 
		return NULL;
	Material* pMaterial=new MohrCoulomb(id,elasticId,c,phi,T);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_Tresca(PyObject *self,PyObject *args)
{
	int id,elasticId;
	double cu;
	double kx=0.,ky=0.,kz=0.;
    if(!PyArg_ParseTuple(args,"iidddd",&id,&elasticId,&cu,&kx,&ky,&kz)) 
		return NULL;
	Material* pMaterial=new Tresca(id,elasticId,cu,kx,ky,kz);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_Creep(PyObject *self,PyObject *args)
{
	int id,elasticId;
	double A,n,k;
    if(!PyArg_ParseTuple(args,"iiddd",&id,&elasticId,&A,&n,&k)) 
		return NULL;
	Material* pMaterial=new Creep(id,elasticId,A,n,k);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_DruckerPrager(PyObject *self,PyObject *args)
{
	int id,elasticId,type;
	double c,phi,psi,T;
    if(!PyArg_ParseTuple(args,"iiidddd",&id,&elasticId,&type,&c,&phi,&psi,&T)) 
		return NULL;
	Material* pMaterial=new DruckerPrager(id,elasticId,type,c,phi,psi,T);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_ModifiedCamClay(PyObject *self,PyObject *args)
{
	int id,elasticId;
	double M,po,kappa,lambda;
    if(!PyArg_ParseTuple(args,"iidddd",&id,&elasticId,&M,&po,&kappa,&lambda)) 
		return NULL;
	Material* pMaterial=new ModifiedCamClay(id,elasticId,M,po,kappa,lambda);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyMaterial_LadeDuncan(PyObject *self,PyObject *args)
{
	int id,elasticId;
	double K;
    if(!PyArg_ParseTuple(args,"iid",&id,&elasticId,&K)) 
		return NULL;
	Material* pMaterial=new LadeDuncan(id,elasticId,K);
	pD->add(pD->getMaterials(),pMaterial);
	createGroupByMaterial(id);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef MaterialMethods[] = 
{
	{"sdof",						pyMaterial_SDof,	
		METH_VARARGS,"Define a single dof material."},
	{"springElastic",				pyMaterial_SpringElastic,	
		METH_VARARGS,"Define a spring elastic material."},
	{"springContact",				pyMaterial_SpringContact,	
		METH_VARARGS,"Define a spring contact material."},
	{"uniElastic",					pyMaterial_UniaxialElastic,	
		METH_VARARGS,"Define a uniaxial elastic material."},
	{"uniElastoPlastic",			pyMaterial_UniaxialElastoPlastic,
		METH_VARARGS,"Define a uniaxial elastoplastic material."},
	{"uniCyclic",					pyMaterial_UniaxialCyclic,
		METH_VARARGS,"Define a cyclic material based on a backbone curve."},
	{"uniGap",						pyMaterial_UniaxialGap,
		METH_VARARGS,"Define a uniaxial elastoplastic material with a gap."},
	{"elastic",						pyMaterial_Elastic,
		METH_VARARGS,"Define a multiaxial elastic material."},
	{"DuncanChang",					pyMaterial_DuncanChang,
		METH_VARARGS,"Define a Duncan-Chang material."},
	{"vonMises",					pyMaterial_VonMises,
		METH_VARARGS,"Define a von-Mises type material."},
	{"MohrCoulomb",					pyMaterial_MohrCoulomb,
		METH_VARARGS,"Define a Mohr-Coulomb type material."},
	{"Tresca",						pyMaterial_Tresca,
		METH_VARARGS,"Define a Tresca type material."},
	{"DruckerPrager",				pyMaterial_DruckerPrager,
		METH_VARARGS,"Define a Drucker-Prager type material."},
	{"modifiedCamClay",				pyMaterial_ModifiedCamClay,
		METH_VARARGS,"Define a Modified Cam-Clay type material."},
	{"LadeDuncan",					pyMaterial_LadeDuncan,
		METH_VARARGS,"Define a Lade-Duncan type material."},
	{"creep",						pyMaterial_Creep,
		METH_VARARGS,"Define a creep material."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Element commands
******************************************************************************/
static PyObject* pyElement_Spring(PyObject *self, PyObject *args)
{
	int id,iNode,jNode,mat;
	double xp1=1.,xp2=0.,xp3=0.,yp1=0.,yp2=1.,yp3=0.;
    if(!PyArg_ParseTuple(args,"iiii|(ddd)(ddd)",&id,&iNode,&jNode,&mat,
		                                    &xp1,&xp2,&xp3,&yp1,&yp2,&yp3))
		return NULL;
	try
	{
		Element* pElement=new Spring(id,iNode,jNode,mat,xp1,xp2,xp3,yp1,yp2,yp3);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Bar2s(PyObject *self, PyObject *args)
{
	int id,iNode,jNode,mat,iSec,jSec=-9999;
    if(!PyArg_ParseTuple(args,"iiiii|i",&id,&iNode,&jNode,&mat,&iSec,&jSec))
		return NULL;
	try
	{	if(jSec==-9999) jSec=iSec;
		Element* pElement=new Bar2s(id,iNode,jNode,mat,iSec,jSec);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Bar2t(PyObject *self, PyObject *args)
{
	int id,iNode,jNode,mat,iSec,jSec=-9999;
    if(!PyArg_ParseTuple(args,"iiiii|i",&id,&iNode,&jNode,&mat,&iSec,&jSec))
		return NULL;
	try
	{
		if(jSec==-9999) jSec=iSec;
		Element* pElement=new Bar2t(id,iNode,jNode,mat,iSec,jSec);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Beam2e(PyObject *self, PyObject *args)
{
	int id,iNode,jNode,mat,sec;
    if(!PyArg_ParseTuple(args,"iiiii",&id,&iNode,&jNode,&mat,&sec))
		return NULL;
	try
	{
		Element* pElement=new Beam2e(id,iNode,jNode,mat,sec);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Beam2t(PyObject *self, PyObject *args)
{
	int id,iNode,jNode,mat,sec,rule=-1;
    if(!PyArg_ParseTuple(args,"iiiii|i",&id,&iNode,&jNode,&mat,&sec,&rule))
		return NULL;
	try
	{
		if(rule==-1) rule=1;
		Element* pElement=new Timoshenko2d(id,iNode,jNode,mat,sec,rule);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Beam3t(PyObject *self, PyObject *args)
{
	int id,iNode,jNode,mNode,mat,sec,rule=-1;
    if(!PyArg_ParseTuple(args,"iiiiii|i",&id,&iNode,&jNode,&mNode,&mat,&sec,&rule))
		return NULL;
	try
	{
		if(rule==-1) rule=2;
		Element* pElement=new Timoshenko2d(id,iNode,jNode,mNode,mat,sec,rule);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Brick8d(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,n5,n6,n7,n8,mat;
	if(!PyArg_ParseTuple(args,"iiiiiiiiii",&id,&n1,&n2,&n3,&n4,&n5,&n6,&n7,&n8,&mat))
		return NULL;
	try
	{
		Element* pElement=new Brick8d(id,n1,n2,n3,n4,n5,n6,n7,n8,mat);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Brick8b(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,n5,n6,n7,n8,mat;
	if(!PyArg_ParseTuple(args,"iiiiiiiiii",&id,&n1,&n2,&n3,&n4,&n5,&n6,&n7,&n8,&mat))
		return NULL;
	try
	{
		Element* pElement=new Brick8b(id,n1,n2,n3,n4,n5,n6,n7,n8,mat);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_SDOF(PyObject *self, PyObject *args)
{
	int id,nodeID,dofID,matID;
    if(!PyArg_ParseTuple(args,"iiii",&id,&nodeID,&dofID,&matID))
		return NULL;
	Element* pElement=new SDofElement(id,nodeID,dofID,matID);
	pD->add(pD->getElements(),pElement);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Quad4d(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,matID;
    if(!PyArg_ParseTuple(args,"iiiiii",
		&id,&n1,&n2,&n3,&n4,&matID))
		return NULL;
	try
	{
		Element* pElement;
		if(pD->getTag()==TAG_DOMAIN_PLANE_STRAIN||pD->getTag()==TAG_DOMAIN_PLANE_STRESS)
			pElement=new Quad4DispPlain(id,n1,n2,n3,n4,matID,2,2);
		else if(pD->getTag()==TAG_DOMAIN_AXISYMMETRIC)
			pElement=new Quad4DispAxisymmetric(id,n1,n2,n3,n4,matID,2,2);
		else throw SException("[nemesis:%d] %s",1110,"A quad4d can be used only in plane strain/stress/axisymmetry.");
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Quad4Test(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,matID;
    if(!PyArg_ParseTuple(args,"iiiiii",
		&id,&n1,&n2,&n3,&n4,&matID))
		return NULL;
	try
	{
		Element* pElement=new Quad4d(id,n1,n2,n3,n4,matID);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Quad4b(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,matID;
    if(!PyArg_ParseTuple(args,"iiiiii",
		&id,&n1,&n2,&n3,&n4,&matID))
		return NULL;
	try
	{
		Element* pElement=new Quad4b(id,n1,n2,n3,n4,matID);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Quad4i(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,matID;
    if(!PyArg_ParseTuple(args,"iiiiii",
		&id,&n1,&n2,&n3,&n4,&matID))
		return NULL;
	try
	{
		Element* pElement=new Quad4i(id,n1,n2,n3,n4,matID);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Quad4e(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,matID;
    if(!PyArg_ParseTuple(args,"iiiiii",
		&id,&n1,&n2,&n3,&n4,&matID))
		return NULL;
	try
	{
		Element* pElement;
		if(pD->getTag()==TAG_DOMAIN_PLANE_STRAIN||pD->getTag()==TAG_DOMAIN_PLANE_STRESS)
			pElement=new Quad4e(id,n1,n2,n3,n4,matID);
		else throw SException("[nemesis:%d] %s",1110,"A quad4e can be used only in plane strain/stress.");
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Triangle3d(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,matID;
    if(!PyArg_ParseTuple(args,"iiiii",&id,&n1,&n2,&n3,&matID))
		return NULL;
	try
	{
		Element* pElement=new Triangle3(id,n1,n2,n3,matID);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Triangle3dXFem(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,matID;
    if(!PyArg_ParseTuple(args,"iiiii",&id,&n1,&n2,&n3,&matID))
		return NULL;
	try
	{
		Element* pElement=new Triangle3XFem(id,n1,n2,n3,matID);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Triangle6d(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,n5,n6,matID;
    if(!PyArg_ParseTuple(args,"iiiiiiii",&id,&n1,&n2,&n3,&n4,&n5,&n6,&matID))
		return NULL;
	try
	{
		Element* pElement=new Triangle6(id,n1,n2,n3,n4,n5,n6,matID);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Tetrahedron4d(PyObject *self, PyObject *args)
{
	int id,n1,n2,n3,n4,matID;
    if(!PyArg_ParseTuple(args,"iiiiii",&id,&n1,&n2,&n3,&n4,&matID))
		return NULL;
	try
	{
		Element* pElement=new Tetrahedron4Disp(id,n1,n2,n3,n4,matID);
		pD->add(pD->getElements(),pElement);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Data(PyObject *self, PyObject *args)
{
	int id;
    if(!PyArg_ParseTuple(args,"i",&id))	return NULL;
 	ostringstream os;
	pD->get<Element>(pD->getElements(),id)->save(os);
	static istringstream is;
	is.str(os.str());
	return buildDict(is);
}
static PyObject* pyElement_Track(PyObject *self, PyObject *args)
{
	int id,index;
    if(!PyArg_ParseTuple(args,"ii",&id,&index))	return NULL;
	pD->get<Element>(pD->getElements(),id)->addTracker(index);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyElement_Path(PyObject *self, PyObject *args)
{
	int id,index;
    if(!PyArg_ParseTuple(args,"ii",&id,&index))	return NULL;
 	ostringstream os;
	pD->get<Element>(pD->getElements(),id)->getTracker(index)->save(os);
	istringstream is(os.str());
	return buildList(is);
}
static PyMethodDef ElementMethods[] = 
{
	{"spring",			pyElement_Spring,	
		METH_VARARGS,"Define a spring element."},
	{"bar",				pyElement_Bar2s,	
		METH_VARARGS,"An alias for bar2s."},
	{"bar2s",			pyElement_Bar2s,	
		METH_VARARGS,"Defines a simple 1d/2d/3d geometrically linear bar."},
	{"bar2t",			pyElement_Bar2t,	
		METH_VARARGS,"Defines a 1d/2d/3d bar based on a Total Lagrangian formulation."},
	{"beam2e",			pyElement_Beam2e,
		METH_VARARGS,"Defines a 2-Node Euler-Bernulli beam."},
	{"beam2t",			pyElement_Beam2t,
		METH_VARARGS,"Defines a 2-Node Timoshenko beam."},
	{"beam3t",			pyElement_Beam3t,
		METH_VARARGS,"Defines a 3-Node Timoshenko beam."},
	{"sdof",			pyElement_SDOF,
		METH_VARARGS,"Define a single dof element."},
	{"quad4d",			pyElement_Quad4d,	
		METH_VARARGS,"Define a 4-Noded standard displacement quad."},
	{"quad4t",			pyElement_Quad4Test,	
		METH_VARARGS,"Define a 4-Noded standard displacement quad."},
	{"quad4b",			pyElement_Quad4b,	
		METH_VARARGS,"Define a 4-Noded B-Bar quad."},
	{"quad4i",			pyElement_Quad4i,	
		METH_VARARGS,"Define a 4-Noded non-conforming quad."},
	{"quad4e",			pyElement_Quad4e,	
		METH_VARARGS,"Define a 4-Noded enhanced assumed strain quad."},
	{"tria3d",				pyElement_Triangle3d,	
		METH_VARARGS,"Define a 3-Noded constant strain triangle."},
	{"tria3x",				pyElement_Triangle3dXFem,	
		METH_VARARGS,"Define a 3-Noded constant strain triangle for XFem."},
	{"tria6d",				pyElement_Triangle6d,	
		METH_VARARGS,"Define a 6-Noded linear strain triangle."},
	{"tetra4d",			pyElement_Tetrahedron4d,	
		METH_VARARGS,"Define a 4-Noded constant strain tetrahedron."},
	{"brick8d",	pyElement_Brick8d,	
		METH_VARARGS,"Define a 8-Noded standard displacement brick."},
	{"brick8b",	pyElement_Brick8b,	
		METH_VARARGS,"Define a 8-Noded Bbar brick."},
	{"data",			pyElement_Data,
		METH_VARARGS,"Access to the element data."},
	{"track",			pyElement_Track,	
		METH_VARARGS,	"Add a Tracker to an element's material."},
	{"path",			pyElement_Path,	
		METH_VARARGS,	"Return the data recorded to the tracker."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Constraint commands
******************************************************************************/
static PyObject* pyConstraint_Set(PyObject *self, PyObject *args)
{
	int nodeID,dofID;
	double a,c;
    if(!PyArg_ParseTuple(args,"d(iid)",&c,&nodeID,&dofID,&a))
		return NULL;
	Constraint* pConstraint=new Constraint;
	pConstraint->setcDof(nodeID,dofID,a);
	pConstraint->setcVal(c);
	pD->add(pD->getConstraints(),pConstraint);
	Py_INCREF(Py_None);
	return Py_None;
}
///@todo make this one more general
static PyObject* pyConstraint_twoDofs(PyObject *self, PyObject *args)
{
	int node1,dof1,node2,dof2;
	double a1,a2,c;
    if(!PyArg_ParseTuple(args,"(iid)(iid)d",&node1,&dof1,&a1,&node2,&dof2,&a2,&c))
		return NULL;
	Constraint* pConstraint=new Constraint;
	pConstraint->setcDof(node1,dof1,a1);
	pConstraint->setcDof(node2,dof2,a2);
	pConstraint->setcVal(c);
	pD->add(pD->getConstraints(),pConstraint);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyConstraint_Linear(PyObject *self, PyObject *args)
{
	double c0;
	PyObject* cList;
    if(!PyArg_ParseTuple(args,"dO",&c0,&cList)) return NULL;
	Constraint* pConstraint=new Constraint;
	pConstraint->setcVal(c0);
	
	// Check if list ok
	if(!PyList_Check(cList))
	{
		PyErr_SetString(PyExc_StandardError, "Expected a list.");
		return NULL;
	}
	// Start reading list
	for(int i=0;i<PyList_Size(cList);i++)
	{
		// Get tuple
		PyObject* cTuple=PyList_GetItem(cList,i);
		// Check if tuple ok
		if(!PyTuple_Check(cTuple) || PyTuple_Size(cTuple)!=3)
		{
			PyErr_SetString(PyExc_StandardError, "Expected a tuple of size 3.");
			return NULL;
		}
		// Get arguments and check if arguments ok
		PyObject* tmpNode=PyTuple_GetItem(cTuple,0);
		PyObject* tmpDof =PyTuple_GetItem(cTuple,1);
		PyObject* tmpCi  =PyTuple_GetItem(cTuple,2);
		if(!PyInt_Check(tmpNode))
		{
			PyErr_SetString(PyExc_StandardError, "Expected (INTEGER,integer,float).");
			return NULL;
		}
		if(!PyInt_Check(tmpDof))
		{
			PyErr_SetString(PyExc_StandardError, "Expected (integer,INTEGER,float).");
			return NULL;
		}
		if(!PyFloat_Check(tmpCi))
		{
			PyErr_SetString(PyExc_StandardError, "Expected (integer,integer,FLOAT).");
			return NULL;
		}
		// Get arguments
		int node =PyInt_AsLong(tmpNode);
		int dof  =PyInt_AsLong(tmpDof);
		double ci=PyFloat_AsDouble(tmpCi);
		// Build triple
		pConstraint->setcDof(node,dof,ci);
	}
	// Add to domain
	pD->add(pD->getConstraints(),pConstraint);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef ConstraintMethods[] = 
{
	{"set",							pyConstraint_Set,	
		METH_VARARGS,"Define a single/multiple dof constraint."},
	{"twoDofs",						pyConstraint_twoDofs,	
		METH_VARARGS,"Define a two dof constraint."},
	{"linear",						pyConstraint_Linear,	
		METH_VARARGS,"Define a linear constraint."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Group commands
******************************************************************************/
static PyObject* pyGroup_Define(PyObject *self, PyObject *args)
{
	int groupId;
	if(!PyArg_ParseTuple(args,"i",&groupId)) return NULL;
	try
	{
		Group* pGroup=new Group(groupId);
		pD->add(pD->getGroups(),pGroup);
		pD->setCurrentGroup(groupId);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyGroup_Set(PyObject *self, PyObject *args)
{
	int groupID,elemID;
	if(!PyArg_ParseTuple(args,"ii",&groupID,&elemID)) return NULL;
	try
	{
		pD->get<Element>(pD->getElements(),elemID)->setGroup(groupID);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyGroup_State(PyObject *self, PyObject *args)
{
	if(currentLC<=0)
	{
		PyErr_SetString(PyExc_StandardError,"No LoadCase yet defined.");
		return NULL;
	}
	int groupID;
	int active;
	double facK=1.;
	double facS=1.;
	double facG=1.;
	double facP=1.;
	if(!PyArg_ParseTuple(args,"ii|ddddd",&groupID,&active,&facK,&facS,&facG,&facP)) return NULL;
	GroupState* pGroupState=new GroupState(groupID,active,facK,facS,facG,facP);
	pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addGroupState(pGroupState);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef GroupMethods[] = 
{
	{"define",	pyGroup_Define,	
		METH_VARARGS,"Define and set current group."},
	{"set",		pyGroup_Set,	
		METH_VARARGS,"Assigns a group to an element."},
	{"state",	pyGroup_State,	
		METH_VARARGS,"Define group state for current laodcase."},
	{NULL,NULL,0,NULL}
};

/******************************************************************************
* Sensitivity commands
******************************************************************************/
static PyObject* pySens_Elem(PyObject *self, PyObject *args)
{
	int elemID,param;
	if(!PyArg_ParseTuple(args,"ii",&elemID,&param)) return NULL;
	try
	{
		ElementSensitivityParameter* pParam=new ElementSensitivityParameter(elemID,param);
		pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addSensitivityParameter(pParam);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef SensitivityMethods[] = 
{
	{"elem",	pySens_Elem,	
		METH_VARARGS,"Define element sensitivity parameter."},
	{NULL,NULL,0,NULL}
};

/******************************************************************************
* Load commands
******************************************************************************/
static PyObject* pyLoad_NodeConstant(PyObject *self, PyObject *args)
{
	if(currentLC<=0)
	{
		PyErr_SetString(PyExc_StandardError,"No LoadCase yet defined.");
		return NULL;
	}
	int nodeID,dofID;
	double val;    
	if(!PyArg_ParseTuple(args,"iid",&nodeID,&dofID,&val))
		return NULL;
	Load* pLoad=new NodalLoadConstant(nodeID,dofID,val);
	pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addLoad(pLoad);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyLoad_Linear(PyObject *self, PyObject *args)
{
	if(currentLC<=0)
	{
		PyErr_SetString(PyExc_StandardError,"No LoadCase yet defined.");
		return NULL;
	}
	int nodeID,dofID;
	double val,grad;    
	if(!PyArg_ParseTuple(args,"iidd",&nodeID,&dofID,&val,&grad))
		return NULL;
	Load* pLoad=new NodalLoadLinear(nodeID,dofID,val,grad);
	pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addLoad(pLoad);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyLoad_Sin(PyObject *self, PyObject *args)
{
	if(currentLC<=0)
	{
		PyErr_SetString(PyExc_StandardError,"No LoadCase yet defined.");
		return NULL;
	}
	int nodeID,dofID;
	double a,omega,phi;    
	if(!PyArg_ParseTuple(args,"iiddd",&nodeID,&dofID,&a,&omega,&phi))
		return NULL;
	Load* pLoad=new NodalLoadSin(nodeID,dofID,a,omega,phi);
	pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addLoad(pLoad);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyLoad_BeamPoint(PyObject *self, PyObject *args)
{
	try
	{
		if(currentLC<=0) ///todo remove exception to lc
			throw SException("[nemesis:%d] %s",1110,"No LoadCase yet defined.");
		int elemID;
		const char* dir;
		double a0,p0;    
		if(!PyArg_ParseTuple(args,"isdd",&elemID,&dir,&a0,&p0)) return NULL;
		Load* pLoad=new BeamLoadPoint(elemID,dir,a0,p0);
		pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addLoad(pLoad);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyLoad_BeamUniform(PyObject *self, PyObject *args)
{
	try
	{
		if(currentLC<=0) ///todo remove exception to lc
			throw SException("[nemesis:%d] %s",1110,"No LoadCase yet defined.");
		int elemID;
		const char* dir;
		double p0;    
		if(!PyArg_ParseTuple(args,"isd",&elemID,&dir,&p0)) return NULL;
		Load* pLoad=new BeamLoadUniform(elemID,dir,p0);
		pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addLoad(pLoad);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef LoadMethods[] = 
{
	{"node",		pyLoad_NodeConstant,	
		METH_VARARGS,"Define a constant nodal load."},
	{"linear",		pyLoad_Linear,	
		METH_VARARGS,"Define a linear in time nodal load."},
	{"sin",			pyLoad_Sin,	
		METH_VARARGS,"Define a sinus nodal load."},
	{"beamPoint",	pyLoad_BeamPoint,	
		METH_VARARGS,"Define a concentrated force load acting on a beam."},
	{"beamUniform",	pyLoad_BeamUniform,	
		METH_VARARGS,"Define a uniform force load acting on a beam."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Ground Motion commands
******************************************************************************/
static PyObject* pyGroundMotion_File(PyObject *self, PyObject *args)
{
	try
	{
		if(currentLC<=0) ///todo remove exception to lc
			throw SException("[nemesis:%d] %s",1110,"No LoadCase yet defined.");
		int dof;
		const char* filename;
		double dt;
		double scale=1.0;
		if(!PyArg_ParseTuple(args,"isd|d",&dof,&filename,&dt,&scale)) return NULL;
		ifstream data(filename);
		Load* pLoad=new GroundMotionFile(dof,data,dt,scale);
		data.close();
		pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addLoad(pLoad);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyGroundMotion_Sin(PyObject *self, PyObject *args)
{
	try
	{
		if(currentLC<=0) ///todo remove exception to lc
			throw SException("[nemesis:%d] %s",1110,"No LoadCase yet defined.");
		int dof;
		double a, omega, phi=0.;
		if(!PyArg_ParseTuple(args,"idd|d",&dof,&a,&omega,&phi)) return NULL;
		Load* pLoad=new GroundMotionSin(dof,a,omega,phi);
		pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addLoad(pLoad);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}	
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef GroundMotionMethods[] = 
{
	{"file",	pyGroundMotion_File,	
		METH_VARARGS,"Define a uniform excitation through a file."},
	{"sin",		pyGroundMotion_Sin,	
		METH_VARARGS,"Define a sinus uniform excitation."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Initial Conditions commands
******************************************************************************/
static PyObject* pyInitialConditions_Disp(PyObject *self, PyObject *args)
{
	if(currentLC<=0)
	{
		PyErr_SetString(PyExc_StandardError,"No LoadCase yet defined.");
		return NULL;
	}
	int nodeID,dofID;
	double u;
	if(!PyArg_ParseTuple(args,"iid",&nodeID,&dofID,&u)) return NULL;
	InitialCondition* pInitial=new InitialDisplacement(nodeID,dofID,u);
	pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addInitialCondition(pInitial);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyInitialConditions_Velc(PyObject *self, PyObject *args)
{
	if(currentLC<=0)
	{
		PyErr_SetString(PyExc_StandardError,"No LoadCase yet defined.");
		return NULL;
	}
	int nodeID,dofID;
	double v;
	if(!PyArg_ParseTuple(args,"iid",&nodeID,&dofID,&v)) return NULL;
	InitialCondition* pInitial=new InitialVelocity(nodeID,dofID,v);
	pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addInitialCondition(pInitial);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyInitialConditions_Stresses(PyObject *self, PyObject *args)
{
	if(currentLC<=0)
	{
		PyErr_SetString(PyExc_StandardError,"No LoadCase yet defined.");
		return NULL;
	}
	int groupID,dir;
	double h1,s1,h2,s2,K0;
	if(!PyArg_ParseTuple(args,"iiddddd",&groupID,&dir,&h1,&s1,&h2,&s2,&K0)) return NULL;
	try
	{
		InitialStresses* pInitial=new InitialStresses(groupID,dir,h1,s1,h2,s2,K0);
		pD->get<LoadCase>(pD->getLoadCases(),currentLC)->addInitialCondition(pInitial);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef InitialConditionsMethods[] = 
{
	{"disp",		pyInitialConditions_Disp,
		METH_VARARGS,	"Set initial displacement to nodal dof."},
	{"velc",		pyInitialConditions_Velc,
		METH_VARARGS,	"Set initial velocity to nodal dof."},
	{"stresses",		pyInitialConditions_Stresses,
		METH_VARARGS,	"Set initial stresses to a group."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Crack commands
******************************************************************************/
static PyObject* pyCrack_OneTip(PyObject *self, PyObject *args)
{
	int id;
	double xS,yS,xT,yT;
	if(!PyArg_ParseTuple(args,"idddd",&id,&xS,&yS,&xT,&yT)) return NULL;
	Crack* pCrack=new Crack(id,xS,yS,xT,yT);
	pD->add<Crack>(pD->getCracks(),pCrack);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef CrackMethods[] = 
{
	{"oneTip",		pyCrack_OneTip,
		METH_VARARGS,	"Define a one tip crack."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* LoadCase commands
******************************************************************************/
static PyObject* pyLC_Define(PyObject *self, PyObject *args)
{
	int id;
	const char* name=0;
	if(!PyArg_ParseTuple(args,"i|s",&id,&name)) return NULL;
	LoadCase* pLC;
	///todo: check this better
	if(name!=0) pLC=new LoadCase(id,name);
	else		pLC=new LoadCase(id,"default");
	pD->add(pD->getLoadCases(),pLC);
	currentLC=id;
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyLC_Info(PyObject *self, PyObject *args)
{
	int id;
    if(!PyArg_ParseTuple(args,"i",&id))	return NULL;
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef LCMethods[] = 
{
	{"define",	pyLC_Define,	
		METH_VARARGS,"Define a loadcase."},
	{"info",	pyLC_Info,
		METH_VARARGS,"Print loadcase info."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Analysis commands
******************************************************************************/
static PyObject* pyAnalysis_Static(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	AnalysisType* pType=new StaticAnalysis();
	pA->setAnalysisType(pType);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAnalysis_Transient(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	AnalysisType* pType=new TransientAnalysis();
	pA->setAnalysisType(pType);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAnalysis_Eigen(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	AnalysisType* pType=new EigenAnalysis();
	pA->setAnalysisType(pType);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAnalysis_Sensitivity(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	AnalysisType* pType=new SensitivityStaticAnalysis();
	try
	{
		pA->setAnalysisType(pType);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAnalysis_XFem(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	AnalysisType* pType=new XFemAnalysis();
	try
	{
		pA->setAnalysisType(pType);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAnalysis_Run(PyObject *self, PyObject *args)
{
	int LC,steps,ret;
	if(!PyArg_ParseTuple(args,"ii",&LC,&steps)) return NULL;
	try
	{
		ret=pA->analyze(LC,steps);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef AnalysisMethods[] = 
{
	{"static",	pyAnalysis_Static,	
		METH_VARARGS,"Define a static analysis."},
	{"transient",	pyAnalysis_Transient,	
		METH_VARARGS,"Define a transient analysis."},
	{"eigen",		pyAnalysis_Eigen,	
		METH_VARARGS,"Define a eigeinvalue analysis."},
	{"sensitivity",	pyAnalysis_Sensitivity,	
		METH_VARARGS,"Define a sensitivity analysis."},
	{"xFem",		pyAnalysis_XFem,	
		METH_VARARGS,"Define an xFem analysis."},
	{"run",	pyAnalysis_Run,	
		METH_VARARGS,"Run given analysis."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Imposer commands
******************************************************************************/
static PyObject* pyImposer_Eliminatation(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	Imposer* pImposer=new EliminationImposer();
	pA->setImposer(pImposer);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyImposer_Lagrange(PyObject *self, PyObject *args)
{
	double a;
	if(!PyArg_ParseTuple(args,"",&a)) return NULL;
	Imposer* pImposer=new LagrangeImposer();
	pA->setImposer(pImposer);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyImposer_Penalty(PyObject *self, PyObject *args)
{
	double a;
	if(!PyArg_ParseTuple(args,"d",&a)) return NULL;
	Imposer* pImposer=new PenaltyImposer(a);
	pA->setImposer(pImposer);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef ImposerMethods[] = 
{
	{"elimination",	pyImposer_Eliminatation,	
		METH_VARARGS,"Use the elimination method to impose constraints."},
	{"lagrange",	pyImposer_Lagrange,	
		METH_VARARGS,"Use the Langange multipliers method to impose constraints."},
	{"penalty",	pyImposer_Penalty,	
		METH_VARARGS,"Use the penalty method to impose constraints."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Control commands
******************************************************************************/
static PyObject* pyControl_Load(PyObject *self, PyObject *args)
{
	double DL0,minDL=1.,maxDL=1.,n=0.5,dt=0.;
	int Id=1;
	if(!PyArg_ParseTuple(args,"d|ddidd",&DL0,&minDL,&maxDL,&Id,&n,&dt)) return NULL;
	Control* pControl=new LoadControl(DL0,minDL,maxDL,Id,n,dt);
	pA->setControl(pControl);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyControl_ArcLength(PyObject *self, PyObject *args)
{
	double DL0,minDL=1.,maxDL=1.,n=0.5,dt=0.;
	int Id=1;
	if(!PyArg_ParseTuple(args,"d|ddidd",&DL0,&minDL,&maxDL,&Id,&n,&dt)) return NULL;
	Control* pControl=new ArcLengthSpherical(DL0,minDL,maxDL,Id,n,dt);
	pA->setControl(pControl);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyControl_ArcLengthUNP(PyObject *self, PyObject *args)
{
	double DL0,minDL=1.,maxDL=1.,n=0.5,dt=0.;
	int Id=1;
	if(!PyArg_ParseTuple(args,"d|ddidd",&DL0,&minDL,&maxDL,&Id,&n,&dt)) return NULL;
	Control* pControl=new ArcLengthUNP(DL0,minDL,maxDL,Id,n,dt);
	pA->setControl(pControl);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyControl_Disp(PyObject *self, PyObject *args)
{
	int nodeID,dofID;
	double Du0,minDu=1.,maxDu=1.,dt=0.;
	int Id=1;
	double n=0.5;
	if(!PyArg_ParseTuple(args,"iid|ddidd",&nodeID,&dofID,&Du0,&minDu,&maxDu,&Id,&n,&dt)) return NULL;
	Control* pControl=new DisplacementControl(nodeID,dofID,Du0,minDu,maxDu,Id,n,dt);
	pA->setControl(pControl);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyControl_Newmark(PyObject *self, PyObject *args)
{
	double beta,gamma,dt;
	if(!PyArg_ParseTuple(args,"ddd",&beta,&gamma,&dt)) return NULL;
	Control* pControl=new Newmark(beta,gamma,dt);
	pA->setControl(pControl);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef ControlMethods[] = 
{
	{"load",			pyControl_Load,	
		METH_VARARGS,"Use load control for the analysis."},
	{"arcLength",		pyControl_ArcLength,	
		METH_VARARGS,"Use arc-length control for the analysis."},
	{"arcLengthUNP",		pyControl_ArcLengthUNP,	
		METH_VARARGS,"Use arc-length control for the analysis."},
	{"disp",			pyControl_Disp,	
		METH_VARARGS,"Use arc-length control for the analysis."},
	{"Newmark",			pyControl_Newmark,	
		METH_VARARGS,"Use Newmark control for the analysis."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Algorithm commands
******************************************************************************/
static PyObject* pyAlgorithm_Linear(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	Algorithm* pAlgorithm=new LinearAlgorithm();
	pA->setAlgorithm(pAlgorithm);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAlgorithm_BFGS(PyObject *self, PyObject *args)
{
	int m;
	int lsearch=0;
	double etaMin=0.1,etaMax=1.0,rTol=0.8;
	int maxIter=10;
	if(!PyArg_ParseTuple(args,"i|idddi",&m,&lsearch,&etaMin,&etaMax,&rTol,&maxIter)) return NULL;
	Algorithm* pAlgorithm;
	if(lsearch==0)	pAlgorithm=new BFGS(m);
	else			pAlgorithm=new BFGS(m,etaMin,etaMax,rTol,maxIter);
	pA->setAlgorithm(pAlgorithm);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAlgorithm_FullNewtonRaphson(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	Algorithm* pAlgorithm=new NewtonRaphsonFull();
	pA->setAlgorithm(pAlgorithm);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAlgorithm_InitialNewtonRaphson(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	Algorithm* pAlgorithm=new NewtonRaphsonInitial();
	pA->setAlgorithm(pAlgorithm);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyAlgorithm_ModifiedNewtonRaphson(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	Algorithm* pAlgorithm=new NewtonRaphsonModified();
	pA->setAlgorithm(pAlgorithm);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef AlgorithmMethods[] = 
{
	{"linear",			pyAlgorithm_Linear,	
		METH_VARARGS,	"Use linear algorithm."},	
	{"BFGS",			pyAlgorithm_BFGS,	
		METH_VARARGS,	"Use BFGS rank-two update algorithm."},
	{"fNR",				pyAlgorithm_FullNewtonRaphson,	
		METH_VARARGS,	"Use full Newton-Raphson algorithm."},
	{"mNR",				pyAlgorithm_ModifiedNewtonRaphson,	
		METH_VARARGS,	"Use modified Newton-Raphson algorithm."},
	{"iNR",				pyAlgorithm_InitialNewtonRaphson,	
		METH_VARARGS,	"Use Newton-Raphson algorithm with initial stiffness."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Convergence commands
******************************************************************************/
static PyObject* pyConvergence_Set(PyObject *self, PyObject *args)
{
	int maxIter;
	double tolRabs;
	double tolRrel=1.0;
	double tolWrel=1.0;
	if(!PyArg_ParseTuple(args,"id|dd",
		&maxIter,&tolRabs,&tolRrel,&tolWrel)) return NULL;
	pA->getConvergenceNorm()->setCheck(maxIter,tolRabs,tolRrel,tolWrel);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef ConvergenceMethods[] = 
{
	{"set",	pyConvergence_Set,	
		METH_VARARGS,"Set convergence check's properties."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* SOE commands
******************************************************************************/
static PyObject* pySOE_FullLinearSOE(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	SOE* pSOE=new FullLinearSOE();
	pA->setSOE(pSOE);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pySOE_SymmLinearSOE(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	SOE* pSOE=new SymmLinearSOE();
	pA->setSOE(pSOE);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pySOE_BandLinearSOE(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	SOE* pSOE=new BandLinearSOE();
	pA->setSOE(pSOE);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pySOE_Info(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	pA->getSOE()->print();
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pySOE_PlotGraph(PyObject *self, PyObject *args)
{
	const char* s;
	if(!PyArg_ParseTuple(args,"s",&s)) return NULL;
	try
	{
		pA->getSOE()->plotGraph(s);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef SOEMethods[] = 
{
	{"full",	pySOE_FullLinearSOE,	
		METH_VARARGS,"Use a full storage scheme to hold the system of equations."},
	{"symm",	pySOE_SymmLinearSOE,	
		METH_VARARGS,"Use a symmetric storage scheme to hold the system of equations."},
	{"band",	pySOE_BandLinearSOE,	
		METH_VARARGS,"Use a band storage scheme to hold the system of equations."},
	{"info",	pySOE_Info,	
		METH_VARARGS,"Plot the Graph under given filename."},
	{"plotGraph",	pySOE_PlotGraph,	
		METH_VARARGS,"Plot the Graph under given filename."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Reorder commands
******************************************************************************/
static PyObject* pyReorder_RCM(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	Reorderer* pReorderer=new ReverseCuthillMckee();
	pA->setReorderer(pReorderer);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyReorder_FCM(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	Reorderer* pReorderer=new ForwardCuthillMckee();
	pA->setReorderer(pReorderer);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyReorder_FSloan(PyObject *self, PyObject *args)
{
	double w1=0.5,w2=0.5;
	if(!PyArg_ParseTuple(args,"|dd",&w1,&w2)) return NULL;
	Reorderer* pReorderer=new ForwardSloan(w1,w2);
	pA->setReorderer(pReorderer);
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyReorder_RSloan(PyObject *self, PyObject *args)
{
	double w1=0.5,w2=0.5;
	if(!PyArg_ParseTuple(args,"|dd",&w1,&w2)) return NULL;
	Reorderer* pReorderer=new ReverseSloan(w1,w2);
	pA->setReorderer(pReorderer);
	Py_INCREF(Py_None);
	return Py_None;
}
/*
static PyObject* pyReorder_King(PyObject *self, PyObject *args)
{
	if(!PyArg_ParseTuple(args,"")) return NULL;
	Reorderer* pReorderer=new King();
	pA->setReorderer(pReorderer);
	Py_INCREF(Py_None);
	return Py_None;
}
*/
static PyMethodDef ReorderMethods[] = 
{
	{"rCM",	pyReorder_RCM,	
		METH_VARARGS,"Use reverse Cuthill-Mckee to reorder nodal numbering."},
	{"fCM",	pyReorder_FCM,	
		METH_VARARGS,"Use forward Cuthill-Mckee to reorder nodal numbering."},
	{"rSloan",	pyReorder_RSloan,	
		METH_VARARGS,"Use reverse Sloan to reorder nodal numbering."},
	{"fSloan",	pyReorder_FSloan,	
		METH_VARARGS,"Use forward Sloan to reorder nodal numbering."},
//	{"mdo",	pyReorder_MinimumDegreeOrdering,	
//		METH_VARARGS,"Use Minimum Degree Ordering to reorder nodal numbering."},
	{NULL,NULL,0,NULL}
};
/******************************************************************************
* Tracker commands
******************************************************************************/
/*
static PyObject* pyTracker_Node(PyObject *self, PyObject *args)
{
	int id,nodeID;
	if(!PyArg_ParseTuple(args,"ii",&id,&nodeID)) return NULL;
	try
	{
		Tracker* pTracker=new Tracker(id,nodeID);
		pD->add(pD->getTrackers(),pTracker);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyObject* pyTracker_Steps(PyObject *self, PyObject *args)
{
	int id;
	int steps=0;
	if(!PyArg_ParseTuple(args,"i",&id)) return NULL;
	try
	{
		steps=pD->get<Tracker>(pD->getTrackers(),id)->getSteps();
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	return Py_BuildValue("i",steps);
}
static PyObject* pyTracker_Time(PyObject *self, PyObject *args)
{
	int id;
	int step;
	double time=0;
	if(!PyArg_ParseTuple(args,"ii",&id,&step)) return NULL;
	try
	{
		time=pD->get<Tracker>(pD->getTrackers(),id)->getTime(step);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	return Py_BuildValue("d",time);
}
static PyObject* pyTracker_Lambda(PyObject *self, PyObject *args)
{
	int id;
	int step;
	double lambda=0;
	if(!PyArg_ParseTuple(args,"ii",&id,&step)) return NULL;
	try
	{
		lambda=pD->get<Tracker>(pD->getTrackers(),id)->getLambda(step);
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	return Py_BuildValue("d",lambda);
}
static PyObject* pyTracker_Data(PyObject *self, PyObject *args)
{
	int id,step;
	if(!PyArg_ParseTuple(args,"ii",&id,&step)) return NULL;
	try
	{
		return buildTuple(pD->get<Tracker>(pD->getTrackers(),id)->getPacket(step));
	}
	catch(SException e)
	{
		PyErr_SetString(PyExc_StandardError,e.what());
		return NULL;
	}
	Py_INCREF(Py_None);
	return Py_None;
}
static PyMethodDef TrackerMethods[] = 
{
	{"node",			pyTracker_Node,	
		METH_VARARGS,	"Add a new tracker to given node."},
	{"steps",			pyTracker_Steps,	
		METH_VARARGS,	"Get number tracker steps."},
	{"time",			pyTracker_Time,	
		METH_VARARGS,	"Get tracker time for a given step."},
	{"Lambda",			pyTracker_Lambda,	
		METH_VARARGS,	"Get tracker lambda for a given step."},
	{"data",			pyTracker_Data,	
		METH_VARARGS,	"Get tracker data for a given step."},
	{NULL,NULL,0,NULL}
};
*/
/******************************************************************************
* Log commands
******************************************************************************/
PyObject* log_CaptureStdout(PyObject* self, PyObject* pArgs)
{
 char* LogStr = NULL;
 if (!PyArg_ParseTuple(pArgs, "s", &LogStr)) return NULL;

 printf("%s", LogStr); 
 // Simply using printf to do the real work. 
 // You could also write it to a .log file or whatever...
 // MessageBox(NULL, LogStr...
 // WriteFile(hFile, LogStr...

 Py_INCREF(Py_None);
 return Py_None;
}

// Notice we have STDERR too.
PyObject* log_CaptureStderr(PyObject* self, PyObject* pArgs)
{
 char* LogStr = NULL;
 if (!PyArg_ParseTuple(pArgs, "s", &LogStr)) return NULL;

 printf("%s", LogStr);

 Py_INCREF(Py_None);
 return Py_None;
}

static PyMethodDef logMethods[] = {
 {"CaptureStdout", log_CaptureStdout, METH_VARARGS, "Logs stdout"},
 {"CaptureStderr", log_CaptureStderr, METH_VARARGS, "Logs stderr"},
 {NULL, NULL, 0, NULL}
};
/******************************************************************************
* PyParser Methods
******************************************************************************/
/**
 * Default constructor.
 */
PyParser::PyParser()
{
	pD=&D;
	pA=&A;
}
/**
 * Destructor.
 */
PyParser::~PyParser()
{
}
/**
* Parse interactively.
*/
int PyParser::parse()
{
	Py_Initialize();
	this->initModules();
	FILE* fp=stdin;
	PyRun_InteractiveLoop(fp,"<Solver-Interactive>"); 
	Py_Finalize();
	return 0;
}
/**
* Parse file.
*/
int PyParser::parse(char* filename)
{
	// Check if the extension is correct
	int len=strlen(filename);
	if((  filename[len-4]!='.' 
		||filename[len-3]!='s' 
		||filename[len-2]!='l'
		||filename[len-1]!='v'))
	{
		cout<<"solver: Files passed with '-p' should have an '.slv' extension."<<endl;
		return 0;
	}
	// Check if file exists
	fstream tmpFile;
	tmpFile.open(filename,ios::binary|ios::in);
	if(tmpFile.fail())
	{
		tmpFile.clear();
		tmpFile.close();
		cout<<"solver: File not found!"<<endl;
		return 0;
	}
	tmpFile.close();
	// Now it is ok to parse the file
	Py_Initialize();
	this->initModules();
	PyObject* PyFileObject = PyFile_FromString(filename,"r");
	PyRun_SimpleFile(PyFile_AsFile(PyFileObject),filename);
	Py_DECREF(PyFileObject);
	Py_Finalize();
	return 0;
}
int PyParser::initModules()
{
	Py_InitModule("db",DatabaseMethods);
    Py_InitModule("domain",DomainMethods);
    Py_InitModule("node",NodeMethods);
    Py_InitModule("section",SectionMethods);
    Py_InitModule("material",MaterialMethods);
    Py_InitModule("element",ElementMethods);
    Py_InitModule("constraint",ConstraintMethods);
    Py_InitModule("load",LoadMethods);
    Py_InitModule("groundMotion",GroundMotionMethods);
    Py_InitModule("initial",InitialConditionsMethods);
    Py_InitModule("lc",LCMethods);
    Py_InitModule("analysis",AnalysisMethods);
    Py_InitModule("imposer",ImposerMethods);
    Py_InitModule("control",ControlMethods);
    Py_InitModule("algorithm",AlgorithmMethods);
    Py_InitModule("convergence",ConvergenceMethods);
	Py_InitModule("soe",SOEMethods);
	Py_InitModule("reorder",ReorderMethods);
	Py_InitModule("group",GroupMethods);
	Py_InitModule("sens",SensitivityMethods);
	Py_InitModule("crack",CrackMethods);
	PyRun_SimpleString("import db");
	PyRun_SimpleString("import domain");
	PyRun_SimpleString("import node");
	PyRun_SimpleString("import section");
	PyRun_SimpleString("import material");
	PyRun_SimpleString("import material");
	PyRun_SimpleString("import element");
	PyRun_SimpleString("import constraint");
	PyRun_SimpleString("import load");
	PyRun_SimpleString("import groundMotion");
	PyRun_SimpleString("import initial");
	PyRun_SimpleString("import lc");
	PyRun_SimpleString("import analysis");
	PyRun_SimpleString("import imposer");
	PyRun_SimpleString("import control");
	PyRun_SimpleString("import algorithm");
	PyRun_SimpleString("import convergence");
	PyRun_SimpleString("import soe");
	PyRun_SimpleString("import reorder");
	PyRun_SimpleString("import group");
	PyRun_SimpleString("import sens");
	PyRun_SimpleString("import crack");
//	PyRun_SimpleString("import tracker");
	PyRun_SimpleString("import nemesis");
	
	Py_InitModule("log",logMethods);
	PyRun_SimpleString(	"import log\n"
						"import sys\n"
						"class StdoutCatcher:\n"
                        "\tdef write(self, str):\n"
                        "\t\tlog.CaptureStdout(str)\n"
                        "class StderrCatcher:\n"
                        "\tdef write(self, str):\n"
                        "\t\tlog.CaptureStderr(str)\n"
                        "sys.stdout = StdoutCatcher()\n"
                        "sys.stderr = StderrCatcher()\n"	);	
	return 0;
}
