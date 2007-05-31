/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
*                                                                             *
*   This program is free software; you can redistribute it and/or modify      *
*   it under the terms of the GNU General Public License version 2, as        *
*   published by the Free Software Foundation.                                *
*                                                                             *
*   This program is distributed in the hope that it will be useful,           *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License along   *
*   with this program; if not, write to the Free Software Foundation, Inc.,   *
*   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.               *
******************************************************************************/

//*****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
//*****************************************************************************

#include <Element.h>
#include <math.h>

Matrix** Element::theStaticMatrices=0;
Vector** Element::theStaticVectors=0;

/**
 * Default Constructor
 */
Element::Element()
{
}
/**
 * Constructor
 * Creates an Element with a given ID.
 */
Element::Element(int ID,int matID)
:DomainObject(ID),myMatID(0),myMaterial(0),activeParameter(0)
{
	if(theStaticMatrices==0)
	{
		theStaticMatrices=new Matrix*[64];
		for(int i=1;i<64;i++) theStaticMatrices[i]=new Matrix(i,i,0.);
		theStaticVectors=new Vector*[64];
		for(int i=1;i<64;i++) theStaticVectors[i]=new Vector(i,0.);
	}
	myMatID=matID;
	if(pD->areGroupsByMaterial())
		myGroup=pD->get<Group>(pD->getGroups(),myMatID);
	else
		myGroup=pD->get<Group>(pD->getGroups(),pD->getCurrentGroup());
}
Element::~Element()
{
	if(theStaticMatrices!=0)
	{
		for(int i=1;i<64;i++) delete theStaticMatrices[i];
		for(int i=1;i<64;i++) delete theStaticVectors[i];
		delete[] theStaticMatrices;
		delete[] theStaticVectors;
		theStaticMatrices=0;
		theStaticVectors=0;
	}
}
const Matrix& Element::getC()
{
	static Matrix C;
	int nDofs=myNodalIDs.size()*myLocalNodalDofs.size();
	C.resize(nDofs,nDofs,0.);
	Vector RayleighFactors=pD->getRayleighFactors();
	if(RayleighFactors.size()==0)
		C.clear();
	else
	{
		C.add_cM(RayleighFactors[0],this->getK());
		C.add_cM(RayleighFactors[1],this->getM());
	}
	return C;
}
int Element::handleCommonInfo()
{
	// Define number of nodes, number of dofs
	int nNodes=myNodalIDs.size();
	int nLocalDofs=myLocalNodalDofs.size();
	int nDofs=nNodes*nLocalDofs;
    myNodes.resize(nNodes);
	// Find own matrix and vector
	myMatrix=theStaticMatrices[nDofs];
	myVector=theStaticVectors[nDofs];
	// Retrieve the node pointers and check if they do exist
	for(int i=0;i<nNodes;i++)
		myNodes[i]=pD->get<Node>(pD->getNodes(),myNodalIDs[i]);
	// Get nodal coordinates
	x.resize(nNodes,3);
	for(int i=0;i<nNodes;i++) 
	{
		x(i,0)=myNodes[i]->getx1();
		x(i,1)=myNodes[i]->getx2();
		x(i,2)=myNodes[i]->getx3();
	}
	// Resize external load vector, self weigth vector
	P.resize(nDofs);
	G.resize(nDofs);
	// Inform each of the nodes that an element is connected to them
	for(int i=0;i<nNodes;i++) myNodes[i]->addEleToNode(this);
	// Inform the nodes that the corresponding Dof's must be activated
	for(int i=0;i<nNodes;i++) 
		for(int j=0;j<nLocalDofs;j++) 
			myNodes[i]->addDofToNode(myLocalNodalDofs[j]);
	// Create load vector
	for(int i=0;i<nDofs;i++) P[i]=0;
	// Retrieve the material pointer
	if(myMatID>0)
	{
		myMaterial=pD->get<Material>(pD->getMaterials(),myMatID);
		// Check if the material can be assigned
		this->checkIfAllows(myMaterial);
		// Self weight
		///@todo make g=9.81 non-constant
		b.resize(3);
		b[0]=9.81*(pD->getGravityDirection()[0])*(myMaterial->getRho());
		b[1]=9.81*(pD->getGravityDirection()[1])*(myMaterial->getRho());
		b[2]=9.81*(pD->getGravityDirection()[2])*(myMaterial->getRho());
	}
	return 0;
}
/**
 * Returns the id's of the Nodes of an element.
 * @return An IDContainer of the Nodes.
 */
const IDContainer& Element::getNodalIDs() const
{
	return myNodalIDs;
}
const IDContainer& Element::getLocalNodalDofs() const
{
	return myLocalNodalDofs;
}
const std::vector<Node*>& Element::getNodes() const
{
	return myNodes;
}
const Vector& Element::getDispTrial()
{
	int nDofs=myLocalNodalDofs.size();
	int nNodes=myNodalIDs.size();
	Vector& disp=*myVector;
	for(int i=0;i<nNodes;i++)
		for(int j=0;j<nDofs;j++)
			disp[i*nDofs+j]=myNodes[i]->getDispTrialAtDof(myLocalNodalDofs[j]);
	return disp;
}
const Vector& Element::getVelcTrial()
{
	int nDofs=myLocalNodalDofs.size();
	int nNodes=myNodalIDs.size();
	Vector& velc=*myVector;
	for(int i=0;i<nNodes;i++)
		for(int j=0;j<nDofs;j++)
			velc[i*nDofs+j]=myNodes[i]->getVelcTrialAtDof(myLocalNodalDofs[j]);
	return velc;
}
const Vector& Element::getAcclTrial()
{
	int nDofs=myLocalNodalDofs.size();
	int nNodes=myNodalIDs.size();
	Vector& accl=*myVector;
	for(int i=0;i<nNodes;i++)
		for(int j=0;j<nDofs;j++)
			accl[i*nDofs+j]=myNodes[i]->getAcclTrialAtDof(myLocalNodalDofs[j]);
	return accl;
}
const Vector& Element::getDispConvg()
{
	int nDofs=myLocalNodalDofs.size();
	int nNodes=myNodalIDs.size();
	Vector& disp=*myVector;
	for(int i=0;i<nNodes;i++)
		for(int j=0;j<nDofs;j++)
			disp[i*nDofs+j]=myNodes[i]->getDispConvgAtDof(myLocalNodalDofs[j]);
	return disp;
}
const Vector& Element::getVelcConvg()
{
	int nDofs=myLocalNodalDofs.size();
	int nNodes=myNodalIDs.size();
	Vector& velc=*myVector;
	for(int i=0;i<nNodes;i++)
		for(int j=0;j<nDofs;j++)
			velc[i*nDofs+j]=myNodes[i]->getVelcConvgAtDof(myLocalNodalDofs[j]);
	return velc;
}
const Vector& Element::getAcclConvg()
{
	int nDofs=myLocalNodalDofs.size();
	int nNodes=myNodalIDs.size();
	Vector& disp=*myVector;
	for(int i=0;i<nNodes;i++)
		for(int j=0;j<nDofs;j++)
			disp[i*nDofs+j]=myNodes[i]->getAcclConvgAtDof(myLocalNodalDofs[j]);
	return disp;
}
const Vector& Element::getDispIncrm()
{
	int nDofs=myLocalNodalDofs.size();
	int nNodes=myNodalIDs.size();
	Vector& disp=*myVector;
	///@todo: CHECK!!!!!
	for(int i=0;i<nNodes;i++)
		for(int j=0;j<nDofs;j++)
			disp[i*nDofs+j]=myNodes[i]->getDispTrialAtDof(myLocalNodalDofs[j])-
			                myNodes[i]->getDispConvgAtDof(myLocalNodalDofs[j]);
	return disp;
}
void Element::zeroLoad()
{
	P.clear();
}
void Element::addLoad(const Vector& val,double fac)
{
	P.add_cV(fac,val);
}
void Element::addInitialStresses(InitialStresses* pInitialStresses)
{
	// needs to be overwritten
}
void Element::addGroundMotion(int dof,double val)
{
	this->getM();
	Matrix& M=*myMatrix;
	int pos=Containers::index_find(myLocalNodalDofs,dof);
	if(pos<0) return;
	int nNodes=myNodalIDs.size();
	int nLocalDofs=myLocalNodalDofs.size();
	int nDofs=nLocalDofs*nNodes;
    for(int i=0;i<nNodes;i++)
		for(int j=0;j<nDofs;j++)
			P[i*nLocalDofs+pos]-=M(j,i*nLocalDofs)*val;
}
bool Element::isActive()
{
	return myGroup->isActive();
}
void Element::setGroup(int groupID)
{
		myGroup=pD->get<Group>(pD->getGroups(),groupID);
}
const int Element::getnPlasticPoints()
{
	int n=0;
//	for(unsigned i=0;i<theMaterialItems.size();i++)
//		if(theMaterialItems[i]->plastified()) n++;
	return n;
}
const Packet& Element::getPacket()
{
	thePacket.zero();
	thePacket.tag=this->getTag();
	thePacket.id=this->getID();
	// Store nodes
	thePacket.intArray[0]=myNodalIDs.size();
	for(unsigned i=0;i<myNodalIDs.size();i++) thePacket.intArray[i+1]=myNodalIDs[i];
	// Store nDofs
	thePacket.intArray[28]=myNodalIDs.size()*myLocalNodalDofs.size();
	// Store plastic points
	thePacket.intArray[29]=myMatID;
	// Store material
	thePacket.intArray[30]=this->getnPlasticPoints();
	// Send data
	return thePacket;
}
void Element::setPacket(const Packet& p)
{
}
void Element::save(std::ostream& s)
{
	s<<"ELEMENT "	<<' ';
	s<<"tag	"	<<1000<<' '<<myTag<<' ';
	s<<"id "	<<1000<<' '<<myID<<' ';
	s<<"mat "	<<1000<<' '<<myMatID<<' ';
	s<<"END "<<' ';
}
