/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]  *
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

#include "elements/Bar.h"

/**
 * Default constructor.
 */
Bar::Bar()	
{
}
/**
 * Constructor.
 * Creates a Bar Element. 
 */
Bar::Bar(int ID,int Node_1,int Node_2,int matID,int iSecID,int jSecID)	
:Element(ID,matID)
{
	// Get dimension
	nDim=pD->getnDim();
	
	// Store the nodes
	myNodalIDs.resize(2);
	myNodalIDs[0]=Node_1;
	myNodalIDs[1]=Node_2;

	// The dofs needed for this element
	myLocalNodalDofs.resize(nDim);
	for(int i=0;i<nDim;i++) myLocalNodalDofs[i]=i;

	// Handle common info
    this->handleCommonInfo();
	
	// Find length
	L0=0.;
	for(int i=0;i<nDim;i++) L0+=(x(1,i)-x(0,i))*(x(1,i)-x(0,i));
	L0=sqrt(L0);
	if(num::tiny(L0)) throw SException("[nemesis:%d] %s",9999,"Zero length bar is not permitted");
	
	// Retrieve the CrossSection pointers and get A0
	iSection=pD->get<CrossSection>(pD->getCrossSections(),iSecID);
	jSection=pD->get<CrossSection>(pD->getCrossSections(),jSecID);
	A0=0.5*(iSection->getA()+jSection->getA());

	// Store material information
	myUniMaterial=static_cast<UniaxialMaterial*>(myMaterial)->getClone();
}
/**
 * Destructor.
 */
Bar::~Bar()
{
	delete myUniMaterial;
}
void Bar::update()
{
	static Vector du(2*nDim);
	du=this->getDispIncrm();
    double dL=0;
	for (int i=0;i<nDim;i++) dL+=(du[i+nDim]-du[i])*cosX[i];
  	double de=dL/L0;
	myUniMaterial->setStrain(de);
}
void Bar::commit()
{
	myUniMaterial->commit();
}
bool Bar::checkIfAllows(FEObject* f)
{
	return true;
}
const Matrix& Bar::getM()
{   
	Matrix& M=*myMatrix;
	M.clear();
	double rho=myUniMaterial->getRho();
	double mass=0.5*L0*A0*rho;
	for(int i=0;i<nDim;i++)
	{
		M(i     ,i     ) = mass;
		M(i+nDim,i+nDim) = mass;
	}
	return M;
}
const Vector& Bar::getReff()
{
	Vector velc=this->getVelcTrial(); ///@todo: problem with memory sharing otherwise
	Vector& Reff=*myVector;
	// +Fint-Fext
	this->getR();
	// +M*aTrial
	double rho=myUniMaterial->getRho();
	double mass=0.5*rho*L0;
	const Vector& a0=myNodes[0]->getAcclTrial();
	const Vector& a1=myNodes[1]->getAcclTrial();
	for (int i=0;i<nDim;i++)
	{
	    Reff[i]      += mass*a0[i];
	    Reff[i+nDim] += mass*a1[i];
	}
	Reff+=(this->getC())*(velc);
	return Reff;
}
void Bar::recoverStresses()
{
	///@todo Stresses from bar to nodes
	static Vector s(6);
	s[0]=myUniMaterial->getStress();
	myNodes[0]->addStress(s);
	myNodes[1]->addStress(s);
}
/**
 * Add a Tracker to the Bar's Material.
 * \a index is checked if is in \a myMatPoints range.
 * @param index The index to the Element's Material.
 */
void Bar::addTracker(int index)
{
	if(index!=1)
		throw SException("[nemesis:%d] %s",9999,"Invalid index.\n");
	myUniMaterial->addTracker();
}
/**
 * Get the Tracker with \a index.
 * \a index is checked if is 1.
 * An exception is thrown if no tracker is set.
 * @param index The index to the Element's Material.
 */
Tracker* Bar::getTracker(int index)
{
	if(index!=1)
		throw SException("[nemesis:%d] %s",9999,"Invalid index.\n");
	if(myUniMaterial->getTracker()==0)
		throw SException("[nemesis:%d] No tracker is set for Element %d, index %d.",9999,myID,index);
	return myUniMaterial->getTracker();
}
/**
 * Add a record to the tracker.
 * For all non null trackers records are added.
 */
void Bar::track()
{
		myUniMaterial->track();
}