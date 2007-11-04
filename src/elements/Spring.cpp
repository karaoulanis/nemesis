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

#include <Spring.h>

/**
 * Default constructor.
 */
Spring::Spring()	
{
}
/**
 * Constructor.
 * Creates a Bar Element. 
 */
Spring::Spring(int ID,int Node_1,int Node_2,int matID,
			   double xp1,double xp2,double xp3,
		       double yp1,double yp2,double yp3)
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
	
	// Store material information
	mySpringMaterial=static_cast<SpringMaterial*>(myMaterial)->getClone();

	// Define transformations
	T.resize(nDim,nDim,0.);
	static Vector xp,yp,zp;
	switch(nDim)
	{
	case 1:
		T(0,0)=1.;
		break;
	case 2:
		xp.resize(2);
		yp.resize(2);
		xp[0]= xp1; xp[1]= xp2;
		yp[0]=-xp2; yp[1]= xp1;
		xp.normalize();
		yp.normalize();
		T.appendRow(xp,0,0);
		T.appendRow(yp,1,0);
		break;
	case 3:
		xp.resize(3);
		yp.resize(3);
		zp.resize(3);
		xp[0]=xp1; xp[1]=xp2; xp[2]=xp3;
		yp[0]=yp1; yp[1]=yp2; yp[2]=yp3;
		zp=cross(xp,yp);
		yp=cross(zp,xp);
		xp.normalize();
		yp.normalize();
		zp.normalize();
		T.appendRow(xp,0,0);
		T.appendRow(yp,1,0);
		T.appendRow(zp,2,0);
		break;
	default:
		break;
	}
}
/**
 * Destructor.
 */
Spring::~Spring()
{
	delete mySpringMaterial;
}
void Spring::update()
{
	static Vector du(2*nDim);
	du=this->getDispIncrm();
	static Vector de(nDim,0);
	de.clear();
	for(int k=0;k<nDim;k++)
		for(int i=0;i<nDim;i++)
			de[i]+=T(i,k)*(du[k+nDim]-du[k]);
	mySpringMaterial->setStrain(de);
}
void Spring::commit()
{
	mySpringMaterial->commit();
}
bool Spring::checkIfAllows(FEObject* f)
{
	return true;
}
const Matrix& Spring::getK()
{
	Matrix& K=*myMatrix;
	K.clear();
	Vector Ct=mySpringMaterial->getC();
	for(int k=0;k<nDim;k++)
	{
		for(int j=0;j<nDim;j++)
		{
			for(int i=0;i<nDim;i++)
			{
				K(i     ,j     )=K(i     ,j     )+T(k,i)*Ct[k]*T(k,j);
				K(i     ,j+nDim)=K(i     ,j+nDim)-T(k,i)*Ct[k]*T(k,j);
				K(i+nDim,j     )=K(i+nDim,j     )-T(k,i)*Ct[k]*T(k,j);
				K(i+nDim,j+nDim)=K(i+nDim,j+nDim)+T(k,i)*Ct[k]*T(k,j);
			}
		}
	}
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	return K;
}
const Matrix& Spring::getM()
{   
	Matrix& M=*myMatrix;
	M.clear();
	return M;
}
const Vector& Spring::getR()
{
	Vector& R=*myVector;
	R.clear();
	// Factors
	if(!(myGroup->isActive()))	return R;
	double facS=myGroup->getFacS();
	double facG=myGroup->getFacG();
	double facP=myGroup->getFacP();
	// R = Fint - Fext 
	Vector sigma=mySpringMaterial->getStress();
	for(int k=0;k<nDim;k++)
	{
		for(int i=0;i<nDim;i++)
		{
			double d=facS*T(k,i)*sigma[k];
			R[i     ]+=-d;
			R[i+nDim]+= d;
		}
	}
	return R;
}
const Vector& Spring::getReff()
{
	///@todo: form Reff
	Vector& Reff=*myVector;
	Reff.clear();
	return Reff;
}
const Vector& Spring::getRgrad()
{
	///@todo: form Rgrad
	Vector& Rgrad=*myVector;
	Rgrad.clear();
	return Rgrad;
}
void Spring::recoverStresses()
{
	///@todo Stresses from bar to nodes
	static Vector s(6);
	s.clear();
	s.append(mySpringMaterial->getStress(),0);
	myNodes[0]->addStress(s);
	myNodes[1]->addStress(s);
}
/**
 * Add a Tracker to the Bar's Material.
 * \a index is checked if is in \a myMatPoints range.
 * @param index The index to the Element's Material.
 */
void Spring::addTracker(int index)
{
	if(index!=1)
		throw SException("[nemesis:%d] %s",9999,"Invalid index.\n");
	mySpringMaterial->addTracker();
}
/**
 * Get the Tracker with \a index.
 * \a index is checked if is 1.
 * An exception is thrown if no tracker is set.
 * @param index The index to the Element's Material.
 */
Tracker* Spring::getTracker(int index)
{
	if(index!=1)
		throw SException("[nemesis:%d] %s",9999,"Invalid index.\n");
	if(mySpringMaterial->getTracker()==0)
		throw SException("[nemesis:%d] No tracker is set for Element %d, index %d.",9999,myID,index);
	return mySpringMaterial->getTracker();
}
/**
 * Add a record to the tracker.
 * For all non null trackers records are added.
 */
void Spring::track()
{
		mySpringMaterial->track();
}

