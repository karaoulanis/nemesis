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

#include <Domain.h>

/**
 * Default Constructor.
 */
Domain::Domain()	
{
	nDim=0;
	myTag=TAG_DOMAIN_NOT_SET;
	myFac=1.0;
	upToDate=false;
	theDatabase=0;
	this->init();
}
/**
 * Destructor.
 */
Domain::~Domain()	
{
	this->clear();
	Containers::map_delete(theGroups);
	if(theDatabase!=0) delete theDatabase;
}
void Domain::init()
{
	// Define group 0
	Group* pGroup=new Group(0);
	pGroup->setDomain(this);
	this->add(this->getGroups(),pGroup);
	// Define that groups are by set by materials
	groupsByMaterial=true;
	currentGroup=0;
	// Initialize time/lambda
	timeCurr=0.;
	timePrev=0.;
	lambdaConvg=0.;
	// Gravity direction vector and default orientation/acceleration
	gravityVect.resize(3,0.);
	gravityVect[1]=-1;
	gravityAccl=9.81;
}
/**
 * Clear the Domain.
 * The containers and the contained objects are deleted.
 */
void Domain::clear()
{
	nDim=0;
	RayleighFactors.resize(0);
	eigenVals.resize(0);
	Containers::map_delete(theNodes);
	Containers::map_delete(theElements);
	Containers::map_delete(theGroups);
	Containers::map_delete(theCrossSections);
	Containers::map_delete(theMaterials);
	Containers::map_delete(theLoadCases);
	Containers::map_delete(theConstraints);
	this->init();
}
/**
 * Set the number of the dimensions of the Domain.
 * @return An integer implying if something is wrong. 
 */
int Domain::setnDim(int nDimensions)
{
	if(nDim!=0) 
		throw SException("[nemesis:%d] %s",9999,"Domain must be cleared before changing dim."); 
	else if((nDimensions<1)||(nDimensions>3))
		throw SException("[nemesis:%d] %s",9999,"Dim value 1,2 or 3 is only allowed.");
	else nDim=nDimensions;
	return 0;
}
/**
 * Access to the number of the dimensions of the Domain.
 * @return The number of the dimensions of the Domain.
 */
int Domain::getnDim() const
{
	return nDim;
}
void Domain::setDatabase(Database* pDB)
{
	if(theDatabase!=0) delete theDatabase;
	theDatabase=pDB;
}
Database* Domain::getDatabase()
{
	return theDatabase;
}
void Domain::zeroNodalStress()
{
	for(NodeIterator nIter=theNodes.begin();
		nIter!=theNodes.end();nIter++) 
			nIter->second->zeroStress();
}
void Domain::state(double facD)
{
	for(NodeIterator nIter=theNodes.begin();
		nIter!=theNodes.end();nIter++) 
			nIter->second->multDisp(facD);
}
void Domain::zeroSensitivityParameters()
{
	for(ElementIterator eIter=theElements.begin();
		eIter!=theElements.end();eIter++)
			eIter->second->activateParameter(0);
}
int Domain::storeState(const char* tableName)
{
	// Begin transaction
	theDatabase->beginTransaction();
	theDatabase->createTable(tableName);
	theDatabase->useTable(tableName);
	// Store nodes
	for(NodeIterator nIter=theNodes.begin();
		nIter!=theNodes.end();nIter++) 
		if(nIter->second->isActive())
			theDatabase->storeData(nIter->second->getPacket());
	// Store elements
	for(ElementIterator eIter=theElements.begin();
		eIter!=theElements.end();eIter++)
		if(eIter->second->isActive())
			theDatabase->storeData(eIter->second->getPacket());
	// Commit transaction
	theDatabase->commitTransaction();
	return 0;
}
int Domain::restoreState(const char* tableName)
{
	// Begin transaction
	theDatabase->beginTransaction();
	theDatabase->useTable(tableName);
	// Restore the nodes
	for(NodeIterator nIter=theNodes.begin();
		nIter!=theNodes.end();nIter++)
			nIter->second->setPacket(theDatabase->retrieveData(nIter->second->getTag(),nIter->second->getID()));
	theDatabase->commitTransaction();
	return 0;
}
// Rayleigh damping
void Domain::setRayleighFactors(const Vector& factors)
{
	RayleighFactors.resize(factors.size());
	RayleighFactors=factors;
}
const Vector& Domain::getRayleighFactors()
{
	return RayleighFactors;
}
// EigenValues
void Domain::setEigenValues(const Vector& vals)
{
	eigenVals.resize(vals.size());
	eigenVals=vals;
}
const Vector& Domain::getEigenValues()
{
	return eigenVals;
}
void Domain::applyLoads(double lambda_,double time_)
{
	for(LoadCaseIterator lcIter=theLoadCases.begin();
		lcIter!=theLoadCases.end();lcIter++) 
			lcIter->second->applyLoads(lambda_,time_);
}
void Domain::zeroLoads()
{
	for(NodeIterator nIter=theNodes.begin();
		nIter!=theNodes.end();nIter++) 
			nIter->second->zeroLoad();
	for(ElementIterator eIter=theElements.begin();
		eIter!=theElements.end();eIter++) 
			eIter->second->zeroLoad();
}
void Domain::zeroGroups()
{
	GroupIterator gi;
	for(gi=theGroups.begin();gi!=theGroups.end();gi++) 
		gi->second->setDefault();
}
/**
 * Sets gravity info.
 * gravityVect Vector giving the direction of gravity. It does not include 
 * special cases for 1D, 2D or 3D cases. It might not be normalized; it is
 * normalized here and remains so thereafter.
 * @param g  The gravity acceleration.
 * @param xG x-coordinate of gravity vector.
 * @param yG y-coordinate of gravity vector.
 * @param zG z-coordinate of gravity vector.
 */
void Domain::setGravity(double g,double xG,double yG,double zG)
{
	gravityAccl=g;
	gravityVect[0]=xG;
	gravityVect[1]=yG;
	gravityVect[2]=zG;
	gravityVect.normalize();
}
/**
 * Returns gravity vector.
 * @return A constant reference to a 3x1 normalized Vector holding gravity
 * directions.
 */
const Vector& Domain::getGravityVect()
{
	return gravityVect;
}
/**
 * Returns gravity acceleration.
 * @return Gravity acceleration
 */
const double Domain::getGravityAccl()
{
	return gravityAccl;
}
