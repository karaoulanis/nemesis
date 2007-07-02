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

#ifndef _ELEMENT_H
#define _ELEMENT_H

// Included Files
#include <Domain.h>
#include <DomainObject.h>
#include <Node.h>
#include <Group.h>
#include <Matrix.h>
#include <Vector.h>
#include <ElementalLoad.h>
#include <Containers.h>
#include <InitialStresses.h>

// Forward Declerations
class Domain;
class Material;
class Node;
class Group;
class ElementalLoad;
class InitialStresses;

/**
 * The Element Class.                                                
 * The Element class is an abstract class, that provides the interface to all
 * kind of elements. One, two or three dimensional elements should be able to 
 * be derived from this class. 
 */
class Element: public DomainObject
{
protected:
	IDContainer myNodalIDs;
	IDContainer myLocalNodalDofs;
	std::vector<Node*> myNodes;
	Group* myGroup;
	
	// Materials
	Material* myMaterial;

	Vector P;
	Vector G;
	Matrix x;
	Vector b;

	Matrix* myMatrix;
	Vector* myVector;
	static Matrix** theStaticMatrices;
	static Vector** theStaticVectors;
	int handleCommonInfo();
	int activeParameter;

public:
	// Constructors and Destructor
	Element();
	Element(int ID,int matID);
	~Element();

	const IDContainer& getNodalIDs() const;
	const IDContainer& getLocalNodalDofs() const;
	const std::vector<Node*>& getNodes() const;
	
	bool isActive();
	virtual const int getnPlasticPoints();

	// Build and return global element matrices
    virtual const Matrix& getK()=0;
    virtual const Matrix& getM()=0;
    virtual const Matrix& getC();
    virtual const Vector& getR()=0;
    virtual const Vector& getReff()		{myVector->clear(); return *myVector;}
	virtual const Vector& getRgrad()	{myVector->clear(); return *myVector;}

	// Handle elemental loads
	void addLoad(const Vector& val,double fac=1.0);
	void zeroLoad();
	void addGroundMotion(int dof,double val);	
	virtual void addInitialStresses(InitialStresses* pInitialStresses);	

	virtual void update()=0;
	virtual void commit()=0;

	virtual const Vector& getDispTrial();
	virtual const Vector& getVelcTrial();
	virtual const Vector& getAcclTrial();
	virtual const Vector& getDispConvg();
	virtual const Vector& getVelcConvg();
	virtual const Vector& getAcclConvg();
	virtual const Vector& getDispIncrm();

	void setGroup(int groupID);
	virtual void recoverStresses() {}
 
	// Send and receive packet
	const Packet& getPacket();
	void setPacket(const Packet& p);
	void save(std::ostream& s);

	void activateParameter(int param) {activeParameter=param;}
};

#endif
