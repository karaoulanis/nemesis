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

#ifndef _MODELNODE_H
#define _MODELNODE_H

#include <ModelObject.h>
#include <Node.h>
#include <Containers.h>

class ModelNode :public ModelObject
{
protected:
	Node* myNode;
	Constraint* myConstraint;
public:
	// Constructors
	ModelNode();
	ModelNode(const IDContainer& FTable,Node* pNode=0,Constraint* pConstraint=0);
	virtual ~ModelNode();
	
	// Access to data members
	Node* getNode()			{return myNode;}

	virtual void add_R(double factor=1.0)=0;
	virtual void add_uTrial(double factor=1.0);
	virtual void add_vTrial(double factor=1.0);

	virtual void incTrialDisp(const Vector& du)=0;
	virtual void incTrialVecs(const Vector& du,const Vector& dv,const Vector& da)=0;
	virtual void setTrialDisp(const Vector& u)=0;
	virtual void setTrialVecs(const Vector& u,const Vector& v,const Vector& a)=0;
	
	virtual void commit()=0;
	virtual void commitSens(const Vector& ds,int param) {}
	virtual void rollback();
};

#endif
