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

#ifndef _LOADCASE_H
#define _LOADCASE_H

// Included files
#include <DomainObject.h>
#include <NodalLoad.h>
#include <Load.h>
#include <ElementalLoad.h>
#include <GroupState.h>
#include <InitialCondition.h>
#include <ElementSensitivityParameter.h>
#include <Containers.h>

// Forward declarations
class NodalLoad;
class InitialCondition;
class Load;
class GroupState;
class Group;
class ElementSensitivityParameter;

// Type definitions 
typedef std::vector<Load*>								LoadVector;
typedef std::vector<InitialCondition*>					InitialConditionVector;
typedef std::vector<GroupState* >						GroupStateVector;
typedef std::vector<ElementSensitivityParameter* >		SensitivityParameterVector;
typedef std::map<int,Group*>							GroupContainer;

class LoadCase: public DomainObject
{
private:
	char myLabel[512];
	double myFac;
	int active;
	LoadVector					myLoads;
	GroupStateVector			myGroupStates;
	InitialConditionVector		myInitialConditions;
	SensitivityParameterVector	mySensitivityParameters;
public:
	// Constructors and destructor
	LoadCase();
	LoadCase(int ID,const char* label);
	~LoadCase();

	// Initialize loadcase
	void init();
	void commit();

	// Add member functions
	void addLoad(Load* pLoad);
	void addGroupState(GroupState* pGroupState);
	void addInitialCondition(InitialCondition* pInitialCondition);
	void addSensitivityParameter(ElementSensitivityParameter* pElementSensitivityParameter);

	// Apply member functions
	void applyLoads(double lambda_,double time_);
	void applySensitivityParameter(int param);
	int getnSensitivityParameters();
};

#endif
