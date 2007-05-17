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

#include <InitialStresses.h>

InitialStresses::InitialStresses()
{
}
InitialStresses::InitialStresses(int groupID_,double h1_,double s1_,double h2_,double s2_,double K0_)
	:InitialCondition()
{
	myTag=TAG_INITIAL_STRESSES;
	theGroupID=groupID_;
	h1=h1_;
	s1=s1_;
	h2=h2_;
	s2=s2_;
	K0=K0_;
	// Check if group exists
	pD->get<Group>(pD->getGroups(),theGroupID);
}
InitialStresses::~InitialStresses()
{
}
int InitialStresses::apply()
{
	ElementContainer& DomainElems=pD->getElements();
	for(ElementIterator ei=DomainElems.begin();ei!=DomainElems.end();ei++)
		ei->second->addInitialStresses(this);
	return 0;
}
