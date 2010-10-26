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

#include "loadcase/group_state.h"

int GroupState::nGroupStates=0;

GroupState::GroupState()
	:DomainObject()
{
}
GroupState::GroupState(int groupID,int active_,
				double facK_,double facS_,double facG_,double facP_)
	:DomainObject(++nGroupStates)
{
	myGroup=pD->get<Group>(pD->getGroups(),groupID);	
	active_==0 ? active=false : active=true;
	facK=facK_;
	facS=facS_;
	facG=facG_;
	facP=facP_;
}
int GroupState::apply()
{
	myGroup->setState(this);
	return 0;
}
