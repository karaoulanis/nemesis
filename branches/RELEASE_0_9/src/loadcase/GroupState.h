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

#ifndef _GROUP_STATE_H
#define _GROUP_STATE_H

#include <DomainObject.h>
#include <Domain.h>
#include <Group.h>

class Group;
class GroupState: public DomainObject
{
private:
	static int nGroupStates;
	Group* myGroup;
	bool active;
	double facK;
	double facS;
	double facG;
	double facP;
public:
	GroupState();
	GroupState(int groupID,int active_,
		double facK_,double facS_,double facG_,double facP_);
	int apply();
	inline bool getActive()				{return active;}
	inline double getFacK()				{return facK;}
	inline double getFacS()				{return facS;}
	inline double getFacG()				{return facG;}
	inline double getFacP()				{return facP;}
};

#endif
