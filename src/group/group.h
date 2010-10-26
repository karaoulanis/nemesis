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

#ifndef _GROUP_H
#define _GROUP_H

// Included files
#include "domain/domain_object.h"
#include "loadcase/group_state.h"

class GroupState;

class Group: public DomainObject
{
private:
	bool active;
	double facK;
	double facS;
	double facG;
	double facP;
public:
	Group();
	~Group();
	Group(int ID);
	void setDefault();
	void setState(GroupState* g);
	inline bool isActive()				{return active;}
	inline double getFacK()				{return facK;}
	inline double getFacS()				{return facS;}
	inline double getFacG()				{return facG;}
	inline double getFacP()				{return facP;}
};

#endif
