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

#ifndef _FEOBJECT_H
#define _FEOBJECT_H

#include <FEObjectTags.h>
#include <Database.h>
#include <Packet.h>
#include <vector>

class FEObject
{
protected:
	FEObjectTag myTag;
	std::vector<FEObjectTag> consistentTags;
	char myLabel[64];
public:
	FEObject();
	FEObject(FEObjectTag tag);
	virtual ~FEObject();
	
	FEObjectTag getTag();
	virtual bool checkIfAllows(FEObject* f);
	virtual void checkConstency(FEObject* f);
};

#endif
