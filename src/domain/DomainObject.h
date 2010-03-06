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

#ifndef _DOMAINOBJECT_H
#define _DOMAINOBJECT_H

#include <FEObject.h>
#include <iostream>

// Forward declarations
class Domain;

/**
 * The DomainObject Class.                                                
 */
class DomainObject: public FEObject
{
protected:
	int myID;
	static Packet thePacket;
	static Domain* pD;
public:
	// Constructors
	DomainObject();		
	DomainObject(int ID);	
	virtual ~DomainObject();

	virtual int getID();

	virtual const Packet& getPacket();
	virtual void setPacket(const Packet& p);
	virtual void save(std::ostream& s)	{}
	virtual void load(std::istream& s)	{}

	void setDomain(Domain* pDomain);
};
#endif
