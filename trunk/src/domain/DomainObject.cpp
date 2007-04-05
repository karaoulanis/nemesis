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

#include <DomainObject.h>

Domain* DomainObject::pD=0;
Packet DomainObject::thePacket;

/**
 * Default Constructor.
 */
DomainObject::DomainObject()
{
	// Does nothing.
}
/**
 * Constructor. 
 * Construct a DomainObject for a given ID.
 * @param ID An integer that initialize a DomainObject.
 */
DomainObject::DomainObject(int ID)
{
	myID=ID;
}
DomainObject::~DomainObject()
{
	// Nothing to destruct here.
}
int DomainObject::getID()					
{
	return myID;
}
const Packet& DomainObject::getPacket()
{
	///@ todo When finished implementing all turn this function into pure.
	std::cout<<"DomainObject::Not implemented yet!"<<std::endl;
	thePacket.zero();
	return thePacket;
}
void DomainObject::setPacket(const Packet& p)
{
	///@ todo When finished implementing all turn this function into pure.
	std::cout<<"DomainObject::Not implemented yet!"<<std::endl;
}
void DomainObject::setDomain(Domain* pDomain)
{
	pD=pDomain;
}
