/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "loadcase/nodal_load.h"

NodalLoad::NodalLoad()
{
}
NodalLoad::NodalLoad(int nodeID,int dofID)
	:Load()
{
	// Retrieve node from the domain and check if exists
	myNode=pD->get<Node>(pD->getNodes(),nodeID);
	// Set the dof and check if activated
	dof=dofID-1;
	if(myNode->getActivatedDof(dof)<0) 
		throw SException("[nemesis:%d] %s",9999,"Dof is not activated.");
}
void NodalLoad::apply(double fact,double time)
{
	myNode->addLoad(dof,this->getValue(time),fact);
}
