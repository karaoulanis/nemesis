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

#include <VonMises.h>
#include <VM.h>

VonMises::VonMises()
{
}
VonMises::VonMises(int ID,int elasticID,double s0)
:MultiaxialElastoPlastic(ID,elasticID)
{
	// Material parameters
	MatParams[0]=s0;
	// Yield/potential surfaces
	fSurfaces.push_back(new VM(s0));
	gSurfaces.push_back(new VM(s0));
	// Material tag
	//myTag=TAG_MATERIAL_MOHR_COULOMB;
}
VonMises::~VonMises()
{
}
MultiaxialMaterial* VonMises::getClone()
{
	// Material parameters
	int myID    = this->getID();
	int elID    = myElastic->getID();
	double s0   = MatParams[ 0];
	// Create clone and return
	VonMises* newClone=new VonMises(myID,elID,s0);
	return newClone;
}
