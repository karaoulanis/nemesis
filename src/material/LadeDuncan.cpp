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

#include <LadeDuncan.h>
#include <LadeDuncanSurface.h>

LadeDuncan::LadeDuncan()
{
}
LadeDuncan::LadeDuncan(int ID,int elasticID,double K)
:MultiaxialElastoPlastic(ID,elasticID)
{
	// Material parameters
	MatParams[0]=K;
	// Yield/potential surfaces
	fSurfaces.push_back(new LadeDuncanSurface(K));
	gSurfaces.push_back(new LadeDuncanSurface(K));
	// Material tag
	//myTag=TAG_MATERIAL_MOHR_COULOMB;
	nHardeningVariables=1;
}
LadeDuncan::~LadeDuncan()
{
}
MultiaxialMaterial* LadeDuncan::getClone()
{
	// Material parameters
	int myID    = this->getID();
	int elID    = myElastic->getID();
	double K    = MatParams[ 0];
	// Create clone and return
	LadeDuncan* newClone=new LadeDuncan(myID,elID,K);
	return newClone;
}