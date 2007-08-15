/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2007 F.E.Karaoulanis [http://www.nemesis-project.org]  *
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

#include <MohrCoulomb.h>
#include <MC_0.h>
#include <MC_1.h>
#include <MC_2.h>
#include <TC.h>

MohrCoulomb::MohrCoulomb()
{
}
MohrCoulomb::MohrCoulomb(int ID,int elasticID,double c,double phi,double T)
:MultiaxialElastoPlastic(ID,elasticID)
{
	// Material parameters
	MatParams[0]=c;
	MatParams[1]=phi;
	MatParams[2]=T;
	// Yield/potential surfaces
	fSurfaces.push_back(new MC_0(c,phi));
	fSurfaces.push_back(new MC_1(c,phi));
	fSurfaces.push_back(new MC_2(c,phi));
	gSurfaces.push_back(new MC_0(c,phi));
	gSurfaces.push_back(new MC_1(c,phi));
	gSurfaces.push_back(new MC_2(c,phi));
	//fSurfaces.push_back(new TC(T));
	//gSurfaces.push_back(new TC(T));
	// Material tag
	myTag=TAG_MATERIAL_MOHR_COULOMB;
}
MohrCoulomb::~MohrCoulomb()
{
}
MultiaxialMaterial* MohrCoulomb::getClone()
{
	// Material parameters
	int myID    = this->getID();
	int elID    = myElastic->getID();
	double c    = MatParams[ 0];
	double phi  = MatParams[ 1];
	double T    = MatParams[ 2];
	// Create clone and return
	MohrCoulomb* newClone=new MohrCoulomb(myID,elID,c,phi,T);
	return newClone;
}
