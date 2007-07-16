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

#include <DruckerPrager.h>
#include <DP_in.h>
#include <DP_out.h>
#include <TC.h>

DruckerPrager::DruckerPrager()
{
}
DruckerPrager::DruckerPrager(int ID,int elasticID,double c,double phi,double T)
:MultiaxialElastoPlastic(ID,elasticID)
{
	// Material parameters
	MatParams[0]=c;
	MatParams[1]=phi;
	MatParams[2]=T;
	// Yield/potential surfaces
/*	fSurfaces.push_back(new MC_0(c,phi));
	fSurfaces.push_back(new MC_1(c,phi));
	fSurfaces.push_back(new MC_2(c,phi));
	gSurfaces.push_back(new MC_0(c,phi));
	gSurfaces.push_back(new MC_1(c,phi));
	gSurfaces.push_back(new MC_2(c,phi));
*/	//fSurfaces.push_back(new TC(T));
	//gSurfaces.push_back(new TC(T));
	// Material tag
	myTag=TAG_MATERIAL_MOHR_COULOMB;
}
DruckerPrager::~DruckerPrager()
{
}
MultiaxialMaterial* DruckerPrager::getClone()
{
	// Material parameters
	int myID    = this->getID();
	int elID    = myElastic->getID();
	double c    = MatParams[ 0];
	double phi  = MatParams[ 1];
	double T    = MatParams[ 2];
	// Create clone and return
	DruckerPrager* newClone=new DruckerPrager(myID,elID,c,phi,T);
	return newClone;
}
