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

#include <ModifiedCamClay.h>
#include <MCC.h>

ModifiedCamClay::ModifiedCamClay()
{
}
ModifiedCamClay::ModifiedCamClay(int ID,int elasticID,double M,double po,double kappa,double lambda)
:MultiaxialElastoPlastic(ID,elasticID)
{
	// Material parameters
	MatParams[0]=M;
	MatParams[1]=po;
	MatParams[2]=kappa;
	MatParams[3]=lambda;
	// Yield/potential surfaces
	fSurfaces.push_back(new MCC(M,po,kappa,lambda));
	gSurfaces.push_back(new MCC(M,po,kappa,lambda));
	// Material tag
	myTag=TAG_NONE;
}
ModifiedCamClay::~ModifiedCamClay()
{
}
MultiaxialMaterial* ModifiedCamClay::getClone()
{
	// Material parameters
	int myID      = this->getID();
	int elID      = myElastic->getID();
	double M      = MatParams[ 0];
	double po     = MatParams[ 1];
	double kappa  = MatParams[ 2];
	double lambda = MatParams[ 3];
	// Create clone and return
	ModifiedCamClay* newClone=new ModifiedCamClay(myID,elID,M,po,kappa,lambda);
	return newClone;
}
