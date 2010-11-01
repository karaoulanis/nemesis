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

#include "material/uniaxial_cyclic.h"

UniaxialCyclic::UniaxialCyclic()
{
}
UniaxialCyclic::UniaxialCyclic(int ID,double E,double nu,double rho,double aT,double tmax,double Gmax)
:UniaxialMaterial(ID,rho,aT)
{
	// Material parameters
	MatParams[0]=E;
	MatParams[1]=nu;
	MatParams[2]=tmax;
	MatParams[3]=Gmax;
	// State variables
	sr=0.;
	er=0.;
	eTrial=0.;
	eConvg=0.;
	// Initial stiffness
	Et=Gmax;
	reversed=false;
	// Tag
	myTag=TAG_MATERIAL_UNIAXIAL_CYCLIC;
}
UniaxialCyclic::~UniaxialCyclic()
{
}
UniaxialMaterial* UniaxialCyclic::getClone()
{
	double E   =MatParams[0];
	double nu  =MatParams[1];
	double tmax=MatParams[2];
	double Gmax=MatParams[3];
	double rho =MatParams[30];
	double aT  =MatParams[31];
	UniaxialMaterial* newClone=new UniaxialCyclic(myID,E,nu,rho,aT,tmax,Gmax);
	return newClone;
}	
void UniaxialCyclic::setStrain(const double De)
{
	double tmax=MatParams[2];
	double Gmax=MatParams[3];
	// At first step check if reversed
	if(((eConvg+De-er))*De<0)
	{
		sr=sConvg;
		er=eConvg;
		reversed=true;
		cout<<er<<'\t'<<sr<<endl;
	}
	eTrial=eConvg+De;
	//cout<<er<<'\t'<<eTrial<<'\t'<<reversed<<endl;
	if(reversed==false)
	{
		sTrial=Gmax*eTrial/(1+(Gmax/tmax)*fabs(eTrial));
		Et=Gmax*tmax*(tmax+Gmax*fabs(eTrial)-Gmax*eTrial*num::sign(eTrial))/((tmax+Gmax*abs(eTrial))*(tmax+Gmax*fabs(eTrial)));
	}
	else
	{
		sTrial=sr+2.0*Gmax*(0.5*(eTrial-er))/(1+(Gmax/tmax)*fabs(0.5*(eTrial-er)));
		Et=0.5*Gmax*tmax*(2*tmax+Gmax*fabs(er-eTrial)-0.5*Gmax*num::sign(er-eTrial)*(er-eTrial))
			   /((tmax+0.5*Gmax*abs(er-eTrial))*(tmax+0.5*Gmax*abs(er-eTrial)));
	}
}
double UniaxialCyclic::getC()
{
	///@todo
	return Et;
}
void UniaxialCyclic::commit()
{
	sConvg=sTrial;
	eConvg=eTrial;
	eTotal=eConvg;
	this->track();
}

