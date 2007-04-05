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

#include <UniaxialCyclic.h>

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
	static bool reversed=false;
	double tmax=MatParams[2];
	double Gmax=MatParams[3];
	double e0=eTrial;
	double s0=sTrial;
	double en=eConvg+De;
	eTrial=en;
	if((en-er)*De<0)
	{
		sr=s0;
		er=e0;
		if(reversed==false) reversed=true;
	}
	if(reversed==false)
		sTrial=Gmax*en/(1+(Gmax/tmax)*fabs(en));
	else
		sTrial=sr+2.0*Gmax*(0.5*(en-er))/(1+(Gmax/tmax)*abs(0.5*(en-er)));
}
const double UniaxialCyclic::getC()
{
	///@todo
	return 0;
}
void UniaxialCyclic::commit()
{
	sTrial=sConvg;
	eTrial=eConvg;
}

