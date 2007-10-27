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

#include <SpringMaterialContact.h>

SpringMaterialContact::SpringMaterialContact()
{
}
SpringMaterialContact::SpringMaterialContact(int ID,
											 double Kn,double sy,double gap,
											 double Ks2,double mue2,
											 double Ks3,double mue3)
:SpringMaterial(ID)
{
	// Material parameters
	MatParams[0]=Kn;
	MatParams[1]=sy;
	MatParams[2]=gap;
	MatParams[3]=Ks2;
	MatParams[4]=mue2;
	MatParams[5]=Ks3;
	MatParams[6]=mue3;
	myTag=TAG_MATERIAL_SPRING;

	// Initialize vectors
	nDim=pD->getnDim();
	sTrial.resize(nDim,0.);
	sConvg.resize(nDim,0.);
	eTrial.resize(nDim,0.);
	eTotal.resize(nDim,0.);
	Ct.resize(3,0.);
	Ct[0]=Kn;
	Ct[1]=Ks2;
	Ct[2]=Ks3;

	// Check
	if(sy>0.) 
		throw SException("[nemesis:%d] Sy must be negative or zero. Instead received %f.",9999,sy);
	if(gap>0.) 
		throw SException("[nemesis:%d] Gap must be negative or zero. Instead received %f.",9999,gap);
}
SpringMaterial* SpringMaterialContact::getClone()
{
	// Material parameters
	double Kn  =MatParams[0];
	double sy  =MatParams[1];
	double gap =MatParams[2];
	double Ks2 =MatParams[3];
	double mue2=MatParams[4];
	double Ks3 =MatParams[5];
	double mue3=MatParams[6];
	// Create clone and return
	SpringMaterial* clone=new SpringMaterialContact(myID,Kn,sy,gap,Ks2,mue2,Ks3,mue3);
	return clone;
}
void SpringMaterialContact::setStrain(const Vector& De)
{
	// Get properties
	double sy   =MatParams[1];
	double gap  =MatParams[2];
	double K[]  ={MatParams[0],MatParams[3],MatParams[5]};
	double mue[]={          0.,MatParams[4],MatParams[6]};
	double eElastMin=gap;
	double eElastMax=gap+sy/K[0];
	
	// Find normal
	eTrial=eTotal+De;
	if(eTrial[0]>eElastMin)
	{
		sTrial[0]=0.;
		Ct[0]=0.;
	}
	else if(eTrial[0]<eElastMax)
	{
		sTrial[0]=sy;
		Ct[0]=0.;
	}
	else
	{
		sTrial[0]=K[0]*(eTrial[0]-eElastMin);
		Ct[0]=K[0];
	}
	// Find Tangents
	for(int i=1;i<nDim;i++)
	{
		// Open
		if(sTrial[0]>0.)
		{
			sTrial[i]=0.;
			Ct[i]=0;
		}
		// Closed and sliding
		else if(fabs(mue[i]*sTrial[0])<fabs(K[i]*eTrial[i]))
		{
			sTrial[i]=0.;
			Ct[i]=0;
		}
		// Closed and stuck
		else
		{
			sTrial[i]=K[i]*eTrial[i];
			Ct[i]=K[i];
		}
	}

}
