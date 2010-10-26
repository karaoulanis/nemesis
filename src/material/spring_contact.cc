/******************************************************************************
*   nemesis. an experimental finite element code.                             *
*   Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]  *
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

#include "material/spring_contact.h"

SpringContact::SpringContact()
{
}
SpringContact::SpringContact(int ID,double Kn,double Ks,double mu,double gap)
:SpringMaterial(ID)
{
	// Material parameters
	MatParams[0]=Kn;
	MatParams[1]=Ks;
	MatParams[2]=mu;
	MatParams[3]=gap;
	myTag=TAG_MATERIAL_SPRING;

	// Initialize Ct
	Ct.clear();
	Ct(0,0)=Kn;
	Ct(1,1)=Ks;
	Ct(2,2)=Ks;

	if(gap!=0.)
		cout<<"Warning  : Gap is not yet implemented."<<endl;
	///@todo gap.
	//eTotal.clear();
	//eTotal[0]=gap;
}
SpringMaterial* SpringContact::getClone()
{
	// Material parameters
	double Kn  =MatParams[0];
	double Ks  =MatParams[1];
	double mu  =MatParams[2];
	double gap =MatParams[3];

	// Create clone and return
	SpringMaterial* clone=new SpringContact(myID,Kn,Ks,mu,gap);
	return clone;
}
void SpringContact::setStrain(const Vector& De)
{
	double Kn  =MatParams[0];
	double Ks  =MatParams[1];
	double mu  =MatParams[2];
	double gap =MatParams[3];

	eTrial=eTotal+De;

	double Fn =Kn*eTrial[0];
	double Fs2=Ks*eTrial[1];
	double Fs3=Ks*eTrial[2];
	double f1=Fn;
	double f2=sqrt(Fs2*Fs2+Fs3*Fs3)+mu*Fn;

	Ct.clear();
	sTrial.clear();
	// Case 1: Closed and sticking
	if(f1<0. && f2<0.)
	{
		cout<<"Spring : "<<myID<<" Closed and sticking."<<endl;
		sTrial[0]=Kn*eTrial[0];
		sTrial[1]=Ks*eTrial[1];
		sTrial[2]=Ks*eTrial[2];
		Ct(0,0)=Kn;
		Ct(1,1)=Ks;
		Ct(2,2)=Ks;	
	}
	// Case 2: Closed and sliding
	else if(f1<0. && f2>=0.)
	{
		cout<<"Spring : "<<myID<<" Closed and sliding."<<endl;
		sTrial[0]=Kn*eTrial[0];
		sTrial[1]=0.;
		sTrial[2]=0.;
		Ct(0,0)=Kn;
		Ct(1,1)=0.;
		Ct(2,2)=0.;	
	}
	// Case 3: Open
	else
	{
		cout<<"Spring : "<<myID<<" Open."<<endl;
		sTrial[0]=0.;
		sTrial[1]=0.;
		sTrial[2]=0.;
		Ct(0,0)=0.;
		Ct(1,1)=0.;
		Ct(2,2)=0.;	
	}

}
