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

#include <DuncanChang.h>

DuncanChang::DuncanChang()
{
}
DuncanChang::DuncanChang(int ID,double E,double nu,double c,double phi,
						 double m,double Rf,double pa,double rho,double aT)
:MultiaxialMaterial(ID,rho,aT)
{
	// Material parameters
	MatParams[0]=E;
	MatParams[1]=nu;
	MatParams[2]=m;
	MatParams[3]=c;
	MatParams[4]=phi;
	MatParams[5]=Rf;
	MatParams[6]=pa;
	// Material tag
	myTag=TAG_MATERIAL_MULTIAXIAL_ELASTIC;
}
MultiaxialMaterial* DuncanChang::getClone()
{
	// Material parameters
	double E   =MatParams[ 0];
	double nu  =MatParams[ 1];
	double m   =MatParams[ 2];
	double c   =MatParams[ 3];
	double phi =MatParams[ 4];
	double Rf  =MatParams[ 5];
	double pa  =MatParams[ 6];
	double rho =MatParams[30];
	double aT  =MatParams[31];
	// Create clone and return
	DuncanChang* newClone=new DuncanChang(myID,E,nu,c,phi,m,Rf,pa,rho,aT);
	return newClone;
}
void DuncanChang::setStrain(const Vector& De)
{
	eTrial=eTotal+De;
	sTrial=(this->getC())*eTrial;
}
const Matrix& DuncanChang::getC()
{
	C.clear();
	const Vector& s=sTrial.eigenvalues();
	double s1=-s[2];
	double s3=-s[0];

	// Material parameters
	double E   =MatParams[ 0];
	double nu  =MatParams[ 1];
	double m   =MatParams[ 2];
	double c   =MatParams[ 3];
	double phi =MatParams[ 4];
	double Rf  =MatParams[ 5];
	double pa  =MatParams[ 6];
	double d=1-Rf*(1-sin(phi)*(s1-s3))/(2*c*cos(phi)+2*s3*sin(phi));
	double Et=d*d*E*pa*pow(s3/pa,m);
	if(num::tiny(Et)) Et=E;
	cout<<d*d<<endl;

	// Find and return C
	double Em=Et/((1.+nu)*(1.-2*nu));
	C(0,0)=Em*(1.-nu);		C(0,1)=Em*nu;		C(0,2)=Em*nu;		
	C(1,0)=Em*nu;			C(1,1)=Em*(1.-nu);	C(1,2)=Em*nu;		
	C(2,0)=Em*nu;			C(2,1)=Em*nu;		C(2,2)=Em*(1.-nu);	
	C(3,3)=Em*0.5*(1.-2*nu);
	C(4,4)=Em*0.5*(1.-2*nu);
	C(5,5)=Em*0.5*(1.-2*nu);
	return C;
}
void DuncanChang::commit()
{
	eTotal=eTrial;
	sConvg=sTrial;
	this->track();
}
/**
 * Add a record to the tracker.
 * If \a myTracker pointer is null (no tracker is added) just return.
 * Otherwise gather info and send them to the tracker.
 * The domain should be already updated!
 */
void DuncanChang::track()
{
	if(myTracker==0) return;
	ostringstream s;
	s<<"DATA "	<<' ';
	s<<"sigm "	<<' '<<sConvg;
	s<<"epst "	<<' '<<eTotal;
//	s<<"epsp "	<<' '<<ePConvg;	///@todo
//	s<<"epsv "	<<1020<<' '<<eTotal[0]+eTotal[1]+eTotal[2]<<' ';
	s<<"p "	    <<1020<<' '<<sConvg.p()<<' ';
	s<<"q "	    <<1020<<' '<<sConvg.q()<<' ';
	s<<"END "<<' ';
	myTracker->track(pD->getLambda(),pD->getTimeCurr(),s.str());
}