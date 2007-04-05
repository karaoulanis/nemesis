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

#include <MultiaxialElastic.h>

MultiaxialElastic::MultiaxialElastic()
{
}
MultiaxialElastic::MultiaxialElastic(int ID,double E,double nu,double rho,double aT)
:MultiaxialMaterial(ID,rho,aT)
{
	// Material parameters
	MatParams[0]=E;
	MatParams[1]=nu;
	// Material tag
	myTag=TAG_MATERIAL_MULTIAXIAL_ELASTIC;
}
MultiaxialMaterial* MultiaxialElastic::getClone()
{
	// Material parameters
	double E   =MatParams[ 0];
	double nu  =MatParams[ 1];
	double rho =MatParams[30];
	double aT  =MatParams[31];
	// Create clone and return
	MultiaxialElastic* newClone=new MultiaxialElastic(myID,E,nu,rho,aT);
	return newClone;
}
void MultiaxialElastic::setStrain(const Vector& De)
{
	sTrial=sConvg+(this->getC())*De;
}
const Matrix& MultiaxialElastic::getC()
{
	C.clear();
	// Material parameters
	double E   =MatParams[ 0];
	double nu  =MatParams[ 1];
	// Find and return C
	double Em=E/((1.+nu)*(1.-2*nu));
	C(0,0)=Em*(1.-nu);		C(0,1)=Em*nu;		C(0,2)=Em*nu;		
	C(1,0)=Em*nu;			C(1,1)=Em*(1.-nu);	C(1,2)=Em*nu;		
	C(2,0)=Em*nu;			C(2,1)=Em*nu;		C(2,2)=Em*(1.-nu);	
	C(3,3)=Em*0.5*(1.-2*nu);
	C(4,4)=Em*0.5*(1.-2*nu);
	C(5,5)=Em*0.5*(1.-2*nu);
	return C;
}
void MultiaxialElastic::commit()
{
	sConvg=sTrial;
}
