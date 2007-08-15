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

#include <UniaxialLoad.h>

/******************************************************************************
* Uniaxial Load
******************************************************************************/
UniaxialLoad::UniaxialLoad()
{
}
UniaxialLoad::UniaxialLoad(int elemID,const char* dir)
:ElementalLoad(elemID)
{
	double xA=myElement->getNodes()[0]->getx1();
	double yA=myElement->getNodes()[0]->getx2();
	double xB=myElement->getNodes()[1]->getx1();
	double yB=myElement->getNodes()[1]->getx2();
	L=sqrt((xB-xA)*(xB-xA)+(yB-yA)*(yB-yA));
	double c=(xB-xA)/L;
	double s=(yB-yA)/L;
	myDirection.resize(3,0.);
	if(!strcmp(dir,"LOCX"))
	{ myDirection[0]= c;  myDirection[1]= s;  myDirection[2]= 0; }
	else if(!strcmp(dir,"LOCY"))
	{ myDirection[0]=-s;  myDirection[1]= c;  myDirection[2]= 0; }
	else if(!strcmp(dir,"POSX"))
	{ myDirection[0]= 1;  myDirection[1]= 0;  myDirection[2]= 0; }
	else if(!strcmp(dir,"POSY"))
	{ myDirection[0]= 0;  myDirection[1]= 1;  myDirection[2]= 0; }
	else
		throw SException("[nemesis:%d] %s",9999,"Unknown/misspelled direction."); 
	projections.resize(2,0.);
	projections[0]=	myDirection[0]*(+c)+
					myDirection[1]*(+s);
	projections[1]=	myDirection[0]*(-s)+
					myDirection[1]*(+c);
}
UniaxialLoad::~UniaxialLoad()
{
}

/******************************************************************************
* Beam load point
******************************************************************************/
BeamLoadPoint::BeamLoadPoint()
{
}
BeamLoadPoint::BeamLoadPoint(int elemID,const char* dir,double a0_,double p0_)
:UniaxialLoad(elemID,dir)
{
	if(a0_<0.||a0_>1.) 
		throw SException("[nemesis:%d] %s",9999,"a must be between 0. and 1.");
	a0=a0_;
	p0=p0_;
}
BeamLoadPoint::~BeamLoadPoint()
{
}
const Vector& BeamLoadPoint::getP()
{
	P.resize(6,0.);
	P[0]=projections[0]*p0*(1-a0);
	P[1]=projections[1]*p0*(1-3*a0*a0+2*a0*a0*a0);
	P[2]=projections[1]*p0*(a0*L-2*a0*a0*L+a0*a0*a0*L);
	P[3]=projections[0]*p0*a0;
	P[4]=projections[1]*p0*(3*a0*a0-2*a0*a0*a0);
	P[5]=projections[1]*p0*(-a0*a0*L+a0*a0*a0*L);
	return P;
}

/******************************************************************************
* Beam load uniform
******************************************************************************/
BeamLoadUniform::BeamLoadUniform()
{
}
BeamLoadUniform::BeamLoadUniform(int elemID,const char* dir,double p0_)
:UniaxialLoad(elemID,dir)
{
	p0=p0_;
}
BeamLoadUniform::~BeamLoadUniform()
{
}
const Vector& BeamLoadUniform::getP()
{
	P.resize(6,0.);
	P[0]= projections[0]*( p0*L/2.);
	P[1]=-projections[1]*(-p0*L/2.);
	P[2]=-projections[1]*(-p0*L*L/12.);
	P[3]= projections[0]*( p0*L/2.);
	P[4]=-projections[1]*(-p0*L/2.);
	P[5]=-projections[1]*(-p0*L*L/12.);
	return P;
}
