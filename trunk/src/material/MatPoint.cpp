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

#include <MatPoint.h>

int MatPoint::IDCounter=0;

MatPoint::MatPoint()
{
	delete myMaterial;
}
MatPoint::MatPoint(MultiaxialMaterial* mat,double xi1_,double xi2_,double xi3_,double w_)
:DomainObject(IDCounter++)
{
	xi1=xi1_;
	xi2=xi2_;
	xi3=xi3_;
	weight=w_;
	myMaterial=mat->getClone();
	myTag=TAG_MATERIAL_POINT;
}
MatPoint::MatPoint(MultiaxialMaterial* mat,int index,int p1)
:DomainObject(IDCounter++)
{
	if(p1>6) exit(-987);
	xi1=GaussCoords[p1][index+1];
	xi2=0;
	xi3=0;
	weight=GaussWeights[p1][index+1];
	myMaterial=mat->getClone();
	myTag=TAG_MATERIAL_POINT;
}
MatPoint::MatPoint(MultiaxialMaterial* mat,int index1,int index2,int p1,int p2)
:DomainObject(IDCounter++)
{
	if(p1>6||p2>6) exit(-987);
	xi1=GaussCoords[p1][index1];
	xi2=GaussCoords[p2][index2];
	xi3=0.;
	weight=GaussWeights[p1][index1]*GaussWeights[p2][index2];
	myMaterial=mat->getClone();
	myTag=TAG_MATERIAL_POINT;
}
MatPoint::MatPoint(MultiaxialMaterial* mat,int index1,int index2,int index3,int p1,int p2,int p3)
:DomainObject(IDCounter++)
{
	if(p1>6||p2>6) exit(-987);
	xi1=GaussCoords[p1][index1];
	xi2=GaussCoords[p2][index2];
	xi3=GaussCoords[p3][index3];
	weight=GaussWeights[p1][index1]*GaussWeights[p2][index2]*GaussWeights[p3][index3];
	myMaterial=mat->getClone();
	myTag=TAG_MATERIAL_POINT;
}
MatPoint::~MatPoint()
{
	delete myMaterial;
}
void MatPoint::setX(double x1_,double x2_,double x3_)
{
	x1=x1_;
	x2=x2_;
	x3=x3_;
}
void MatPoint::setInitialStresses(InitialStresses* pInitialStresses)
{
	static Vector s0(6);
	s0.clear();
	double H1=pInitialStresses->getH1();
	double H2=pInitialStresses->getH2();
	double s1=pInitialStresses->getS1();
	double s2=pInitialStresses->getS2();
	if(x2<H1 && x2>H2)
	{
        s0[1]=s1+(s2-s1)*(x2-H1)/(H2-H1);
		s0[0]=pInitialStresses->getK0()*s0[1];
		s0[2]=pInitialStresses->getK0()*s0[1];
		myMaterial->addStress(s0);
	}
}

const double MatPoint::GaussCoords[7][7]=	   {{ 0.000000000000000,  // Rule 0
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000}, //_______
 												{ 0.000000000000000,  // Rule 1
                                                  0.000000000000000,  // [1][1] 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 2
                                                 -0.577350269189626,  // [2][1]
                                                 +0.577350269189626,  // [2][2]
                                                  0.000000000000000,
                                                  0.000000000000000,
                                                  0.000000000000000,
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 3
                                                 -0.774596669241483,  // [3][1]
                                                  0.000000000000000,  // [3][2]
                                                 +0.774596669241483,  // [3][3]
                                                  0.000000000000000,
                                                  0.000000000000000,
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 4
                                                 -0.861136311594053,  // [4][1]
                                                 -0.339981043584856,  // [4][2]
                                                 +0.339981043584856,  // [4][3]
                                                 +0.861136311594053,  // [4][4]
                                                  0.000000000000000,
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 5
                                                 -0.906179845938664,  // [5][1]
                                                 -0.538469310105683,  // [5][2]
                                                  0.000000000000000,  // [5][3]
                                                 +0.538469310105683,  // [5][4]
                                                 +0.906179845938664,  // [5][5]
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 6
                                                 -0.932469514203152,  // [6][1]
                                                 -0.661209386466265,  // [6][2]
                                                 -0.238619186083197,  // [6][3]
                                                  0.238619186083197,  // [6][4]
                                                  0.661209386466265,  // [1][5]
                                                  0.932469514203152}};// [6][6]

const double MatPoint::GaussWeights[7][7]=     {{ 0.000000000000000,  // Rule 0
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 1
                                                 +2.000000000000000,  // [1][1] 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 2
                                                 +1.000000000000000,  // [2][1]
                                                 +1.000000000000000,  // [2][1]
                                                  0.000000000000000,
                                                  0.000000000000000,
                                                  0.000000000000000,
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 3
                                                 +0.555555555555556,  // [3][1]
                                                 +0.888888888888889,  // [3][2]
                                                 +0.555555555555556,  // [3][3]
                                                  0.000000000000000,
                                                  0.000000000000000,
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 4
                                                 +0.347854845137454,  // [4][1]
                                                 +0.652145154862546,  // [4][2]
                                                 +0.652145154862546,  // [4][3]
                                                 +0.347854845137454,  // [4][4]
                                                  0.000000000000000,
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 5
                                                 +0.236926885056189,  // [5][1]
                                                 +0.478628670499366,  // [5][2]
                                                 +0.568888888888889,  // [5][3]
                                                 +0.478628670499366,  // [5][4]
                                                 +0.236926885056189,  // [5][5]
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 6
                                                 +0.171324492379170,  // [6][1]
                                                 +0.360761573048139,  // [6][2]
                                                 +0.467913934572691,  // [6][3]
                                                 +0.467913934572691,  // [6][4]
                                                 +0.360761573048139,  // [1][5]
                                                 +0.171324492379170}};// [6][6]
