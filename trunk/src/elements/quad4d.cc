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

#include "elements/quad4d.h"
#include "elements/shape_functions.h"
#include "main/nemesis_debug.h"

double Quad4d::detJ[4];
double Quad4d::shp[4][3][4];
std::vector<int> Quad4d::perm(4);

Quad4d::Quad4d()
{
}
Quad4d::Quad4d(int ID,int Node_1,int Node_2,int Node_3,int Node_4,int MatID)
:Quad4(ID,Node_1,Node_2,Node_3,Node_4,MatID,2,2)
{
	perm[0]=0;
	perm[1]=1;
	perm[2]=2;
	perm[3]=3;
}
Quad4d::~Quad4d()
{
}
const Matrix& Quad4d::getK()
{
	Matrix &K=*myMatrix;
	static Matrix Ba;
	static Matrix Bb;

	this->shapeFunctions();
	K.clear();
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		const Matrix& C=myMatPoints[k]->getMaterial()->getC();
		for(unsigned a=0;a<myNodes.size();a++)
		{
			this->getB(Ba,a,k);
			for(unsigned b=0;b<myNodes.size();b++)
			{
				this->getB(Bb,b,k);
				double dV=(pD->getFac())*detJ[k];
				K.add_BTCB(2*a,2*b,&perm[0],Ba,C,Bb,dV,1.0);
			}
		}
	}
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	return K;
}
const Matrix& Quad4d::getM()
{
	Matrix &M=*myMatrix;
	M.clear();
	double rho=myMaterial->getRho();
	double volume=0.;
	
	this->shapeFunctions();
	for(unsigned k=0;k<myMatPoints.size();k++)
		volume+=detJ[k]*(pD->getFac())*(myMatPoints[k]->get_w()); 
	double mass=rho*volume;
	for(int i=0;i<8;i++) M(i,i)=0.25*mass;
	return M;
}
const Vector& Quad4d::getR()
{
	Vector& R=*myVector;
	static Vector sigma(6);
	static Matrix Ba;

	R.clear();
	if(!(myGroup->isActive()))	return R;
	///@todo check facS, facG
	//double facS=myGroup->getFacS();
	//double facG=myGroup->getFacG();
	double facP=myGroup->getFacP();
	
	this->shapeFunctions();
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		sigma=myMatPoints[k]->getMaterial()->getStress();
		for(unsigned a=0;a<myNodes.size();a++)
		{
			this->getB(Ba,a,k);
			double dV=(pD->getFac())*detJ[k];
			add_BTv(R,2*a,&perm[0],Ba,sigma,dV,1.0);
		}
	}
	R-=facP*P;
	return R;
}
void Quad4d::update()
{
	static Vector u(8);
	static Vector epsilon(6);
	static Matrix Ba;

	if(!(myGroup->isActive()))	return;
	u=this->getDispIncrm();

	this->shapeFunctions();
	// For each material point
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		epsilon.clear();
		for(unsigned a=0;a<myNodes.size();a++)
		{
			this->getB(Ba,a,k);
			///@todo check dV
			//double dV=(pD->getFac())*detJ[k];
			add_Bv(epsilon,2*a,&perm[0],Ba,u,1.0,1.0);
		}
		myMatPoints[k]->getMaterial()->setStrain(epsilon);
	}
}
void Quad4d::shapeFunctions()
{
	shape4(x,shp,detJ);
	// Axisymmetry
	if(pD->getTag()==TAG_DOMAIN_AXISYMMETRIC)
	{
		for(int k=0;k<4;k++)		// matpoints
		{
			double r=0.;
			for(int i=0;i<4;i++)	// nodes
			    r+=x(i,0)*shp[i][0][k];
			detJ[k]*=r;
		}
	}
}
void Quad4d::getB(Matrix& B,int node,int gPoint)
{
	///@todo: check gCrds
	//const double dsq3=0.577350269189626;
	//static const double gCrds[4]={-dsq3,+dsq3,-dsq3,+dsq3};

	B.resize(4,2);
	double Bb1=0.,Bb2=0.,vol=0.;
	for(int k=0;k<4;k++)  //matpoints
	{
		double dV=detJ[k]*(pD->getFac());
		Bb1+=shp[node][1][k]*dV;
		Bb2+=shp[node][2][k]*dV;
		vol+=dV;
	}
	Bb1/=vol;
	Bb2/=vol;
	
    double B0=0.,Bb0=0.;
	// Axisymmetry
	double r=0.;
	if(pD->getTag()==TAG_DOMAIN_AXISYMMETRIC)
	{
		for(int i=0;i<4;i++)		// nodes
			    r+=x(i,0)*shp[i][0][gPoint];
		B0=shp[node][0][gPoint]/r;
		for(int k=0;k<4;k++)		// matpoints
		{
			r=0.;
			for(int i=0;i<4;i++)	// nodes
			    r+=x(i,0)*shp[i][0][k];
			Bb0+=shp[node][0][k]*detJ[k]*(pD->getFac())/vol/r;
		}
	}

	// B-factors
	double B1=shp[node][1][gPoint];
	double B2=shp[node][2][gPoint];
	double B4=(Bb1-B1)/3.;
	double B6=(Bb2-B2)/3.;
	double B7=B2+B6;
	double B10=B4+(Bb0-B0)/3.;
	double B11=B0+B10;
	double B12=B1+B10;

	// B-matrix
	B(0,0)=B12;	B(0,1)=B6;
	B(1,0)=B10;	B(1,1)=B7;
	B(2,0)=B11;	B(2,1)=B6;
	B(3,0)=B2;	B(3,1)=B1;
}
