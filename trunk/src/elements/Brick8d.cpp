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

#include <Brick8d.h>
#include <ShapeFunctions.h>
#include <NemesisDebug.h>

double Brick8d::detJ[8];
double Brick8d::shp[8][4][8];

Brick8d::Brick8d()
{
}
Brick8d::Brick8d(int ID,
					   int Node_1,int Node_2,int Node_3,int Node_4,
					   int Node_5,int Node_6,int Node_7,int Node_8,
					   int matID)
:Element(ID,matID)
{
	myTag=TAG_ELEM_BRICK_8_DISP;
	// Get nodal data
	myNodalIDs.resize(8);
	myNodalIDs[0]=Node_1;
	myNodalIDs[1]=Node_2;
	myNodalIDs[2]=Node_3;
	myNodalIDs[3]=Node_4;
	myNodalIDs[4]=Node_5;
	myNodalIDs[5]=Node_6;
	myNodalIDs[6]=Node_7;
	myNodalIDs[7]=Node_8;
	// Set local nodal dofs
	myLocalNodalDofs.resize(3);
	myLocalNodalDofs[0]=0;
	myLocalNodalDofs[1]=1;	
	myLocalNodalDofs[2]=2;	
	// Handle common info
	this->handleCommonInfo();

	// Materials
	myMatPoints.resize(8);
	MultiaxialMaterial* pMat=static_cast<MultiaxialMaterial*>(myMaterial); 
	myMatPoints[0]=new MatPoint(pMat,1,1,1,2,2,2);
	myMatPoints[1]=new MatPoint(pMat,2,1,1,2,2,2);
	myMatPoints[2]=new MatPoint(pMat,2,2,1,2,2,2);
	myMatPoints[3]=new MatPoint(pMat,1,2,1,2,2,2);
	myMatPoints[4]=new MatPoint(pMat,1,1,2,2,2,2);
	myMatPoints[5]=new MatPoint(pMat,2,1,2,2,2,2);
	myMatPoints[6]=new MatPoint(pMat,2,2,2,2,2,2);
	myMatPoints[7]=new MatPoint(pMat,1,2,2,2,2,2);
	
	shape8(x,shp,detJ);
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		double xG=0,yG=0,zG=0;
		for(unsigned a=0;a<myNodes.size();a++)
		{
			xG+=shp[a][0][k]*x(a,0);
			yG+=shp[a][0][k]*x(a,1);
			zG+=shp[a][0][k]*x(a,2);
		}
		myMatPoints[k]->setX(xG,yG,zG);
	}
}
Brick8d::~Brick8d()
{
	Containers::vector_delete(myMatPoints);
}
const Matrix& Brick8d::getK()
{
	Matrix &K=*myMatrix;
	K.clear();
	shape8(x,shp,detJ);
	static Matrix Ba(6,3);
	static Matrix Bb(6,3);
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		const Matrix& C=myMatPoints[k]->getMaterial()->getC();
		for(unsigned a=0;a<myNodes.size();a++)
		{
			this->getB(Ba,a,k);
			for(unsigned b=0;b<myNodes.size();b++)
			{
				this->getB(Bb,b,k);
				add(K,3*a,3*b,Ba,C,Bb,detJ[k],1.0);
			}
		}
	}
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	return K;
}
const Matrix& Brick8d::getM()
{
	Matrix &M=*myMatrix;
	M.clear();
	return M;
}
const Vector& Brick8d::getR()
{
	static Vector sigma(6);
	static Matrix B(6,3);
	Vector& R=*myVector;
	R.clear();
	if(!(myGroup->isActive()))	return R;
	double facS=myGroup->getFacS();
	double facG=myGroup->getFacG();
	double facP=myGroup->getFacP();
	
	shape8(x,shp,detJ);
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		sigma=myMatPoints[k]->getMaterial()->getStress();
		for(unsigned a=0;a<myNodes.size();a++)
		{
			this->getB(B,a,k);
			add(R,3*a,B,sigma,detJ[k],1.0);
		}
	}
	R-=facP*P;
	return R;
}
void Brick8d::update()
{
	if(!(myGroup->isActive()))	return;
	static Vector u(24);
	u=this->getDispIncrm();
	shape8(x,shp,detJ);

	// For each material point
	static Vector epsilon(6);
	static Matrix B(6,3);
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		epsilon.clear();
		for(unsigned a=0;a<myNodes.size();a++)
		{
			this->getB(B,a,k);
			add2(epsilon,3*a,B,u,1.0,1.0);
		}
		myMatPoints[k]->getMaterial()->setStrain(epsilon);
	}
}
void Brick8d::getB(Matrix& B,int node,int gPoint)
{
	double B1=shp[node][1][gPoint];
	double B2=shp[node][2][gPoint];
	double B3=shp[node][3][gPoint];

	B(0,0)=B1;	B(0,1)=0.;	B(0,2)=0.;
	B(1,0)=0.;	B(1,1)=B2;	B(1,2)=0.;
	B(2,0)=0.;	B(2,1)=0.;	B(2,2)=B3;
	B(3,0)=B2;	B(3,1)=B1;	B(3,2)=0.;
	B(4,0)=0.;	B(4,1)=B3;	B(4,2)=B2;
	B(5,0)=B3;	B(5,1)=0.;	B(5,2)=B1;
}
void Brick8d::commit()
{
	for(unsigned int i=0;i<myMatPoints.size();i++) 
		myMatPoints[i]->getMaterial()->commit();
}
bool Brick8d::checkIfAllows(FEObject* f)
{
	return true;
}
void Brick8d::addInitialStresses(InitialStresses* pInitialStresses)
{
	if(myGroup->isActive()&&pInitialStresses->getGroupID()==myGroup->getID())
		for(unsigned i=0;i<myMatPoints.size();i++)
			myMatPoints[i]->setInitialStresses(pInitialStresses);
}
void Brick8d::recoverStresses()
{
 	static Vector sigma(6);
 	static Matrix E(8,8);
	const double d=0.125;
	const double a=1+num::sq3;
	const double b=1-num::sq3;

	E(0,0)=a*a*a*d; E(0,1)=b*a*a*d; E(0,2)=b*b*a*d; E(0,3)=b*a*a*d;
	E(1,0)=b*a*a*d; E(1,1)=a*a*a*d; E(1,2)=b*a*a*d; E(1,3)=b*b*a*d;
	E(2,0)=b*b*a*d; E(2,1)=b*a*a*d; E(2,2)=a*a*a*d; E(2,3)=b*a*a*d;
	E(3,0)=b*a*a*d; E(3,1)=b*b*a*d; E(3,2)=b*a*a*d; E(3,3)=a*a*a*d;
	E(4,0)=b*a*a*d; E(4,1)=b*b*a*d; E(4,2)=b*b*b*d; E(4,3)=b*b*a*d;
	E(5,0)=b*b*a*d; E(5,1)=b*a*a*d; E(5,2)=b*b*a*d; E(5,3)=b*b*b*d;
	E(6,0)=b*b*b*d; E(6,1)=b*b*a*d; E(6,2)=b*a*a*d; E(6,3)=b*b*a*d;
	E(7,0)=b*b*a*d; E(7,1)=b*b*b*d; E(7,2)=b*b*a*d; E(7,3)=b*a*a*d;
	E(0,4)=b*a*a*d; E(0,5)=b*b*a*d; E(0,6)=b*b*b*d; E(0,7)=b*b*a*d;
	E(1,4)=b*b*a*d; E(1,5)=b*a*a*d; E(1,6)=b*b*a*d; E(1,7)=b*b*b*d;
	E(2,4)=b*b*b*d; E(2,5)=b*b*a*d; E(2,6)=b*a*a*d; E(2,7)=b*b*a*d;
	E(3,4)=b*b*a*d; E(3,5)=b*b*b*d; E(3,6)=b*b*a*d; E(3,7)=b*a*a*d;
	E(4,4)=a*a*a*d; E(4,5)=b*a*a*d; E(4,6)=b*b*a*d; E(4,7)=b*a*a*d;
	E(5,4)=b*a*a*d; E(5,5)=a*a*a*d; E(5,6)=b*a*a*d; E(5,7)=b*b*a*d;
	E(6,4)=b*b*a*d; E(6,5)=b*a*a*d; E(6,6)=a*a*a*d; E(6,7)=b*a*a*d;
	E(7,4)=b*a*a*d; E(7,5)=b*b*a*d; E(7,6)=b*a*a*d; E(7,7)=a*a*a*d;

	for(unsigned i=0;i<8;i++)			// nodes
	{
		sigma.clear();
		for(unsigned j=0;j<6;j++)		// sigma
			for(unsigned k=0;k<8;k++)	// material points
				sigma[j]+=E(i,k)*(myMatPoints[k]->getMaterial()->getStress())[j];
		myNodes[i]->addStress(sigma);
	}
}
const int Brick8d::getnPlasticPoints()
{
	int n=0;
	for(unsigned i=0;i<myMatPoints.size();i++)
		if(myMatPoints[i]->getMaterial()->isPlastic()) n+=1;
	return n;
}
