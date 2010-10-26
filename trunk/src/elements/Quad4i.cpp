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

#include "elements/Quad4i.h"
#include "elements/ShapeFunctions.h"
#include "main/NemesisDebug.h"

double Quad4i::detJ[4];
double Quad4i::shpStd[4][3][4];
double Quad4i::shpInc[2][3][4];
std::vector<int> Quad4i::perm(3);

Quad4i::Quad4i()
{
}
Quad4i::Quad4i(int ID,int Node_1,int Node_2,int Node_3,int Node_4,int MatID)
:Quad4(ID,Node_1,Node_2,Node_3,Node_4,MatID,2,2)
{
	perm[0]=0;
	perm[1]=1;
	perm[2]=3;
	aTrial.resize(4,0.);
	aConvg.resize(4,0.);
}
Quad4i::~Quad4i()
{
}
const Matrix& Quad4i::getK()
{
	// Get a reference to myMatrix as K
	Matrix &K=*myMatrix;
	// Define local static matrices
	static Matrix Ba(3,2),Bb(3,2);
	static Matrix Kdd(8,8),Kda(8,4),Kaa(4,4);
	// Form shape functions
	this->shapeFunctions();
	// Get Kdd,Kda,Kaa Matrices
	this->getKdd(Kdd);
	this->getKda(Kda);
	this->getKaa(Kaa);
	// Form K
	K=Kdd-Kda*Inverse(Kaa)*Transpose(Kda);
	// Get group factors
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	// Return K
	return K;
}
const Matrix& Quad4i::getM()
{
	// Get a reference to myMatrix as K
	Matrix &M=*myMatrix;
	M.clear();
	// Find total mass
	double rho=myMaterial->getRho();
	double volume=0.;
	this->shapeFunctions();
	for(unsigned k=0;k<myMatPoints.size();k++)
		volume+=detJ[k]*(pD->getFac())*(myMatPoints[k]->get_w()); 
	double mass=rho*volume;
	// Set corresponding mass to diagonal terms
	for(int i=0;i<8;i++) 
		M(i,i)=0.25*mass;
	// Return M
	return M;
}
/**
 * Element residual vector.
 */
const Vector& Quad4i::getR()
{
	// Get a reference to myVector as R
	Vector& R=*myVector;
	R.clear();
	// Static vectors and matrices
	static Vector sigma(6);
	static Matrix Ba(3,2);
	// Factors
	if(!(myGroup->isActive()))	return R;
	double facS=myGroup->getFacS();
	double facG=myGroup->getFacG();
	double facP=myGroup->getFacP();
	
	// Find shape functions for all GaussPoints
	this->shapeFunctions();

	// R = facS*Fint - facG*SelfWeigth - facP*ElementalLoads
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		sigma=myMatPoints[k]->getMaterial()->getStress();
		double dV=(pD->getFac())*detJ[k];
		for(unsigned a=0;a<myNodes.size();a++)
		{
			// +facS*Fint
			this->getBStd(Ba,a,k);
			add_BTv(R,2*a,&perm[0],Ba,sigma,facS*dV,1.0);
			// -facG*SelfWeigth
			for(int i=0;i<2;i++)
				R[2*a+i]-=facG*shpStd[a][0][k]*b[i]*dV;
		}
	}
	// -facP*ElementalLoads
	R-=facP*P;

	// Return R
	return R;
}
/**
 * Element update.
 */
void Quad4i::update()
{
	// Check for a quick return
	if(!(myGroup->isActive()))	return;
	// Static vectors and matrices
	static Vector Du(8),Da(4),epsilon(6);
	static Matrix Ba(3,2);
	static Matrix Kda(8,4),Kaa(4,4);
	// Form shape functions
	this->shapeFunctions();
	// Get incremental displacements
	Du=this->getDispIncrm();
	// Get incremental alphas
	this->getKda(Kda);
	this->getKaa(Kaa);
	Da=-Inverse(Kaa)*Transpose(Kda)*Du;
	aTrial=aConvg+Da;
	// For each material point
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		epsilon.clear();
		for(unsigned a=0;a<4;a++)
		{
			this->getBStd(Ba,a,k);
			double dV=(pD->getFac())*detJ[k];
			add_Bv(epsilon,2*a,&perm[0],Ba,Du,1.0,1.0);
		}
		for(unsigned a=0;a<2;a++)
		{
			this->getBInc(Ba,a,k);
			double dV=(pD->getFac())*detJ[k];
			add_Bv(epsilon,2*a,&perm[0],Ba,Da,1.0,1.0);
		}
		myMatPoints[k]->getMaterial()->setStrain(epsilon);
	}
}
/**
 * Element commit.
 * Overwrites base function, because incompatible modes must also
 * be commited. This takes place in element level.
 */
void Quad4i::commit()
{
	for(unsigned int i=0;i<myMatPoints.size();i++) 
		myMatPoints[i]->getMaterial()->commit();
	aConvg=aTrial;
}
/**
 * Call shape functions and fill in corresponding arrays.
 */
void Quad4i::shapeFunctions()
{
	shape4(x,shpStd,detJ);
	shapeQM6(x,shpInc);
}
/**
 * Get standard displacement B-matrix.
 * Shape function array shpStd and detJs must be defined.
 * @param B B-Matrix.
 * @param node The corresponding node.
 * @param gPoint The corresponding Gauss Point.
 */
void Quad4i::getBStd(Matrix& B,int node,int gPoint)
{
	// B-factors
	double B1=shpStd[node][1][gPoint];
	double B2=shpStd[node][2][gPoint];

	// B-matrix
	B(0,0)=B1;	B(0,1)=0.;
	B(1,0)=0.;	B(1,1)=B2;
	B(2,0)=B2;	B(2,1)=B1;
}
/**
 * Get B-matrix of incompatible modes.
 * Shape function array shpInc and detJs must be defined.
 * @param B B-Matrix.
 * @param node The corresponding node.
 * @param gPoint The corresponding Gauss Point.
 */
void Quad4i::getBInc(Matrix& B,int node,int gPoint)
{
	// B-factors
	double B1=shpInc[node][1][gPoint]/detJ[gPoint];
	double B2=shpInc[node][2][gPoint]/detJ[gPoint];

	// B-matrix
	B(0,0)=B1;	B(0,1)=0.;
	B(1,0)=0.;	B(1,1)=B2;
	B(2,0)=B2;	B(2,1)=B1;
}
/**
 * Get Kdd Matrix.
 * This is the usual displacement matrix.
 * Formed as the integral of Bstd^T.C.Bstd on dOmega.
 * @todo Increase performance.
 * @param K The matrix to be filled.
 */
void Quad4i::getKdd(Matrix& K)
{
	K.clear();
	static Matrix Ba(3,2),Bb(3,2);
	// For all Gauss points
	for(unsigned k=0;k<4;k++)
	{
		const Matrix& C=myMatPoints[k]->getMaterial()->getC();
		for(unsigned a=0;a<4;a++)
		{
			getBStd(Ba,a,k);
			for(unsigned b=0;b<4;b++)
			{
				getBStd(Bb,b,k);
				double dV=(pD->getFac())*detJ[k];
				K.add_BTCB(2*a,2*b,&perm[0],Ba,C,Bb,dV,1.0);
			}
		}
	}

}
/**
 * Get Kda Matrix.
 * Formed as the integral of Bstd^T.C.Binc on dOmega.
 * @todo Increase performance.
 * @param K The matrix to be filled.
 */
void Quad4i::getKda(Matrix& K)
{
	K.clear();
	static Matrix Ba(3,2),Bb(3,2);
	// For all Gauss points
	for(unsigned k=0;k<4;k++)
	{
		const Matrix& C=myMatPoints[k]->getMaterial()->getC();
		for(unsigned a=0;a<4;a++)
		{
			getBStd(Ba,a,k);
			for(unsigned b=0;b<2;b++)
			{
				getBInc(Bb,b,k);
				double dV=1.0*detJ[k];
				K.add_BTCB(2*a,2*b,&perm[0],Ba,C,Bb,dV,1.0);
			}
		}
	}
}
/**
 * Get Kda Matrix.
 * Formed as the integral of Binc^T.C.Binc on dOmega.
 * @todo Increase performance.
 * @param K The matrix to be filled.
 */
void Quad4i::getKaa(Matrix& K)
{
	K.clear();
	static Matrix Ba(3,2),Bb(3,2);
	// For all Gauss points
	for(unsigned k=0;k<4;k++)
	{
		const Matrix& C=myMatPoints[k]->getMaterial()->getC();
		for(unsigned a=0;a<2;a++)
		{
			getBInc(Ba,a,k);
			for(unsigned b=0;b<2;b++)
			{
				getBInc(Bb,b,k);
				double dV=1.0*detJ[k];
				K.add_BTCB(2*a,2*b,&perm[0],Ba,C,Bb,dV,1.0);
			}
		}
	}
}
