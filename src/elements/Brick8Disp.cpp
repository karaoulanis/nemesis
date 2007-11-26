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

#include <Brick8Disp.h>
#include <NemesisDebug.h>

Matrix Brick8Disp::N(8,4);
double Brick8Disp::detJ;

Brick8Disp::Brick8Disp()
{
}
Brick8Disp::Brick8Disp(int ID,
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
	for(unsigned i=0;i<8;i++)
	{
		this->findShapeFunctionsAt(myMatPoints[i]);
		double xG=0,yG=0,zG=0;
		for(unsigned j=0;j<8;j++)
		{
			xG+=N(j,0)*x(j,0);
			yG+=N(j,0)*x(j,1);
			zG+=N(j,0)*x(j,2);
		}
		myMatPoints[i]->setX(xG,yG,zG);
	}
}
Brick8Disp::~Brick8Disp()
{
	Containers::vector_delete(myMatPoints);
}
int Brick8Disp::findShapeFunctionsAt(MatPoint* pMatPoint)
{
	double xi=  pMatPoint->get_r();
	double eta= pMatPoint->get_s();
	double zeta=pMatPoint->get_t();

	N(0,0)=0.125*(1-xi)*(1-eta)*(1-zeta);			// N1
	N(1,0)=0.125*(1+xi)*(1-eta)*(1-zeta);			// N2
	N(2,0)=0.125*(1+xi)*(1+eta)*(1-zeta);			// N3
	N(3,0)=0.125*(1-xi)*(1+eta)*(1-zeta);			// N4
	N(4,0)=0.125*(1-xi)*(1-eta)*(1+zeta);			// N5
	N(5,0)=0.125*(1+xi)*(1-eta)*(1+zeta);			// N6
	N(6,0)=0.125*(1+xi)*(1+eta)*(1+zeta);			// N7
	N(7,0)=0.125*(1-xi)*(1+eta)*(1+zeta);			// N8

	static Matrix J(3,3);
	J.clear();
	J(0,0)=0.125*(-x(7,0)+x(6,0)*eta*zeta-x(3,0)*eta-x(5,0)*eta*zeta-x(7,0)*eta*zeta+x(3,0)*eta*zeta+x(2,0)*eta+x(4,0)*eta*zeta-x(1,0)*eta-x(0,0)+x(1,0)-x(7,0)*eta+x(5,0)*zeta-x(2,0)*zeta-x(1,0)*zeta+x(2,0)-x(3,0)+x(3,0)*zeta-x(4,0)*zeta+x(6,0)*eta+x(6,0)*zeta+x(1,0)*eta*zeta-x(5,0)*eta-x(2,0)*eta*zeta+x(4,0)*eta+x(0,0)*eta+x(0,0)*zeta-x(4,0)+x(5,0)+x(6,0)-x(7,0)*zeta-x(0,0)*eta*zeta);
	J(0,1)=0.125*(+x(1,1)*eta*zeta+x(2,1)*eta-x(7,1)*eta*zeta+x(6,1)*eta*zeta-x(0,1)+x(1,1)+x(3,1)*zeta+x(3,1)*eta*zeta-x(3,1)*eta-x(4,1)*zeta+x(4,1)*eta+x(5,1)*zeta-x(5,1)*eta-x(2,1)*eta*zeta+x(6,1)*eta+x(6,1)*zeta-x(7,1)*eta-x(7,1)*zeta+x(0,1)*eta-x(1,1)*eta+x(0,1)*zeta-x(1,1)*zeta+x(5,1)+x(4,1)*eta*zeta-x(5,1)*eta*zeta-x(0,1)*eta*zeta-x(2,1)*zeta+x(6,1)-x(7,1)+x(2,1)-x(3,1)-x(4,1));
	J(0,2)=0.125*(-x(1,2)*zeta+x(0,2)*zeta-x(5,2)*eta*zeta-x(0,2)+x(1,2)+x(2,2)-x(3,2)-x(4,2)+x(6,2)-x(7,2)+x(3,2)*eta*zeta+x(6,2)*eta*zeta+x(1,2)*eta*zeta-x(7,2)*eta*zeta-x(0,2)*eta*zeta-x(1,2)*eta+x(4,2)*eta*zeta+x(2,2)*eta-x(2,2)*zeta+x(3,2)*zeta-x(3,2)*eta+x(4,2)*eta-x(4,2)*zeta+x(5,2)*zeta-x(5,2)*eta+x(6,2)*zeta+x(6,2)*eta-x(7,2)*eta-x(7,2)*zeta-x(2,2)*eta*zeta+x(5,2)+x(0,2)*eta);
	J(1,0)=0.125*(+x(7,0)-x(5,0)*xi*zeta-x(7,0)*xi*zeta+x(3,0)*xi*zeta-x(1,0)*xi+x(4,0)*xi*zeta-x(2,0)*xi*zeta-x(7,0)*xi-x(0,0)-x(1,0)-x(0,0)*xi*zeta-x(5,0)*zeta-x(2,0)*zeta+x(2,0)*xi+x(1,0)*zeta+x(2,0)+x(3,0)+x(6,0)*xi*zeta-x(3,0)*zeta-x(4,0)*zeta+x(6,0)*xi+x(6,0)*zeta+x(0,0)*xi-x(5,0)*xi+x(4,0)*xi+x(0,0)*zeta+x(1,0)*xi*zeta-x(3,0)*xi-x(4,0)-x(5,0)+x(6,0)+x(7,0)*zeta);
	J(1,1)=0.125*(-x(5,1)*xi*zeta+x(3,1)*xi*zeta-x(7,1)*xi*zeta+x(0,1)*xi+x(6,1)*xi*zeta-x(0,1)-x(1,1)-x(1,1)*xi-x(3,1)*zeta-x(3,1)*xi-x(4,1)*zeta+x(4,1)*xi-x(5,1)*zeta-x(5,1)*xi+x(6,1)*xi+x(6,1)*zeta+x(7,1)*zeta+x(4,1)*xi*zeta-x(7,1)*xi+x(0,1)*zeta-x(2,1)*xi*zeta-x(0,1)*xi*zeta+x(1,1)*zeta+x(2,1)*xi-x(5,1)+x(1,1)*xi*zeta-x(2,1)*zeta+x(6,1)+x(7,1)+x(2,1)+x(3,1)-x(4,1));
	J(1,2)=0.125*(+x(3,2)*xi*zeta+x(1,2)*zeta+x(0,2)*zeta-x(0,2)-x(1,2)+x(2,2)+x(3,2)-x(4,2)+x(6,2)+x(7,2)+x(6,2)*xi*zeta+x(4,2)*xi*zeta-x(7,2)*xi*zeta-x(1,2)*xi+x(2,2)*xi-x(2,2)*zeta-x(3,2)*zeta-x(4,2)*zeta-x(3,2)*xi+x(4,2)*xi-x(5,2)*zeta-x(5,2)*xi+x(6,2)*zeta+x(6,2)*xi-x(7,2)*xi+x(7,2)*zeta-x(2,2)*xi*zeta-x(5,2)*xi*zeta-x(5,2)+x(0,2)*xi-x(0,2)*xi*zeta+x(1,2)*xi*zeta);
	J(2,0)=0.125*(+x(7,0)+x(6,0)*xi*eta-x(5,0)*xi*eta-x(3,0)*eta-x(7,0)*xi*eta+x(1,0)*xi*eta-x(2,0)*eta-x(1,0)*xi-x(0,0)*xi*eta+x(4,0)*xi*eta-x(7,0)*xi+x(1,0)*eta-x(0,0)-x(1,0)+x(7,0)*eta-x(2,0)*xi*eta+x(3,0)*xi*eta-x(2,0)*xi-x(2,0)-x(3,0)+x(6,0)*xi+x(6,0)*eta+x(0,0)*xi+x(5,0)*xi-x(5,0)*eta-x(4,0)*xi-x(4,0)*eta+x(0,0)*eta+x(3,0)*xi+x(4,0)+x(5,0)+x(6,0));
	J(2,1)=0.125*(-x(2,1)*eta+x(4,1)*xi*eta+x(3,1)*xi*eta+x(0,1)*xi+x(6,1)*xi*eta-x(5,1)*xi*eta-x(7,1)*xi*eta-x(0,1)-x(1,1)-x(1,1)*xi+x(3,1)*xi-x(3,1)*eta-x(4,1)*xi-x(4,1)*eta-x(5,1)*eta+x(5,1)*xi-x(2,1)*xi*eta+x(6,1)*xi+x(6,1)*eta+x(7,1)*eta-x(7,1)*xi+x(0,1)*eta+x(1,1)*eta-x(2,1)*xi+x(5,1)-x(0,1)*xi*eta+x(1,1)*xi*eta+x(6,1)+x(7,1)-x(2,1)-x(3,1)+x(4,1));
	J(2,2)=0.125*(-x(5,2)*xi*eta-x(0,2)-x(1,2)-x(2,2)-x(3,2)+x(4,2)+x(6,2)+x(7,2)+x(6,2)*xi*eta+x(4,2)*xi*eta-x(7,2)*xi*eta+x(1,2)*eta-x(2,2)*xi*eta-x(1,2)*xi-x(2,2)*xi-x(2,2)*eta-x(3,2)*eta-x(4,2)*eta+x(3,2)*xi-x(4,2)*xi+x(5,2)*xi-x(5,2)*eta+x(3,2)*xi*eta+x(6,2)*xi+x(6,2)*eta-x(7,2)*xi+x(7,2)*eta+x(1,2)*xi*eta-x(0,2)*xi*eta+x(5,2)+x(0,2)*xi+x(0,2)*eta);

	detJ=det(J);
	if(detJ<0.) cout<<"Determinant negative!"<<endl;
	static Matrix invJ(3,3);
	invJ=Inverse(J);
	double dxidx   = invJ(0,0);	  double detadx  = invJ(0,1);   double dzetadx = invJ(0,2);
	double dxidy   = invJ(1,0);	  double detady  = invJ(1,1);   double dzetady = invJ(1,2);
	double dxidz   = invJ(2,0);	  double detadz  = invJ(2,1);   double dzetadz = invJ(2,2);

	N(0,1)=0.125*(-(1-eta)*(1-zeta)*dxidx-(1-xi)*(1-zeta)*detadx-(1-xi)*(1-eta)*dzetadx);
	N(0,2)=0.125*(-(1-eta)*(1-zeta)*dxidy-(1-xi)*(1-zeta)*detady-(1-xi)*(1-eta)*dzetady);
	N(0,3)=0.125*(-(1-eta)*(1-zeta)*dxidz-(1-xi)*(1-zeta)*detadz-(1-xi)*(1-eta)*dzetadz);
	N(1,1)=0.125*(+(1-eta)*(1-zeta)*dxidx-(1+xi)*(1-zeta)*detadx-(1+xi)*(1-eta)*dzetadx);
	N(1,2)=0.125*(+(1-eta)*(1-zeta)*dxidy-(1+xi)*(1-zeta)*detady-(1+xi)*(1-eta)*dzetady);
	N(1,3)=0.125*(+(1-eta)*(1-zeta)*dxidz-(1+xi)*(1-zeta)*detadz-(1+xi)*(1-eta)*dzetadz);
	N(2,1)=0.125*(+(1+eta)*(1-zeta)*dxidx+(1+xi)*(1-zeta)*detadx-(1+xi)*(1+eta)*dzetadx);
	N(2,2)=0.125*(+(1+eta)*(1-zeta)*dxidy+(1+xi)*(1-zeta)*detady-(1+xi)*(1+eta)*dzetady);
	N(2,3)=0.125*(+(1+eta)*(1-zeta)*dxidz+(1+xi)*(1-zeta)*detadz-(1+xi)*(1+eta)*dzetadz);
	N(3,1)=0.125*(-(1+eta)*(1-zeta)*dxidx+(1-xi)*(1-zeta)*detadx-(1-xi)*(1+eta)*dzetadx);
	N(3,2)=0.125*(-(1+eta)*(1-zeta)*dxidy+(1-xi)*(1-zeta)*detady-(1-xi)*(1+eta)*dzetady);
	N(3,3)=0.125*(-(1+eta)*(1-zeta)*dxidz+(1-xi)*(1-zeta)*detadz-(1-xi)*(1+eta)*dzetadz);
	N(4,1)=0.125*(-(1-eta)*(1+zeta)*dxidx-(1-xi)*(1+zeta)*detadx+(1-xi)*(1-eta)*dzetadx);
	N(4,2)=0.125*(-(1-eta)*(1+zeta)*dxidy-(1-xi)*(1+zeta)*detady+(1-xi)*(1-eta)*dzetady);
	N(4,3)=0.125*(-(1-eta)*(1+zeta)*dxidz-(1-xi)*(1+zeta)*detadz+(1-xi)*(1-eta)*dzetadz);
	N(5,1)=0.125*(+(1-eta)*(1+zeta)*dxidx-(1+xi)*(1+zeta)*detadx+(1+xi)*(1-eta)*dzetadx);
	N(5,2)=0.125*(+(1-eta)*(1+zeta)*dxidy-(1+xi)*(1+zeta)*detady+(1+xi)*(1-eta)*dzetady);
	N(5,3)=0.125*(+(1-eta)*(1+zeta)*dxidz-(1+xi)*(1+zeta)*detadz+(1+xi)*(1-eta)*dzetadz);
	N(6,1)=0.125*(+(1+eta)*(1+zeta)*dxidx+(1+xi)*(1+zeta)*detadx+(1+xi)*(1+eta)*dzetadx);
	N(6,2)=0.125*(+(1+eta)*(1+zeta)*dxidy+(1+xi)*(1+zeta)*detady+(1+xi)*(1+eta)*dzetady);
	N(6,3)=0.125*(+(1+eta)*(1+zeta)*dxidz+(1+xi)*(1+zeta)*detadz+(1+xi)*(1+eta)*dzetadz);
	N(7,1)=0.125*(-(1+eta)*(1+zeta)*dxidx+(1-xi)*(1+zeta)*detadx+(1-xi)*(1+eta)*dzetadx);
	N(7,2)=0.125*(-(1+eta)*(1+zeta)*dxidy+(1-xi)*(1+zeta)*detady+(1-xi)*(1+eta)*dzetady);
	N(7,3)=0.125*(-(1+eta)*(1+zeta)*dxidz+(1-xi)*(1+zeta)*detadz+(1-xi)*(1+eta)*dzetadz);

	return 0;
}
const Matrix& Brick8Disp::getK()
{
	Matrix &K=*myMatrix;
	K.clear();
	for(unsigned int k=0;k<myMatPoints.size();k++)
	{
		this->findShapeFunctionsAt(myMatPoints[k]);
		double dV=detJ*(pD->getFac())*(myMatPoints[k]->get_w());
		const Matrix& C=myMatPoints[k]->getMaterial()->getC();
		int ii=0;
		for(int i=0;i<8;i++) 
		{
			int jj=0;
			for(int j=0;j<8;j++) 
			{
				static Matrix CB(6,3);
				CB(0,0)=C(0,0)*N(j,1)+C(0,3)*N(j,2)+C(0,5)*N(j,3);
				CB(1,0)=C(1,0)*N(j,1)+C(1,3)*N(j,2)+C(1,5)*N(j,3);
				CB(2,0)=C(2,0)*N(j,1)+C(2,3)*N(j,2)+C(2,5)*N(j,3);
				CB(3,0)=C(3,0)*N(j,1)+C(3,3)*N(j,2)+C(3,5)*N(j,3);
				CB(4,0)=C(4,0)*N(j,1)+C(4,3)*N(j,2)+C(4,5)*N(j,3);
				CB(5,0)=C(5,0)*N(j,1)+C(5,3)*N(j,2)+C(5,5)*N(j,3);
				CB(0,1)=C(0,1)*N(j,2)+C(0,3)*N(j,1)+C(0,4)*N(j,3);
				CB(1,1)=C(1,1)*N(j,2)+C(1,3)*N(j,1)+C(1,4)*N(j,3);
				CB(2,1)=C(2,1)*N(j,2)+C(2,3)*N(j,1)+C(2,4)*N(j,3);
				CB(3,1)=C(3,1)*N(j,2)+C(3,3)*N(j,1)+C(3,4)*N(j,3);
				CB(4,1)=C(4,1)*N(j,2)+C(4,3)*N(j,1)+C(4,4)*N(j,3);
				CB(5,1)=C(5,1)*N(j,2)+C(5,3)*N(j,1)+C(5,4)*N(j,3);
				CB(0,2)=C(0,2)*N(j,3)+C(0,4)*N(j,2)+C(0,5)*N(j,1);
				CB(1,2)=C(1,2)*N(j,3)+C(1,4)*N(j,2)+C(1,5)*N(j,1);
				CB(2,2)=C(2,2)*N(j,3)+C(2,4)*N(j,2)+C(2,5)*N(j,1);
				CB(3,2)=C(3,2)*N(j,3)+C(3,4)*N(j,2)+C(3,5)*N(j,1);
				CB(4,2)=C(4,2)*N(j,3)+C(4,4)*N(j,2)+C(4,5)*N(j,1);
				CB(5,2)=C(5,2)*N(j,3)+C(5,4)*N(j,2)+C(5,5)*N(j,1);
				K(ii  ,jj  )+=(N(i,1)*CB(0,0)+N(i,2)*CB(3,0)+N(i,3)*CB(5,0))*dV;
				K(ii  ,jj+1)+=(N(i,1)*CB(0,1)+N(i,2)*CB(3,1)+N(i,3)*CB(5,1))*dV;
				K(ii  ,jj+2)+=(N(i,1)*CB(0,2)+N(i,2)*CB(3,2)+N(i,3)*CB(5,2))*dV;
				K(ii+1,jj  )+=(N(i,2)*CB(1,0)+N(i,1)*CB(3,0)+N(i,3)*CB(4,0))*dV;
				K(ii+1,jj+1)+=(N(i,2)*CB(1,1)+N(i,1)*CB(3,1)+N(i,3)*CB(4,1))*dV;
				K(ii+1,jj+2)+=(N(i,2)*CB(1,2)+N(i,1)*CB(3,2)+N(i,3)*CB(4,2))*dV;
				K(ii+2,jj  )+=(N(i,3)*CB(2,0)+N(i,2)*CB(4,0)+N(i,1)*CB(5,0))*dV;
				K(ii+2,jj+1)+=(N(i,3)*CB(2,1)+N(i,2)*CB(4,1)+N(i,1)*CB(5,1))*dV;
				K(ii+2,jj+2)+=(N(i,3)*CB(2,2)+N(i,2)*CB(4,2)+N(i,1)*CB(5,2))*dV;
				jj+=3;
			}
			ii+=3;
		}
	}
	//cout<<K;
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	return K;
}
const Matrix& Brick8Disp::getM()
{
	Matrix &M=*myMatrix;
	M.clear();
	return M;
}
const Vector& Brick8Disp::getR()
{
	static Vector sigma(6);
	Vector& R=*myVector;
	R.clear();
	if(!(myGroup->isActive()))	return R;
	double facS=myGroup->getFacS();
	double facG=myGroup->getFacG();
	double facP=myGroup->getFacP();
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		sigma=myMatPoints[k]->getMaterial()->getStress();
		this->findShapeFunctionsAt(myMatPoints[k]);
		int ii=0;
		double dV=detJ*(pD->getFac())*(myMatPoints[k]->get_w());
		for(int i=0;i<8;i++)
		{
			R[ii  ]+=facS*(N(i,1)*sigma[0]+N(i,2)*sigma[3]+N(i,3)*sigma[5])*dV;
			R[ii+1]+=facS*(N(i,2)*sigma[1]+N(i,1)*sigma[3]+N(i,3)*sigma[4])*dV;
			R[ii+2]+=facS*(N(i,3)*sigma[2]+N(i,2)*sigma[4]+N(i,1)*sigma[5])*dV;
			///@todo check
			R[ii  ]-=facG*(N(i,0)*b[0]*dV);
			R[ii+1]-=facG*(N(i,0)*b[1]*dV);
			R[ii+2]-=facG*(N(i,0)*b[2]*dV);
			ii+=3;
		}
	}
	R-=facP*P;
	return R;
}
void Brick8Disp::update()
{
	if(!(myGroup->isActive()))	return;
	static Vector u(24);
	u=this->getDispIncrm();
	// For each material point
	for(unsigned int i=0;i<myMatPoints.size();i++)
	{
		this->findShapeFunctionsAt(myMatPoints[i]);
		// Determine the strain
		static Vector epsilon(6);
  		epsilon.clear();
		epsilon[0]=N(0,1)*u[ 0]+N(1,1)*u[ 3]+N(2,1)*u[ 6]+N(3,1)*u[ 9]
				  +N(4,1)*u[12]+N(5,1)*u[15]+N(6,1)*u[18]+N(7,1)*u[21];
		epsilon[1]=N(0,2)*u[ 1]+N(1,2)*u[ 4]+N(2,2)*u[ 7]+N(3,2)*u[10]
				  +N(4,2)*u[13]+N(5,2)*u[16]+N(6,2)*u[19]+N(7,2)*u[22];
		epsilon[2]=N(0,3)*u[ 2]+N(1,3)*u[ 5]+N(2,3)*u[ 8]+N(3,3)*u[11]
				  +N(4,3)*u[14]+N(5,3)*u[17]+N(6,3)*u[20]+N(7,3)*u[23];
		epsilon[3]=N(0,2)*u[ 0]+N(0,1)*u[ 1]+N(1,2)*u[ 3]+N(1,1)*u[ 4]
				  +N(2,2)*u[ 6]+N(2,1)*u[ 7]+N(3,2)*u[ 9]+N(3,1)*u[10]
				  +N(4,2)*u[12]+N(4,1)*u[13]+N(5,2)*u[15]+N(5,1)*u[16]
				  +N(6,2)*u[18]+N(6,1)*u[19]+N(7,2)*u[21]+N(7,1)*u[22];
		epsilon[4]=N(0,3)*u[ 1]+N(0,2)*u[ 2]+N(1,3)*u[ 4]+N(1,2)*u[ 5]
				  +N(2,3)*u[ 7]+N(2,2)*u[ 8]+N(3,3)*u[10]+N(3,2)*u[11]
				  +N(4,3)*u[13]+N(4,2)*u[14]+N(5,3)*u[16]+N(5,2)*u[17]
				  +N(6,3)*u[19]+N(6,2)*u[20]+N(7,3)*u[22]+N(7,2)*u[23];
		epsilon[5]=N(0,3)*u[ 0]+N(0,1)*u[ 2]+N(1,3)*u[ 3]+N(1,1)*u[ 5]
                  +N(2,3)*u[ 6]+N(2,1)*u[ 8]+N(3,3)*u[ 9]+N(3,1)*u[11]
				  +N(4,3)*u[12]+N(4,1)*u[14]+N(5,3)*u[15]+N(5,1)*u[17]
				  +N(6,3)*u[18]+N(6,1)*u[20]+N(7,3)*u[21]+N(7,1)*u[23];
		// And send it to the material point
		myMatPoints[i]->getMaterial()->setStrain(epsilon);
	}
}

void Brick8Disp::commit()
{
	for(unsigned int i=0;i<myMatPoints.size();i++) 
		myMatPoints[i]->getMaterial()->commit();
}
bool Brick8Disp::checkIfAllows(FEObject* f)
{
	return true;
}
void Brick8Disp::addInitialStresses(InitialStresses* pInitialStresses)
{
	if(myGroup->isActive()&&pInitialStresses->getGroupID()==myGroup->getID())
		for(unsigned i=0;i<myMatPoints.size();i++)
			myMatPoints[i]->setInitialStresses(pInitialStresses);
}
void Brick8Disp::recoverStresses()
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
const int Brick8Disp::getnPlasticPoints()
{
	int n=0;
	for(unsigned i=0;i<myMatPoints.size();i++)
		if(myMatPoints[i]->getMaterial()->isPlastic()) n+=1;
	return n;
}
