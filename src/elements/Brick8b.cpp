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

#include <Brick8b.h>
#include <NemesisDebug.h>

void add(Matrix& K,int row,int col,const Matrix& B1,const Matrix& C,const Matrix B2,double c1,double c0=0.)
{
	int m=B1.rows();
	int n=B2.cols();
	int pos=row*K.cols()+col;
	int colK=K.cols();
	double* pK =K.data();
	double* pB1=B1.data();
	double* pB2=B2.data();
	double* pC=C.data();
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			pK[pos+i*colK+j]*=c0;
	for(int i=0;i<n;i++)
		for(int k=0;k<m;k++)
			for(int l=0;l<m;l++)
				for(int j=0;j<n;j++)
					pK[pos+i*colK+j]+=c1*pB1[k*n+i]*pC[k*m+l]*pB2[l*n+j];
}
void add(Vector& R,int row,const Matrix& BT,const Vector& V,double c1,double c0=0.)
{
	int m=BT.rows();
	int n=BT.cols();
	double* pV =V.data();
	double* pR =R.data();
	double* pBT=BT.data();
	for(int i=0;i<n;i++)
		pR[row+i]*=c0;
	for(int i=0;i<n;i++)
		for(int j=0;j<m;j++)
			pR[row+i]+=c1*pBT[j*n+i]*pV[j];
}
void add2(Vector& R,int row,const Matrix& BT,const Vector& V,double c1,double c0=0.)
{
	int m=BT.rows();
	int n=BT.cols();
	double* pV =V.data();
	double* pR =R.data();
	double* pBT=BT.data();
	for(int i=0;i<n;i++)
		pR[row+i]*=c0;
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
			pR[i]+=c1*pBT[i*n+j]*pV[row+j];
}
Vector Brick8b::N(8);
Matrix Brick8b::B[8];
double Brick8b::detJ;
bool Brick8b::first=true;

Brick8b::Brick8b()
{
}
Brick8b::Brick8b(int ID,
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
	
	if(first)
	{
		for(int i=0;i<8;i++) 
			B[i].resize(6,3);
		first=false;
	}
	
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

	Bb1.resize(8,0.);
	Bb2.resize(8,0.);
	Bb3.resize(8,0.);
	this->formBb(Bb1,Bb2,Bb3);	

	for(unsigned i=0;i<myMatPoints.size();i++)
	{
		this->shapeFunctions(myMatPoints[i]);
		double xG=0,yG=0,zG=0;
		for(unsigned j=0;j<myNodes.size();j++)
		{
			xG+=N[j]*x(j,0);
			yG+=N[j]*x(j,1);
			zG+=N[j]*x(j,2);
		}
		myMatPoints[i]->setX(xG,yG,zG);
	}

}
Brick8b::~Brick8b()
{
	Containers::vector_delete(myMatPoints);
}
void Brick8b::shapeFunctions(MatPoint* pMatPoint)
{
	double xi=  pMatPoint->get_r();
	double eta= pMatPoint->get_s();
	double zeta=pMatPoint->get_t();

	N[0]=0.125*(1-xi)*(1-eta)*(1-zeta);			// N1
	N[1]=0.125*(1+xi)*(1-eta)*(1-zeta);			// N2
	N[2]=0.125*(1+xi)*(1+eta)*(1-zeta);			// N3
	N[3]=0.125*(1-xi)*(1+eta)*(1-zeta);			// N4
	N[4]=0.125*(1-xi)*(1-eta)*(1+zeta);			// N5
	N[5]=0.125*(1+xi)*(1-eta)*(1+zeta);			// N6
	N[6]=0.125*(1+xi)*(1+eta)*(1+zeta);			// N7
	N[7]=0.125*(1-xi)*(1+eta)*(1+zeta);			// N8

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

	static double dNdX[8][3];
	dNdX[0][0]=0.125*(-(1-eta)*(1-zeta)*dxidx-(1-xi)*(1-zeta)*detadx-(1-xi)*(1-eta)*dzetadx);
	dNdX[0][1]=0.125*(-(1-eta)*(1-zeta)*dxidy-(1-xi)*(1-zeta)*detady-(1-xi)*(1-eta)*dzetady);
	dNdX[0][2]=0.125*(-(1-eta)*(1-zeta)*dxidz-(1-xi)*(1-zeta)*detadz-(1-xi)*(1-eta)*dzetadz);

	dNdX[1][0]=0.125*(+(1-eta)*(1-zeta)*dxidx-(1+xi)*(1-zeta)*detadx-(1+xi)*(1-eta)*dzetadx);
	dNdX[1][1]=0.125*(+(1-eta)*(1-zeta)*dxidy-(1+xi)*(1-zeta)*detady-(1+xi)*(1-eta)*dzetady);
	dNdX[1][2]=0.125*(+(1-eta)*(1-zeta)*dxidz-(1+xi)*(1-zeta)*detadz-(1+xi)*(1-eta)*dzetadz);

	dNdX[2][0]=0.125*(+(1+eta)*(1-zeta)*dxidx+(1+xi)*(1-zeta)*detadx-(1+xi)*(1+eta)*dzetadx);
	dNdX[2][1]=0.125*(+(1+eta)*(1-zeta)*dxidy+(1+xi)*(1-zeta)*detady-(1+xi)*(1+eta)*dzetady);
	dNdX[2][2]=0.125*(+(1+eta)*(1-zeta)*dxidz+(1+xi)*(1-zeta)*detadz-(1+xi)*(1+eta)*dzetadz);

	dNdX[3][0]=0.125*(-(1+eta)*(1-zeta)*dxidx+(1-xi)*(1-zeta)*detadx-(1-xi)*(1+eta)*dzetadx);
	dNdX[3][1]=0.125*(-(1+eta)*(1-zeta)*dxidy+(1-xi)*(1-zeta)*detady-(1-xi)*(1+eta)*dzetady);
	dNdX[3][2]=0.125*(-(1+eta)*(1-zeta)*dxidz+(1-xi)*(1-zeta)*detadz-(1-xi)*(1+eta)*dzetadz);
	
	dNdX[4][0]=0.125*(-(1-eta)*(1+zeta)*dxidx-(1-xi)*(1+zeta)*detadx+(1-xi)*(1-eta)*dzetadx);
	dNdX[4][1]=0.125*(-(1-eta)*(1+zeta)*dxidy-(1-xi)*(1+zeta)*detady+(1-xi)*(1-eta)*dzetady);
	dNdX[4][2]=0.125*(-(1-eta)*(1+zeta)*dxidz-(1-xi)*(1+zeta)*detadz+(1-xi)*(1-eta)*dzetadz);
	
	dNdX[5][0]=0.125*(+(1-eta)*(1+zeta)*dxidx-(1+xi)*(1+zeta)*detadx+(1+xi)*(1-eta)*dzetadx);
	dNdX[5][1]=0.125*(+(1-eta)*(1+zeta)*dxidy-(1+xi)*(1+zeta)*detady+(1+xi)*(1-eta)*dzetady);
	dNdX[5][2]=0.125*(+(1-eta)*(1+zeta)*dxidz-(1+xi)*(1+zeta)*detadz+(1+xi)*(1-eta)*dzetadz);
	
	dNdX[6][0]=0.125*(+(1+eta)*(1+zeta)*dxidx+(1+xi)*(1+zeta)*detadx+(1+xi)*(1+eta)*dzetadx);
	dNdX[6][1]=0.125*(+(1+eta)*(1+zeta)*dxidy+(1+xi)*(1+zeta)*detady+(1+xi)*(1+eta)*dzetady);
	dNdX[6][2]=0.125*(+(1+eta)*(1+zeta)*dxidz+(1+xi)*(1+zeta)*detadz+(1+xi)*(1+eta)*dzetadz);
	
	dNdX[7][0]=0.125*(-(1+eta)*(1+zeta)*dxidx+(1-xi)*(1+zeta)*detadx+(1-xi)*(1+eta)*dzetadx);
	dNdX[7][1]=0.125*(-(1+eta)*(1+zeta)*dxidy+(1-xi)*(1+zeta)*detady+(1-xi)*(1+eta)*dzetady);
	dNdX[7][2]=0.125*(-(1+eta)*(1+zeta)*dxidz+(1-xi)*(1+zeta)*detadz+(1-xi)*(1+eta)*dzetadz);

	for(int i=0;i<8;i++)
	{
		double B1=dNdX[i][0];
		double B2=dNdX[i][1];
		double B3=dNdX[i][2];
		double B4=(Bb1[i]-B1)/3.;
		double B5=B1+B4;
		double B6=(Bb2[i]-B2)/3.;
		double B7=B2+B6;
		double B8=(Bb3[i]-B3)/3.;
		double B9=B3+B8;
		B[i](0,0)=B5;	B[i](0,1)=B6;	B[i](0,2)=B8;
		B[i](1,0)=B4;	B[i](1,1)=B7;	B[i](1,2)=B8;
		B[i](2,0)=B4;	B[i](2,1)=B6;	B[i](2,2)=B9;
		B[i](3,0)=B2;	B[i](3,1)=B1;	B[i](3,2)=0.;
		B[i](4,0)=0.;	B[i](4,1)=B3;	B[i](4,2)=B2;
		B[i](5,0)=B3;	B[i](5,1)=0.;	B[i](5,2)=B1;
		//report(B[i]);
	}

}
const Matrix& Brick8b::getK()
{
	Matrix &K=*myMatrix;
	K.clear();
	for(unsigned int k=0;k<myMatPoints.size();k++)
	{
		this->shapeFunctions(myMatPoints[k]);
		double dV=detJ*(pD->getFac())*(myMatPoints[k]->get_w());
		const Matrix& C=myMatPoints[k]->getMaterial()->getC();
		for(unsigned a=0;a<myNodes.size();a++)
				for(unsigned b=0;b<myNodes.size();b++)
					add(K,3*a,3*b,B[a],C,B[b],dV,1.0);
	}
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	return K;
}
const Matrix& Brick8b::getM()
{
	Matrix &M=*myMatrix;
	M.clear();
	return M;
}
const Vector& Brick8b::getR()
{
	static Vector sigma(6);
	Vector& R=*myVector;
	R.clear();
	double facS=myGroup->getFacS();
	double facG=myGroup->getFacG();
	double facP=myGroup->getFacP();
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		sigma=myMatPoints[k]->getMaterial()->getStress();
		this->shapeFunctions(myMatPoints[k]);
		double dV=detJ*(pD->getFac())*(myMatPoints[k]->get_w());
		for(unsigned a=0;a<myNodes.size();a++)
			add(R,3*a,B[a],sigma,dV,1.);
	}
	//R.add_cV(-facP,P);
	return R;
}
void Brick8b::update()
{
	if(!(myGroup->isActive()))	return;
	static Vector u(24);
	u=this->getDispIncrm();
	//report(u,"dsp");
	static Vector epsilon(6);
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		epsilon.clear();
		this->shapeFunctions(myMatPoints[k]);
		for(unsigned a=0;a<myNodes.size();a++)
			add2(epsilon,3*a,B[a],u,1.0,1.0);
		myMatPoints[k]->getMaterial()->setStrain(epsilon);
		//report(epsilon,"eps",20,18);
	}
}

void Brick8b::commit()
{
	for(unsigned int i=0;i<myMatPoints.size();i++) 
		myMatPoints[i]->getMaterial()->commit();
}
bool Brick8b::checkIfAllows(FEObject* f)
{
	return true;
}
void Brick8b::addInitialStresses(InitialStresses* pInitialStresses)
{
	if(myGroup->isActive()&&pInitialStresses->getGroupID()==myGroup->getID())
		for(unsigned i=0;i<myMatPoints.size();i++)
			myMatPoints[i]->setInitialStresses(pInitialStresses);
}
void Brick8b::recoverStresses()
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
const int Brick8b::getnPlasticPoints()
{
	int n=0;
	for(unsigned i=0;i<myMatPoints.size();i++)
		if(myMatPoints[i]->getMaterial()->isPlastic()) n+=1;
	return n;
}
void Brick8b::formBb(Vector& Bb1,Vector& Bb2,Vector& Bb3)
{
	// find volume
	double vol=0.;
	static Matrix J(3,3);
	static Matrix invJ(3,3);
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		double xi=  myMatPoints[k]->get_r();
		double eta= myMatPoints[k]->get_s();
		double zeta=myMatPoints[k]->get_t();
		double w=myMatPoints[k]->get_w();

		N[0]=0.125*(1-xi)*(1-eta)*(1-zeta);			// N1
		N[1]=0.125*(1+xi)*(1-eta)*(1-zeta);			// N2
		N[2]=0.125*(1+xi)*(1+eta)*(1-zeta);			// N3
		N[3]=0.125*(1-xi)*(1+eta)*(1-zeta);			// N4
		N[4]=0.125*(1-xi)*(1-eta)*(1+zeta);			// N5
		N[5]=0.125*(1+xi)*(1-eta)*(1+zeta);			// N6
		N[6]=0.125*(1+xi)*(1+eta)*(1+zeta);			// N7
		N[7]=0.125*(1-xi)*(1+eta)*(1+zeta);			// N8

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

		vol+=det(J);
	}
	// find Bb's
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		double xi=  myMatPoints[k]->get_r();
		double eta= myMatPoints[k]->get_s();
		double zeta=myMatPoints[k]->get_t();
		double w=myMatPoints[k]->get_w();

		N[0]=0.125*(1-xi)*(1-eta)*(1-zeta);			// N1
		N[1]=0.125*(1+xi)*(1-eta)*(1-zeta);			// N2
		N[2]=0.125*(1+xi)*(1+eta)*(1-zeta);			// N3
		N[3]=0.125*(1-xi)*(1+eta)*(1-zeta);			// N4
		N[4]=0.125*(1-xi)*(1-eta)*(1+zeta);			// N5
		N[5]=0.125*(1+xi)*(1-eta)*(1+zeta);			// N6
		N[6]=0.125*(1+xi)*(1+eta)*(1+zeta);			// N7
		N[7]=0.125*(1-xi)*(1+eta)*(1+zeta);			// N8

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
		invJ=Inverse(J);

		double dxidx   = invJ(0,0);	  double detadx  = invJ(0,1);   double dzetadx = invJ(0,2);
		double dxidy   = invJ(1,0);	  double detady  = invJ(1,1);   double dzetady = invJ(1,2);
		double dxidz   = invJ(2,0);	  double detadz  = invJ(2,1);   double dzetadz = invJ(2,2);

		static double dNdX[8][3];
		dNdX[0][0]=0.125*(-(1-eta)*(1-zeta)*dxidx-(1-xi)*(1-zeta)*detadx-(1-xi)*(1-eta)*dzetadx);
		dNdX[0][1]=0.125*(-(1-eta)*(1-zeta)*dxidy-(1-xi)*(1-zeta)*detady-(1-xi)*(1-eta)*dzetady);
		dNdX[0][2]=0.125*(-(1-eta)*(1-zeta)*dxidz-(1-xi)*(1-zeta)*detadz-(1-xi)*(1-eta)*dzetadz);

		dNdX[1][0]=0.125*(+(1-eta)*(1-zeta)*dxidx-(1+xi)*(1-zeta)*detadx-(1+xi)*(1-eta)*dzetadx);
		dNdX[1][1]=0.125*(+(1-eta)*(1-zeta)*dxidy-(1+xi)*(1-zeta)*detady-(1+xi)*(1-eta)*dzetady);
		dNdX[1][2]=0.125*(+(1-eta)*(1-zeta)*dxidz-(1+xi)*(1-zeta)*detadz-(1+xi)*(1-eta)*dzetadz);

		dNdX[2][0]=0.125*(+(1+eta)*(1-zeta)*dxidx+(1+xi)*(1-zeta)*detadx-(1+xi)*(1+eta)*dzetadx);
		dNdX[2][1]=0.125*(+(1+eta)*(1-zeta)*dxidy+(1+xi)*(1-zeta)*detady-(1+xi)*(1+eta)*dzetady);
		dNdX[2][2]=0.125*(+(1+eta)*(1-zeta)*dxidz+(1+xi)*(1-zeta)*detadz-(1+xi)*(1+eta)*dzetadz);

		dNdX[3][0]=0.125*(-(1+eta)*(1-zeta)*dxidx+(1-xi)*(1-zeta)*detadx-(1-xi)*(1+eta)*dzetadx);
		dNdX[3][1]=0.125*(-(1+eta)*(1-zeta)*dxidy+(1-xi)*(1-zeta)*detady-(1-xi)*(1+eta)*dzetady);
		dNdX[3][2]=0.125*(-(1+eta)*(1-zeta)*dxidz+(1-xi)*(1-zeta)*detadz-(1-xi)*(1+eta)*dzetadz);
		
		dNdX[4][0]=0.125*(-(1-eta)*(1+zeta)*dxidx-(1-xi)*(1+zeta)*detadx+(1-xi)*(1-eta)*dzetadx);
		dNdX[4][1]=0.125*(-(1-eta)*(1+zeta)*dxidy-(1-xi)*(1+zeta)*detady+(1-xi)*(1-eta)*dzetady);
		dNdX[4][2]=0.125*(-(1-eta)*(1+zeta)*dxidz-(1-xi)*(1+zeta)*detadz+(1-xi)*(1-eta)*dzetadz);
		
		dNdX[5][0]=0.125*(+(1-eta)*(1+zeta)*dxidx-(1+xi)*(1+zeta)*detadx+(1+xi)*(1-eta)*dzetadx);
		dNdX[5][1]=0.125*(+(1-eta)*(1+zeta)*dxidy-(1+xi)*(1+zeta)*detady+(1+xi)*(1-eta)*dzetady);
		dNdX[5][2]=0.125*(+(1-eta)*(1+zeta)*dxidz-(1+xi)*(1+zeta)*detadz+(1+xi)*(1-eta)*dzetadz);
		
		dNdX[6][0]=0.125*(+(1+eta)*(1+zeta)*dxidx+(1+xi)*(1+zeta)*detadx+(1+xi)*(1+eta)*dzetadx);
		dNdX[6][1]=0.125*(+(1+eta)*(1+zeta)*dxidy+(1+xi)*(1+zeta)*detady+(1+xi)*(1+eta)*dzetady);
		dNdX[6][2]=0.125*(+(1+eta)*(1+zeta)*dxidz+(1+xi)*(1+zeta)*detadz+(1+xi)*(1+eta)*dzetadz);
		
		dNdX[7][0]=0.125*(-(1+eta)*(1+zeta)*dxidx+(1-xi)*(1+zeta)*detadx+(1-xi)*(1+eta)*dzetadx);
		dNdX[7][1]=0.125*(-(1+eta)*(1+zeta)*dxidy+(1-xi)*(1+zeta)*detady+(1-xi)*(1+eta)*dzetady);
		dNdX[7][2]=0.125*(-(1+eta)*(1+zeta)*dxidz+(1-xi)*(1+zeta)*detadz+(1-xi)*(1+eta)*dzetadz);

		for(int i=0;i<8;i++)
		{
			double dNdx=dNdX[i][0];
			double dNdy=dNdX[i][1];
			double dNdz=dNdX[i][2];
			Bb1[i]+=dNdx*w*detJ/vol;
			Bb2[i]+=dNdy*w*detJ/vol;
			Bb3[i]+=dNdz*w*detJ/vol;
		}
	}
}