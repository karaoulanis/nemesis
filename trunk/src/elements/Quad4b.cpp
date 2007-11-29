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

#include <Quad4b.h>
#include <NemesisDebug.h>

Quad4b::Quad4b()
{
}
Quad4b::Quad4b(int ID,int Node_1,int Node_2,int Node_3,int Node_4,int MatID)
:Quad4(ID,Node_1,Node_2,Node_3,Node_4,MatID,2,2)
{
	static Matrix N(2,8);
	static Matrix B(4,8);
	double detJ;
	double dV=0.;
	volume=1.0;
	Bb1.resize(4,0.);
	Bb2.resize(4,0.);

	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		double xi =myMatPoints[k]->get_r();
		double eta=myMatPoints[k]->get_s();
		this->shape(xi,eta,N,B,detJ);
		dV+=detJ*(pD->getFac())*(myMatPoints[k]->get_w());
	}
	volume=dV;
	for(int i=0;i<4;i++) //nodes
	{
		for(unsigned k=0;k<myMatPoints.size();k++)
		{
			double xi =myMatPoints[k]->get_r();
			double eta=myMatPoints[k]->get_s();
			this->shape(xi,eta,N,B,detJ);
			Bb1[i]+=B(2,2*i+1)*detJ*(pD->getFac())*(myMatPoints[k]->get_w())/volume;
			Bb2[i]+=B(2,2*i  )*detJ*(pD->getFac())*(myMatPoints[k]->get_w())/volume;
		}
	}
}
Quad4b::~Quad4b()
{
}
const Matrix& Quad4b::getK()
{
	Matrix &K=*myMatrix;
	static Matrix C(4,4);
	static Matrix B(4,8);
	static Matrix N(2,8);
	double detJ,dV;

	K.clear();
	for(unsigned int k=0;k<myMatPoints.size();k++)
	{
		double xi =myMatPoints[k]->get_r();
		double eta=myMatPoints[k]->get_s();
		this->shape(xi,eta,N,B,detJ);
		dV=detJ*(pD->getFac())*(myMatPoints[k]->get_w());
		const Matrix& C0=myMatPoints[k]->getMaterial()->getC();
		C(0,0)=C0(0,0);	C(0,1)=C0(0,1);	C(0,2)=C0(0,3); C(0,3)=C0(0,2);  
		C(1,0)=C0(1,0);	C(1,1)=C0(1,1);	C(1,2)=C0(1,3); C(1,3)=C0(1,2);  
		C(2,0)=C0(3,0);	C(2,1)=C0(3,1);	C(2,2)=C0(3,3); C(2,3)=C0(3,2);
		C(3,0)=C0(2,0);	C(3,1)=C0(2,1);	C(3,2)=C0(2,3); C(3,3)=C0(2,2);
		K+=dV*(Transpose(B)*C*B);
	}
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	return K;
}
const Matrix& Quad4b::getM()
{
	Matrix &M=*myMatrix;
	M.clear();
	return M;
}
const Vector& Quad4b::getR()
{
 	static Vector sigma(6);
 	static Vector sigma4(4);
	static Matrix B(4,8);
	static Matrix N(2,8);
	double detJ,dV;
	static Vector b2(2);
	b2[0]=b[0];
	b2[1]=b[1];

	Vector& R=*myVector;
	R.clear();

	if(!(myGroup->isActive()))	return R;
	double facS=myGroup->getFacS();
	double facG=myGroup->getFacG();
	double facP=myGroup->getFacP();

	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		sigma=myMatPoints[k]->getMaterial()->getStress();
		sigma4[0]=sigma[0];
		sigma4[1]=sigma[1];
		sigma4[2]=sigma[3];
		sigma4[3]=sigma[2];

		double xi =myMatPoints[k]->get_r();
		double eta=myMatPoints[k]->get_s();
		this->shape(xi,eta,N,B,detJ);
		dV=detJ*(pD->getFac())*(myMatPoints[k]->get_w());

		R+=dV*(Transpose(B)*sigma4+Transpose(N)*b2);
	}
	R-=facP*P;
	//report(R);
	return R;
}
void Quad4b::update()
{
	static Vector Du(8);
	static Matrix B(4,8);
	static Matrix N(2,8);
	double detJ;

	if(!(myGroup->isActive()))	return;
	Du=this->getDispIncrm();
	
	// For each material point
	for(unsigned int k=0;k<myMatPoints.size();k++)
	{
		static Vector De(6),De4(4);
		double xi =myMatPoints[k]->get_r();
		double eta=myMatPoints[k]->get_s();
		// Determine the strain
		this->shape(xi,eta,N,B,detJ);
		De4=B*Du;
		De.clear();
		De[0]=De4[0];
		De[1]=De4[1];
		De[2]=De4[3];
		De[3]=De4[2];
		// And send it to the material point
		myMatPoints[k]->getMaterial()->setStrain(De);
//		report(De);
	}
}
void Quad4b::shape(double xi,double eta,Matrix& N,Matrix& B,double& detJ)
{
	double N1=0.25*(1-xi)*(1-eta);						// N1
	double N2=0.25*(1+xi)*(1-eta);						// N2
	double N3=0.25*(1+xi)*(1+eta);						// N3
	double N4=0.25*(1-xi)*(1+eta);						// N4
	N(0,0)=N1;	N(0,1)=0.;	N(0,2)=N2;	N(0,3)=0.;	N(0,4)=N3;	N(0,5)=0.;	N(0,6)=N4;	N(0,7)=0.;	
	N(1,0)=0.;	N(1,1)=N1;	N(1,2)=0.;	N(1,3)=N2;	N(1,4)=0.;	N(1,5)=N3;	N(1,6)=0.;	N(1,7)=N4;	

	static Matrix J(2,2);
	J(0,0)=0.25*(-x(0,0)*(1-eta) +x(1,0)*(1-eta) +x(2,0)*(1+eta) -x(3,0)*(1+eta));
	J(0,1)=0.25*(-x(0,0)*(1-xi)  -x(1,0)*(1+xi)  +x(2,0)*(1+xi)  +x(3,0)*(1-xi));
	J(1,0)=0.25*(-x(0,1)*(1-eta) +x(1,1)*(1-eta) +x(2,1)*(1+eta) -x(3,1)*(1+eta));
	J(1,1)=0.25*(-x(0,1)*(1-xi)  -x(1,1)*(1+xi)  +x(2,1)*(1+xi)  +x(3,1)*(1-xi));

	detJ=J(0,0)*J(1,1)-J(0,1)*J(1,0);			// detJ
	
	double dxidx  = J(1,1)/detJ;
	double detadx =-J(1,0)/detJ;
	double dxidy  =-J(0,1)/detJ;
	double detady = J(0,0)/detJ;

    double N1x=-0.25*(1-eta)*dxidx -0.25*(1-xi)*detadx;	// N1,1
    double N2x=+0.25*(1-eta)*dxidx -0.25*(1+xi)*detadx;	// N2,1
    double N3x=+0.25*(1+eta)*dxidx +0.25*(1+xi)*detadx;	// N3,1
    double N4x=-0.25*(1+eta)*dxidx +0.25*(1-xi)*detadx;	// N4,1

    double N1y=-0.25*(1-eta)*dxidy -0.25*(1-xi)*detady;	// N1,2
    double N2y=+0.25*(1-eta)*dxidy -0.25*(1+xi)*detady;	// N2,2
    double N3y=+0.25*(1+eta)*dxidy +0.25*(1+xi)*detady;	// N3,2
    double N4y=-0.25*(1+eta)*dxidy +0.25*(1-xi)*detady;	// N4,2

	double B1,B2,B4,B6,B7,B10,B11,B12;

	// Node 1
	B1=N1x;
	B2=N1y;
	B4=(Bb1[0]-B1)/3.;
	B6=(Bb2[0]-B2)/3.;
	B7=B2+B6;
	B10=B4;
	B11=B4;
	B12=B1+B4;
	B(0,0)=B12;	B(0,1)=B6;
	B(1,0)=B10;	B(1,1)=B7;
	B(2,0)=B2;	B(2,1)=B1;
	B(3,0)=B11;	B(3,1)=B6;

	// Node 2
	B1=N2x;
	B2=N2y;
	B4=(Bb1[1]-B1)/3.;
	B6=(Bb2[1]-B2)/3.;
	B7=B2+B6;
	B10=B4;
	B11=B4;
	B12=B1+B4;
	B(0,2)=B12;	B(0,3)=B6;
	B(1,2)=B10;	B(1,3)=B7;
	B(2,2)=B2;	B(2,3)=B1;
	B(3,2)=B11;	B(3,3)=B6;

	// Node 3
	B1=N3x;
	B2=N3y;
	B4=(Bb1[2]-B1)/3.;
	B6=(Bb2[2]-B2)/3.;
	B7=B2+B6;
	B10=B4;
	B11=B4;
	B12=B1+B4;
	B(0,4)=B12;	B(0,5)=B6;
	B(1,4)=B10;	B(1,5)=B7;
	B(2,4)=B2;	B(2,5)=B1;
	B(3,4)=B11;	B(3,5)=B6;

	// Node 4
	B1=N4x;
	B2=N4y;
	B4=(Bb1[3]-B1)/3.;
	B6=(Bb2[3]-B2)/3.;
	B7=B2+B6;
	B10=B4;
	B11=B4;
	B12=B1+B4;
	B(0,6)=B12;	B(0,7)=B6;
	B(1,6)=B10;	B(1,7)=B7;
	B(2,6)=B2;	B(2,7)=B1;
	B(3,6)=B11;	B(3,7)=B6;

}
