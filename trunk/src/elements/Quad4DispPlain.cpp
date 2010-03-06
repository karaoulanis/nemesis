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

#include <Quad4DispPlain.h>

Quad4DispPlain::Quad4DispPlain()
{
}
Quad4DispPlain::Quad4DispPlain(int ID,int Node_1,int Node_2,int Node_3,int Node_4,int MatID,
								int integrationRuleXi,int integrationRuleEta)
:Quad4(ID,Node_1,Node_2,Node_3,Node_4,MatID,integrationRuleXi,integrationRuleEta)
{
}
Quad4DispPlain::~Quad4DispPlain()
{
}
const Matrix& Quad4DispPlain::getK()
{
	Matrix &K=*myMatrix;
	K.clear();
	for(unsigned int k=0;k<myMatPoints.size();k++)
	{
		this->findShapeFunctionsAt(myMatPoints[k]);
		double dV=detJ*(pD->getFac())*(myMatPoints[k]->get_w());
		const Matrix& C=myMatPoints[k]->getMaterial()->getC();
		int ii=0;
		for(int i=0;i<4;i++) 
		{
			int jj=0;
			for(int j=0;j<4;j++) 
			{
				static Matrix CB(3,2);
				CB(0,0)=C(0,0)*N(1,j) + C(0,3)*N(2,j);
				CB(1,0)=C(1,0)*N(1,j) + C(1,3)*N(2,j);
				CB(2,0)=C(3,0)*N(1,j) + C(3,3)*N(2,j);
				CB(0,1)=C(0,1)*N(2,j) + C(0,3)*N(1,j);
				CB(1,1)=C(1,1)*N(2,j) + C(1,3)*N(1,j);
				CB(2,1)=C(3,1)*N(2,j) + C(3,3)*N(1,j);

				K(ii  ,jj  ) += (N(1,i)*CB(0,0) + N(2,i)*CB(2,0))*dV;
				K(ii  ,jj+1) += (N(1,i)*CB(0,1) + N(2,i)*CB(2,1))*dV;
				K(ii+1,jj  ) += (N(2,i)*CB(1,0) + N(1,i)*CB(2,0))*dV;
				K(ii+1,jj+1) += (N(2,i)*CB(1,1) + N(1,i)*CB(2,1))*dV;
				jj+=2;
			}
			ii+=2;
		}
	}
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	return K;
}
const Matrix& Quad4DispPlain::getM()
{
	Matrix &M=*myMatrix;
	M.clear();
	double rho=myMaterial->getRho();
	double volume=0.;
	for(unsigned k=0;k<myMatPoints.size();k++)
	{
		this->findShapeFunctionsAt(myMatPoints[k]);
		volume+=detJ*(pD->getFac())*(myMatPoints[k]->get_w()); 
	}
	double mass=rho*volume;
	for(int i=0;i<8;i++) M(i,i)=0.25*mass;
	return M;
}
const Vector& Quad4DispPlain::getR()
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
		for(int i=0;i<4;i++)
		{
			R[ii  ]+=facS*(N(1,i)*sigma[0]+N(2,i)*sigma[3])*dV;
			R[ii+1]+=facS*(N(2,i)*sigma[1]+N(1,i)*sigma[3])*dV;
			R[ii  ]-=facG*(N(0,i)*b[0]*dV);
			R[ii+1]-=facG*(N(0,i)*b[1]*dV);
			ii+=2;
		}
	}
	R-=facP*P;
	return R;
}
void Quad4DispPlain::update()
{
	if(!(myGroup->isActive()))	return;
	Vector& u=*myVector;
	///todo: Change this with incrm
	static Vector u1(8);
	static Vector u2(8);
	u2=this->getDispConvg();
	u1=this->getDispTrial();
	u=u1-u2;
	// For each material point
	for(unsigned int i=0;i<myMatPoints.size();i++)
	{
		this->findShapeFunctionsAt(myMatPoints[i]);
		// Determine the strain
		static Vector epsilon(6);
  		epsilon.clear();
		epsilon[0]=N(1,0)*u[0]+N(1,1)*u[2]+N(1,2)*u[4]+N(1,3)*u[6];
		epsilon[1]=N(2,0)*u[1]+N(2,1)*u[3]+N(2,2)*u[5]+N(2,3)*u[7];
		epsilon[3]=N(2,0)*u[0]+N(1,0)*u[1]+N(2,1)*u[2]+N(1,1)*u[3]+
				   N(2,2)*u[4]+N(1,2)*u[5]+N(2,3)*u[6]+N(1,3)*u[7];
		// And send it to the material point
		myMatPoints[i]->getMaterial()->setStrain(epsilon);
	}
}
