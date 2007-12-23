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

#include <Triangle3XFem.h>

Triangle3XFem::Triangle3XFem()
{
}
Triangle3XFem::Triangle3XFem(int ID,int Node_1,int Node_2,int Node_3,int matID)
:Element(ID,matID)
{
	myTag=TAG_ELEM_TRIANGLE_3_PRESSURE;
	// Get nodal data
	myNodalIDs.resize(3);
	myNodalIDs[0]=Node_1;
	myNodalIDs[1]=Node_2;
	myNodalIDs[2]=Node_3;
	// Set local nodal dofs
	myLocalNodalDofs.resize(2);
	myLocalNodalDofs[0]=0;
	myLocalNodalDofs[1]=1;	
	// Handle common info
	this->handleCommonInfo();

	// Find geometrical relations
	a1=x(1,0)*x(2,1)-x(2,0)*x(1,1);		b1=x(1,1)-x(2,1);	c1=x(2,0)-x(1,0);
	a2=x(2,0)*x(0,1)-x(0,0)*x(2,1);		b2=x(2,1)-x(0,1);	c2=x(0,0)-x(2,0);
	a3=x(0,0)*x(1,1)-x(1,0)*x(0,1);		b3=x(0,1)-x(1,1);	c3=x(1,0)-x(0,0);
	A=0.5*(x(1,0)*x(2,1)-x(2,0)*x(1,1)+x(2,0)*x(0,1)-x(0,0)*x(2,1)+x(0,0)*x(1,1)-x(1,0)*x(0,1));

	// Material
	myMatPoints.resize(1);
	MultiaxialMaterial* pMat=static_cast<MultiaxialMaterial*>(myMaterial); 
	myMatPoints[0]=new MatPoint(pMat,1,1,1,1);
	myMatPoints[0]->setX(num::d13*(x(0,0)+x(1,0)+x(2,0)),num::d13*(x(0,1)+x(1,1)+x(2,1)));
}
Triangle3XFem::~Triangle3XFem()
{
}
void Triangle3XFem::addInitialStresses(InitialStresses* pInitialStresses)
{
	if(myGroup->isActive()&&pInitialStresses->getGroupID()==myGroup->getID())
		for(unsigned i=0;i<myMatPoints.size();i++)
			myMatPoints[i]->setInitialStresses(pInitialStresses);
}
const Matrix& Triangle3XFem::getK()
{
	Matrix& K=*myMatrix;
	const Matrix& C=myMatPoints[0]->getMaterial()->getC();
	double coeff=pD->getFac()*0.25/A;
	K(0,0)=coeff*((b1*C(0,0)+c1*C(3,0))*b1+(b1*C(0,3)+c1*C(3,3))*c1);	K(0,1)=coeff*((b1*C(0,1)+c1*C(3,1))*c1+(b1*C(0,3)+c1*C(3,3))*b1);	K(0,2)=coeff*((b1*C(0,0)+c1*C(3,0))*b2+(b1*C(0,3)+c1*C(3,3))*c2);	K(0,3)=coeff*((b1*C(0,1)+c1*C(3,1))*c2+(b1*C(0,3)+c1*C(3,3))*b2);	K(0,4)=coeff*((b1*C(0,0)+c1*C(3,0))*b3+(b1*C(0,3)+c1*C(3,3))*c3);	K(0,5)=coeff*((b1*C(0,1)+c1*C(3,1))*c3+(b1*C(0,3)+c1*C(3,3))*b3);
	K(1,0)=coeff*((c1*C(1,0)+b1*C(3,0))*b1+(c1*C(1,3)+b1*C(3,3))*c1);	K(1,1)=coeff*((c1*C(1,1)+b1*C(3,1))*c1+(c1*C(1,3)+b1*C(3,3))*b1);	K(1,2)=coeff*((c1*C(1,0)+b1*C(3,0))*b2+(c1*C(1,3)+b1*C(3,3))*c2);	K(1,3)=coeff*((c1*C(1,1)+b1*C(3,1))*c2+(c1*C(1,3)+b1*C(3,3))*b2);	K(1,4)=coeff*((c1*C(1,0)+b1*C(3,0))*b3+(c1*C(1,3)+b1*C(3,3))*c3);	K(1,5)=coeff*((c1*C(1,1)+b1*C(3,1))*c3+(c1*C(1,3)+b1*C(3,3))*b3);
	K(2,0)=coeff*((b2*C(0,0)+c2*C(3,0))*b1+(b2*C(0,3)+c2*C(3,3))*c1);	K(2,1)=coeff*((b2*C(0,1)+c2*C(3,1))*c1+(b2*C(0,3)+c2*C(3,3))*b1);	K(2,2)=coeff*((b2*C(0,0)+c2*C(3,0))*b2+(b2*C(0,3)+c2*C(3,3))*c2);	K(2,3)=coeff*((b2*C(0,1)+c2*C(3,1))*c2+(b2*C(0,3)+c2*C(3,3))*b2);	K(2,4)=coeff*((b2*C(0,0)+c2*C(3,0))*b3+(b2*C(0,3)+c2*C(3,3))*c3);	K(2,5)=coeff*((b2*C(0,1)+c2*C(3,1))*c3+(b2*C(0,3)+c2*C(3,3))*b3);
	K(3,0)=coeff*((c2*C(1,0)+b2*C(3,0))*b1+(c2*C(1,3)+b2*C(3,3))*c1);	K(3,1)=coeff*((c2*C(1,1)+b2*C(3,1))*c1+(c2*C(1,3)+b2*C(3,3))*b1);	K(3,2)=coeff*((c2*C(1,0)+b2*C(3,0))*b2+(c2*C(1,3)+b2*C(3,3))*c2);	K(3,3)=coeff*((c2*C(1,1)+b2*C(3,1))*c2+(c2*C(1,3)+b2*C(3,3))*b2);	K(3,4)=coeff*((c2*C(1,0)+b2*C(3,0))*b3+(c2*C(1,3)+b2*C(3,3))*c3);	K(3,5)=coeff*((c2*C(1,1)+b2*C(3,1))*c3+(c2*C(1,3)+b2*C(3,3))*b3);
	K(4,0)=coeff*((b3*C(0,0)+c3*C(3,0))*b1+(b3*C(0,3)+c3*C(3,3))*c1);	K(4,1)=coeff*((b3*C(0,1)+c3*C(3,1))*c1+(b3*C(0,3)+c3*C(3,3))*b1);	K(4,2)=coeff*((b3*C(0,0)+c3*C(3,0))*b2+(b3*C(0,3)+c3*C(3,3))*c2);	K(4,3)=coeff*((b3*C(0,1)+c3*C(3,1))*c2+(b3*C(0,3)+c3*C(3,3))*b2);	K(4,4)=coeff*((b3*C(0,0)+c3*C(3,0))*b3+(b3*C(0,3)+c3*C(3,3))*c3);	K(4,5)=coeff*((b3*C(0,1)+c3*C(3,1))*c3+(b3*C(0,3)+c3*C(3,3))*b3);
	K(5,0)=coeff*((c3*C(1,0)+b3*C(3,0))*b1+(c3*C(1,3)+b3*C(3,3))*c1);	K(5,1)=coeff*((c3*C(1,1)+b3*C(3,1))*c1+(c3*C(1,3)+b3*C(3,3))*b1);	K(5,2)=coeff*((c3*C(1,0)+b3*C(3,0))*b2+(c3*C(1,3)+b3*C(3,3))*c2);	K(5,3)=coeff*((c3*C(1,1)+b3*C(3,1))*c2+(c3*C(1,3)+b3*C(3,3))*b2);	K(5,4)=coeff*((c3*C(1,0)+b3*C(3,0))*b3+(c3*C(1,3)+b3*C(3,3))*c3);	K(5,5)=coeff*((c3*C(1,1)+b3*C(3,1))*c3+(c3*C(1,3)+b3*C(3,3))*b3);
	double facK=1e-7;
	if(myGroup->isActive()) facK=myGroup->getFacK();
	K*=facK;
	return K;
}
const Matrix& Triangle3XFem::getM()
{
	return *myMatrix;
}
const Vector& Triangle3XFem::getR()
{
	Vector& R=*myVector;
	R.clear();
	if(!(myGroup->isActive()))	return R;
	double facS=myGroup->getFacS();
	double facG=myGroup->getFacG();
	double facP=myGroup->getFacP();
	double s1=(myMatPoints[0]->getMaterial()->getStress())[0];
	double s2=(myMatPoints[0]->getMaterial()->getStress())[1];
	double s3=(myMatPoints[0]->getMaterial()->getStress())[3];
	double fac=facS*0.5*(pD->getFac());
	double facb0=facG*(pD->getFac())*A*num::d13*b[0];
	double facb1=facG*(pD->getFac())*A*num::d13*b[1];
	R[0]=fac*(b1*s1+c1*s3)-facb0;
	R[1]=fac*(c1*s2+b1*s3)-facb1;
	R[2]=fac*(b2*s1+c2*s3)-facb0;
	R[3]=fac*(c2*s2+b2*s3)-facb1;
	R[4]=fac*(b3*s1+c3*s3)-facb0;
	R[5]=fac*(c3*s2+b3*s3)-facb1;
	R-=facP*P;
	return R;
}
void Triangle3XFem::update()
{
	if(!(myGroup->isActive()))	return;
	Vector& u=*myVector;
	u=this->getDispIncrm();
	// Determine the strain
	static Vector epsilon(6);
  	epsilon.clear();
	epsilon[0]=(0.5/A)*(b1*u[0]+b2*u[2]+b3*u[4]);
	epsilon[1]=(0.5/A)*(c1*u[1]+c2*u[3]+c3*u[5]);
	epsilon[3]=(0.5/A)*(c1*u[0]+b1*u[1]+c2*u[2]+b2*u[3]+c3*u[4]+b3*u[5]);
	// And send it to the material point
	myMatPoints[0]->getMaterial()->setStrain(epsilon);
}
void Triangle3XFem::commit()
{
	myMatPoints[0]->getMaterial()->commit();
}
void Triangle3XFem::recoverStresses()
{
	static Vector sigma(6);
 	sigma=myMatPoints[0]->getMaterial()->getStress();
	for(int i=0;i<3;i++) myNodes[i]->addStress(sigma);
}
bool Triangle3XFem::checkIfAllows(FEObject* f)
{
	return true;
}
const int Triangle3XFem::getnPlasticPoints()
{
	int n=0;
	if(myMatPoints[0]->getMaterial()->isPlastic()) n=1;
	return n;
}
void Triangle3XFem::enrich()
{
}