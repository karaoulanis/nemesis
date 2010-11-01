/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2010 F.E.Karaoulanis [http://www.nemesis-project.org]     *
*                                                                              *
* This program is free software; you can redistribute it and/or modify         *
* it under the terms of the GNU General Public License version 3, as           *
* published by the Free Software Foundation.                                   *
*                                                                              *
* This program is distributed in the hope that it will be useful,              *
* but WITHOUT ANY WARRANTY; without even the implied warranty of               *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
* GNU General Public License for more details.                                 *
*                                                                              *
* You should have received a copy of the GNU General Public License            *
* along with this program.  If not, see < http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "elements/timoshenko2d.h"

Timoshenko2d::Timoshenko2d() {
}
Timoshenko2d::Timoshenko2d(int ID, int Node_1, int Node_2, int matID, int secID, int rule)
:Element(ID, matID) {
  myTag = TAG_ELEM_BEAM_2D_TIMOSHENKO_2;
  myNodalIDs.resize(2);
  myNodalIDs[0]=Node_1;
  myNodalIDs[1]=Node_2;
  myLocalNodalDofs.resize(3);
  myLocalNodalDofs[0]=0;
  myLocalNodalDofs[1]=1;
  myLocalNodalDofs[2]=5;
  mySecID = secID;
  // Handle common info
  this->handleCommonInfo();
  mySecID = secID;
  mySection = pD->get < CrossSection>(pD->getCrossSections(), mySecID);
  L = sqrt((x(1, 1)-x(0, 1))*(x(1, 1)-x(0, 1))+(x(1, 0)-x(0, 0))*(x(1, 0)-x(0, 0)));
  myUniMaterial = static_cast < UniaxialMaterial*>(myMaterial);
  cosX[0]=(x(1, 0)-x(0, 0))/L;
  cosX[1]=(x(1, 1)-x(0, 1))/L;
  // Self weight - Transform vector b to local system
  ///@todo: check this
  double A = mySection->getA();
  double b0A = b[0]*A;
  double b1A = b[1]*A;
  b[0]= cosX[0]*b0A+cosX[1]*b1A;
  b[1]=-cosX[1]*b0A+cosX[0]*b1A;
  gPoints = rule;
  if (rule != 1&&rule != 2) throw SException("[nemesis:%d] %s", 9999, "Integration rule should be 1 or 2.");
}
Timoshenko2d::Timoshenko2d(int ID, int Node_1, int Node_2, int Node_3, int matID, int secID, int rule)
:Element(ID, matID) {
  myTag = TAG_ELEM_BEAM_2D_TIMOSHENKO_3;
  myNodalIDs.resize(3);
  myNodalIDs[0]=Node_1;
  myNodalIDs[1]=Node_2;
  myNodalIDs[2]=Node_3;
  myLocalNodalDofs.resize(3);
  myLocalNodalDofs[0]=0;
  myLocalNodalDofs[1]=1;
  myLocalNodalDofs[2]=5;
  mySecID = secID;
  // Handle common info
  this->handleCommonInfo();
  mySecID = secID;
  mySection = pD->get < CrossSection>(pD->getCrossSections(), mySecID);
  L = sqrt((x(1, 1)-x(0, 1))*(x(1, 1)-x(0, 1))+(x(1, 0)-x(0, 0))*(x(1, 0)-x(0, 0)));
  myUniMaterial = static_cast < UniaxialMaterial*>(myMaterial);
  cosX[0]=(x(1, 0)-x(0, 0))/L;
  cosX[1]=(x(1, 1)-x(0, 1))/L;
  // Self weight - Transform vector b to local system
  ///@todo: check this
  double A = mySection->getA();
  double b0A = b[0]*A;
  double b1A = b[1]*A;
  b[0]= cosX[0]*b0A+cosX[1]*b1A;
  b[1]=-cosX[1]*b0A+cosX[0]*b1A;
  gPoints = rule;
  if (rule != 2&&rule != 3) throw SException("[nemesis:%d] %s", 9999, "Integration rule should be 2 or 3.");
}
/**
 * Destructor.
 */
Timoshenko2d::~Timoshenko2d() {
}
void Timoshenko2d::shapeFunctions(int n, double xi, double &N, double &dN) {
  switch(myNodes.size())
  {
    case 2:
      if (n == 0)    {N = 0.5*(1-xi);  dN=-1/L;}
      else      {N = 0.5*(1+xi);  dN=+1/L;}
      break;
    case 3:
      if (n == 0)    {N =-0.5*xi*(1-xi); dN=(-1+2*xi)/L;}
      else if (n == 1) {N = 0.5*xi*(1+xi); dN=(+1+2*xi)/L;}
      else      {N = (1+xi)*(1-xi); dN=-4*xi/L;}
      break;
    default:
      break;
  }
}
const Matrix& Timoshenko2d::getK() {
  Matrix& K=*myMatrix;
  K.clear();
  double E =myUniMaterial->getParam(0);
  double nu = myUniMaterial->getParam(1);
  double A =mySection->getA();
  double J =mySection->getJ3();
  double a =5./6.; ///@todo

  double c = cosX[0];
  double s = cosX[1];
  double C1 = E*A;
  double C2 = E*J;
  double C3 = a*E/(2*(1+nu))*A;

  for (int k = 0; k < gPoints; k++) {
    double xi = GaussCoords[gPoints][k+1];
    double dx = GaussWeights[gPoints][k+1]*0.5*L;
    for (unsigned i = 0; i < myNodes.size(); i++) {
      for (unsigned j = 0; j < myNodes.size(); j++) {
        double Ni, Nj, dNi, dNj;
        this->shapeFunctions(i, xi, Ni, dNi);
        this->shapeFunctions(j, xi, Nj, dNj);
        K(3*i+0, 3*j+0)+=(dNi*dNj*(C1*c*c+C3-c*c*C3)  )*dx;
        K(3*i+0, 3*j+1)+=(c*dNi*dNj*s*(-C3+C1)    )*dx;
        K(3*i+0, 3*j+2)+=(s*dNi*C3*Nj         )*dx;
        K(3*i+1, 3*j+0)+=(c*dNi*dNj*s*(-C3+C1)    )*dx;
        K(3*i+1, 3*j+1)+=(-dNi*dNj*(-c*c*C3-C1+C1*c*c))*dx;
        K(3*i+1, 3*j+2)+=(-c*dNi*C3*Nj        )*dx;
        K(3*i+2, 3*j+0)+=(Ni*C3*dNj*s         )*dx;
        K(3*i+2, 3*j+1)+=(-Ni*C3*dNj*c        )*dx;
        K(3*i+2, 3*j+2)+=(dNi*C2*dNj+Ni*C3*Nj     )*dx;
      }
    }
  }
  if (myGroup->isActive()) K*=myGroup->getFacK();
  else          K*=1e-12;
  return K;
}
const Matrix& Timoshenko2d::getM() {
  Matrix& M=*myMatrix;
  M.clear();
  return M;
}
const Vector& Timoshenko2d::getRgrad() {
  myVector->clear();
  return *myVector;
}
const Vector& Timoshenko2d::getR() {
  Vector& R=*myVector;
  R.clear();
  if (!(myGroup->isActive()))  return R;
  double facS = myGroup->getFacS();
  double facG = myGroup->getFacG();
  double facP = myGroup->getFacP();

  double E =myUniMaterial->getParam(0);
  double nu = myUniMaterial->getParam(1);
  double A =mySection->getA();
  double J =mySection->getJ3();
  double a =5./6.; ///@todo
  double c = cosX[0];
  double s = cosX[1];
  double Ni, dNi;
  
  Vector u(3*myNodes.size());
  static Vector epsilon(3);
  u = this->getDispTrial();
  R.clear();  // getDispTrial() and R use the same static matrix
  for (int k = 0; k < gPoints; k++) {
    // Find epsilon 
    ///@todo This should happen to update() for elastoplastic computations
    double xi = GaussCoords[gPoints][k+1];
    double dx = GaussWeights[gPoints][k+1]*0.5*L;
    epsilon.clear();
    for (unsigned i = 0; i < myNodes.size(); i++) {
      this->shapeFunctions(i, xi, Ni, dNi);
      epsilon[0]+=dNi*(c*u[0+3*i]+s*u[1+3*i]);
      epsilon[1]+=dNi*u[2+3*i];
      epsilon[2]+=dNi*(-s*u[0+3*i]+c*u[1+3*i])-Ni*u[2+3*i];
    }
    // Find sigma 
    ///@todo This should come from material for elastoplastic computations
    static Vector sigma(3);
    sigma[0]=E*A*epsilon[0];
    sigma[1]=E*J*epsilon[1];
    sigma[2]=a*E/(2*(1+nu))*A*epsilon[2];
    // R = Fint
    for (unsigned i = 0; i < myNodes.size(); i++) {
      this->shapeFunctions(i, xi, Ni, dNi);
      R[3*i+0]+=facS*((dNi*sigma[0])        )*dx;
      R[3*i+1]+=facS*((dNi*sigma[2])        )*dx;
      R[3*i+2]+=facS*((dNi*sigma[1]-Ni*sigma[2]))*dx;
    }
  }
  // R = R - Fext(Gravity)
  switch(myNodes.size())
  {
    case 2:
      R[0]-=facG*(0.5*b[0]*L);
      R[1]-=facG*(0.5*b[1]*L);
      R[3]-=facG*(0.5*b[0]*L);
      R[4]-=facG*(0.5*b[1]*L);
      break;
    case 3:
      R[0]-=facG*(num::d16*0.5*b[0]*L);
      R[1]-=facG*(num::d16*0.5*b[1]*L);
      R[3]-=facG*(num::d16*0.5*b[0]*L);
      R[4]-=facG*(num::d16*0.5*b[1]*L);
      R[6]-=facG*(num::d13*0.5*b[0]*L);
      R[7]-=facG*(num::d13*0.5*b[1]*L);
      break;
    default:
      break;
  }
  // R = R - Fext(Loads)
  R-=facP*P;
  
  // Local R to global R
  for (unsigned i = 0; i < myNodes.size(); i++) {
    double d1 = c*R[3*i+0]-s*R[3*i+1];
    double d2 = s*R[3*i+0]+c*R[3*i+1];
    R[3*i+0]=d1;
    R[3*i+1]=d2;
  }
  return R;
}
void Timoshenko2d::recoverStresses() {
  ///@todo This might affect nodal stresses.
  return;
}
const double Timoshenko2d::GaussCoords[4][4]=  {{ 0.000000000000000,  // Rule 0
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                          0.000000000000000}, //_______
                        { 0.000000000000000,  // Rule 1
                                                  0.000000000000000,  // [1][1] 
                                                  0.000000000000000, 
                          0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 2
                                                 -0.577350269189626,  // [2][1]
                                                 +0.577350269189626,  // [2][2]
                         0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 3
                                                 -0.774596669241483,  // [3][1]
                                                  0.000000000000000,  // [3][2]
                         +0.774596669241483}};// [3][3]
 
const double Timoshenko2d::GaussWeights[4][4]= {{ 0.000000000000000,  // Rule 0
                                                  0.000000000000000, 
                                                  0.000000000000000, 
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 1
                                                 +2.000000000000000,  // [1][1] 
                                                  0.000000000000000, 
                          0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 2
                                                 +1.000000000000000,  // [2][1]
                                                 +1.000000000000000,  // [2][1]
                                                  0.000000000000000}, //_______
                                                { 0.000000000000000,  // Rule 3
                                                 +0.555555555555556,  // [3][1]
                                                 +0.888888888888889,  // [3][2]
                         +0.555555555555556}};// [3][3]
