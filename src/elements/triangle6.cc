/*******************************************************************************
* nemesis. an experimental finite element code.                                *
* Copyright (C) 2004-2011 F.E.Karaoulanis [http://www.nemesis-project.org]     *
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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "elements/triangle6.h"
#include "group/group_data.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"
#include "node/node.h"

Triangle6::Triangle6() {
}
Triangle6::Triangle6(int ID,
             int Node_1, int Node_2, int Node_3,
             int Node_4, int Node_5, int Node_6,
             int matID)
:Element(ID, matID) {
  myTag = TAG_ELEM_TRIANGLE_6;
  // Get nodal data
  myNodalIDs.resize(6);
  myNodalIDs[0]=Node_1;
  myNodalIDs[1]=Node_2;
  myNodalIDs[2]=Node_3;
  myNodalIDs[3]=Node_4;
  myNodalIDs[4]=Node_5;
  myNodalIDs[5]=Node_6;
  // Set local nodal dofs
  myLocalNodalDofs.resize(2);
  myLocalNodalDofs[0]=0;
  myLocalNodalDofs[1]=1;
  // Handle common info
  this->handleCommonInfo();

  // Materials
  myMatPoints.resize(3);
  MultiaxialMaterial* pMat = static_cast < MultiaxialMaterial*>(myMaterial);
  myMatPoints[0]=new MatPoint(pMat, num::d23, num::d16, num::d16, num::d13);
  myMatPoints[1]=new MatPoint(pMat, num::d16, num::d23, num::d16, num::d13);
  myMatPoints[2]=new MatPoint(pMat, num::d16, num::d16, num::d23, num::d13);
  static Matrix N(6, 3);
  static double detJ;
  for (unsigned i = 0; i < 3; i++) {
    this->get_shape_functions(myMatPoints[i], N, detJ);
    double xG = N(0, 0)*x(0, 0)+N(1, 0)*x(2, 0)+N(1, 0)*x(2, 0)
           +N(3, 0)*x(3, 0)+N(4, 0)*x(4, 0)+N(5, 0)*x(5, 0);
    double yG = N(0, 0)*x(0, 1)+N(1, 0)*x(1, 1)+N(2, 0)*x(2, 1)
           +N(3, 0)*x(3, 1)+N(4, 0)*x(4, 1)+N(5, 0)*x(5, 1);
    myMatPoints[i]->set_X(xG, yG);
  }
}
Triangle6::~Triangle6() {
  Containers::vector_delete(myMatPoints);
}
void Triangle6::get_shape_functions(MatPoint* pMatPoint,
                                    Matrix& N, double& detJ) {
  double z1 = pMatPoint->get_r();
  double z2 = pMatPoint->get_s();
  double z3 = pMatPoint->get_t();
  double x21 = x(1, 0)-x(0, 0);
  double x32 = x(2, 0)-x(1, 0);
  double x13 = x(0, 0)-x(2, 0);
  double y12 = x(0, 1)-x(1, 1);
  double y23 = x(1, 1)-x(2, 1);
  double y31 = x(2, 1)-x(0, 1);
  double Dx4 = x(3, 0)-0.5*(x(0, 0)+x(1, 0));
  double Dx5 = x(4, 0)-0.5*(x(1, 0)+x(2, 0));
  double Dx6 = x(5, 0)-0.5*(x(2, 0)+x(0, 0));
  double Dy4 = x(3, 1)-0.5*(x(0, 1)+x(1, 1));
  double Dy5 = x(4, 1)-0.5*(x(1, 1)+x(2, 1));
  double Dy6 = x(5, 1)-0.5*(x(2, 1)+x(0, 1));
  double Jx21 = x21+4*(Dx4*(z1-z2)+(Dx5-Dx6)*z3);
  double Jx32 = x32+4*(Dx5*(z2-z3)+(Dx6-Dx4)*z1);
  double Jx13 = x13+4*(Dx6*(z3-z1)+(Dx4-Dx5)*z2);
  double Jy12 = y12+4*(Dy4*(z2-z1)+(Dy6-Dy5)*z3);
  double Jy23 = y23+4*(Dy5*(z3-z2)+(Dy4-Dy6)*z1);
  double Jy31 = y31+4*(Dy6*(z1-z3)+(Dy5-Dy4)*z2);
  detJ = 0.5*(Jx21*Jy31-Jy12*Jx13);
  double J = 0.5/detJ;

  N(0, 0)=z1*(2*z1-1);
  N(1, 0)=z2*(2*z2-1);
  N(2, 0)=z3*(2*z3-1);
  N(3, 0)=4*z1*z2;
  N(4, 0)=4*z2*z3;
  N(5, 0)=4*z3*z1;
  N(0, 1)=J*(4*z1-1)*Jy23;
  N(1, 1)=J*(4*z2-1)*Jy31;
  N(2, 1)=J*(4*z3-1)*Jy12;
  N(3, 1)=J*4*(z2*Jy23+z1*Jy31);
  N(4, 1)=J*4*(z3*Jy31+z2*Jy12);
  N(5, 1)=J*4*(z1*Jy12+z3*Jy23);
  N(0, 2)=J*(4*z1-1)*Jx32;
  N(1, 2)=J*(4*z2-1)*Jx13;
  N(2, 2)=J*(4*z3-1)*Jx21;
  N(3, 2)=J*4*(z2*Jx32+z1*Jx13);
  N(4, 2)=J*4*(z3*Jx13+z2*Jx21);
  N(5, 2)=J*4*(z1*Jx21+z3*Jx32);
}
const Matrix& Triangle6::get_K() {
  Matrix &K=*myMatrix;
  K.Clear();
  for (unsigned int k = 0; k < myMatPoints.size(); k++) {
    static Matrix N(6, 3);
    static double detJ;
    this->get_shape_functions(myMatPoints[k], N, detJ);
    double dV = detJ*(pD->get_fac())*(myMatPoints[k]->get_w());
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    int ii = 0;
    for (int i = 0; i < 6; i++) {
      int jj = 0;
      for (int j = 0; j < 6; j++) {
        static Matrix CB(3, 2);
        CB(0, 0)=C(0, 0)*N(j, 1) + C(0, 3)*N(j, 2);
        CB(1, 0)=C(1, 0)*N(j, 1) + C(1, 3)*N(j, 2);
        CB(2, 0)=C(3, 0)*N(j, 1) + C(3, 3)*N(j, 2);
        CB(0, 1)=C(0, 1)*N(j, 2) + C(0, 3)*N(j, 1);
        CB(1, 1)=C(1, 1)*N(j, 2) + C(1, 3)*N(j, 1);
        CB(2, 1)=C(3, 1)*N(j, 2) + C(3, 3)*N(j, 1);

        K(ii  , jj)   += (N(i, 1)*CB(0, 0) + N(i, 2)*CB(2, 0))*dV;
        K(ii  , jj+1) += (N(i, 1)*CB(0, 1) + N(i, 2)*CB(2, 1))*dV;
        K(ii+1, jj)   += (N(i, 2)*CB(1, 0) + N(i, 1)*CB(2, 0))*dV;
        K(ii+1, jj+1) += (N(i, 2)*CB(1, 1) + N(i, 1)*CB(2, 1))*dV;
        jj+=2;
      }
      ii+=2;
    }
  }
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  return K;
}
const Matrix& Triangle6::get_M() {
  Matrix &M=*myMatrix;
  M.Clear();
  return M;
}
const Vector& Triangle6::get_R() {
  Vector& R=*myVector;
  R.clear();
  // Quick return if incactive
  if (!(groupdata_->active)) {
    return R;
  }
  // Get group factors
  double facS = groupdata_->factor_S;
  double facG = groupdata_->factor_G;
  double facP = groupdata_->factor_P;
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    // Get B-matrix
    static Matrix N(6, 3);
    static double detJ;
    this->get_shape_functions(myMatPoints[k], N, detJ);
    double dV = detJ*(pD->get_fac())*(myMatPoints[k]->get_w());
    // Get stress vector
    static Vector sigma(6);
    sigma = myMatPoints[k]->get_material()->get_stress();
    int ii = 0;
    for (int i = 0; i < 6; i++) {
      R[ii  ]+=facS*(N(i, 1)*sigma[0]+N(i, 2)*sigma[3])*dV;
      R[ii+1]+=facS*(N(i, 2)*sigma[1]+N(i, 1)*sigma[3])*dV;
      R[ii  ]-=facG*(N(i, 0)*b[0]*dV);
      R[ii+1]-=facG*(N(i, 0)*b[1]*dV);
      ii+=2;
    }
  }
  R-=facP*P;
  return R;
}
void Triangle6::update() {
  // Quick return if inactive
  if (!(groupdata_->active)) {
    return;
  }
  // Get incremental displacement vector
  Vector& u=*myVector;
  u = this->get_disp_incrm();
  // For each material point
  for (unsigned int i = 0; i < myMatPoints.size(); i++) {
    // Get B-matrix
    static Matrix N(6, 3);
    static double detJ;
    this->get_shape_functions(myMatPoints[i], N, detJ);
    // Determine the strain
    static Vector epsilon(6);
      epsilon.clear();
    epsilon[0] = N(0, 1)*u[ 0]+N(1, 1)*u[ 2]+N(2, 1)*u[ 4]
                +N(3, 1)*u[ 6]+N(4, 1)*u[ 8]+N(5, 1)*u[10];
    epsilon[1] = N(0, 2)*u[ 1]+N(1, 2)*u[ 3]+N(2, 2)*u[ 5]
                +N(3, 2)*u[ 7]+N(4, 2)*u[ 9]+N(5, 2)*u[11];
    epsilon[3] = N(0, 2)*u[ 0]+N(0, 1)*u[ 1]+N(1, 2)*u[ 2]
                +N(1, 1)*u[ 3]+N(2, 2)*u[ 4]+N(2, 1)*u[ 5]
                +N(3, 2)*u[ 6]+N(3, 1)*u[ 7]+N(4, 2)*u[ 8]
                +N(4, 1)*u[ 9]+N(5, 2)*u[10]+N(5, 1)*u[11];
    // And send it to the material point
    myMatPoints[i]->get_material()->set_strain(epsilon);
  }
}
void Triangle6::commit() {
  for (unsigned int i = 0;i < myMatPoints.size();i++)
    myMatPoints[i]->get_material()->commit();
}
bool Triangle6::checkIfAllows(FEObject* /*f*/) {
  return true;
}
void Triangle6::recoverStresses() {
  /// @todo check Gauss extrapolation
  static Vector sigma(6);
  static Matrix xi(6, 3);
  const double d1 = 1.0+num::d23;
  const double d2=    num::d13;
  const double d3=    num::d23;

  xi(0, 0) =  d1;
  xi(0, 1) = -d2;
  xi(0, 2) = -d2;
  xi(1, 0) = -d2;
  xi(1, 1) =  d1;
  xi(1, 2) = -d2;
  xi(2, 0) = -d2;
  xi(2, 1) = -d2;
  xi(2, 2) =  d1;
  xi(3, 0) =  d3;
  xi(3, 1) =  d3;
  xi(3, 2) = -d2;
  xi(4, 0) = -d2;
  xi(4, 1) =  d3;
  xi(4, 2) =  d3;
  xi(5, 0) =  d3;
  xi(5, 1) = -d2;
  xi(5, 2) =  d3;

  for (unsigned i = 0; i < 6; i++) {
    sigma.clear();
    for (unsigned j = 0; j < 3; j++) {
      sigma+=xi(i, j)*(myMatPoints[j]->get_material()->get_stress());
    }
    nodes_[i]->addStress(sigma);
  }
}

/**
 * Initial stresses.
 */
void Triangle6::AddInitialStresses(int direction,
                                 double h1, double s1,
                                 double h2, double s2, double K0) {
  for (unsigned i = 0;i < myMatPoints.size();i++) {
    myMatPoints[i]->AddInitialStresses(direction, h1, s1, h2, s2, K0);
  }
}

int Triangle6::get_num_plastic_points() {
  int n = 0;
  for (unsigned int i = 0;i < myMatPoints.size();i++)
    if (myMatPoints[i]->get_material()->isPlastic()) n++;
  return n;
}
