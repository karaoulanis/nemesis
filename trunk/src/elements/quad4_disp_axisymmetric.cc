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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.       *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "elements/quad4_disp_axisymmetric.h"
#include "loadcase/group_data.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"

Quad4DispAxisymmetric::Quad4DispAxisymmetric() {
}
Quad4DispAxisymmetric::Quad4DispAxisymmetric(int ID, int Node_1, int Node_2,
                                             int Node_3, int Node_4, int MatID,
                                             int integrationRuleXi,
                                             int integrationRuleEta)
:Quad4(ID, Node_1, Node_2, Node_3, Node_4, MatID, integrationRuleXi, integrationRuleEta) {
}
Quad4DispAxisymmetric::~Quad4DispAxisymmetric() {
}
const Matrix& Quad4DispAxisymmetric::get_K() {
  Matrix &K=*myMatrix;
  K.clear();
  for (unsigned int k = 0; k < myMatPoints.size(); k++) {
    this->findShapeFunctionsAt(myMatPoints[k]);
    double r = x(0, 0)*N(0, 0)+x(1, 0)*N(0, 1)+x(2, 0)*N(0, 2)+x(3, 0)*N(0, 3);
    double dV = detJ*r*(pD->get_fac())*(myMatPoints[k]->get_w());
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    int ii = 0;
    for (int i = 0; i < 4; i++) {
      int jj = 0;
      for (int j = 0; j < 4; j++) {
        static Matrix CB(4, 2);

        CB(0, 0)=C(0, 0)*N(1, j) + C(0, 2)*N(0, j)/r + C(0, 3)*N(2, j);
        CB(1, 0)=C(1, 0)*N(1, j) + C(1, 2)*N(0, j)/r + C(1, 3)*N(2, j);
        CB(2, 0)=C(2, 0)*N(1, j) + C(2, 2)*N(0, j)/r + C(2, 3)*N(2, j);
        CB(3, 0)=C(3, 0)*N(1, j) + C(3, 2)*N(0, j)/r + C(3, 3)*N(2, j);

        CB(0, 1)=C(0, 1)*N(2, j)                   + C(0, 3)*N(1, j);
        CB(1, 1)=C(1, 1)*N(2, j)                   + C(1, 3)*N(1, j);
        CB(2, 1)=C(2, 1)*N(2, j)                   + C(2, 3)*N(1, j);
        CB(3, 1)=C(3, 1)*N(2, j)                   + C(3, 3)*N(1, j);

        K(ii  , jj  ) += (N(1, i)*CB(0, 0) + N(0, i)/r*CB(2, 0) + N(2, i)*CB(3, 0))*dV;
        K(ii  , jj+1) += (N(1, i)*CB(0, 1) + N(0, i)/r*CB(2, 1) + N(2, i)*CB(3, 1))*dV;
        K(ii+1, jj  ) += (N(2, i)*CB(1, 0)                      + N(1, i)*CB(3, 0))*dV;
        K(ii+1, jj+1) += (N(2, i)*CB(1, 1)                      + N(1, i)*CB(3, 1))*dV;
        jj+=2;
      }
      ii+=2;
    }
  }
  double facK = groupdata_->active_ ? groupdata_->factor_K_: 1e-7;
  K*=facK;
  return K;
}
/// @todo
const Matrix& Quad4DispAxisymmetric::get_M() {
  Matrix &M=*myMatrix;
  M.clear();
  double rho = myMaterial->get_rho();
  double volume = 0.;
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    this->findShapeFunctionsAt(myMatPoints[k]);
    volume+=detJ*(pD->get_fac())*(myMatPoints[k]->get_w());
  }
  double mass = rho*volume;
  for (int i = 0;i < 8;i++) M(i, i)=0.25*mass;
  return M;
}

const Vector& Quad4DispAxisymmetric::get_R() {
  static Vector sigma(6);
  Vector& R=*myVector;
  R.clear();
  // Factors
  if (!(groupdata_->active_))  return R;
  double facS = groupdata_->factor_S_;
  double facG = groupdata_->factor_G_;
  double facP = groupdata_->factor_P_;
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    sigma = myMatPoints[k]->get_material()->get_stress();
    this->findShapeFunctionsAt(myMatPoints[k]);
    double r = x(0, 0)*N(0, 0)+x(1, 0)*N(0, 1)+x(2, 0)*N(0, 2)+x(3, 0)*N(0, 3);
    double dV = detJ*r*(pD->get_fac())*(myMatPoints[k]->get_w());
    int ii = 0;
    for (int i = 0; i < 4; i++) {
      R[ii  ]+=facS*(N(1, i)*sigma[0] + N(0, i)/r*sigma[2]  + N(2, i)*sigma[3])*dV;
      R[ii+1]+=facS*(N(2, i)*sigma[1]                       + N(1, i)*sigma[3])*dV;
      R[ii  ]-=facG*(N(0, i)*b[0]*dV);
      R[ii+1]-=facG*(N(0, i)*b[1]*dV);
      ii+=2;
    }
  }
  R-=facP*P;
  return R;
}
void Quad4DispAxisymmetric::update() {
  if (!(groupdata_->active_))  return;
  Vector& u=*myVector;
  static Vector u1(8);
  static Vector u2(8);
  u2 = this->get_disp_convg();
  u1 = this->get_disp_trial();
  u = u1-u2;
  // For each material point
  for (unsigned int i = 0; i < myMatPoints.size(); i++) {
    this->findShapeFunctionsAt(myMatPoints[i]);
    double r = x(0, 0)*N(0, 0)+x(1, 0)*N(0, 1)+x(2, 0)*N(0, 2)+x(3, 0)*N(0, 3);
    // Determine the strain
    static Vector epsilon(6);
      epsilon.clear();
    epsilon[0]=N(1, 0)*u[0]   + N(1, 1)*u[2]   + N(1, 2)*u[4]   + N(1, 3)*u[6];
    epsilon[1]=N(2, 0)*u[1]   + N(2, 1)*u[3]   + N(2, 2)*u[5]   + N(2, 3)*u[7];
    epsilon[2]=N(0, 0)/r*u[0] + N(0, 1)/r*u[2] + N(0, 2)/r*u[4] + N(0, 3)/r*u[6];
    epsilon[3]=N(2, 0)*u[0]   + N(1, 0)*u[1]   + N(2, 1)*u[2]   + N(1, 1)*u[3] +
           N(2, 2)*u[4]   + N(1, 2)*u[5]   + N(2, 3)*u[6]   + N(1, 3)*u[7];
    // And send it to the material point
    myMatPoints[i]->get_material()->set_strain(epsilon);
  }
}
