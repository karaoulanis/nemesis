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

#include "elements/triangle3.h"
#include "group/group_data.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"
#include "node/node.h"

Triangle3::Triangle3()
    : thickness_(0),
      myMatPoints(0),
      a1(0.),
      a2(0.),
      a3(0.),
      b1(0.),
      b2(0.),
      b3(0.),
      c1(0.),
      c2(0.),
      c3(0.),
      A(0.) {
}

Triangle3::Triangle3(int id, std::vector<Node*> nodes,
                     MultiaxialMaterial* material, double thickness)
    : Element(id, nodes),
      thickness_(thickness) {
  // Get nodal data
  myNodalIDs.resize(3);
  myNodalIDs[0] = nodes_[0]->get_id();
  myNodalIDs[1] = nodes_[1]->get_id();
  myNodalIDs[2] = nodes_[2]->get_id();
  // Set local nodal dofs
  myLocalNodalDofs.resize(2);
  myLocalNodalDofs[0]=0;
  myLocalNodalDofs[1]=1;
  // Handle common info: Start -------------------------------------------------
  // Find own matrix and vector
  myMatrix = theStaticMatrices[6];
  myVector = theStaticVectors[6];
  // Get nodal coordinates
  /// @todo replace with const references
  x.Resize(3, 3);
  for (int i = 0; i < 3; i++) {
    x(i, 0) = nodes_[i]->get_x1();
    x(i, 1) = nodes_[i]->get_x2();
    x(i, 2) = nodes_[i]->get_x3();
  }
  // Inform the nodes that the corresponding Dof's must be activated
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 2; j++)
      nodes_[i]->addDofToNode(myLocalNodalDofs[j]);
  // Load vector
  P.Resize(6, 0.);
  // Self weight
  G.Resize(6, 0.);
  myMaterial = material;
  this->AssignGravityLoads();
  // Handle common info: End ---------------------------------------------------

  // Find geometrical relations
  a1 = x(1, 0)*x(2, 1)-x(2, 0)*x(1, 1);
  a2 = x(2, 0)*x(0, 1)-x(0, 0)*x(2, 1);
  a3 = x(0, 0)*x(1, 1)-x(1, 0)*x(0, 1);
  b1 = x(1, 1)-x(2, 1);
  b2 = x(2, 1)-x(0, 1);
  b3 = x(0, 1)-x(1, 1);
  c1 = x(2, 0)-x(1, 0);
  c2 = x(0, 0)-x(2, 0);
  c3 = x(1, 0)-x(0, 0);
  A = 0.5*(x(1, 0)*x(2, 1)-x(2, 0)*x(1, 1)
          +x(2, 0)*x(0, 1)-x(0, 0)*x(2, 1)
          +x(0, 0)*x(1, 1)-x(1, 0)*x(0, 1));

  // Material
  myMatPoints.resize(1);
  MultiaxialMaterial* pMat = static_cast<MultiaxialMaterial*>(myMaterial);
  myMatPoints[0]=new MatPoint(pMat, 1, 1, 1, 1);
  myMatPoints[0]->set_X(num::d13*(x(0, 0)+x(1, 0)+x(2, 0)),
                       num::d13*(x(0, 1)+x(1, 1)+x(2, 1)));
}
Triangle3::~Triangle3() {
}

/**
 * Initial stresses.
 */
void Triangle3::AddInitialStresses(int direction,
                                 double h1, double s1,
                                 double h2, double s2, double K0) {
  for (unsigned i = 0; i < myMatPoints.size(); i++) {
    myMatPoints[i]->AddInitialStresses(direction, h1, s1, h2, s2, K0);
  }
}

const Matrix& Triangle3::get_K() {
  Matrix& K=*myMatrix;
  const Matrix& C = myMatPoints[0]->get_material()->get_C();
  double coeff = thickness_*0.25/A;
  K(0, 0) = coeff*((b1*C(0, 0)+c1*C(3, 0))*b1+(b1*C(0, 3)+c1*C(3, 3))*c1);
  K(0, 1) = coeff*((b1*C(0, 1)+c1*C(3, 1))*c1+(b1*C(0, 3)+c1*C(3, 3))*b1);
  K(0, 2) = coeff*((b1*C(0, 0)+c1*C(3, 0))*b2+(b1*C(0, 3)+c1*C(3, 3))*c2);
  K(0, 3) = coeff*((b1*C(0, 1)+c1*C(3, 1))*c2+(b1*C(0, 3)+c1*C(3, 3))*b2);
  K(0, 4) = coeff*((b1*C(0, 0)+c1*C(3, 0))*b3+(b1*C(0, 3)+c1*C(3, 3))*c3);
  K(0, 5) = coeff*((b1*C(0, 1)+c1*C(3, 1))*c3+(b1*C(0, 3)+c1*C(3, 3))*b3);
  K(1, 0) = coeff*((c1*C(1, 0)+b1*C(3, 0))*b1+(c1*C(1, 3)+b1*C(3, 3))*c1);
  K(1, 1) = coeff*((c1*C(1, 1)+b1*C(3, 1))*c1+(c1*C(1, 3)+b1*C(3, 3))*b1);
  K(1, 2) = coeff*((c1*C(1, 0)+b1*C(3, 0))*b2+(c1*C(1, 3)+b1*C(3, 3))*c2);
  K(1, 3) = coeff*((c1*C(1, 1)+b1*C(3, 1))*c2+(c1*C(1, 3)+b1*C(3, 3))*b2);
  K(1, 4) = coeff*((c1*C(1, 0)+b1*C(3, 0))*b3+(c1*C(1, 3)+b1*C(3, 3))*c3);
  K(1, 5) = coeff*((c1*C(1, 1)+b1*C(3, 1))*c3+(c1*C(1, 3)+b1*C(3, 3))*b3);
  K(2, 0) = coeff*((b2*C(0, 0)+c2*C(3, 0))*b1+(b2*C(0, 3)+c2*C(3, 3))*c1);
  K(2, 1) = coeff*((b2*C(0, 1)+c2*C(3, 1))*c1+(b2*C(0, 3)+c2*C(3, 3))*b1);
  K(2, 2) = coeff*((b2*C(0, 0)+c2*C(3, 0))*b2+(b2*C(0, 3)+c2*C(3, 3))*c2);
  K(2, 3) = coeff*((b2*C(0, 1)+c2*C(3, 1))*c2+(b2*C(0, 3)+c2*C(3, 3))*b2);
  K(2, 4) = coeff*((b2*C(0, 0)+c2*C(3, 0))*b3+(b2*C(0, 3)+c2*C(3, 3))*c3);
  K(2, 5) = coeff*((b2*C(0, 1)+c2*C(3, 1))*c3+(b2*C(0, 3)+c2*C(3, 3))*b3);
  K(3, 0) = coeff*((c2*C(1, 0)+b2*C(3, 0))*b1+(c2*C(1, 3)+b2*C(3, 3))*c1);
  K(3, 1) = coeff*((c2*C(1, 1)+b2*C(3, 1))*c1+(c2*C(1, 3)+b2*C(3, 3))*b1);
  K(3, 2) = coeff*((c2*C(1, 0)+b2*C(3, 0))*b2+(c2*C(1, 3)+b2*C(3, 3))*c2);
  K(3, 3) = coeff*((c2*C(1, 1)+b2*C(3, 1))*c2+(c2*C(1, 3)+b2*C(3, 3))*b2);
  K(3, 4) = coeff*((c2*C(1, 0)+b2*C(3, 0))*b3+(c2*C(1, 3)+b2*C(3, 3))*c3);
  K(3, 5) = coeff*((c2*C(1, 1)+b2*C(3, 1))*c3+(c2*C(1, 3)+b2*C(3, 3))*b3);
  K(4, 0) = coeff*((b3*C(0, 0)+c3*C(3, 0))*b1+(b3*C(0, 3)+c3*C(3, 3))*c1);
  K(4, 1) = coeff*((b3*C(0, 1)+c3*C(3, 1))*c1+(b3*C(0, 3)+c3*C(3, 3))*b1);
  K(4, 2) = coeff*((b3*C(0, 0)+c3*C(3, 0))*b2+(b3*C(0, 3)+c3*C(3, 3))*c2);
  K(4, 3) = coeff*((b3*C(0, 1)+c3*C(3, 1))*c2+(b3*C(0, 3)+c3*C(3, 3))*b2);
  K(4, 4) = coeff*((b3*C(0, 0)+c3*C(3, 0))*b3+(b3*C(0, 3)+c3*C(3, 3))*c3);
  K(4, 5) = coeff*((b3*C(0, 1)+c3*C(3, 1))*c3+(b3*C(0, 3)+c3*C(3, 3))*b3);
  K(5, 0) = coeff*((c3*C(1, 0)+b3*C(3, 0))*b1+(c3*C(1, 3)+b3*C(3, 3))*c1);
  K(5, 1) = coeff*((c3*C(1, 1)+b3*C(3, 1))*c1+(c3*C(1, 3)+b3*C(3, 3))*b1);
  K(5, 2) = coeff*((c3*C(1, 0)+b3*C(3, 0))*b2+(c3*C(1, 3)+b3*C(3, 3))*c2);
  K(5, 3) = coeff*((c3*C(1, 1)+b3*C(3, 1))*c2+(c3*C(1, 3)+b3*C(3, 3))*b2);
  K(5, 4) = coeff*((c3*C(1, 0)+b3*C(3, 0))*b3+(c3*C(1, 3)+b3*C(3, 3))*c3);
  K(5, 5) = coeff*((c3*C(1, 1)+b3*C(3, 1))*c3+(c3*C(1, 3)+b3*C(3, 3))*b3);
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  return K;
}
const Matrix& Triangle3::get_M() {
  return *myMatrix;
}
const Vector& Triangle3::get_R() {
  Vector& R=*myVector;
  R.Clear();
  // Quick return if inactive
  if (!(groupdata_->active)) {
    return R;
  }
  // Get factors
  double facS = groupdata_->factor_S;
  double facG = groupdata_->factor_G;
  double facP = groupdata_->factor_P;
  double s1=(myMatPoints[0]->get_material()->get_stress())[0];
  double s2=(myMatPoints[0]->get_material()->get_stress())[1];
  double s3=(myMatPoints[0]->get_material()->get_stress())[3];
  double fac = facS*0.5*thickness_;
  double facb0 = facG*thickness_*A*num::d13*b[0];
  double facb1 = facG*thickness_*A*num::d13*b[1];
  R[0]=fac*(b1*s1+c1*s3)-facb0;
  R[1]=fac*(c1*s2+b1*s3)-facb1;
  R[2]=fac*(b2*s1+c2*s3)-facb0;
  R[3]=fac*(c2*s2+b2*s3)-facb1;
  R[4]=fac*(b3*s1+c3*s3)-facb0;
  R[5]=fac*(c3*s2+b3*s3)-facb1;
  R-=facP*P;
  return R;
}

void Triangle3::Update(const double Dt) {
  // Quick return if inactive
  if (!(groupdata_->active)) {
    return;
  }
  Vector& u=*myVector;
  u = this->get_disp_incrm();
  // Determine the strain
  static Vector epsilon(6);
  epsilon.Clear();
  epsilon[0]=(0.5/A)*(b1*u[0]+b2*u[2]+b3*u[4]);
  epsilon[1]=(0.5/A)*(c1*u[1]+c2*u[3]+c3*u[5]);
  epsilon[3]=(0.5/A)*(c1*u[0]+b1*u[1]+c2*u[2]+b2*u[3]+c3*u[4]+b3*u[5]);
  // And send it to the material point
  myMatPoints[0]->get_material()->set_strain(epsilon);
}


void Triangle3::Commit() {
  myMatPoints[0]->get_material()->Commit();
}


void Triangle3::recoverStresses() {
  static Vector sigma(6);
  sigma = myMatPoints[0]->get_material()->get_stress();
  for (int i = 0; i < 3; i++) nodes_[i]->addStress(sigma);
}


int Triangle3::get_num_plastic_points() {
  int n = 0;
  if (myMatPoints[0]->get_material()->isPlastic()) n = 1;
  return n;
}
