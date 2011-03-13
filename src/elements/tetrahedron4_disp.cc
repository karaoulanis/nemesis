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

#include "elements/tetrahedron4_disp.h"
#include "loadcase/group_data.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"
#include "node/node.h"

Matrix Tetrahedron4Disp::N(4, 4);
double Tetrahedron4Disp::V = 0;

Tetrahedron4Disp::Tetrahedron4Disp() {
}
Tetrahedron4Disp::Tetrahedron4Disp(int ID,
             int Node_1, int Node_2, int Node_3, int Node_4,
             int matID)
:Element(ID, matID) {
  myTag = TAG_ELEM_TETRHEDRON_4_DISP;
  // Get nodal data
  myNodalIDs.resize(4);
  myNodalIDs[0]=Node_1;
  myNodalIDs[1]=Node_2;
  myNodalIDs[2]=Node_3;
  myNodalIDs[3]=Node_4;
  // Set local nodal dofs
  myLocalNodalDofs.resize(3);
  myLocalNodalDofs[0]=0;
  myLocalNodalDofs[1]=1;
  myLocalNodalDofs[2]=2;
  // Handle common info
  this->handleCommonInfo();

  // Materials
  myMatPoints.resize(1);
  MultiaxialMaterial* pMat = static_cast < MultiaxialMaterial*>(myMaterial);
  myMatPoints[0]=new MatPoint(pMat, 1, 1, 1, 1, 1, 1);
}
Tetrahedron4Disp::~Tetrahedron4Disp() {
  Containers::vector_delete(myMatPoints);
}
void Tetrahedron4Disp::findShapeFunctions() {
  V = 1./6.*((x(1, 0)-x(0, 0))*((x(2, 1)-x(0, 1))*(x(3, 2)-x(0, 2))-(x(3, 1)-x(0, 1))*(x(2, 2)-x(0, 2)))
        +(x(1, 1)-x(0, 1))*((x(3, 0)-x(0, 0))*(x(2, 2)-x(0, 2))-(x(2, 0)-x(0, 0))*(x(3, 2)-x(0, 2)))
      +(x(1, 2)-x(0, 2))*((x(2, 0)-x(0, 0))*(x(3, 1)-x(0, 1))-(x(3, 0)-x(0, 0))*(x(2, 1)-x(0, 1))));
  double V6 = 1/(6.*V);
  N(0, 0)=0.;
  N(1, 0)=0.;
  N(2, 0)=0.;
  N(3, 0)=0.;
  N(0, 1)=( x(1, 1)*(x(3, 2)-x(2, 2))-x(2, 1)*(x(3, 2)-x(1, 2))+x(3, 1)*(x(2, 2)-x(1, 2)))*V6;
  N(1, 1)=(-x(0, 1)*(x(3, 2)-x(2, 2))+x(2, 1)*(x(3, 2)-x(0, 2))-x(3, 1)*(x(2, 2)-x(0, 2)))*V6;
  N(2, 1)=( x(0, 1)*(x(3, 2)-x(1, 2))-x(1, 1)*(x(3, 2)-x(0, 2))+x(3, 1)*(x(1, 2)-x(0, 2)))*V6;
  N(3, 1)=(-x(0, 1)*(x(2, 2)-x(1, 2))+x(1, 1)*(x(2, 2)-x(0, 2))-x(2, 1)*(x(1, 2)-x(0, 2)))*V6;
  N(0, 2)=(-x(1, 0)*(x(3, 2)-x(2, 2))+x(2, 0)*(x(3, 2)-x(1, 2))-x(3, 0)*(x(2, 2)-x(1, 2)))*V6;
  N(1, 2)=( x(0, 0)*(x(3, 2)-x(2, 2))-x(2, 0)*(x(3, 2)-x(0, 2))+x(3, 0)*(x(2, 2)-x(0, 2)))*V6;
  N(2, 2)=(-x(0, 0)*(x(3, 2)-x(1, 2))+x(1, 0)*(x(3, 2)-x(0, 2))-x(3, 0)*(x(1, 2)-x(0, 2)))*V6;
  N(3, 2)=( x(0, 0)*(x(2, 2)-x(1, 2))-x(1, 0)*(x(2, 2)-x(0, 2))+x(2, 0)*(x(1, 2)-x(0, 2)))*V6;
  N(0, 3)=( x(1, 0)*(x(3, 1)-x(2, 1))-x(2, 0)*(x(3, 1)-x(1, 1))+x(3, 0)*(x(2, 1)-x(1, 1)))*V6;
  N(1, 3)=(-x(0, 0)*(x(3, 1)-x(2, 1))+x(2, 0)*(x(3, 1)-x(0, 1))-x(3, 0)*(x(2, 1)-x(0, 1)))*V6;
  N(2, 3)=( x(0, 0)*(x(3, 1)-x(1, 1))-x(1, 0)*(x(3, 1)-x(0, 1))+x(3, 0)*(x(1, 1)-x(0, 1)))*V6;
  N(3, 3)=(-x(0, 0)*(x(2, 1)-x(1, 1))+x(1, 0)*(x(2, 1)-x(0, 1))-x(2, 0)*(x(1, 1)-x(0, 1)))*V6;
}
const Matrix& Tetrahedron4Disp::get_K() {
  Matrix &K=*myMatrix;
  K.clear();
  this->findShapeFunctions();
  const Matrix& C = myMatPoints[0]->get_material()->get_C();
  static Matrix CB(6, 3);
  int ii = 0;
  for (int i = 0; i < 4; i++) {
    int jj = 0;
    for (int j = 0; j < 4; j++) {
      CB(0, 0)=C(0, 0)*N(j, 1)+C(0, 3)*N(j, 2)+C(0, 5)*N(j, 3);
      CB(1, 0)=C(1, 0)*N(j, 1)+C(1, 3)*N(j, 2)+C(1, 5)*N(j, 3);
      CB(2, 0)=C(2, 0)*N(j, 1)+C(2, 3)*N(j, 2)+C(2, 5)*N(j, 3);
      CB(3, 0)=C(3, 0)*N(j, 1)+C(3, 3)*N(j, 2)+C(3, 5)*N(j, 3);
      CB(4, 0)=C(4, 0)*N(j, 1)+C(4, 3)*N(j, 2)+C(4, 5)*N(j, 3);
      CB(5, 0)=C(5, 0)*N(j, 1)+C(5, 3)*N(j, 2)+C(5, 5)*N(j, 3);
      CB(0, 1)=C(0, 1)*N(j, 2)+C(0, 3)*N(j, 1)+C(0, 4)*N(j, 3);
      CB(1, 1)=C(1, 1)*N(j, 2)+C(1, 3)*N(j, 1)+C(1, 4)*N(j, 3);
      CB(2, 1)=C(2, 1)*N(j, 2)+C(2, 3)*N(j, 1)+C(2, 4)*N(j, 3);
      CB(3, 1)=C(3, 1)*N(j, 2)+C(3, 3)*N(j, 1)+C(3, 4)*N(j, 3);
      CB(4, 1)=C(4, 1)*N(j, 2)+C(4, 3)*N(j, 1)+C(4, 4)*N(j, 3);
      CB(5, 1)=C(5, 1)*N(j, 2)+C(5, 3)*N(j, 1)+C(5, 4)*N(j, 3);
      CB(0, 2)=C(0, 2)*N(j, 3)+C(0, 4)*N(j, 2)+C(0, 5)*N(j, 1);
      CB(1, 2)=C(1, 2)*N(j, 3)+C(1, 4)*N(j, 2)+C(1, 5)*N(j, 1);
      CB(2, 2)=C(2, 2)*N(j, 3)+C(2, 4)*N(j, 2)+C(2, 5)*N(j, 1);
      CB(3, 2)=C(3, 2)*N(j, 3)+C(3, 4)*N(j, 2)+C(3, 5)*N(j, 1);
      CB(4, 2)=C(4, 2)*N(j, 3)+C(4, 4)*N(j, 2)+C(4, 5)*N(j, 1);
      CB(5, 2)=C(5, 2)*N(j, 3)+C(5, 4)*N(j, 2)+C(5, 5)*N(j, 1);
      K(ii  , jj  )+=(N(i, 1)*CB(0, 0)+N(i, 2)*CB(3, 0)+N(i, 3)*CB(5, 0))*V;
      K(ii  , jj+1)+=(N(i, 1)*CB(0, 1)+N(i, 2)*CB(3, 1)+N(i, 3)*CB(5, 1))*V;
      K(ii  , jj+2)+=(N(i, 1)*CB(0, 2)+N(i, 2)*CB(3, 2)+N(i, 3)*CB(5, 2))*V;
      K(ii+1, jj  )+=(N(i, 2)*CB(1, 0)+N(i, 1)*CB(3, 0)+N(i, 3)*CB(4, 0))*V;
      K(ii+1, jj+1)+=(N(i, 2)*CB(1, 1)+N(i, 1)*CB(3, 1)+N(i, 3)*CB(4, 1))*V;
      K(ii+1, jj+2)+=(N(i, 2)*CB(1, 2)+N(i, 1)*CB(3, 2)+N(i, 3)*CB(4, 2))*V;
      K(ii+2, jj  )+=(N(i, 3)*CB(2, 0)+N(i, 2)*CB(4, 0)+N(i, 1)*CB(5, 0))*V;
      K(ii+2, jj+1)+=(N(i, 3)*CB(2, 1)+N(i, 2)*CB(4, 1)+N(i, 1)*CB(5, 1))*V;
      K(ii+2, jj+2)+=(N(i, 3)*CB(2, 2)+N(i, 2)*CB(4, 2)+N(i, 1)*CB(5, 2))*V;
      jj+=3;
    }
    ii+=3;
  }
  double facK = groupdata_->active_ ? groupdata_->factor_K_: 1e-7;
  K*=facK;
  return K;
}
const Matrix& Tetrahedron4Disp::get_M() {
  Matrix &M=*myMatrix;
  M.clear();
  return M;
}
const Vector& Tetrahedron4Disp::get_R() {
  static Vector sigma(6);
  Vector& R=*myVector;
  R.clear();

  if (!(groupdata_->active_))  return R;
  double facS = groupdata_->factor_S_;
  // double facG = groupdata_->factor_G_;
  double facP = groupdata_->factor_P_;

  sigma = myMatPoints[0]->get_material()->get_stress();
  this->findShapeFunctions();
  int ii = 0;
  for (int i = 0; i < 4; i++) {
    R[ii  ]+=facS*(N(i, 1)*sigma[0]+N(i, 2)*sigma[3]+N(i, 3)*sigma[5])*V;
    R[ii+1]+=facS*(N(i, 2)*sigma[1]+N(i, 1)*sigma[3]+N(i, 3)*sigma[4])*V;
    R[ii+2]+=facS*(N(i, 3)*sigma[2]+N(i, 2)*sigma[4]+N(i, 1)*sigma[5])*V;
    /// @todo check
    // R[ii  ]-=facG*(N(0, i)*b[0]*dV);
    // R[ii+1]-=facG*(N(0, i)*b[1]*dV);
    // R[ii+1]-=facG*(N(0, i)*b[2]*dV);
    ii+=3;
  }

  R-=facP*P;
  return R;
}
void Tetrahedron4Disp::update() {
  if (!(groupdata_->active_))  return;
  static Vector u(12);
  u = this->get_disp_incrm();
  // For each material point
  this->findShapeFunctions();
  // Determine the strain
  static Vector epsilon(6);
  epsilon[0]=N(0, 1)*u[0]+N(1, 1)*u[3]+N(2, 1)*u[6]+N(3, 1)*u[ 9];
  epsilon[1]=N(0, 2)*u[1]+N(1, 2)*u[4]+N(2, 2)*u[7]+N(3, 2)*u[10];
  epsilon[2]=N(0, 3)*u[2]+N(1, 3)*u[5]+N(2, 3)*u[8]+N(3, 3)*u[11];
  epsilon[3]=N(0, 2)*u[0]+N(0, 1)*u[1]+N(1, 2)*u[3]+N(1, 1)*u[4]+
        N(2, 2)*u[6]+N(2, 1)*u[7]+N(3, 2)*u[9]+N(3, 1)*u[10];
  epsilon[4]=N(0, 3)*u[1]+N(0, 2)*u[2]+N(1, 3)*u[4]+N(1, 2)*u[5]+
        N(2, 3)*u[7]+N(2, 2)*u[8]+N(3, 3)*u[10]+N(3, 2)*u[11];
  epsilon[5]=N(0, 3)*u[0]+N(0, 1)*u[2]+N(1, 3)*u[3]+N(1, 1)*u[5]+
        N(2, 3)*u[6]+N(2, 1)*u[8]+N(3, 3)*u[9]+N(3, 1)*u[11];
  // And send it to the material point
  myMatPoints[0]->get_material()->set_strain(epsilon);
}
void Tetrahedron4Disp::commit() {
  for (unsigned int i = 0;i < myMatPoints.size();i++)
    myMatPoints[i]->get_material()->commit();
}
bool Tetrahedron4Disp::checkIfAllows(FEObject* /*f*/) {
  return true;
}
void Tetrahedron4Disp::recoverStresses() {
  static Vector sigma(6);
  sigma = myMatPoints[0]->get_material()->get_stress();
  for (int i = 0;i < 4;i++) myNodes[i]->addStress(sigma);
}
