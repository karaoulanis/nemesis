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

#include "elements/brick8.h"
#include "elements/shape_functions.h"
#include "group/group_data.h"
#include "main/nemesis_debug.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"
#include "node/node.h"

double Brick8::detJ[8];
double Brick8::shp[8][4][8];
std::vector<int> Brick8::perm(8);

/**
 * Default constructor.
 */
Brick8::Brick8()
    : myMatPoints(0) {
}

/**
 * Constructor.
 */
Brick8::Brick8(int id, std::vector<Node*> nodes,
                     MultiaxialMaterial* material)
    : Element(id, nodes),
      myMatPoints(8) {
  // Get nodal data
  myNodalIDs.resize(8);
  myNodalIDs[0] = nodes_[0]->get_id();
  myNodalIDs[1] = nodes_[1]->get_id();
  myNodalIDs[2] = nodes_[2]->get_id();
  myNodalIDs[3] = nodes_[3]->get_id();
  myNodalIDs[4] = nodes_[4]->get_id();
  myNodalIDs[5] = nodes_[5]->get_id();
  myNodalIDs[6] = nodes_[6]->get_id();
  myNodalIDs[7] = nodes_[7]->get_id();

  // Set local nodal dofs
  myLocalNodalDofs.resize(3);
  myLocalNodalDofs[0]=0;
  myLocalNodalDofs[1]=1;
  myLocalNodalDofs[2]=2;

  // Handle common info: Start -------------------------------------------------
  // Find own matrix and vector
  myMatrix = theStaticMatrices[24];
  myVector = theStaticVectors[24];
  // Get nodal coordinates
  /// @todo replace with const references
  x.Resize(8, 3);
  for (int i = 0; i < 8; i++) {
    x(i, 0) = nodes_[i]->get_x1();
    x(i, 1) = nodes_[i]->get_x2();
    x(i, 2) = nodes_[i]->get_x3();
  }
  // Inform the nodes that the corresponding Dof's must be activated
  for (int i = 0; i < 8; i++)
    for (int j = 0; j < 2; j++)
      nodes_[i]->addDofToNode(myLocalNodalDofs[j]);
  // Load vector
  P.Resize(24, 0.);
  // Self weight
  G.Resize(24, 0.);
  myMaterial = material;
  this->AssignGravityLoads();
  // Handle common info: End ---------------------------------------------------

  // Materials
  MultiaxialMaterial* pMat = static_cast<MultiaxialMaterial*>(myMaterial);
  myMatPoints[0]=new MatPoint(pMat, 1, 1, 1, 2, 2, 2);
  myMatPoints[1]=new MatPoint(pMat, 2, 1, 1, 2, 2, 2);
  myMatPoints[2]=new MatPoint(pMat, 2, 2, 1, 2, 2, 2);
  myMatPoints[3]=new MatPoint(pMat, 1, 2, 1, 2, 2, 2);
  myMatPoints[4]=new MatPoint(pMat, 1, 1, 2, 2, 2, 2);
  myMatPoints[5]=new MatPoint(pMat, 2, 1, 2, 2, 2, 2);
  myMatPoints[6]=new MatPoint(pMat, 2, 2, 2, 2, 2, 2);
  myMatPoints[7]=new MatPoint(pMat, 1, 2, 2, 2, 2, 2);

  // Find shape functions for all GaussPoints (double shp[node][N, i][GPoint])
  this->shapeFunctions();

  // Materials coordinates
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    double xG = 0, yG = 0, zG = 0;
    for (unsigned a = 0; a < nodes_.size(); a++) {
      xG+=shp[a][0][k]*x(a, 0);
      yG+=shp[a][0][k]*x(a, 1);
      zG+=shp[a][0][k]*x(a, 2);
    }
    myMatPoints[k]->set_X(xG, yG, zG);
  }

  // Permutation index
  for (int i = 0;i < 8;i++) perm[i]=i;
}
/**
 * Destructor.
 */
Brick8::~Brick8() {
  Containers::vector_delete(myMatPoints);
}
/**
 * Element stiffness matrix.
 */
const Matrix& Brick8::get_K() {
  // Static variables and references
  static Matrix Ba;
  static Matrix Bb;
  Matrix &K=*myMatrix;
  // Find shape functions for all GaussPoints (double shp[node][N, i][GPoint])
  this->shapeFunctions();
  // Stiffness matrix K+=B[a]T.C.B[b].dV (for each GaussPoint)
  K.Clear();
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    double dV = detJ[k];
    for (unsigned a = 0; a < nodes_.size(); a++) {
      this->get_B(&Ba, a, k);
      for (unsigned b = 0; b < nodes_.size(); b++) {
        this->get_B(&Bb, b, k);
        K.Add_BTCB(3*a, 3*b, &perm[0], Ba, C, Bb, dV, 1.0);
      }
    }
  }
  // Multiply by facK (active/deactivated)
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  return K;
}
/**
 * Element mass matrix.
 */
const Matrix& Brick8::get_M() {
  Matrix &M=*myMatrix;
  M.Clear();
  return M;
}
/**
 * Element residual vector.
 */
const Vector& Brick8::get_R() {
  // Static variables and references
  static Vector sigma(6);
  static Matrix Ba(6, 3);
  Vector& R=*myVector;
  R.Clear();
  // Quick return if not active
  if (!(groupdata_->active)) {
    return R;
  }
  // Get factors
  double facS = groupdata_->factor_S;
  double facG = groupdata_->factor_G;
  double facP = groupdata_->factor_P;
  // Find shape functions for all GaussPoints (double shp[node][N, i][GPoint])
  this->shapeFunctions();
  // R = facS*Fint - facG*SelfWeigth - facP*ElementalLoads
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    sigma = myMatPoints[k]->get_material()->get_stress();
    double dV = detJ[k];
    for (unsigned a = 0; a < nodes_.size(); a++) {
      // +facS*Fint
      this->get_B(&Ba, a, k);
      add_BTv(R, 3*a, &perm[0], Ba, sigma, facS*dV, 1.0);
      // -facG*SelfWeigth
      for (int i = 0;i < 3;i++)
        R[3*a+i]-=facG*shp[a][0][k]*b[i]*dV;
    }
  }
  // -facP*ElementalLoads
  R-=facP*P;
  // Return
  return R;
}
/**
 * Element update.
 */
void Brick8::update() {
  // Static variables and references
  static Vector u(24);
  static Vector epsilon(6);
  static Matrix B;
  // Quick return if inactive
  if (!(groupdata_->active)) {
    return;
  }
  // Get incremental displacements
  u = this->get_disp_incrm();
  // Find shape functions for all GaussPoints (double shp[node][N, i][GPoint])
  this->shapeFunctions();
  // Incremental strains: De+=B[a].Du  (for each GaussPoint)
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    epsilon.Clear();
    for (unsigned a = 0; a < nodes_.size(); a++) {
      this->get_B(&B, a, k);
      add2(epsilon, 3*a, B, u, 1.0, 1.0);
    }
    myMatPoints[k]->get_material()->set_strain(epsilon);
  }
}
/**
 * Element commit.
 */
void Brick8::commit() {
  for (unsigned int i = 0;i < myMatPoints.size();i++)
    myMatPoints[i]->get_material()->commit();
}
/**
 * Element shape functions.
 */
void Brick8::shapeFunctions() {
  shape8(x, shp, detJ);
}

/**
 * Initial stresses.
 */
void Brick8::AddInitialStresses(int direction,
                                 double h1, double s1,
                                 double h2, double s2, double K0) {
  for (unsigned i = 0;i < myMatPoints.size();i++)
    myMatPoints[i]->AddInitialStresses(direction, h1, s1, h2, s2, K0);
}

/**
 * Element stress recovery.
 */
void Brick8::recoverStresses() {
  static Vector sigma(6);
  static Matrix E(8, 8);
  const double d = 0.125;
  const double a = 1+num::sq3;
  const double b = 1-num::sq3;

  E(0, 0)=a*a*a*d;
  E(0, 1)=b*a*a*d;
  E(0, 2)=b*b*a*d;
  E(0, 3)=b*a*a*d;
  E(1, 0)=b*a*a*d;
  E(1, 1)=a*a*a*d;
  E(1, 2)=b*a*a*d;
  E(1, 3)=b*b*a*d;
  E(2, 0)=b*b*a*d;
  E(2, 1)=b*a*a*d;
  E(2, 2)=a*a*a*d;
  E(2, 3)=b*a*a*d;
  E(3, 0)=b*a*a*d;
  E(3, 1)=b*b*a*d;
  E(3, 2)=b*a*a*d;
  E(3, 3)=a*a*a*d;
  E(4, 0)=b*a*a*d;
  E(4, 1)=b*b*a*d;
  E(4, 2)=b*b*b*d;
  E(4, 3)=b*b*a*d;
  E(5, 0)=b*b*a*d;
  E(5, 1)=b*a*a*d;
  E(5, 2)=b*b*a*d;
  E(5, 3)=b*b*b*d;
  E(6, 0)=b*b*b*d;
  E(6, 1)=b*b*a*d;
  E(6, 2)=b*a*a*d;
  E(6, 3)=b*b*a*d;
  E(7, 0)=b*b*a*d;
  E(7, 1)=b*b*b*d;
  E(7, 2)=b*b*a*d;
  E(7, 3)=b*a*a*d;
  E(0, 4)=b*a*a*d;
  E(0, 5)=b*b*a*d;
  E(0, 6)=b*b*b*d;
  E(0, 7)=b*b*a*d;
  E(1, 4)=b*b*a*d;
  E(1, 5)=b*a*a*d;
  E(1, 6)=b*b*a*d;
  E(1, 7)=b*b*b*d;
  E(2, 4)=b*b*b*d;
  E(2, 5)=b*b*a*d;
  E(2, 6)=b*a*a*d;
  E(2, 7)=b*b*a*d;
  E(3, 4)=b*b*a*d;
  E(3, 5)=b*b*b*d;
  E(3, 6)=b*b*a*d;
  E(3, 7)=b*a*a*d;
  E(4, 4)=a*a*a*d;
  E(4, 5)=b*a*a*d;
  E(4, 6)=b*b*a*d;
  E(4, 7)=b*a*a*d;
  E(5, 4)=b*a*a*d;
  E(5, 5)=a*a*a*d;
  E(5, 6)=b*a*a*d;
  E(5, 7)=b*b*a*d;
  E(6, 4)=b*b*a*d;
  E(6, 5)=b*a*a*d;
  E(6, 6)=a*a*a*d;
  E(6, 7)=b*a*a*d;
  E(7, 4)=b*a*a*d;
  E(7, 5)=b*b*a*d;
  E(7, 6)=b*a*a*d;
  E(7, 7)=a*a*a*d;

  for (unsigned i = 0;i < 8;i++) {      // nodes
    sigma.Clear();
    for (unsigned j = 0;j < 6;j++) {    // sigma
      for (unsigned k = 0;k < 8;k++) {  // material points
        sigma[j]+=E(i, k)*(myMatPoints[k]->get_material()->get_stress())[j];
      }
    }
    nodes_[i]->addStress(sigma);
  }
}
/**
 * Element non-linear points.
 */
int Brick8::get_num_plastic_points() {
  int n = 0;
  for (unsigned i = 0;i < myMatPoints.size();i++)
    if (myMatPoints[i]->get_material()->isPlastic()) n+=1;
  return n;
}
