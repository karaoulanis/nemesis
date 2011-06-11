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

#include "elements/brick8i.h"
#include <vector>
#include "elements/shape_functions.h"
#include "group/group_data.h"
#include "main/nemesis_debug.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"

double Brick8i::detJ[8];
double Brick8i::shpStd[8][4][8];
double Brick8i::shpInc[3][4][8];

Brick8i::Brick8i() {
}


Brick8i::Brick8i(int id, std::vector<Node*> nodes,
                     MultiaxialMaterial* material)
    : Brick8(id, nodes, material) {
  aTrial.resize(9, 0.);
  aConvg.resize(9, 0.);
}


Brick8i::~Brick8i() {
}


const Matrix& Brick8i::get_K() {
  // Get a reference to myMatrix as K
  Matrix &K=*myMatrix;
  // Define local static matrices
  static Matrix Ba(6, 3), Bb(6, 3);
  static Matrix Kdd(24, 24), Kda(24, 9), Kaa(9, 9);
  // Form shape functions
  this->shapeFunctions();
  // Get Kdd, Kda, Kaa Matrices
  this->get_Kdd(Kdd);
  this->get_Kda(Kda);
  this->get_Kaa(Kaa);
  // Form K
  K = Kdd-Kda*Inverse(Kaa)*Transpose(Kda);
  // Get group factors
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  // Return K
  return K;
}


const Matrix& Brick8i::get_M() {
  // Get a reference to myMatrix as K
  Matrix &M=*myMatrix;
  M.Clear();
  // Find total mass
  double rho = myMaterial->get_rho();
  double volume = 0.;
  this->shapeFunctions();
  for (unsigned k = 0;k < myMatPoints.size();k++) {
    volume+=detJ[k]*(myMatPoints[k]->get_w());
  }
  double mass = rho*volume;
  // Set corresponding mass to diagonal terms
  for (int i = 0;i < 24;i++) {
    M(i, i)=0.25*mass;
  }
  // Return M
  return M;
}


const Vector& Brick8i::get_R() {
  // Static vectors and matrices
  static Vector sigma(6);
  static Matrix Ba(6, 3);
  // Get a reference to myVector as R
  Vector& R=*myVector;
  R.clear();
  // Quick return if inactive
  if (!(groupdata_->active)) {
    return R;
  }
  // Get factors
  double facS = groupdata_->factor_S;
  double facG = groupdata_->factor_G;
  double facP = groupdata_->factor_P;

  // Find shape functions for all GaussPoints
  this->shapeFunctions();

  // R = facS*Fint - facG*SelfWeigth - facP*ElementalLoads
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    sigma = myMatPoints[k]->get_material()->get_stress();
    double dV = detJ[k];
    for (unsigned a = 0; a < nodes_.size(); a++) {
      // +facS*Fint
      this->get_Bstd(Ba, a, k);
      add_BTv(R, 3*a, &perm[0], Ba, sigma, facS*dV, 1.0);
      // -facG*SelfWeigth
      for (int i = 0;i < 3;i++)
        R[3*a+i]-=facG*shpStd[a][0][k]*b[i]*dV;
    }
  }
  // -facP*ElementalLoads
  R-=facP*P;

  // Return R
  return R;
}


void Brick8i::update() {
  // Check for a quick return
  if (!(groupdata_->active)) {
    return;
  }
  // Static vectors and matrices
  static Vector Du(24), Da(9), epsilon(6);
  static Matrix Ba(6, 3);
  static Matrix Kda(24, 9), Kaa(9, 9);
  // Form shape functions
  this->shapeFunctions();
  // Get incremental displacements
  Du = this->get_disp_incrm();
  // Get incremental alphas
  this->get_Kda(Kda);
  this->get_Kaa(Kaa);
  Da=-Inverse(Kaa)*Transpose(Kda)*Du;
  aTrial = aConvg+Da;
  // For each material point
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    epsilon.clear();
    for (unsigned a = 0; a < 8; a++) {
      this->get_Bstd(Ba, a, k);
      add_Bv(epsilon, 3*a, &perm[0], Ba, Du, 1.0, 1.0);
    }
    for (unsigned a = 0; a < 3; a++) {
      this->get_BInc(Ba, a, k);
      add_Bv(epsilon, 3*a, &perm[0], Ba, Da, 1.0, 1.0);
    }
    myMatPoints[k]->get_material()->set_strain(epsilon);
  }
}


void Brick8i::commit() {
  for (unsigned int i = 0;i < myMatPoints.size();i++) {
    myMatPoints[i]->get_material()->commit();
  }
  aConvg = aTrial;
}


void Brick8i::shapeFunctions() {
  shape8(x, shpStd, detJ);
  shapeQM9(x, shpInc);
}


void Brick8i::get_Bstd(Matrix& B, int node, int gPoint) {
  // B-factors
  double B1 = shpStd[node][1][gPoint];
  double B2 = shpStd[node][2][gPoint];
  double B3 = shpStd[node][3][gPoint];

  // B-matrix
  B(0, 0) = B1;
  B(0, 1) = 0.;
  B(0, 2) = 0.;
  B(1, 0) = 0.;
  B(1, 1) = B2;
  B(1, 2) = 0.;
  B(2, 0) = 0.;
  B(2, 1) = 0.;
  B(2, 2) = B3;
  B(3, 0) = B2;
  B(3, 1) = B1;
  B(3, 2) = 0.;
  B(4, 0) = 0.;
  B(4, 1) = B3;
  B(4, 2) = B2;
  B(5, 0) = B3;
  B(5, 1) = 0.;
  B(5, 2) = B1;
}


void Brick8i::get_BInc(Matrix& B, int node, int gPoint) {
  // B-factors
  double B1 = shpInc[node][1][gPoint]/detJ[gPoint];
  double B2 = shpInc[node][2][gPoint]/detJ[gPoint];
  double B3 = shpInc[node][3][gPoint]/detJ[gPoint];

  // B-matrix
  B(0, 0) = B1;
  B(0, 1) = 0.;
  B(0, 2) = 0.;
  B(1, 0) = 0.;
  B(1, 1) = B2;
  B(1, 2) = 0.;
  B(2, 0) = 0.;
  B(2, 1) = 0.;
  B(2, 2) = B3;
  B(3, 0) = B2;
  B(3, 1) = B1;
  B(3, 2) = 0.;
  B(4, 0) = 0.;
  B(4, 1) = B3;
  B(4, 2) = B2;
  B(5, 0) = B3;
  B(5, 1) = 0.;
  B(5, 2) = B1;
}


void Brick8i::get_Kdd(Matrix& K) {
  K.Clear();
  static Matrix Ba(6, 3), Bb(6, 3);
  // For all Gauss points
  for (unsigned k = 0; k < 8; k++) {
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    for (unsigned a = 0; a < 8; a++) {
      get_Bstd(Ba, a, k);
      for (unsigned b = 0; b < 8; b++) {
        get_Bstd(Bb, b, k);
        double dV = detJ[k];
        K.add_BTCB(3*a, 3*b, &perm[0], Ba, C, Bb, dV, 1.0);
      }
    }
  }
}


void Brick8i::get_Kda(Matrix& K) {
  K.Clear();
  static Matrix Ba(6, 3), Bb(6, 3);
  // For all Gauss points
  for (unsigned k = 0; k < 8; k++) {
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    for (unsigned a = 0; a < 8; a++) {
      get_Bstd(Ba, a, k);
      for (unsigned b = 0; b < 3; b++) {
        get_BInc(Bb, b, k);
        double dV = 1.0*detJ[k];
        K.add_BTCB(3*a, 3*b, &perm[0], Ba, C, Bb, dV, 1.0);
      }
    }
  }
}


void Brick8i::get_Kaa(Matrix& K) {
  K.Clear();
  static Matrix Ba(6, 3), Bb(6, 3);
  // For all Gauss points
  for (unsigned k = 0; k < 8; k++) {
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    for (unsigned a = 0; a < 3; a++) {
      get_BInc(Ba, a, k);
      for (unsigned b = 0; b < 3; b++) {
        get_BInc(Bb, b, k);
        double dV = 1.0*detJ[k];
        K.add_BTCB(3*a, 3*b, &perm[0], Ba, C, Bb, dV, 1.0);
      }
    }
  }
}
