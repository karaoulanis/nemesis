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

#include "elements/quad4i.h"
#include "elements/shape_functions.h"
#include "group/group_data.h"
#include "main/nemesis_debug.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"

double Quad4i::detJ[4];
double Quad4i::shpStd[4][3][4];
double Quad4i::shpInc[2][3][4];
std::vector<int> Quad4i::perm(3);

Quad4i::Quad4i()
  : aTrial(4, 0.),
    aConvg(4, 0.) {
}

Quad4i::Quad4i(int id, std::vector<Node*> nodes, MultiaxialMaterial* material,
               double thickness)
               : Quad4(id, nodes, material, thickness),
                 aTrial(4, 0.),
                 aConvg(4, 0.) {
  perm[0]=0;
  perm[1]=1;
  perm[2]=3;
}

Quad4i::~Quad4i() {
}

const Matrix& Quad4i::get_K() {
  // Get a reference to myMatrix as K
  Matrix &K=*myMatrix;
  // Define local static matrices
  static Matrix Ba(3, 2), Bb(3, 2);
  static Matrix Kdd(8, 8), Kda(8, 4), Kaa(4, 4);
  // Form shape functions
  this->shapeFunctions();
  // Get Kdd, Kda, Kaa Matrices
  this->get_Kdd(&Kdd);
  this->get_Kda(&Kda);
  this->get_Kaa(&Kaa);
  // Form K
  K = Kdd-Kda*Inverse(Kaa)*Transpose(Kda);
  // Get group factor
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  // Return K
  return K;
}
const Matrix& Quad4i::get_M() {
  // Get a reference to myMatrix as K
  Matrix &M=*myMatrix;
  M.Clear();
  // Find total mass
  double rho = myMaterial->get_rho();
  double volume = 0.;
  this->shapeFunctions();
  for (unsigned k = 0;k < myMatPoints.size();k++)
    volume += detJ[k]*thickness_*(myMatPoints[k]->get_w());
  double mass = rho*volume;
  // Set corresponding mass to diagonal terms
  for (int i = 0;i < 8;i++)
    M(i, i)=0.25*mass;
  // Return M
  return M;
}
/**
 * Element residual vector.
 */
const Vector& Quad4i::get_R() {
  // Get a reference to myVector as R
  Vector& R=*myVector;
  R.Clear();
  // Static vectors and matrices
  static Vector sigma(6);
  static Matrix Ba(3, 2);

  // Quick return if inactive
  if (!(groupdata_->active)) {
    return R;
  }

  // Factors
  double facS = groupdata_->factor_S;
  double facG = groupdata_->factor_G;
  double facP = groupdata_->factor_P;

  // Find shape functions for all GaussPoints
  this->shapeFunctions();

  // R = facS*Fint - facG*SelfWeigth - facP*ElementalLoads
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    sigma = myMatPoints[k]->get_material()->get_stress();
    double dV = thickness_*detJ[k];
    for (unsigned a = 0; a < nodes_.size(); a++) {
      // +facS*Fint
      this->get_Bstd(&Ba, a, k);
      add_BTv(&R, 2*a, &perm[0], Ba, sigma, facS*dV, 1.0);
      // -facG*SelfWeigth
      for (int i = 0;i < 2;i++)
        R[2*a+i]-=facG*shpStd[a][0][k]*b[i]*dV;
    }
  }
  // -facP*ElementalLoads
  R-=facP*P;

  // Return R
  return R;
}
/**
 * Element update.
 */
void Quad4i::update() {
  // Check for a quick return
  if (!(groupdata_->active)) {
    return;
  }
  // Static vectors and matrices
  static Vector Du(8), Da(4), epsilon(6);
  static Matrix Ba(3, 2);
  static Matrix Kda(8, 4), Kaa(4, 4);
  // Form shape functions
  this->shapeFunctions();
  // Get incremental displacements
  Du = this->get_disp_incrm();
  // Get incremental alphas
  this->get_Kda(&Kda);
  this->get_Kaa(&Kaa);
  Da=-Inverse(Kaa)*Transpose(Kda)*Du;
  aTrial = aConvg+Da;
  // For each material point
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    epsilon.Clear();
    for (unsigned a = 0; a < 4; a++) {
      this->get_Bstd(&Ba, a, k);
      /// @todo check
      // double dV = thickness_*detJ[k];
      add_Bv(&epsilon, 2*a, &perm[0], Ba, Du, 1.0, 1.0);
    }
    for (unsigned a = 0; a < 2; a++) {
      this->get_BInc(&Ba, a, k);
      /// @todo check
      // double dV = thickness_*detJ[k];
      add_Bv(&epsilon, 2*a, &perm[0], Ba, Da, 1.0, 1.0);
    }
    myMatPoints[k]->get_material()->set_strain(epsilon);
  }
}
/**
 * Element commit.
 * Overwrites base function, because incompatible modes must also
 * be commited. This takes place in element level.
 */
void Quad4i::commit() {
  for (unsigned int i = 0;i < myMatPoints.size();i++)
    myMatPoints[i]->get_material()->commit();
  aConvg = aTrial;
}
/**
 * Call shape functions and fill in corresponding arrays.
 */
void Quad4i::shapeFunctions() {
  shape4(x, shpStd, detJ);
  shapeQM6(x, shpInc);
}
/**
 * Get standard displacement B-matrix.
 * Shape function array shpStd and detJs must be defined.
 * @param B B-Matrix.
 * @param node The corresponding node.
 * @param gPoint The corresponding Gauss Point.
 */
void Quad4i::get_Bstd(Matrix* B, int node, int gPoint) {
  // B-factors
  double B1 = shpStd[node][1][gPoint];
  double B2 = shpStd[node][2][gPoint];

  // B-matrix
  (*B)(0, 0) = B1;
  (*B)(0, 1) = 0.;
  (*B)(1, 0) = 0.;
  (*B)(1, 1) = B2;
  (*B)(2, 0) = B2;
  (*B)(2, 1) = B1;
}

/**
 * Get B-matrix of incompatible modes.
 * Shape function array shpInc and detJs must be defined.
 * @param B B-Matrix.
 * @param node The corresponding node.
 * @param gPoint The corresponding Gauss Point.
 */
void Quad4i::get_BInc(Matrix* B, int node, int gPoint) {
  // B-factors
  double B1 = shpInc[node][1][gPoint]/detJ[gPoint];
  double B2 = shpInc[node][2][gPoint]/detJ[gPoint];

  // B-matrix
  (*B)(0, 0) = B1;
  (*B)(0, 1) = 0.;
  (*B)(1, 0) = 0.;
  (*B)(1, 1) = B2;
  (*B)(2, 0) = B2;
  (*B)(2, 1) = B1;
}

/**
 * Get Kdd Matrix.
 * This is the usual displacement matrix.
 * Formed as the integral of Bstd^T.C.Bstd on dOmega.
 * @todo Increase performance.
 * @param K The matrix to be filled.
 */
void Quad4i::get_Kdd(Matrix* K) {
  K->Clear();
  static Matrix Ba(3, 2), Bb(3, 2);
  // For all Gauss points
  for (unsigned k = 0; k < 4; k++) {
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    for (unsigned a = 0; a < 4; a++) {
      get_Bstd(&Ba, a, k);
      for (unsigned b = 0; b < 4; b++) {
        get_Bstd(&Bb, b, k);
        double dV = thickness_*detJ[k];
        K->Add_BTCB(2*a, 2*b, &perm[0], Ba, C, Bb, dV, 1.0);
      }
    }
  }
}

/**
 * Get Kda Matrix.
 * Formed as the integral of Bstd^T.C.Binc on dOmega.
 * @todo Increase performance.
 * @param K The matrix to be filled.
 */
void Quad4i::get_Kda(Matrix* K) {
  K->Clear();
  static Matrix Ba(3, 2), Bb(3, 2);
  // For all Gauss points
  for (unsigned k = 0; k < 4; k++) {
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    for (unsigned a = 0; a < 4; a++) {
      get_Bstd(&Ba, a, k);
      for (unsigned b = 0; b < 2; b++) {
        get_BInc(&Bb, b, k);
        double dV = 1.0*detJ[k];
        K->Add_BTCB(2*a, 2*b, &perm[0], Ba, C, Bb, dV, 1.0);
      }
    }
  }
}
/**
 * Get Kda Matrix.
 * Formed as the integral of Binc^T.C.Binc on dOmega.
 * @todo Increase performance.
 * @param K The matrix to be filled.
 */
void Quad4i::get_Kaa(Matrix* K) {
  K->Clear();
  static Matrix Ba(3, 2), Bb(3, 2);
  // For all Gauss points
  for (unsigned k = 0; k < 4; k++) {
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    for (unsigned a = 0; a < 2; a++) {
      get_BInc(&Ba, a, k);
      for (unsigned b = 0; b < 2; b++) {
        get_BInc(&Bb, b, k);
        double dV = 1.0*detJ[k];
        K->Add_BTCB(2*a, 2*b, &perm[0], Ba, C, Bb, dV, 1.0);
      }
    }
  }
}
