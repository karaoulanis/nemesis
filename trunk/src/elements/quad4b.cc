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

#include "elements/quad4b.h"
#include "elements/shape_functions.h"
#include "group/group_data.h"
#include "main/nemesis_debug.h"
#include "material/matpoint.h"
#include "material/multiaxial_material.h"

double Quad4b::detJ[4];
double Quad4b::shp[4][3][4];
std::vector<int> Quad4b::perm(4);

Quad4b::Quad4b()
    : axisymmetric_(false) {
}

Quad4b::Quad4b(int id, std::vector<Node*> nodes, MultiaxialMaterial* material,
               double thickness, bool axisymmetric)
    : Quad4(id, nodes, material, thickness),
      axisymmetric_(axisymmetric) {
  perm[0]=0;
  perm[1]=1;
  perm[2]=2;
  perm[3]=3;
}

Quad4b::~Quad4b() {
}

const Matrix& Quad4b::get_K() {
  Matrix &K=*myMatrix;
  static Matrix Ba;
  static Matrix Bb;

  this->shapeFunctions();
  K.Clear();
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    const Matrix& C = myMatPoints[k]->get_material()->get_C();
    for (unsigned a = 0; a < nodes_.size(); a++) {
      this->get_B(Ba, a, k);
      for (unsigned b = 0; b < nodes_.size(); b++) {
        this->get_B(Bb, b, k);
        double dV = thickness_*detJ[k];
        K.Add_BTCB(2*a, 2*b, &perm[0], Ba, C, Bb, dV, 1.0);
      }
    }
  }
  double facK = groupdata_->active ? groupdata_->factor_K : 1e-7;
  K*=facK;
  return K;
}

const Matrix& Quad4b::get_M() {
  Matrix &M=*myMatrix;
  M.Clear();
  double rho = myMaterial->get_rho();
  double volume = 0.;

  this->shapeFunctions();
  for (unsigned k = 0;k < myMatPoints.size();k++)
    volume+=detJ[k]*thickness_*(myMatPoints[k]->get_w());
  double mass = rho*volume;
  for (int i = 0;i < 8;i++) M(i, i)=0.25*mass;
  return M;
}

/**
 * Element residual vector.
 */
const Vector& Quad4b::get_R() {
  // Static variables and references
  static Vector sigma(6);
  static Matrix Ba;
  Vector& R=*myVector;
  R.Clear();
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
      this->get_B(Ba, a, k);
      add_BTv(R, 2*a, &perm[0], Ba, sigma, facS*dV, 1.0);
      // -facG*SelfWeigth
      for (int i = 0;i < 2;i++)
        R[2*a+i]-=facG*shp[a][0][k]*b[i]*dV;
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
void Quad4b::update() {
  static Vector u(8);
  static Vector epsilon(6);
  static Matrix Ba;

  // Quick return if inactive
  if (!(groupdata_->active)) {
    return;
  }
  u = this->get_disp_incrm();

  this->shapeFunctions();
  // For each material point
  for (unsigned k = 0; k < myMatPoints.size(); k++) {
    epsilon.Clear();
    for (unsigned a = 0; a < nodes_.size(); a++) {
      this->get_B(Ba, a, k);
      /// @todo check dV
      // double dV = thickness_*detJ[k];
      add_Bv(epsilon, 2*a, &perm[0], Ba, u, 1.0, 1.0);
    }
    myMatPoints[k]->get_material()->set_strain(epsilon);
  }
}

void Quad4b::shapeFunctions() {
  shape4(x, shp, detJ);
  // Axisymmetry
  if (axisymmetric_) {
    for (int k = 0; k < 4; k++) {    // matpoints
      double r = 0.;
      for (int i = 0; i < 4; i++) {  // nodes
          r+=x(i, 0)*shp[i][0][k];
      }
      detJ[k]*=r;
    }
  }
}

void Quad4b::get_B(Matrix& B, int node, int gPoint) {
  B.Resize(4, 2);
  double Bb1 = 0., Bb2 = 0., vol = 0.;
  for (int k = 0;k < 4;k++) {  // matpoints
    double dV = thickness_*detJ[k];
    Bb1+=shp[node][1][k]*dV;
    Bb2+=shp[node][2][k]*dV;
    vol+=dV;
  }
  Bb1/=vol;
  Bb2/=vol;

  double B0 = 0., Bb0 = 0.;
  // Axisymmetry
  if (axisymmetric_) {
    double r = 0.;
    for (int i = 0;i < 4;i++) {
          r+=x(i, 0)*shp[i][0][gPoint];
    }
    B0 = shp[node][0][gPoint]/r;
    for (int k = 0;k < 4;k++) {
      r = 0.;
      for (int i = 0;i < 4;i++) {
          r+=x(i, 0)*shp[i][0][k];
      }
      Bb0+=shp[node][0][k]*detJ[k]*thickness_/vol/r;
    }
  }

  // B-factors
  double B1  = shp[node][1][gPoint];
  double B2  = shp[node][2][gPoint];
  double B4  = (Bb1-B1)/3.;
  double B6  = (Bb2-B2)/3.;
  double B7  = B2+B6;
  double B10 = B4+(Bb0-B0)/3.;
  double B11 = B0+B10;
  double B12 = B1+B10;

  // B-matrix
  B(0, 0) = B12;
  B(0, 1) = B6;
  B(1, 0) = B10;
  B(1, 1) = B7;
  B(2, 0) = B11;
  B(2, 1) = B6;
  B(3, 0) = B2;
  B(3, 1) = B1;
}
