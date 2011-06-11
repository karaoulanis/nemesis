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

#include "material/matpoint.h"
#include "exception/sexception.h"
#include "material/multiaxial_material.h"

int MatPoint::IDCounter = 0;

MatPoint::MatPoint()
    : myMaterial(0),
      r(0.),
      s(0.),
      t(0.),
      w(0.),
      x(0.),
      y(0.),
      z(0.) {
}


MatPoint::MatPoint(MultiaxialMaterial* mat, int index,
                   int p1)
    : DomainObject(IDCounter++),
      r(GaussCrds[p1][index+1]),
      s(0.),
      t(0.),
      w(GaussWght[p1][index+1]),
      x(0.),
      y(0.),
      z(0.) {
  /// @todo: Chech rule.
  if (p1>6) {
    throw SException("[nemesis:%d] %s", 9999, "Rule does not exist.");
  }
  myMaterial = mat->get_clone();
}


MatPoint::MatPoint(MultiaxialMaterial* mat, int index1, int index2,
                   int p1, int p2)
    : DomainObject(IDCounter++),
      r(GaussCrds[p1][index1]),
      s(GaussCrds[p2][index2]),
      t(0.),
      w(GaussWght[p1][index1]*GaussWght[p2][index2]),
      x(0.),
      y(0.),
      z(0.) {
  /// @todo: Chech rule.
  if (p1>6||p2>6) {
    throw SException("[nemesis:%d] %s", 9999, "Integration rule error.");
  }
  myMaterial = mat->get_clone();
}


MatPoint::MatPoint(MultiaxialMaterial* mat, int index1, int index2, int index3,
                   int p1, int p2, int p3)
    : DomainObject(IDCounter++),
      r(GaussCrds[p1][index1]),
      s(GaussCrds[p2][index2]),
      t(GaussCrds[p3][index3]),
      w(GaussWght[p1][index1]*GaussWght[p2][index2]*GaussWght[p3][index3]),
      x(0.),
      y(0.),
      z(0.) {
  /// @todo: Chech rule.
  if (p1>6||p2>6) {
    throw SException("[nemesis:%d] %s", 9999, "Integration rule error.");
  }
  myMaterial = mat->get_clone();
}


MatPoint::MatPoint(MultiaxialMaterial* mat, double r_, double s_, double t_,
                   double w_)
    : DomainObject(IDCounter++),
      r(r_),
      s(s_),
      t(t_),
      w(w_),
      x(0.),
      y(0.),
      z(0.) {
  myMaterial = mat->get_clone();
}


MatPoint::~MatPoint() {
  delete myMaterial;
}


void MatPoint::set_X(double x_, double y_, double z_) {
  x = x_;
  y = y_;
  z = z_;
  myMaterial->set_X(x, y, z);
}


bool MatPoint::isPlastic() {
  return myMaterial->isPlastic();
}

/**
 * Initial stresses.
 */
void MatPoint::AddInitialStresses(int direction,
                                 double h1, double s1,
                                 double h2, double s2, double K0) {
  static Vector s0(6);
  s0.clear();
  switch (direction) {
    case 1:
      if (x < h1 && x > h2) {
      s0[0] = s1+(s2-s1)*(x-h1)/(h2-h1);
      s0[1] = K0*s0[0];
      s0[2] = K0*s0[0];
      myMaterial->addStress(s0);
      }
      break;
  case 2:
    if (y < h1 && y > h2) {
      s0[1] = s1+(s2-s1)*(y-h1)/(h2-h1);
      s0[0] = K0*s0[1];
      s0[2] = K0*s0[1];
      myMaterial->addStress(s0);
    }
    break;
  case 3:
    if (z < h1 && z>h2) {
      s0[2] = s1+(s2-s1)*(z-h1)/(h2-h1);
      s0[0] = K0*s0[2];
      s0[1] = K0*s0[2];
      myMaterial->addStress(s0);
    }
    break;
  default:
    break;
  }
}


void MatPoint::Save(std::ostream* /*os*/) {
  /// @todo Implement this method.
}


const double MatPoint::GaussCrds[7][7]=   {  { 0.000000000000000,     // Rule 0
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 1
                                               0.000000000000000,     // [1][1]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 2
                                              -0.577350269189626,     // [2][1]
                                              +0.577350269189626,     // [2][2]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 3
                                              -0.774596669241483,     // [3][1]
                                               0.000000000000000,     // [3][2]
                                              +0.774596669241483,     // [3][3]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 4
                                              -0.861136311594053,     // [4][1]
                                              -0.339981043584856,     // [4][2]
                                              +0.339981043584856,     // [4][3]
                                              +0.861136311594053,     // [4][4]
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 5
                                              -0.906179845938664,     // [5][1]
                                              -0.538469310105683,     // [5][2]
                                               0.000000000000000,     // [5][3]
                                              +0.538469310105683,     // [5][4]
                                              +0.906179845938664,     // [5][5]
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 6
                                              -0.932469514203152,     // [6][1]
                                              -0.661209386466265,     // [6][2]
                                              -0.238619186083197,     // [6][3]
                                               0.238619186083197,     // [6][4]
                                               0.661209386466265,     // [1][5]
                                               0.932469514203152}};   // [6][6]

const double MatPoint::GaussWght[7][7]=   {  { 0.000000000000000,     // Rule 0
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 1
                                              +2.000000000000000,     // [1][1]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 2
                                              +1.000000000000000,     // [2][1]
                                              +1.000000000000000,     // [2][1]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 3
                                              +0.555555555555556,     // [3][1]
                                              +0.888888888888889,     // [3][2]
                                              +0.555555555555556,     // [3][3]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 4
                                              +0.347854845137454,     // [4][1]
                                              +0.652145154862546,     // [4][2]
                                              +0.652145154862546,     // [4][3]
                                              +0.347854845137454,     // [4][4]
                                               0.000000000000000,
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 5
                                              +0.236926885056189,     // [5][1]
                                              +0.478628670499366,     // [5][2]
                                              +0.568888888888889,     // [5][3]
                                              +0.478628670499366,     // [5][4]
                                              +0.236926885056189,     // [5][5]
                                               0.000000000000000},    // ______
                                             { 0.000000000000000,     // Rule 6
                                              +0.171324492379170,     // [6][1]
                                              +0.360761573048139,     // [6][2]
                                              +0.467913934572691,     // [6][3]
                                              +0.467913934572691,     // [6][4]
                                              +0.360761573048139,     // [1][5]
                                              +0.171324492379170}};   // [6][6]
