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

#include "material/matpoint.h"

int MatPoint::IDCounter = 0;

MatPoint::MatPoint() {
  delete myMaterial;
}
MatPoint::MatPoint(MultiaxialMaterial* mat, int index, int p1)
:DomainObject(IDCounter++) {
  if (p1>6) exit(-987);
  r = GaussCrds[p1][index+1];
  s = 0;
  t = 0;
  w = GaussWght[p1][index+1];
  myMaterial = mat->getClone();
  myTag = TAG_MATERIAL_POINT;
}
MatPoint::MatPoint(MultiaxialMaterial* mat, int index1, int index2, int p1, int p2)
:DomainObject(IDCounter++) {
  if (p1>6||p2>6) exit(-987);
  r = GaussCrds[p1][index1];
  s = GaussCrds[p2][index2];
  t = 0.;
  w = GaussWght[p1][index1]*GaussWght[p2][index2];
  myMaterial = mat->getClone();
  myTag = TAG_MATERIAL_POINT;
}
MatPoint::MatPoint(MultiaxialMaterial* mat, int index1, int index2, int index3, int p1, int p2, int p3)
:DomainObject(IDCounter++) {
  if (p1>6||p2>6) exit(-987);
  r = GaussCrds[p1][index1];
  s = GaussCrds[p2][index2];
  t = GaussCrds[p3][index3];
  w = GaussWght[p1][index1]*GaussWght[p2][index2]*GaussWght[p3][index3];
  myMaterial = mat->getClone();
  myTag = TAG_MATERIAL_POINT;
}
MatPoint::MatPoint(MultiaxialMaterial* mat, double r_, double s_, double t_, double w_)
:DomainObject(IDCounter++) {
  r = r_;
  s = s_;
  t = t_;
  w = w_;
  myMaterial = mat->getClone();
  myTag = TAG_MATERIAL_POINT;
}
MatPoint::~MatPoint() {
  delete myMaterial;
}
void MatPoint::setX(double x_, double y_, double z_) {
  x = x_;
  y = y_;
  z = z_;
  myMaterial->setX(x, y, z);
}
void MatPoint::setInitialStresses(InitialStresses* pInitialStresses) {
  static Vector s0(6);
  s0.clear();
  int dir = pInitialStresses->getDir();
  double H1 = pInitialStresses->getH1();
  double H2 = pInitialStresses->getH2();
  double s1 = pInitialStresses->getS1();
  double s2 = pInitialStresses->getS2();
  switch(dir)
  {
  case 1:
    if (x < H1 && x>H2)
    {
      s0[0]=s1+(s2-s1)*(x-H1)/(H2-H1);
      s0[1]=pInitialStresses->getK0()*s0[0];
      s0[2]=pInitialStresses->getK0()*s0[0];
      myMaterial->addStress(s0);
    }
    break;
  case 2:
    if (y < H1 && y>H2)
    {
      s0[1]=s1+(s2-s1)*(y-H1)/(H2-H1);
      s0[0]=pInitialStresses->getK0()*s0[1];
      s0[2]=pInitialStresses->getK0()*s0[1];
      myMaterial->addStress(s0);
    }
    break;
  case 3:
    if (z < H1 && z>H2)
    {
      s0[2]=s1+(s2-s1)*(z-H1)/(H2-H1);
      s0[0]=pInitialStresses->getK0()*s0[2];
      s0[1]=pInitialStresses->getK0()*s0[2];
      myMaterial->addStress(s0);
    }
    break;
  default:
    break;
  }
}

const double MatPoint::GaussCrds[7][7]=   {{ 0.000000000000000,  // Rule 0
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000}, //_______
                       { 0.000000000000000,  // Rule 1
                                               0.000000000000000,  // [1][1] 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 2
                                              -0.577350269189626,  // [2][1]
                                              +0.577350269189626,  // [2][2]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 3
                                              -0.774596669241483,  // [3][1]
                                               0.000000000000000,  // [3][2]
                                              +0.774596669241483,  // [3][3]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 4
                                              -0.861136311594053,  // [4][1]
                                              -0.339981043584856,  // [4][2]
                                              +0.339981043584856,  // [4][3]
                                              +0.861136311594053,  // [4][4]
                                               0.000000000000000,
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 5
                                              -0.906179845938664,  // [5][1]
                                              -0.538469310105683,  // [5][2]
                                               0.000000000000000,  // [5][3]
                                              +0.538469310105683,  // [5][4]
                                              +0.906179845938664,  // [5][5]
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 6
                                              -0.932469514203152,  // [6][1]
                                              -0.661209386466265,  // [6][2]
                                              -0.238619186083197,  // [6][3]
                                               0.238619186083197,  // [6][4]
                                               0.661209386466265,  // [1][5]
                                               0.932469514203152}};// [6][6]

const double MatPoint::GaussWght[7][7]=   {{ 0.000000000000000,  // Rule 0
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 1
                                              +2.000000000000000,  // [1][1] 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000, 
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 2
                                              +1.000000000000000,  // [2][1]
                                              +1.000000000000000,  // [2][1]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 3
                                              +0.555555555555556,  // [3][1]
                                              +0.888888888888889,  // [3][2]
                                              +0.555555555555556,  // [3][3]
                                               0.000000000000000,
                                               0.000000000000000,
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 4
                                              +0.347854845137454,  // [4][1]
                                              +0.652145154862546,  // [4][2]
                                              +0.652145154862546,  // [4][3]
                                              +0.347854845137454,  // [4][4]
                                               0.000000000000000,
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 5
                                              +0.236926885056189,  // [5][1]
                                              +0.478628670499366,  // [5][2]
                                              +0.568888888888889,  // [5][3]
                                              +0.478628670499366,  // [5][4]
                                              +0.236926885056189,  // [5][5]
                                               0.000000000000000}, //_______
                                             { 0.000000000000000,  // Rule 6
                                              +0.171324492379170,  // [6][1]
                                              +0.360761573048139,  // [6][2]
                                              +0.467913934572691,  // [6][3]
                                              +0.467913934572691,  // [6][4]
                                              +0.360761573048139,  // [1][5]
                                              +0.171324492379170}};// [6][6]
