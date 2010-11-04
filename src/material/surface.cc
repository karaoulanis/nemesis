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

#include "material/surface.h"

Vector Surface::myVector(6, 0.);
Matrix Surface::myMatrix(6, 6, 0.);
Vector Surface::a1(6, 0.);
Vector Surface::a2(6, 0.);
Vector Surface::a3(6, 0.);
Matrix Surface::da2(6, 6, 0.);
Matrix Surface::da3(6, 6, 0.);
double Surface::C1 = 0.;
double Surface::C2 = 0.;
double Surface::C3 = 0.;
double Surface::C4 = 0.;
double Surface::C11 = 0.;
double Surface::C22 = 0.;
double Surface::C23 = 0.;
double Surface::C32 = 0.;
double Surface::C33 = 0.;

Surface::Surface()
:active(false) {
}

Surface::~Surface() {
}
void Surface::find_a(const Vector& s) {
  double sb = s.sb();
  double J2 = s.J2();
    double sx = s[0]-sb;
  double sy = s[1]-sb;
  double sz = s[2]-sb;
  double txy = s[3];
  double tyz = s[4];
  double tzx = s[5];

  // Flow vectors a1, a2, a3 [Crisfield II, p.105, (14.20), (14.21), (14.22)].
  a1[0]=1.0;        a2[0]=sx;         a3[0]=sy*sz-tyz*tyz+J2/3.;
  a1[1]=1.0;        a2[1]=sy;         a3[1]=sz*sx-tzx*tzx+J2/3.;
  a1[2]=1.0;        a2[2]=sz;         a3[2]=sx*sy-txy*txy+J2/3.;
  a1[3]=0.0;        a2[3]=2.*txy;     a3[3]=2.*(tyz*tzx-sz*txy);
  a1[4]=0.0;        a2[4]=2.*tyz;     a3[4]=2.*(tzx*txy-sx*tyz);
  a1[5]=0.0;        a2[5]=2.*tzx;     a3[5]=2.*(txy*tyz-sy*tzx);
}
void Surface::find_da(const Vector& s) {
  double sb = s.sb();
  double sx = s[0]-sb;
  double sy = s[1]-sb;
  double sz = s[2]-sb;
  double txy = s[3];
  double tyz = s[4];
  double tzx = s[5];
  // Matrix da2 [Crisfield II, p.105, (14.24)].
  da2.clear();
  da2(0, 0)= 2;  da2(0, 1)=-1;  da2(0, 2)=-1;
  da2(1, 0)=-1;  da2(1, 1)= 2;  da2(1, 2)=-1;
  da2(2, 0)=-1;  da2(2, 1)=-1;  da2(2, 2)= 2;
  da2(3, 3)= 6;
  da2(4, 4)= 6;
  da2(5, 5)= 6;
  da2*=num::d13;

  // Matrix da3 [Crisfield II, p.105, (14.24)].
  da3(0, 0)= sx;      da3(0, 1)= sz;      da3(0, 2)= sy;      da3(0, 3)= txy;    da3(0, 4)=-2*tyz;  da3(0, 5)= tzx;
  da3(1, 0)= sz;      da3(1, 1)= sy;      da3(1, 2)= sx;      da3(1, 3)= txy;    da3(1, 4)= tyz;    da3(1, 5)= 2*tzx;
  da3(2, 0)= sy;      da3(2, 1)= sx;      da3(2, 2)= sz;      da3(2, 3)=-2*txy;  da3(2, 4)= tyz;    da3(2, 5)= tzx;
  da3(3, 0)= txy;     da3(3, 1)= txy;     da3(3, 2)=-2*txy;   da3(3, 3)=-3*sz;   da3(3, 4)= 3*tzx;  da3(3, 5)= 3*tyz;
  da3(4, 0)=-2*tyz;   da3(4, 1)= tyz;     da3(4, 2)= tyz;     da3(4, 3)= 3*tzx;  da3(4, 4)=-3*sx;   da3(4, 5)= 3*txy;
  da3(5, 0)= tzx;     da3(5, 1)= 2*tzx;   da3(5, 2)= tzx;     da3(5, 3)= 3*tyz;  da3(5, 4)= 3*txy;  da3(5, 5)=-3*sy;
  da3*=num::d23;
}
void Surface::find_C(const Vector& /*s*/, const double /*a*/) {
  C11 = 0.; C1 = 0; C2 = 0; C3 = 0;
  C23 = 0;  C32 = 0;  C22 = 0;  C33 = 0;
}
const Vector& Surface::get_dfds(const Vector& s, const double a) {
  this->find_C(s, a);
  this->find_a(s);
  myVector = C1*a1+C2*a2+C3*a3;
  return myVector;
}
double Surface::get_dfda(const Vector& /*s*/, const double /*a*/) {
  return 0;
}
const Matrix& Surface::get_df2dss(const Vector& s, const double a) {
  this->find_C(s, a);
  this->find_a(s);
  this->find_da(s);
  myMatrix = C2*da2+C3*da3+C11*VVT(a1, a1)
              +C22*VVT(a2, a2)
              +C23*VVT(a2, a3)
              +C32*VVT(a3, a2)
              +C33*VVT(a3, a3);
  return myMatrix;
}
const Vector& Surface::get_df2dsa(const Vector& /*s*/, const double /*a*/) {
  myVector.clear();
  return myVector;
}
double  Surface::get_df2daa(const Vector& /*s*/, const double /*a*/) {
  return 0;
}
