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

#include "material/ys.h"

YS::YS() {
  a.resize(3);
  a1.resize(3);
  a2.resize(3);
  da.resize(3, 3);
  da2.resize(3, 3);
  da11.resize(3, 3);
  da22.resize(3, 3);
  active = false;
}
void YS::set_sigma(const Vector& s) {
  // Deviatoric stresses
  double s0=(s[0]+s[1]+s[2])/3.;
  s1 = s[0]-s0;
  s2 = s[1]-s0;
  s3 = s[2]-s0;
  // Stress invariants
  I1 = s[0]+s[1]+s[2];
  J2=((s1-s2)*(s1-s2)+(s2-s3)*(s2-s3)+(s3-s1)*(s3-s1))/6.;
  // Vector dI/dsigma
  a1[0]=1.;
  a1[1]=1.;
  a1[2]=1.;
  // Vector dJ/dsigma
  a2[0]=s1;
  a2[1]=s2;
  a2[2]=s3;
  // Matrix d2J/dsigma2
  da2(0, 0)=+2.;  da2(0, 1)=-1.;  da2(0, 2)=-1.;
  da2(1, 0)=-1.;  da2(1, 1)=+2.;  da2(1, 2)=-1.;
  da2(2, 0)=-1.;  da2(2, 1)=-1.;  da2(2, 2)=+2.;
  da2*=1/3.;
  // Matrix dI1/dsigma * dI1/dsigma
  da11(0, 0)=1.;  da11(0, 1)=1.;  da11(0, 2)=1.;
  da11(1, 0)=1.;  da11(1, 1)=1.;  da11(1, 2)=1.;
  da11(2, 0)=1.;  da11(2, 1)=1.;  da11(2, 2)=1.;
  // Matrix dJ/dsigma * dJ/dsigma
  da22(0, 0)=s1*s1;  da22(0, 1)=s1*s2;  da22(0, 2)=s1*s3;
  da22(1, 0)=s2*s1;  da22(1, 1)=s2*s2;  da22(1, 2)=s2*s3;
  da22(2, 0)=s3*s1;  da22(2, 1)=s3*s2;  da22(2, 2)=s3*s3;
}
