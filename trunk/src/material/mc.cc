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

#include "material/mc.h"

MC::MC()
:c(0.), phi(0.) {
}
MC::MC(double c_, double phi_)
:c(c_), phi(phi_) {
}
MC::~MC() {
}
void MC::find_C(const Vector& s, const double /*a*/) {
  double theta = s.theta();
  double J2 = s.J2();
  double sqJ2 = sqrt(J2);
  double A, dA, d2A;
  this->find_A(s, A, dA, d2A);

  // Constants Ci, Cii [Crisfield II, p.106, (Tab. 14.3), p.106, (Tab. 14.2)]
  if (fabs(theta)*180./num::pi < 29.99) {
    C1 = sin(phi)*num::d13;
    C2 = 0.5/sqJ2*(A-tan(3.*theta)*dA);
    C3 =-sqrt(3.)*dA/(2.*J2*cos(3.*theta));
    C4 = d2A+3.*tan(3.*theta)*dA;
    C11 = 0.;
    C22 =-(A-tan(3.*theta)*tan(3.*theta)*C4-3.*tan(3.*theta)*dA)/(4.*J2*sqJ2);
    C23 = (0.5*tan(3*theta)*C4+dA)*sqrt(3.)/(2.*J2*J2*cos(3*theta));
    C32 = C23;
    C33 = 3.*C4/(4.*J2*J2*sqJ2*cos(3.*theta)*cos(3.*theta));
  } else {
    C1  = 2*sin(phi)/(num::sq3*(3+sin(phi)));
    C2  = 0.5/sqrt(s.J2());
    C3  = 0.;
    C11 = 0.;
    C22 =-0.25*pow(s.J2(), -1.5);
    C23 = 0.;
    C32 = 0.;
    C33 = 0.;
    // C1 = 2*sin(phi)/(num::sq3*(3-sin(phi)));
    // C2 = 0.5/sqrt(s.J2());
    // C3 = 0.;
    // C11 = 0.;
    // C22=-0.25*pow(s.J2(), -1.5);
    // C23 = 0.; C32 = 0.; C33 = 0.;
  }
}
double MC::get_f(const Vector& s, const double /*q*/) {
  double A, dA, d2A;
  this->find_A(s, A, dA, d2A);
  double I1 = s.I1();
  double J2 = s.J2();
  return I1*sin(phi)/3.+sqrt(J2)*A-c*cos(phi);
}
