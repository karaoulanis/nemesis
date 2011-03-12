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

#include "material/vm.h"

VM::VM() {
}
VM::VM(double s0_, double K_) {
  s0 = s0_;
  K = K_;
}
VM::~VM() {
}
double VM::get_f(const Vector& s, const double kappa) {
  return num::sq3*sqrt(s.J2())-s0-K*kappa;
}
void VM::find_C(const Vector& s, const double /*a*/) {
  C1 = 0.;
  C2 = sqrt(0.75)/sqrt(s.J2());
  C3 = 0.;
  C11 = 0.;
  C22=-0.25*num::sq3*pow(s.J2(), -1.5);
  C23 = 0.;
  C32 = 0.;
  C33 = 0.;
}
double VM::get_dfda(const Vector& /*s*/, const double /*a*/) {
  return -K;
}
const Vector& VM::get_df2dsa(const Vector& /*s*/, const double /*a*/) {
  static Vector ret(6, 0.);
  ret.clear();
  return ret;
}
double VM::get_df2daa(const Vector& /*s*/, const double /*a*/) {
  return 0.;
}
