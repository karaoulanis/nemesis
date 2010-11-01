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

#include "material/mcc.h"

MCC::MCC() {
}
MCC::MCC(double M_, double po_, double kappa_, double lambda_) {
  M = M_;
  po = po_;
  kappa = kappa_;
  lambda = lambda_;
}
MCC::~MCC() {
}
double MCC::get_f(const Vector& s, const double /*q_*/) {
  double q = s.q();
  double p = s.p();
  double pc = 0.;//po*exp(e.I1()/(kappa+lambda));
  return q*q+M*M*p*(p-pc);
}
void MCC::find_C(const Vector& s, const double /*a*/) {
  C1 = 1/3.*M*M*(2/3.*s.I1()+po);
  C2 = 3.;
  C3 = 0.;
  C11 = 2./9.*M*M;
  C22 = 3.;
  C23 = 0.; C32 = 0.; C33 = 0.;
}
double MCC::get_dfda(const Vector& /*s*/, const double /*a*/) {
//  double dfda=(po/(lambda-kappa))*exp(e.I1()/(lambda-kappa));
//  return num::d13*M*M*s.I1()*dfda;
  return 0;
}
const Vector& MCC::get_df2dsa(const Vector& /*s*/, const double /*a*/) {
//  double dfda=(po/(lambda-kappa))*exp(e.I1()/(lambda-kappa));
//  double c = num::d13*M*M*s.I1()*dfda;
  static Vector ret(6, 0.);
  ret.clear();
  ret[0]=1.0; ret[1]=1.0; ret[2]=1.0;
//  ret*=c;
    return ret;
}
double MCC::get_df2daa(const Vector& /*s*/, const double /*a*/) {
//  double dfda=(po/(lambda-kappa))*exp(e.I1()/(lambda-kappa));
//  double d2fda2=(po/(lambda-kappa))*dfda;
//  return num::d13*M*M*s.I1()*d2fda2;
  return 0;
}

