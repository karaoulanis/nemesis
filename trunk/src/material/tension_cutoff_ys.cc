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

#include "material/tension_cutoff_ys.h"


TensionCutOffYS::TensionCutOffYS(double T_) {
  T = T_;
}
double TensionCutOffYS::getf(const Vector& sigma, const double /*kappa*/) {
  this->setSigma(sigma);
  return I1-T;
}
const Vector& TensionCutOffYS::getdfds(const Vector& sigma, const double /*kappa*/) {
  this->setSigma(sigma);
  a = a1;
  return a;
}
const Matrix& TensionCutOffYS::getd2fdsds(const Vector& /*sigma*/, const double /*kappa*/) {
  da.clear();
  return da;
}
double TensionCutOffYS::getdfdk(const Vector& /*sigma*/, const double /*kappa*/) {
  return 0;
}
const Vector& TensionCutOffYS::getf2dkds(const Vector& /*sigma*/, const double /*kappa*/) {
  a.clear();
  return a;
}

