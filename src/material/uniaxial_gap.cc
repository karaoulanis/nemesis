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

#include "material/uniaxial_gap.h"

UniaxialGap::UniaxialGap()
    : eTrial(0.),
      eElastMin(0.),
      eElastMax(0.),
      sy(0.),
      gap(0.),
      Et(0.) {
}


UniaxialGap::UniaxialGap(int ID, double E, double nu, double rho, double aT,
                         double sy_, double gap_)
    : UniaxialMaterial(ID, rho, aT),
      eTrial(0.),
      eElastMin(gap_),
      eElastMax(gap+sy/E),
      sy(sy_),
      gap(gap_),
      Et(E) {
  // Material parameters
  MatParams[0] = E;
  MatParams[1] = nu;
  MatParams[2] = sy;
  MatParams[3] = gap;
}


UniaxialMaterial* UniaxialGap::get_clone() {
  // Material parameters
  double E   = MatParams[ 0];
  double nu  = MatParams[ 1];
  double sy  = MatParams[ 2];
  double gap = MatParams[ 3];
  double rho = MatParams[30];
  double aT  = MatParams[31];
  // Create clone and return
  UniaxialMaterial* clone = new UniaxialGap(id_, E, nu, rho, aT, sy, gap);
  return clone;
}

void UniaxialGap::set_strain(const double De) {
  double E = MatParams[ 0];
  eTrial = eTotal+De;
  if (sy >= 0) {
    if (eTrial < eElastMin) {
      sTrial = 0.;
      Et = 0.;
    } else if (eTrial>eElastMax) {
      sTrial = sy;
      Et = 0.;
    } else {
      sTrial = E*(eTrial-eElastMin);
      Et = E;
    }
  } else {
    if (eTrial>eElastMin) {
      sTrial = 0.;
      Et = 0.;
    } else if (eTrial < eElastMax) {
      sTrial = sy;
      Et = 0.;
    } else {
      sTrial = E*(eTrial-eElastMin);
      Et = E;
    }
  }
}

double UniaxialGap::get_C() {
  return Et;
}

void UniaxialGap::Commit() {
  sConvg = sTrial;
  eTotal = eTrial;
}
