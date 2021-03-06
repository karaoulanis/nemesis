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

#include "material/uniaxial_cyclic.h"
#include <cmath>
#include <iostream>

UniaxialCyclic::UniaxialCyclic()
    : sr(0.),
      er(0.),
      eTrial(0.),
      eConvg(0.),
      Et(0.),
      reversed(false) {
}

UniaxialCyclic::UniaxialCyclic(int id, double E, double nu, double rho,
                               double aT, double tmax, double Gmax)
    : UniaxialMaterial(id, rho, aT),
      sr(0.),
      er(0.),
      eTrial(0.),
      eConvg(0.),
      Et(Gmax),
      reversed(false) {
  // Material parameters
  MatParams[0] = E;
  MatParams[1] = nu;
  MatParams[2] = tmax;
  MatParams[3] = Gmax;
}

UniaxialCyclic::~UniaxialCyclic() {
}

UniaxialMaterial* UniaxialCyclic::get_clone() {
  double E    = MatParams[0];
  double nu   = MatParams[1];
  double tmax = MatParams[2];
  double Gmax = MatParams[3];
  double rho  = MatParams[30];
  double aT   = MatParams[31];
  UniaxialMaterial* newClone = new UniaxialCyclic(id_, E, nu, rho, aT, tmax,
                                                  Gmax);
  return newClone;
}

void UniaxialCyclic::set_strain(const double De) {
  double tmax = MatParams[2];
  double Gmax = MatParams[3];
  // At first step check if reversed
  if (((eConvg+De-er))*De < 0) {
    sr = sConvg;
    er = eConvg;
    reversed = true;
    std::cout << er << '\t' << sr << std::endl;
  }
  eTrial = eConvg+De;
  // cout << er << '\t '<< eTrial << '\t' << reversed << endl;
  if (reversed == false) {
    sTrial = Gmax*eTrial/(1+(Gmax/tmax)*std::fabs(eTrial));
    Et = Gmax*tmax*
        (tmax+Gmax*std::fabs(eTrial)-Gmax*eTrial*num::sign(eTrial))/
        ((tmax+Gmax*std::fabs(eTrial))*(tmax+Gmax*fabs(eTrial)));
  } else {
    sTrial = sr+2.0*Gmax*(0.5*(eTrial-er))/
            (1+(Gmax/tmax)*fabs(0.5*(eTrial-er)));
    Et = 0.5*Gmax*tmax*
        (2*tmax+Gmax*std::fabs(er-eTrial)
        -0.5*Gmax*num::sign(er-eTrial)*(er-eTrial))/
        ((tmax+0.5*Gmax*std::fabs(er-eTrial))*
         (tmax+0.5*Gmax*std::fabs(er-eTrial)));
  }
}

double UniaxialCyclic::get_C() {
  /// @todo
  return Et;
}

void UniaxialCyclic::Commit() {
  sConvg = sTrial;
  eConvg = eTrial;
  eTotal = eConvg;
}
