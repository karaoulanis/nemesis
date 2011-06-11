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

#include "material/uniaxial_elastic_plastic.h"

UniaxialElastoPlastic::UniaxialElastoPlastic()
    : fTrial(0.),
      aTrial(0.),
      aConvg(0.),
      qTrial(0.),
      qConvg(0.),
      ePTrial(0.),
      ePConvg(0.),
      eTrial(0.) {
}

UniaxialElastoPlastic::UniaxialElastoPlastic(int ID, double E, double nu,
                                             double rho, double aT, double sy,
                                             double Hiso, double Hkin,
                                             double eta)
    : UniaxialMaterial(ID, rho, aT),
      fTrial(0.),
      aTrial(0.),
      aConvg(0.),
      qTrial(0.),
      qConvg(0.),
      ePTrial(0.),
      ePConvg(0.),
      eTrial(0.) {
  // Material parameters
  MatParams[0] = E;
  MatParams[1] = nu;
  MatParams[2] = sy;
  MatParams[3] = Hiso;
  MatParams[4] = Hkin;
  MatParams[5] = eta;
  // Tag
  myTag = TAG_MATERIAL_UNIAXIAL_ELASTIC;
}
UniaxialMaterial* UniaxialElastoPlastic::get_clone() {
  // Material parameters
  double E   =MatParams[ 0];
  double nu  =MatParams[ 1];
  double sy  =MatParams[ 2];
  double Hiso = MatParams[ 3];
  double Hkin = MatParams[ 4];
  double eta =MatParams[ 5];
  double rho =MatParams[30];
  double aT  =MatParams[31];
  // Create clone and return
  UniaxialElastoPlastic* clone =
    new UniaxialElastoPlastic(myID, E, nu, rho, aT, sy, Hiso, Hkin, eta);
  return clone;
}
void UniaxialElastoPlastic::set_strain(const double De) {
  // Material parameters
  eTrial = eTotal+De;
  double E    = MatParams[ 0];
  double sy   = MatParams[ 2];
  double Hiso = MatParams[ 3];
  // double Hkin = MatParams[ 4];
  // double eta =MatParams[ 5];

  // Return mapping
  /// @todo: implement kinematic hardening+viscoplasticity
  sTrial = sConvg+E*De;
  double xi = sTrial;  // -Hkin*q;
  // double dt = pD->get_time_incr();
  fTrial = fabs(xi)-(Hiso*aConvg+sy);
  aTrial = aConvg;
  ePTrial = ePConvg;
  if (fTrial >= 0) {
    // double dt = pD->get_time_incr();
    // dt = 1.0;
    // double dg = fTrial/(E+Hiso+eta/dt);
    double dg = fTrial/(E+Hiso);
    sTrial -=dg*E*num::sign(xi);
    ePTrial+=dg*num::sign(xi);
    aTrial +=dg;
  }
}
double UniaxialElastoPlastic::get_C() {
  // Material parameters
  double E    = MatParams[ 0];
  double Hiso = MatParams[ 3];
  double Hkin = MatParams[ 4];
  double eta  = MatParams[ 5];
  // double dt = pD->get_time_incr();
  double dt = 1.;
  if (fTrial>0) {
    E = E*(Hkin+Hiso+eta/dt)/(E+Hkin+Hiso+eta/dt);
  }
  return E;
}
void UniaxialElastoPlastic::commit() {
  sConvg = sTrial;
  aConvg = aTrial;
  qConvg = qTrial;
  eTotal = eTrial;
  ePConvg = ePTrial;
}
