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

#include "material/drucker_prager_ys.h"

DruckerPragerYS::DruckerPragerYS(double c_, double phi_, double Kc_,
                                 double Kphi_) {
  c0 = c_;
  phi0 = phi_;
  Kc = Kc_;
  Kphi = Kphi_;
}
double DruckerPragerYS::getf(const Vector& sigma, const double kappa) {
  this->setSigma(sigma);
  double c = c0+Kc*kappa;
  double phi = phi0+Kphi*kappa;
  double rho = 2*sin(phi)/(sqrt(3.)*(3-sin(phi)));
  double k = 6*c*cos(phi)/(sqrt(3.)*(3-sin(phi)));
  return rho*I1+sqrt(J2)-k;
}
const Vector& DruckerPragerYS::getdfds(const Vector& sigma,
                                       const double kappa) {
  this->setSigma(sigma);
  double phi = phi0+Kphi*kappa;
  double rho = 2*sin(phi)/(sqrt(3.)*(3-sin(phi)));

  double C1 = rho;
  double C2 = 1./(2.*sqrt(J2));
  a = a1*C1+a2*C2;
  return a;
}
const Matrix& DruckerPragerYS::getd2fdsds(const Vector& sigma,
                                          const double /*kappa*/) {
  this->setSigma(sigma);

  double C2 = 1./(2.*sqrt(J2));
  double C22=-1./(4.*J2*sqrt(J2));
  da = C2*da2+C22*da22;
  return da;
}
double DruckerPragerYS::getdfdk(const Vector& sigma, const double kappa) {
  this->setSigma(sigma);
  double c = c0+Kc*kappa;
  double phi = phi0+Kphi*kappa;

  double denom = -10.+6.*sin(phi+Kphi*kappa)
                 +cos(phi+Kphi*kappa)*cos(phi+Kphi*kappa);
  double drhodkappa = -2.*cos(phi+Kphi*kappa)*Kphi*sqrt(3.)/denom;
  double dkdkappa = -2*cos(phi)*sqrt(3.)*((3.*Kc)
                        -Kc*sin(phi+Kphi*kappa)
                        +Kphi*cos(phi+Kphi*kappa)*c
                        +Kphi*cos(phi+Kphi*kappa)*Kc*kappa)/denom;
  return I1*drhodkappa-dkdkappa;
}
const Vector& DruckerPragerYS::getf2dkds(const Vector& sigma,
                                         const double kappa) {
  this->setSigma(sigma);
  double phi = phi0+Kphi*kappa;

  double denom = -10.+6.*sin(phi+Kphi*kappa)
                 +cos(phi+Kphi*kappa)*cos(phi+Kphi*kappa);
  double drhodkappa=-2.*cos(phi+Kphi*kappa)*Kphi*sqrt(3.)/denom;

  double C1 = drhodkappa;
  a = a1*C1;
  return a;
}

