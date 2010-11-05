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

#ifndef SRC_MATERIAL_DRUCKER_PRAGER_YS_H_
#define SRC_MATERIAL_DRUCKER_PRAGER_YS_H_

#include "material/ys.h"

class DruckerPragerYS: public YS {
  private:
  double c0, phi0, Kc, Kphi;
  public:
  DruckerPragerYS(double c_, double phi_, double Kc_, double Kphi_);
  double get_f(const Vector& sigma, const double kappa);
  const Vector& get_dfds(const Vector& sigma, const double kappa);
  const Matrix& get_d2fdsds(const Vector& sigma, const double kappa);
  double get_dfdk(const Vector& sigma, const double kappa);
  const Vector& get_f2dkds(const Vector& sigma, const double kappa);
};

#endif  // SRC_MATERIAL_DRUCKER_PRAGER_YS_H_
