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

#include "material/drucker_prager.h"
#include "material/dp_in.h"
#include "material/dp_out.h"
#include "material/tc.h"

DruckerPrager::DruckerPrager()
    : type(0) {
}

DruckerPrager::DruckerPrager(int id, MultiaxialMaterial* elastic, int type_,
                             double c, double phi, double psi, double T)
    : MultiaxialElastoPlastic(id, elastic),
      type(type_) {
  // Material parameters
  MatParams[0] = c;
  MatParams[1] = phi;
  MatParams[2] = psi;
  MatParams[3] = T;
  // Yield/potential surfaces
  switch (type) {
    case 1:  // Inner no Tension cut-off
      fSurfaces.push_back(new DP_in(c, phi));
      gSurfaces.push_back(new DP_in(c, psi));
      break;
    case 2:  // Inner with Tension cut-off
      fSurfaces.push_back(new DP_in(c, phi));
      gSurfaces.push_back(new DP_in(c, psi));
      fSurfaces.push_back(new TC(T));
      gSurfaces.push_back(new TC(T));
      break;
    case 3:  // Outer no Tension cut-off
      fSurfaces.push_back(new DP_in(c, phi));
      gSurfaces.push_back(new DP_in(c, psi));
      break;
    case 4:  // Outer with Tension cut-off
      fSurfaces.push_back(new DP_in(c, phi));
      gSurfaces.push_back(new DP_in(c, psi));
      fSurfaces.push_back(new TC(T));
      gSurfaces.push_back(new TC(T));
      break;
    default:
      throw SException("[nemesis:%d] %s", 9999, "Invalid index %d.\n", type);
  }
  nHardeningVariables = 1;
}

DruckerPrager::~DruckerPrager() {
}

MultiaxialMaterial* DruckerPrager::get_clone() {
  // Material parameters
  int id    = this->get_id();
  double c    = MatParams[ 0];
  double phi  = MatParams[ 1];
  double psi  = MatParams[ 2];
  double T    = MatParams[ 3];
  // Create clone and return
  DruckerPrager* clone = new DruckerPrager(id, myElastic, type, c, phi, psi, T);
  return clone;
}
