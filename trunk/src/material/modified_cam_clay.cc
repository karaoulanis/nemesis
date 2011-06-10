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

#include "material/modified_cam_clay.h"
#include "material/mcc.h"

ModifiedCamClay::ModifiedCamClay() {
}

ModifiedCamClay::ModifiedCamClay(int id, MultiaxialMaterial* elastic, double M,
                                 double po, double kappa, double lambda)
:MultiaxialElastoPlastic(id, elastic) {
  // Material parameters
  MatParams[0] = M;
  MatParams[1] = po;
  MatParams[2] = kappa;
  MatParams[3] = lambda;
  // Yield/potential surfaces
  fSurfaces.push_back(new MCC(M, po, kappa, lambda));
  gSurfaces.push_back(new MCC(M, po, kappa, lambda));
  // Material tag
  myTag = TAG_NONE;
  nHardeningVariables = 1;
}

ModifiedCamClay::~ModifiedCamClay() {
}

MultiaxialMaterial* ModifiedCamClay::get_clone() {
  // Material parameters
  int myID      = this->get_id();
  double M      = MatParams[ 0];
  double po     = MatParams[ 1];
  double kappa  = MatParams[ 2];
  double lambda = MatParams[ 3];
  // Create clone and return
  ModifiedCamClay* newClone = new ModifiedCamClay(myID, myElastic, M, po,
                                                  kappa, lambda);
  return newClone;
}
