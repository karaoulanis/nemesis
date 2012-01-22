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

#include "material/sdof_material.h"

SDofMaterial::SDofMaterial() {
}


SDofMaterial::SDofMaterial(int ID, double E, double rho)
    : Material(ID, rho, 0.) {
  // Material parameters
  MatParams[0] = E;
}


SDofMaterial* SDofMaterial::get_clone() {
  // Material parameters
  double E   =MatParams[ 0];
  double rho =MatParams[30];
  // Create clone and return
  SDofMaterial* clone = new SDofMaterial(id_, E, rho);
  return clone;
}


void SDofMaterial::Commit() {
}


void SDofMaterial::Save(std::ostream* /*os*/) {
  /// @todo Implement this method.
}

