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

#include "material/multiaxial_material.h"

Matrix MultiaxialMaterial::C(6, 6);

MultiaxialMaterial::MultiaxialMaterial() {
}
MultiaxialMaterial::MultiaxialMaterial(int ID, double rho, double aT)
:Material(ID, rho, aT) {
  sTrial.Resize(6, 0.);
  sConvg.Resize(6, 0.);
  eTrial.Resize(6, 0.);
  eTotal.Resize(6, 0.);
}

MultiaxialMaterial::~MultiaxialMaterial() {
}

void MultiaxialMaterial::set_stress(const Vector& s) {
  sTrial =s;
  sConvg =s;
}

void MultiaxialMaterial::addStress(const Vector& s) {
  sTrial+=s;
  sConvg+=s;
}

const Vector& MultiaxialMaterial::get_stress() {
  return sTrial;
}

bool MultiaxialMaterial::isPlastic() {
  return false;
}
