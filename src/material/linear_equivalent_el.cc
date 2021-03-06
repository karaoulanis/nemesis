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

#include "material/linear_equivalent_el.h"

LinearEquivalentEL::LinearEquivalentEL() {
}
LinearEquivalentEL::~LinearEquivalentEL() {
}
double LinearEquivalentEL::get_h(const Vector& v) {
  double eq = 2./3.*(v[0]*v[0]+    v[1]*v[1]+    v[2]*v[2]
          +0.5*v[3]*v[3]+0.5*v[4]*v[4]+0.5*v[5]*v[5]);
  return eq;
}
double LinearEquivalentEL::get_dhds(const Vector& /*sTrial*/,
                                    const Vector& /*ePTrial*/) {
  return 0;
}
double LinearEquivalentEL::get_dhda(const Vector& /*sTrial*/,
                                    const Vector& /*ePTrial*/) {
  return 0;
}
