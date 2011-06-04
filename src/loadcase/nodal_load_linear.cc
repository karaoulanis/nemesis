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

#include "loadcase/nodal_load_linear.h"

NodalLoadLinear::NodalLoadLinear()
    : initial_value_(0.),
      gradient_(0.) {
}

NodalLoadLinear::NodalLoadLinear(Node* node, int dof, double initial_value,
                                 double gradient)
    : NodalLoad(node, dof),
      initial_value_(initial_value),
      gradient_(gradient) {
}

double NodalLoadLinear::GetValue(double time) {
  return initial_value_ + (gradient_ * time);
}
