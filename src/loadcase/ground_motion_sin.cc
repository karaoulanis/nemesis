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

#include "loadcase/ground_motion_sin.h"
#include "elements/element.h"

GroundMotionSin::GroundMotionSin() {
}

GroundMotionSin::GroundMotionSin(const std::map<int, Element*>* elements,
                                 int dof, double a, double omega, double phi)
                                 :Load() {
  elements_ = elements;
  dof_ = dof - 1;
  a_ = a;
  omega_ = omega;
  phi_ = phi;
}
void GroundMotionSin::Apply(double /*factor*/, double time) {
  std::map<int, Element*>::const_iterator i;
  for (i = elements_->begin(); i != elements_->end(); i++)
    i->second->addGroundMotion(dof_, a_*std::sin(omega_*time+phi_));
}
