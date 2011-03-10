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

#include "loadcase/initial_stresses.h"
#include "elements/element.h"

InitialStresses::InitialStresses() {
}
InitialStresses::InitialStresses(const std::map<int, Element*>* elements,
                                 int group_id, int dir,
                                 double h1, double s1,
                                 double h2, double s2, double K0)
                                 :InitialCondition() {
  elements_ = elements;
  dir_ = dir;
  group_id_ = group_id;
  h1_ = h1;
  s1_ = s1;
  h2_ = h2;
  s2_ = s2;
  K0_ = K0;
  // Check direction
  ///@todo Give SException the right id and not 9999.
  if (dir < 1 || dir>3)
    throw SException("[nemesis:%d] Direction should be 1, 2 or 3.", 9999);
}
InitialStresses::~InitialStresses() {
}
int InitialStresses::Apply() {
  for (ElementIterator ei = elements_->begin(); ei != elements_->end(); ei++)
    ei->second->addInitialStresses(this);
  return 0;
}
