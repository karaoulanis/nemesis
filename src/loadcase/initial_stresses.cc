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

InitialStresses::InitialStresses() {
}
InitialStresses::InitialStresses(int groupID_, int dir_, double h1_, double s1_,
                                 double h2_, double s2_, double K0_)
  :InitialCondition() {
  myTag = TAG_INITIAL_STRESSES;
  theGroupID = groupID_;
  dir = dir_;
  h1 = h1_;
  s1 = s1_;
  h2 = h2_;
  s2 = s2_;
  K0 = K0_;
  // Check if group exists
  pD->get<Group>(pD->get_groups(), theGroupID);
  // Check direction
  if (dir < 1 || dir>3)
    throw SException("[nemesis:%d] Direction should be 1, 2 or 3.", 9999);
}
InitialStresses::~InitialStresses() {
}
int InitialStresses::apply() {
  ElementContainer& DomainElems = pD->get_elements();
  for (ElementIterator ei = DomainElems.begin(); ei != DomainElems.end(); ei++)
    ei->second->addInitialStresses(this);
  return 0;
}
