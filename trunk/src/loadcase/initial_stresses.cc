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
#include "group/group.h"

InitialStresses::InitialStresses() {
}
InitialStresses::InitialStresses(Group* group, int direction,
                                 double h1, double s1,
                                 double h2, double s2, double K0)
                                 :InitialCondition() {
  group_ = group;
  direction_ = direction;
  h1_ = h1;
  s1_ = s1;
  h2_ = h2;
  s2_ = s2;
  K0_ = K0;
  // Check direction
  /// @todo Give SException the right id and not 9999.
  if (direction_ < 1 || direction_ > 3) {
    throw SException("[nemesis:%d] Direction should be 1, 2 or 3.", 9999);
  }
}
InitialStresses::~InitialStresses() {
}

void InitialStresses::Apply() {
  for (unsigned i = 0; i < group_->get_elements().size(); i++) {
    group_->get_elements()[i]->AddInitialStresses(direction_,
                                                  h1_, s1_, h2_, s2_, K0_);
  }
}
