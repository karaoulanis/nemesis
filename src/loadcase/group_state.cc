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

#include "loadcase/group_state.h"
#include "elements/element.h"
#include "group/group.h"

GroupState::GroupState() {
}

GroupState::GroupState(Group* group, int active,
                       double factor_K, double factor_S, double factor_G,
                       double factor_P) {
  group_    = group;
  groupdata_.active_   == 0 ? active = false : active = true;
  groupdata_.factor_K_ = factor_K;
  groupdata_.factor_S_ = factor_S;
  groupdata_.factor_G_ = factor_G;
  groupdata_.factor_P_ = factor_P;
}

void GroupState::Apply() {
  for (unsigned i = 0; i < group_->get_elements().size(); i++) {
    group_->get_elements()[i]->SetGroupData(&groupdata_);
  }
}
