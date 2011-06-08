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

#include "group/group.h"
#include "group/group_data.h"
#include "elements/element.h"

Group::Group()
    : groupdata_(0),
      elements_(0) {
}

Group::Group(int id)
    : DomainObject(id),
      groupdata_(new GroupData()),
      elements_(0) {
  groupdata_->Reset();
}

Group::~Group() {
  delete groupdata_;
}

void Group::SetGroupData(const GroupData* groupdata) {
  // Just copy data.
  groupdata_->active   = groupdata->active;
  groupdata_->factor_G = groupdata->factor_G;
  groupdata_->factor_K = groupdata->factor_K;
  groupdata_->factor_P = groupdata->factor_P;
  groupdata_->factor_S = groupdata->factor_S;
}

void Group::AddElement(Element* element) {
  elements_.push_back(element);
  element->SetGroupData(groupdata_);
}
