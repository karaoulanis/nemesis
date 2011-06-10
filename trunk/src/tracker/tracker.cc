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

#include "tracker/tracker.h"
#include <iostream>
#include <sstream>

Tracker::Tracker()
    : records_(0) {
}

Tracker::Tracker(int id)
    : DomainObject(id),
      records_(0) {
}


Tracker::~Tracker() {
}

int Tracker::get_steps() {
  return records_.size();
}

void Tracker::Save(std::ostream* os) {
  // start saving
  (*os) << "[";
  // save data
  for (unsigned i = 0; i < records_.size(); i++) {
    if (i > 0) (*os) <<",";
    (*os) << records_[i];
  }
  // finalize
  (*os) << "]";
}
