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

#include "tracker/tracker.h"
// C++ files
#include <sstream>
#include <iostream>

Tracker::Tracker() {
}
Tracker::~Tracker() {
}
int Tracker::get_steps() {
  return records_.size();
}
/**
 * /// @todo: check with const char*
 */
void Tracker::track(string record) {
  records_.push_back(record);
}
const char* Tracker::data() {
  // define an output string stream
  std::ostringstream s;
  // start saving
  s << "[";
  // save data
  for (int i = 0; i < records_.size(); i++) {
    if (i > 0) s <<",";
    s << records_[i];
  }
  // finalize
  s << "]";
  // convert to c style string and return
  // needs to be converted to a static string before
  /// @todo: check for refactoring
  static string tmp;
  tmp = s.str();
  return tmp.c_str();
}
