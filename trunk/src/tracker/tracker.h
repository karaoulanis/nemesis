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

#ifndef SRC_TRACKER_TRACKER_H_
#define SRC_TRACKER_TRACKER_H_

// C++ system files
#include <string>
#include <vector>
#include "domain/domain_object.h"

// namespaces
using std::string;
using std::vector;

class Tracker: public DomainObject {
 public:
  Tracker();
  explicit Tracker(int id);
  virtual ~Tracker();

  int get_steps();
  virtual void Track(double lambda, double time)=0;
  void Save(std::ostream* os);

 protected:
  vector<string> records_;
};
#endif  // SRC_TRACKER_TRACKER_H_
