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

#include "loadcase/ground_motion_file.h"
#include <cmath>
#include "elements/element.h"
#include "exception/sexception.h"

GroundMotionFile::GroundMotionFile() {
}

GroundMotionFile::GroundMotionFile(const std::map<int, Element*>* elements,
                                   int dof, std::istream& s, double dt,
                                   double scale)
                                   :Load() {
  elements_ = elements;
  dof_ = dof - 1;
  dt_ = dt;
  scale_ = scale;
  while (!s.eof()) {
    double d;
    s >> d;
    data_.push_back(d);
  }
}

void GroundMotionFile::Apply(double factor, double time) {
  if (time < 0. || dt_ < 0.) {
    throw SException("[nemesis:%d] %s", 9999, "Time/dt must pe positive.");
  }
  unsigned n = static_cast<int>(std::floor(time/dt_));
  if (n >= data_.size()) {
    throw SException("[nemesis:%d] %s", 9999, "Time exceeds given values.");
  }
  double t1 = n*dt_;
  double t2 = (n+1)*dt_;
  double d1 = data_[n];
  double d2 = data_[n+1];
  double d = d1+(time-t1)*(d2-d1)/(t2-t1);
  std::map<int, Element*>::const_iterator i;
  for (i = elements_->begin(); i != elements_->end(); i++) {
    i->second->addGroundMotion(dof_, factor*scale_*d);
  }
}
