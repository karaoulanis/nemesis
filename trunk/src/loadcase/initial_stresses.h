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

#ifndef SRC_LOADCASE_INITIAL_STRESSES_H_
#define SRC_LOADCASE_INITIAL_STRESSES_H_

#include "loadcase/initial_condition.h"
#include <map>

class Group;

class InitialStresses: public InitialCondition {
 public:
  InitialStresses();
  InitialStresses(Group* group, int direction,
                  double h1, double s1, double h2, double s2, double K0);
  ~InitialStresses();
  void Apply();

 private:
  Group* group_;
  int direction_;
  double h1_;
  double s1_;
  double h2_;
  double s2_;
  double K0_;
};
#endif  // SRC_LOADCASE_INITIAL_STRESSES_H_
