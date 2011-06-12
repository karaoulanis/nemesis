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

#ifndef SRC_IMPOSER_PENALTY_IMPOSER_H_
#define SRC_IMPOSER_PENALTY_IMPOSER_H_

#include "imposer/imposer.h"

class PenaltyImposer: public Imposer {
 public:
  PenaltyImposer();
  PenaltyImposer(double a,
                 const std::map<int, Node*>& nodes,
                 const std::map<int, Element*>& elements,
                 const std::map<int, Constraint*>& constraints);
  int impose(Model* model);
 private:
  double a_;
};
#endif  // SRC_IMPOSER_PENALTY_IMPOSER_H_
