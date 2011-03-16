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

#ifndef SRC_CONTROL_NEWMARK_H_
#define SRC_CONTROL_NEWMARK_H_

#include "control/transient_control.h"

/**
 * The Newmark Class.
  */
class Newmark :public TransientControl {
  private:
  double beta;
  double gamma;
  double dt;
  public:
  // Constructor and destructor
  Newmark(double beta_, double gamma_, double dt_);
  ~Newmark();

  // Methods for incremental/iterative algorithms
  void predict();
  void correct();
};

#endif  // SRC_CONTROL_NEWMARK_H_
