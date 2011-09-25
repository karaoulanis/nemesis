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

#ifndef SRC_CONTROL_LOAD_CONTROL_H_
#define SRC_CONTROL_LOAD_CONTROL_H_

#include "control/static_control.h"

/**
 * The LoadControl Class.
 * LoadControl is a static control that increases \f$\lambda_n\f$ according
 * to the relation \f$\lambda_n=\lambda_o+\Delta\lambda\f$, where 'o' stands
 * for old (previous) step and 'n' stands for new (current step).
 * \f$\Delta\lambda\f$ may vary between different steps (auto-incrementation)
 * but is constant within each step, i.e. \f$\delta\lambda\f$=0.
 * \f$\lambda_n\f$ is in fact a trial \f$\lambda_n\f$, and becomes
 * \f$\lambda_o\f$ for the next step if and only if this step converges.\n
 */
class LoadControl :public StaticControl {
  private:
  public:
  // Constructor and destructor
  LoadControl(double DL0, double minDL, double maxDL, int IterDesired, double n,
              double DeltaTime);
  ~LoadControl();

  // Methods for incremental/iterative algorithms
  void Predict();
  void Correct();
};

#endif  // SRC_CONTROL_LOAD_CONTROL_H_
