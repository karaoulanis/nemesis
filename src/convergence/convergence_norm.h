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

#ifndef SRC_CONVERGENCE_CONVERGENCE_NORM_H_
#define SRC_CONVERGENCE_CONVERGENCE_NORM_H_

#include "analysis/analysis_object.h"
#include "numeric/vector.h"

class ConvergenceNorm: public AnalysisObject {
 protected:
  Vector tol;
  int maxIter;
  int iter;
  double ro, uo, wo;
  int LC, nSteps, step;
  public:
  ConvergenceNorm();
  ~ConvergenceNorm();
  void set_check(int maxIterations, double tolRabs, double tolRrel,
                double tolWrel);
  void init(int LCid, int steps);
  void newStep();
  int update();
};
#endif  // SRC_CONVERGENCE_CONVERGENCE_NORM_H_
