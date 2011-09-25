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

#include "algorithm/linear_algorithm.h"
#include "analysis/analysis.h"
#include "control/control.h"
#include "convergence/convergence_norm.h"

LinearAlgorithm::LinearAlgorithm() {
  myTag = TAG_ALGORITHM_LINEAR;
}
LinearAlgorithm::~LinearAlgorithm() {
}
int LinearAlgorithm::SolveStep(int /*n*/) {
  pA->get_control()->formTangent();
  pA->get_convergence_norm()->newStep();
  pA->get_control()->predict();
  pA->get_control()->formResidual(pA->get_control()->get_lambda());
  pA->get_convergence_norm()->update();
  return 0;
}
