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

#include "algorithm/newton_raphson_full.h"
#include "analysis/analysis.h"
#include "control/control.h"
#include "convergence/convergence_norm.h"
#include "soe/soe.h"

NewtonRaphsonFull::NewtonRaphsonFull() {
  myTag = TAG_ALGORITHM_NEWTON_RAPHSON_FULL;
}
NewtonRaphsonFull::~NewtonRaphsonFull() {
}
int NewtonRaphsonFull::solveStep(int /*n*/) {
  // Predictor phase
  pA->get_control()->formTangent();
  pA->get_control()->predict();
  pA->get_convergence_norm()->newStep();
  pA->get_control()->formResidual(pA->get_control()->get_lambda());

  // Corrector phase
  int check;
  while ((check = pA->get_convergence_norm()->update())>0) {
    pA->get_control()->formTangent();
    pA->get_soe()->solve();
    pA->get_control()->correct();
    pA->get_control()->formResidual(pA->get_control()->get_lambda());
  }
  return check;
}
