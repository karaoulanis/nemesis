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

#include "analysis/xfem_analysis.h"

XFemAnalysis::XFemAnalysis()
  :AnalysisType() {
  myTag = TAG_ANALYSIS_STATIC;
  // defaults
  pA->setImposer(new EliminationImposer());
  pA->setControl(new LoadControl(1., 1., 1., 1, 0.5, 0.));
  pA->setAlgorithm(new LinearAlgorithm());
  pA->setSOE(new FullLinearSOE());
}
bool XFemAnalysis::checkIfAllows(FEObject* /*f*/) {
  return true;
}
int XFemAnalysis::run(int nLC, int nLoadSteps) {
  ///@todo Checks

  // Create model by applying the constraints
  pA->getImposer()->impose();

  // Now that model is complete, reorder the model
  if (pA->getReorderer() !=0 ) pA->getReorderer()->reorder();

  // Now that model is complete, the SOE can be initialized
  pA->getSOE()->setTheSize();

  // Initialize
  pA->getDomain()->get<LoadCase>(pA->getDomain()->getLoadCases(), nLC)->init();
  pA->getControl()->init();
  pA->getConvergenceNorm()->init(nLC, nLoadSteps);

  // Call algorithm to solve step
  // int check = pA->getAlgorithm()->solveStep(0); // gives warning for check
  pA->getAlgorithm()->solveStep(0);  ///@todo:check

  // Commit
  pA->getControl()->commit();

  // Finalize
  pA->getDomain()->get<LoadCase>(pA->getDomain()->getLoadCases(), nLC)->commit();
  pA->getModel()->setNodalStress();
  pA->getDomain()->commit(); ///@todo For commiting time; find more elegant way.
  return 0;
}
