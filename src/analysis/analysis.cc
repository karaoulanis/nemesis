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

#include "analysis/analysis.h"
#include "algorithm/algorithm.h"
#include "analysis/analysis_type.h"
#include "control/control.h"
#include "convergence/convergence_norm.h"
#include "domain/domain.h"
#include "imposer/imposer.h"
#include "reorderer/reorderer.h"
#include "soe/soe.h"

Analysis::Analysis(Domain* pDomain)
:M(pDomain), theDomain(pDomain) {
  M.set_analysis(this);
  theNorm = new ConvergenceNorm();
  theAnalysisType = 0;
  theImposer = 0;
  theControl = 0;
  theAlgorithm = 0;
  theSOE = 0;
  theReorderer = 0;
}
Analysis::~Analysis() {
  delete theNorm;
  this->clear();
}
int Analysis::analyze(LoadCase* loadcase, int num_loadsteps) {
  if (theAnalysisType == 0)
    throw SException("[nemesis:%d] %s", 9999, "No analysis type set.");
  // Run the analysis
  return theAnalysisType->run(loadcase, num_loadsteps);
}
void Analysis::clear() {
  M.clear();
  if (theAnalysisType != 0) {
    delete theAnalysisType;
    theAnalysisType = 0;
  }
  if (theImposer != 0) {
    delete theImposer;
    theImposer = 0;
  }
  if (theControl != 0) {
    delete theControl;
    theControl = 0;
  }
  if (theAlgorithm != 0) {
    delete theAlgorithm;
    theAlgorithm = 0;
  }
  if (theSOE != 0) {
    delete theSOE;
    theSOE = 0;
  }
  if (theReorderer != 0) {
    delete theReorderer;
    theReorderer = 0;
  }
}
void Analysis::set_analysis_type(AnalysisType* p) {
  if (theAnalysisType != 0) delete theAnalysisType;
  theAnalysisType = p;
}
void Analysis::set_algorithm(Algorithm* p) {
  if (theAlgorithm != 0) delete theAlgorithm;
  theAlgorithm = p;
}
void Analysis::set_control(Control* p) {
  if (theControl != 0) delete theControl;
  theControl = p;
}
void Analysis::set_imposer(Imposer* p) {
  if (theImposer != 0) delete theImposer;
  theImposer = p;
}
void Analysis::set_reorderer(Reorderer* p) {
  if (theReorderer != 0) delete theReorderer;
  theReorderer = p;
}
void Analysis::set_soe(SOE* p) {
  if (theSOE != 0) delete theSOE;
  theSOE = p;
}
