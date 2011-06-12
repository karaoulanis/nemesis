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

#include "analysis/eigen_analysis.h"
#include "analysis/analysis.h"
#include "control/eigen_control.h"
#include "domain/domain.h"
#include "imposer/elimination_imposer.h"
#include "reorderer/reorderer.h"
#include "soe/eigen_soe.h"

EigenAnalysis::EigenAnalysis()
  :AnalysisType() {
  myTag = TAG_NONE;
  // defaults
  pA->set_imposer(new EliminationImposer());
  pA->set_control(new EigenControl());
  pA->set_soe(new EigenSOE(pA->get_model()));
}
EigenAnalysis::~EigenAnalysis() {
}
bool EigenAnalysis::checkIfAllows(FEObject* /*f*/) {
  return false;
}
int EigenAnalysis::run(LoadCase* /*loadcase*/, int /*num_loadsteps*/) {
  // Create model by applying the constraints
  pA->get_imposer()->impose(pA->get_model());

  // Now that model is complete, reorder the model
  if (pA->get_reorderer() !=0 ) pA->get_reorderer()->reorder();

  // Now that model is complete, the SOE can be initialized
  pA->get_soe()->set_size();

  // Initialize the control
  pA->get_control()->formTangent();
  //  pA->get_soe()->print();
  pA->get_soe()->solve();
  pA->get_domain()->set_eigenvalues(pA->get_soe()->get_X());

  return 0;
}
