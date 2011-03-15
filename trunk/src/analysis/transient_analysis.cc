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

#include "analysis/transient_analysis.h"
#include "algorithm/linear_algorithm.h"
#include "control/load_control.h"
#include "convergence/convergence_norm.h"
#include "domain/domain.h"
#include "imposer/elimination_imposer.h"
#include "loadcase/loadcase.h"
#include "soe/full_linear_soe.h"
#include "reorderer/reorderer.h"

TransientAnalysis::TransientAnalysis()
  :AnalysisType() {
  myTag = TAG_ANALYSIS_TRANSIENT;
}
bool TransientAnalysis::checkIfAllows(FEObject* /*f*/) {
/*  if ( f->get_tag()==TAG_CONTROL_GENERALIZED_A        ||
    f->get_tag()==TAG_ALGORITHM_LINEAR         ||
    f->get_tag()==TAG_ALGORITHM_NEWTON_RAPHSON_FULL    ||
    f->get_tag()==TAG_ALGORITHM_NEWTON_RAPHSON_MODIFED ||
    f->get_tag()==TAG_ALGORITHM_NEWTON_RAPHSON_INITIAL ||
    f->get_tag()==TAG_ALGORITHM_NEWTON_RAPHSON_PERIODIC  ||
    f->get_tag()==TAG_IMPOSER_ELIMINATION        ||
    f->get_tag()==TAG_IMPOSER_PENALTY          ||
    f->get_tag()==TAG_IMPOSER_LAGRANGE         ||
    f->get_tag()==TAG_SOE_FULL_GENERIC_POSITIVE_DEFINE ||
    f->get_tag()==TAG_SOE_LINEAR_FULL          ||
    f->get_tag()==TAG_SOE_LINEAR_SYMM          ||
    f->get_tag()==TAG_SOE_LINEAR_BAND          ||
    f->get_tag()==TAG_CONVERGENCE_NORM_EUCLIDEAN     ||
    f->get_tag()==TAG_CONVERGENCE_NORM_MAXIMUM     ||
    f->get_tag()==TAG_REORDERER_NONE           ||
    f->get_tag()==TAG_REORDERER_FORWARD_CUTHILL_MCKEE  ||
    f->get_tag()==TAG_REORDERER_REVERSE_CUTHILL_MCKEE  ||
    f->get_tag()==TAG_REORDERER_FORWARD_SLOAN      ||
    f->get_tag()==TAG_REORDERER_REVERSE_SLOAN       )
      return true;
  return false;*/
  return true;
}
int TransientAnalysis::run(LoadCase* loadcase, int num_loadsteps) {
  // Check the imposer
  if (pA->get_imposer() == 0)
    throw SException("[nemesis:%d] %s", 9999, "No imposer has been set.");
  if (!this->checkIfAllows(pA->get_imposer()))
    throw SException("[nemesis:%d] %s", 9999, "Imposer type is incorrect.");
  // Check the control
  if (pA->get_control() == 0)
    throw SException("[nemesis:%d] %s", 9999, "No control has been set.");
  if (!this->checkIfAllows(pA->get_control()))
    throw SException("[nemesis:%d] %s", 9999, "Control type is incorrect.");
  // Check the algorithm
  if (pA->get_algorithm() == 0)
    throw SException("[nemesis:%d] %s", 9999, "No algorithm has been set.");
  if (!this->checkIfAllows(pA->get_algorithm()))
    throw SException("[nemesis:%d] %s", 9999, "Algorithm type is incorrect.");
  // Check the SOE
  if (pA->get_soe() == 0)
    throw SException("[nemesis:%d] %s", 9999, "No soe has been set.");
  if (!this->checkIfAllows(pA->get_soe()))
    throw SException("[nemesis:%d] %s", 9999, "Soe type is incorrect.");
  // Check the Norm
  if (pA->get_convergence_norm() == 0)
    throw SException("[nemesis:%d] %s", 9999, "No norm has been set.");
  if (!this->checkIfAllows(pA->get_convergence_norm()))
    throw SException("[nemesis:%d] %s", 9999, "Norm type is incorrect.");
  // Check the Reorderer
  if ((pA->get_reorderer() != 0) && (!this->checkIfAllows(pA->get_reorderer())))
    throw SException("[nemesis:%d] %s", 9999, "Reorderer type is incorrect.");

  // Create model by applying the constraints
  pA->get_imposer()->impose();

  // Now that model is complete, reorder the model
  if (pA->get_reorderer() != 0) pA->get_reorderer()->reorder();

  // Now that model is complete, the SOE can be initialized
  pA->get_soe()->set_size();

  // Initialize
  loadcase->Initialize();
  pA->get_control()->init();
  pA->get_convergence_norm()->init(loadcase->get_id(), num_loadsteps);
  /// @todo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // pA->get_domain()->keepTrack(pA->get_control()->get_ambda(),
  // pA->get_control()->get_time());

  for (int i = 0; i < num_loadsteps; i++) {
    // Call algorithm to solve step
    int check = pA->get_algorithm()->solveStep(i);

    // Algorithm failed
    if (check < 0) {
      if (check == -1)
        cout << "Warning  : Solution is diverging." << endl;
      else if (check == -2)
        cout << "Warning  : Maximum number of iteration was exceeded." << endl;
      pA->get_control()->rollback();
      break;
    }
    // Algorithm succeeded
    pA->get_control()->commit();
    /// @todo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // pA->get_domain()->keepTrack(pA->get_control()->get_lambda(),
    // pA->get_control()->get_time());
  }
  loadcase->Commit();
  return 0;
}
