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

#include "analysis/static_analysis.h"

StaticAnalysis::StaticAnalysis()
  :AnalysisType() {
  myTag = TAG_ANALYSIS_STATIC;
  // defaults
  pA->set_imposer(new EliminationImposer());
  pA->set_control(new LoadControl(1., 1., 1., 1, 0.5, 0.));
  pA->set_algorithm(new LinearAlgorithm());
  pA->set_soe(new FullLinearSOE());
}
bool StaticAnalysis::checkIfAllows(FEObject* f) {
  if (f->get_tag() == TAG_CONTROL_LOAD                      ||
      f->get_tag() == TAG_CONTROL_DISPLACEMENT              ||
      f->get_tag() == TAG_CONTROL_ARC_LENGTH_SPHERICAL      ||
      f->get_tag() == TAG_CONTROL_ARC_LENGTH_UNP            ||
      f->get_tag() == TAG_ALGORITHM_LINEAR                  ||
      f->get_tag() == TAG_ALGORITHM_NEWTON_RAPHSON_FULL     ||
      f->get_tag() == TAG_ALGORITHM_NEWTON_RAPHSON_MODIFED  ||
      f->get_tag() == TAG_ALGORITHM_NEWTON_RAPHSON_INITIAL  ||
      f->get_tag() == TAG_ALGORITHM_NEWTON_RAPHSON_PERIODIC ||
      f->get_tag() == TAG_IMPOSER_ELIMINATION               ||
      f->get_tag() == TAG_IMPOSER_PENALTY                   ||
      f->get_tag() == TAG_IMPOSER_LAGRANGE                  ||
      f->get_tag() == TAG_SOE_FULL_GENERIC_POSITIVE_DEFINE  ||
      f->get_tag() == TAG_SOE_LINEAR_FULL                   ||
      f->get_tag() == TAG_SOE_LINEAR_SYMM                   ||
      f->get_tag() == TAG_SOE_LINEAR_BAND                   ||
      f->get_tag() == TAG_CONVERGENCE_NORM_EUCLIDEAN        ||
      f->get_tag() == TAG_CONVERGENCE_NORM_MAXIMUM          ||
      f->get_tag() == TAG_REORDERER_NONE                    ||
      f->get_tag() == TAG_REORDERER_FORWARD_CUTHILL_MCKEE   ||
      f->get_tag() == TAG_REORDERER_REVERSE_CUTHILL_MCKEE   ||
      f->get_tag() == TAG_REORDERER_FORWARD_SLOAN           ||
      f->get_tag() == TAG_REORDERER_REVERSE_SLOAN           )
      return true;
  return false;
}
int StaticAnalysis::run(LoadCase* loadcase, int num_loadsteps) {
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
  // Check the Reorderer
  if ((pA->get_reorderer() != 0) && (!this->checkIfAllows(pA->get_reorderer())))
    throw SException("[nemesis:%d] %s", 9999, "Reorderer type is incorrect.");

  // Create model by applying the constraints
  pA->get_imposer()->impose();

  // Now that model is complete, reorder the model
  if (pA->get_reorderer() !=0 ) pA->get_reorderer()->reorder();

  // Now that model is complete, the SOE can be initialized
  pA->get_soe()->set_size();

  // Initialize
  loadcase->Initialize();
  pA->get_control()->init();
  pA->get_convergence_norm()->init(loadcase->get_id(), num_loadsteps);

  int ret = 0;
  for (int i = 0; i < num_loadsteps; i++) {
    // Call algorithm to solve step
    int check = pA->get_algorithm()->solveStep(i);
    // Algorithm failed
    if (check < 0) {
      if (check == -1)
        cout << "Warning  : Solution is diverging." << endl;
      else if (check == -2)
        cout << "Warning  : Maximum number of iteration was exceeded." << endl;
      // pA->get_control()->returnToConverged();
      // pA->get_control()->commit();
      // break;
      ret = -1;
      // return -1;
    }
    // Algorithm succeeded
    pA->get_control()->commit();
  }
  // Finalize
  loadcase->Commit();
  pA->get_model()->set_nodal_stress();
  pA->get_domain()->commit();  /// @todo Commiting time; find more elegant way.
  return ret;
}
