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

#include "analysis/sensitivity_static_analysis.h"
#include "algorithm/linear_algorithm.h"
#include "control/load_control.h"
#include "control/sensitivity_control.h"
#include "convergence/convergence_norm.h"
#include "imposer/elimination_imposer.h"
#include "loadcase/loadcase.h"
#include "reorderer/reorderer.h"
#include "soe/full_linear_soe.h"

SensitivityStaticAnalysis::SensitivityStaticAnalysis()
  :AnalysisType() {
  myTag = TAG_ANALYSIS_STATIC;
  // defaults
  pA->set_imposer(new EliminationImposer());
  pA->set_control(new LoadControl(1., 1., 1., 1, 0.5, 0.));
  pA->set_algorithm(new LinearAlgorithm());
  pA->set_soe(new FullLinearSOE());
  theSensitivityControl = new SensitivityControl;
}

SensitivityStaticAnalysis::~SensitivityStaticAnalysis() {
  delete theSensitivityControl;
}

bool SensitivityStaticAnalysis::checkIfAllows(FEObject* /*f*/) {
  return false;
}

int SensitivityStaticAnalysis::run(LoadCase* loadcase, int num_loadsteps) {
  // Create model by applying the constraints
  pA->get_imposer()->impose();

  // Now that model is complete, reorder the model
  if (pA->get_reorderer() != 0) pA->get_reorderer()->reorder();

  // Now that model is complete, the SOE can be initialized
  pA->get_soe()->set_size();

  // Apply the loads carried by the given loadcase
  //  pA->get_control()->set_loadcase(pA->get_domain()->get<LoadCase>(
  //    pA->get_domain()->get_loadcases(), nLC));
  //  theSensitivityControl->set_loadcase(pA->get_domain()->get<LoadCase>(
  //    pA->get_domain()->get_loadcases(), nLC));

  // Initialize control
  pA->get_control()->init();
  theSensitivityControl->init();

  // Initialize the convergence check
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
      ret=-1;
    }
    // Algorithm succeeded
    pA->get_control()->commit();

    // Sensitivity
    int nParams =loadcase->GetNumSensitivityParameters();
    for (int j = 0; j < nParams; j++) {
      theSensitivityControl->formTangent();
      theSensitivityControl->formResidual(0.);
      // pA->get_soe()->print();
      pA->get_soe()->solve();
      // cout << pA->get_soe()->get_X();
      theSensitivityControl->commit();
    }
  }

  // Finalize loadcase
  pA->get_model()->set_nodal_stress();
  return ret;
}
