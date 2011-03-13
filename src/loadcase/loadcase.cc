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

#include "loadcase/loadcase.h"
#include <stdio.h>
#include "loadcase/element_sensitivity_parameter.h"
#include "loadcase/group_state.h"
#include "loadcase/initial_condition.h"
#include "loadcase/load.h"

/**
 * Default constructor.
 */
LoadCase::LoadCase() {
}

/**
 * Constructor.
 */

LoadCase::LoadCase(int id, const char* label)
:DomainObject(id) {
  if (!strcmp(label, "default")) {
    sprintf(myLabel, "LC_%04d", id);
  } else {
    sprintf(myLabel, "%s", label);
  }
  applied_ = false;
  active_  = false;
  factor_  = 0.;
}

/**
 * Destructor.
 */
LoadCase::~LoadCase() {
  Containers::vector_delete(loads_);
  Containers::vector_delete(group_states_);
  Containers::vector_delete(initial_conditions_);
  Containers::vector_delete(sensitivity_parameters_);
}

/**
 * Add a Load to the Loadcase.
 * The Load object should be first added to the Domain. Then a
 * pointer is passed to the corresponding LoadCase.
 */
void LoadCase::AddLoad(Load* load) {
  loads_.push_back(load);
}

/**
 * Add an Initial Condition to the Loadcase.
 * The InitialCondition object should be first added to the Domain. Then a
 * pointer is passed to the corresponding LoadCase.
 */
void LoadCase::AddInitialCondition(InitialCondition* initial_condition) {
  initial_conditions_.push_back(initial_condition);
}

/**
 * Add a Group State to the Loadcase.
 * The GroupState object should be first added to the Domain. Then a
 * pointer is passed to the corresponding LoadCase.
 */
void LoadCase::AddGroupState(GroupState* group_state) {
  group_states_.push_back(group_state);
}

void LoadCase::Initialize() {
  // Quick return if loadcase is already applied or is now active
  /// @todo Check if this conditional is really necessary.
  if (active_ || applied_) return;
  // First time encountered so this is now active
  active_ = true;
  for (unsigned i = 0; i < group_states_.size(); i++) {
    group_states_[i]->Apply();
  }
  // Apply initial conditions
  for (unsigned int i = 0;i < initial_conditions_.size();i++) {
    initial_conditions_[i]->Apply();
  }
  // Increase domain time
  /// @todo
}

void LoadCase::ApplyLoads(double lambda, double time) {
  if (lambda  ==  0.) return;
  if (!active_ && !applied_)  return;
  else if (active_ == 0) factor_ = lambda;
  // Apply nodal and elemental loads
  for (unsigned int i = 0; i < loads_.size(); i++) {
    loads_[i]->Apply(factor_, time);
  }
}

void LoadCase::Commit() {
  applied_ = true;
  active_  = false;
}

void LoadCase::ApplySensitivityParameter(int param) {
  sensitivity_parameters_[param]->apply();
}
int LoadCase::GetNumSensitivityParameters() {
  return sensitivity_parameters_.size();
}
void LoadCase::AddSensitivityParameter(ElementSensitivityParameter* pElementSensitivityParameter) {
  sensitivity_parameters_.push_back(pElementSensitivityParameter);
}
