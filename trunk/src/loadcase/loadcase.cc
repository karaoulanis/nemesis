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

#include "loadcase/loadcase.h"
#include <stdio.h>
#include <string.h>
#include "containers/containers.h"
#include "loadcase/element_sensitivity_parameter.h"
#include "loadcase/group_state.h"
#include "loadcase/initial_condition.h"
#include "loadcase/load.h"

/**
 * Default constructor.
 */
LoadCase::LoadCase()
    : factor_(0.),
      applied_(false),
      active_(false),
      loads_(0),
      groupstates_(0),
      initialconditions_(0),
      sensitivityparameters_(0) {
}

/**
 * Constructor.
 */

LoadCase::LoadCase(int id, const char* label)
    : DomainObject(id),
      factor_(0.),
      applied_(false),
      active_(false),
      loads_(0),
      groupstates_(0),
      initialconditions_(0),
      sensitivityparameters_(0) {
  if (!strcmp(label, "default")) {
    sprintf(label_, "LC_%04d", id);
  } else {
    sprintf(label_, "%s", label);
  }
}

/**
 * Destructor.
 */
LoadCase::~LoadCase() {
  Containers::vector_delete(loads_);
  Containers::vector_delete(groupstates_);
  Containers::vector_delete(initialconditions_);
  Containers::vector_delete(sensitivityparameters_);
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
void LoadCase::AddInitialCondition(InitialCondition* initialcondition) {
  initialconditions_.push_back(initialcondition);
}

/**
 * Add a Group State to the Loadcase.
 * The GroupState object should be first added to the Domain. Then a
 * pointer is passed to the corresponding LoadCase.
 */
void LoadCase::AddGroupState(GroupState* groupstate) {
  groupstates_.push_back(groupstate);
}

void LoadCase::Initialize() {
  // Quick return if loadcase is already applied or is now active
  /// @todo Check if this conditional is really necessary.
  if (active_ || applied_) return;
  // First time encountered so this is now active
  active_ = true;
  for (unsigned i = 0; i < groupstates_.size(); i++) {
    groupstates_[i]->Apply();
  }
  // Apply initial conditions
  for (unsigned int i = 0;i < initialconditions_.size();i++) {
    initialconditions_[i]->Apply();
  }
  // Increase domain time
  /// @todo
}

void LoadCase::ApplyLoads(double lambda, double time) {
  if (lambda  ==  0.) return;
  if (!active_ && !applied_)  return;
  else if (active_ == true) factor_ = lambda;
  // Apply nodal and elemental loads
  for (unsigned int i = 0; i < loads_.size(); i++) {
    loads_[i]->Apply(factor_, time);
  }
}

void LoadCase::Finalize() {
  applied_ = true;
  active_  = false;
}

void LoadCase::ApplySensitivityParameter(int param) {
  sensitivityparameters_[param]->apply();
}
int LoadCase::GetNumSensitivityParameters() {
  return sensitivityparameters_.size();
}
void LoadCase::AddSensitivityParameter(ElementSensitivityParameter*
                                       pElementSensitivityParameter) {
  sensitivityparameters_.push_back(pElementSensitivityParameter);
}

void LoadCase::Save(std::ostream* /*os*/) {
  /// @todo Implement this method.
}

