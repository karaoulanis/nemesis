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

/**
 * Default constructor.
 */
LoadCase::LoadCase() {
}
/**
 * Constructor.
 */
LoadCase::LoadCase(int ID, const char* label)
:DomainObject(ID) {
  if (!strcmp(label, "default")) {
    sprintf(myLabel, "LC_%04d", ID);
  } else {
    sprintf(myLabel, "%s", label);
  }
  active=-1;
  myFac = 0;
}
/**
 * Destructor.
 */
LoadCase::~LoadCase() {
  Containers::vector_delete(myLoads);
  Containers::vector_delete(myGroupStates);
  Containers::vector_delete(myInitialConditions);
  Containers::vector_delete(mySensitivityParameters);
}
void LoadCase::init() {
  if (active == 1) return;
  active = 0;
  pD->zeroGroups();
  // Apply group states
  for (unsigned i = 0;i < myGroupStates.size();i++) myGroupStates[i]->apply();
  // Apply initial conditions
  for (unsigned int i = 0;i < myInitialConditions.size();i++)
    myInitialConditions[i]->apply();
  // Increase domain time
}
/**
 * Add a Load to the Loadcase.
 * The Load object should be first added to the Domain. Then a 
 * pointer is passed to the corresponding LoadCase.
 */
void LoadCase::addLoad(Load* pLoad) {
  myLoads.push_back(pLoad);
}
/**
 * Add an Initial Condition to the Loadcase.
 * The InitialCondition object should be first added to the Domain. Then a 
 * pointer is passed to the corresponding LoadCase.
 */
void LoadCase::addInitialCondition(InitialCondition* pInitialCondition) {
  myInitialConditions.push_back(pInitialCondition);
}
/**
 * Add a Group State to the Loadcase.
 * The GroupState object should be first added to the Domain. Then a 
 * pointer is passed to the corresponding LoadCase.
 */
void LoadCase::addGroupState(GroupState* pGroupState) {
  myGroupStates.push_back(pGroupState);
}
/**
 * Apply nodal and elemental loads.
 */
void LoadCase::applyLoads(double lambda_, double time_) {
  if (lambda_ == 0.) return;
  if (active  ==-1)  return;
  else if (active == 0) myFac = lambda_;
  // Apply nodal and elemental loads
  for (unsigned int i = 0;i < myLoads.size();i++)
    myLoads[i]->apply(myFac, time_);
}
/**
 * Commit.
 */
void LoadCase::commit() {
  active = 1;
  pD->commit();
}
void LoadCase::applySensitivityParameter(int param) {
  mySensitivityParameters[param]->apply();
}
int LoadCase::getnSensitivityParameters() {
  return mySensitivityParameters.size();
}
void LoadCase::addSensitivityParameter(ElementSensitivityParameter* pElementSensitivityParameter) {
  mySensitivityParameters.push_back(pElementSensitivityParameter);
}
