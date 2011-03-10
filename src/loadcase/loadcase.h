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

#ifndef SRC_LOADCASE_LOADCASE_H_
#define SRC_LOADCASE_LOADCASE_H_

#include <vector>
#include <map>

// Project files (alphabetically)
#include "containers/containers.h"
#include "domain/domain_object.h"
#include "loadcase/elemental_load.h"
#include "loadcase/element_sensitivity_parameter.h"
#include "loadcase/group_state.h"
#include "loadcase/initial_condition.h"
#include "loadcase/load.h"
#include "loadcase/nodal_load.h"

// Forward declarations
class ElementSensitivityParameter;
class GroupState;
class Group;
class InitialCondition;
class Load;
class NodalLoad;

// Type definitions
typedef std::vector<Load*>                        LoadVector;
typedef std::vector<InitialCondition*>            InitialConditionVector;
typedef std::vector<GroupState*>                  GroupStateVector;
typedef std::vector<ElementSensitivityParameter*> SensitivityParameterVector;
typedef std::map<int, Group*>                     GroupContainer;

class LoadCase: public DomainObject {
 public:
  // Constructors and destructor
  LoadCase();
  LoadCase(int id, const char* label);
  ~LoadCase();
  // Initialize loadcase
  void Initialize();
  void Commit();
  // Add member functions
  void AddLoad(Load* pLoad);
  void AddGroupState(GroupState* pGroupState);
  void AddInitialCondition(InitialCondition* pInitialCondition);
  void AddSensitivityParameter(ElementSensitivityParameter* pElementSensitivityParameter);
  // Apply member functions
  void ApplyLoads(double lambda_, double time_);
  void ApplySensitivityParameter(int param);
  int GetNumSensitivityParameters();

 private:
  char label_[512];
  double factor_;
  int active_;
  LoadVector                  loads_;
  GroupStateVector            group_states_;
  InitialConditionVector      initial_conditions_;
  SensitivityParameterVector  sensitivity_parameters_;
};
#endif  // SRC_LOADCASE_LOADCASE_H_
