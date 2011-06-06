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

#ifndef SRC_LOADCASE_LOADCASE_H_
#define SRC_LOADCASE_LOADCASE_H_

#include <map>
#include <vector>
#include "domain/domain_object.h"

// Forward declarations
class ElementSensitivityParameter;
class GroupState;
class InitialCondition;
class Load;
class NodalLoad;

class LoadCase: public DomainObject {
 public:
  // Constructors and destructor
  LoadCase();
  LoadCase(int id, const char* label);
  ~LoadCase();

  // Do stuff at the begining/end of a loadcase.
  void Initialize();
  void Finalize();

  // Add member functions
  void AddLoad(Load* pLoad);
  void AddGroupState(GroupState* pGroupState);
  void AddInitialCondition(InitialCondition* pInitialCondition);
  void AddSensitivityParameter(ElementSensitivityParameter*
                                  pElementSensitivityParameter);

  // Apply member functions
  void ApplyLoads(double lambda_, double time_);
  void ApplySensitivityParameter(int param);
  int GetNumSensitivityParameters();

 private:
  char label_[512];
  double factor_;
  bool applied_;
  bool active_;
  std::vector<Load*> loads_;
  std::vector<GroupState*> groupstates_;
  std::vector<InitialCondition*> initialconditions_;
  std::vector<ElementSensitivityParameter*> sensitivityparameters_;
};
#endif  // SRC_LOADCASE_LOADCASE_H_
