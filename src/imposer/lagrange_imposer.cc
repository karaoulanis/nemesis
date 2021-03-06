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

#include "imposer/lagrange_imposer.h"
#include <map>
#include <utility>
#include <vector>
#include "constraints/constraint.h"
#include "model/lagrange_model_element.h"
#include "model/lagrange_model_node.h"
#include "model/model.h"
#include "model/standard_model_element.h"
#include "model/standard_model_node.h"
#include "node/node.h"

LagrangeImposer::LagrangeImposer()
    : Imposer() {
}


LagrangeImposer::LagrangeImposer(const std::map<int, Node*>& nodes,
                                 const std::map<int, Element*>& elements,
                                 const std::map<int, Constraint*>& constraints)
    : Imposer(nodes, elements, constraints) {
}


int LagrangeImposer::impose(Model* model) {
  // Check if constrains are already imposed
  if (model->isConstrained()) {
    return 0;
  }
  model->clear();

  // Create a map to hold node ids and freedom tables.
  // e.g. {1:[0, 1], 2:[2, 3], ...}
  std::map<int, std::vector<int> > ftables;

  // Step 1.: Run through nodes and collect node ids and active dofs.
  //          Number (global) these dofs.
  //          Store node freedom tables.
  //          Set number of equations.
  int num_dofs = 0;
  for (std::map<int, Node*>::const_iterator ni = nodes_->begin();
                                            ni != nodes_->end();
                                            ++ni) {
    Node* node = ni->second;
    // Get vector<int> of active dofs. Inactive dofs are equal to -1.
    // e.g. [0, 1, -1, -1, -1, ...]
    std::vector<int> ftable = node->get_activated_dofs();
    // Now assign global dofs,
    // e.g. [1, 2, -1, -1, -1, ...]
    for (unsigned i = 0; i < ftable.size(); i++) {
      if (ftable[i] >= 0) {
        ftable[i] = num_dofs++;
      }
    }
    // Insert to ftables map for later reference (for Model Elements).
    // e.g. {node_id:[1, 2, -1, -1, -1, ...]}
    int node_id = node->get_id();
    ftables.insert(std::pair<int, std::vector<int> >(node_id, ftable));
  }

  // Step 2.: Check and apply constraints.
  //          Create Lagrange model nodes and elements.
  //          Set the number of equations.
  for (std::map<int, Constraint*>::const_iterator ci = constraints_->begin();
                                                  ci != constraints_->end();
                                                  ++ci) {
    Constraint* constraint = ci->second;
    // Check constraints
    if (constraint->get_num_cdofs() == 0) {
      continue;
    }
    // Create model nodes
    std::vector<int> node_ftable = std::vector<int>(1);
    node_ftable[0] = num_dofs;
    LagrangeModelNode* model_node =
        new LagrangeModelNode(node_ftable, 0, constraint);
    model->addModelNode(model_node);
    // Create model elements
    int num_constrained_dofs  = constraint->get_num_cdofs();
    std::vector<int> elem_ftable = std::vector<int>(num_constrained_dofs+1);
    for (int i = 0; i < num_constrained_dofs; i++) {
      int constrained_node = constraint->get_cdof(i).node->get_id();
      int constrained_dof  = constraint->get_cdof(i).dof;
      elem_ftable[i] = ftables[constrained_node][constrained_dof];
    }
    elem_ftable[num_constrained_dofs] = num_dofs;
    LagrangeModelElement* model_element =
        new LagrangeModelElement(elem_ftable, constraint);
    model->addModelElement(model_element);
    // Increase the number of equations
    num_dofs++;
  }
  // Now the number of equations can be set.
  model->set_equations(num_dofs);

  // Step 3.: Remove inactive dofs.
  //          Create standard model nodes
  for (std::map<int, Node*>::const_iterator ni = nodes_->begin();
                                            ni != nodes_->end();
                                            ++ni) {
    Node* node = ni->second;
    std::vector<int>& ftable = ftables[node->get_id()];
    // Remove inactive dofs.
    /// @todo Check which is more efficient.
    // table.erase(std::remove(table.begin(), table.end(), 0), table.end());
    ftable.erase(std::remove_if(ftable.begin(),
                                 ftable.end(), num::is_negative), ftable.end());
    // Create model nodes.
    StandardModelNode* model_node = new StandardModelNode(ftable, node);
    model->addModelNode(model_node);
  }

  // Step 4.: Create standard model elements.
  for (std::map<int, Element*>::const_iterator ei = elements_->begin();
                                               ei != elements_->end();
                                               ++ei) {
    Element* elem = ei->second;
    const std::vector<Node*> nodes = elem->get_nodes();
    /*// Copy ftable from the first node.
    std::vector<int> ftable = ftables[nodes[0]->get_id()];
    for (unsigned i = 1; i < nodes.size(); i++) {
      // Append ftables from the rest of the nodes
      const std::vector<int> next_ftable = ftables[nodes[i]->get_id()];
      ftable.insert(ftable.end(), next_ftable.begin(), next_ftable.end());
    }*/
    std::vector<int> ftable;
    const std::vector<int> local_dofs = elem->get_local_nodal_dofs();
    for (unsigned i = 0; i < nodes.size(); i++) {
      for (unsigned j = 0; j < local_dofs.size(); j++) {
        ftable.push_back(ftables[nodes[i]->get_id()][local_dofs[j]]);
      }
    }
    // Create model element
    StandardModelElement* model_elem = new StandardModelElement(ftable, elem);
    // Add it to the model
    model->addModelElement(model_elem);
  }

  // Set the model as constrained
  model->set_constrained(true);
  model->print();
  return 0;
}
