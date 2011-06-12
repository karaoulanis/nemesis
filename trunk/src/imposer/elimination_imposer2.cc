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

#include "imposer/elimination_imposer2.h"
#include <map>
#include <utility>
#include <vector>
#include "analysis/analysis.h"
#include "constraints/constraint.h"
#include "domain/domain.h"
#include "model/elimination_model_element.h"
#include "model/elimination_model_node.h"
#include "model/model.h"
#include "model/standard_model_element.h"
#include "model/standard_model_node.h"
#include "node/node.h"

// #include <iostream>
// using std::cout;
// using std::endl;

EliminationImposer2::EliminationImposer2()
    : Imposer(),
      theNewDofs(0) {
  myTag = TAG_IMPOSER_ELIMINATION;
}


EliminationImposer2::~EliminationImposer2() {
}


int EliminationImposer2::impose() {
  // Check if constrains are already imposed
  if (theModel->isConstrained()) {
    return 0;
  }
  theModel->clear();

  // Create a map to hold node ids and freedom tables.
  // e.g. {1:[0, 1], 2:[2, 3], ...}
  std::map<int, std::vector<int> > ftables;

  // Step 1.: Run through nodes and collect node ids and active dofs.
  for (NodeIterator ni = pA->get_domain()->get_nodes().begin();
            ni != pA->get_domain()->get_nodes().end(); ni++) {
    Node* node = ni->second;
    // Get vector<int> of active dofs. Inactive dofs are equal to -1.
    // e.g. [0, 1, -1, -1, -1, ...]
    std::vector<int> active_dofs = node->get_activated_dofs();
    // Increment all members by one. Inactive dofs should be 0 now.
    // This helps later with the constraints being the only negative ones.
    // e.g. [1, 2, 0, 0, 0, ...]
    for (unsigned i = 0; i < active_dofs.size(); i++) {
      active_dofs[i] += 1;
    }
    // Insert to node_tables map.
    // e.g. {node_id:[1, 2, 0, 0, 0, ...]}
    int node_id = node->get_id();
    ftables.insert(std::pair<int, std::vector<int> >(node_id, active_dofs));
  }

  // Step 2.: Check and apply constraints.
  for (ConstraintIterator ci = theConstraints->begin();
            ci != theConstraints->end(); ci++) {
    Constraint* constraint = ci->second;
    // Check constraints
    if (constraint->get_num_cdofs() == 0) {
      continue;
    }
    if (constraint->get_num_cdofs() > 1) {
      throw SException("[nemesis:%d] %s", 9999,
        "Elimination method cannot be used for MultiFreedom constraints.");
    }
    if (fabs(constraint->get_val())>1e-18) {
      throw SException("[nemesis:%d] %s", 9999,
        "Elimination method cannot be used for non homogeneous constraints.");
    }
    // Apply constraints.
    // Constraints should be the only negatives ones.
    // Inactive dofs should be zero.
    // Active dofs should be greater than zero.
    // e.g. {node_id:[1, -1, 0, 0, 0, ...]}
    int constrained_node = constraint->get_cdof(0).node->get_id();
    int constrained_dof  = constraint->get_cdof(0).dof;
    ftables[constrained_node][constrained_dof] = -constraint->get_id();
  }

  // Step 3.: Renumber free dofs.
  // Again constraints should be the only negatives ones.
  // Inactive dofs should be zero.
  // Active dofs should be greater than zero.
  // That is why num_dof is preincremented (i.e. starts from 1).
  int num_dofs = 0;
  for (std::map<int, std::vector<int> >::iterator ti = ftables.begin();
        ti != ftables.end(); ti++) {
    std::vector<int>& ftable = ti->second;
    for (unsigned i = 0; i < ftable.size(); i++) {
      if (ftable[i] > 0) {
        ftable[i] = ++num_dofs;
      }
    }
    // Remove inactive dofs.
    /// @todo Check which is more efficient.
    // table.erase(std::remove(table.begin(), table.end(), 0), table.end());
    ftable.erase(
      std::remove_if(ftable.begin(), ftable.end(), num::is_zero), ftable.end());
    // Decrement all (positive) members by one.
    // Dof numbering should be start from 0 now.
    // e.g. [0, 1, -1]
    for (unsigned i = 0; i < ftable.size(); i++) {
      if (ftable[i] >0) {
        ftable[i] -= 1;
      }
    }
  }
  // Now the number of equations can be set.
  theModel->set_equations(num_dofs);

  // Step 4.: Create model nodes
  for (NodeIterator ni = pA->get_domain()->get_nodes().begin();
            ni != pA->get_domain()->get_nodes().end(); ni++) {
    Node* node = ni->second;
    ModelNode* model_node;
    const std::vector<int>& ftable = ftables[node->get_id()];
    // If constraints are found within the FTable
     if (Containers::all_positive(ftable)) {
       // Create StandardModelNode
      model_node = new StandardModelNode(ftable, node);
    } else {
      // Create EliminationModelNode
      model_node = new EliminationModelNode(ftable, node);
    }
    theModel->addModelNode(model_node);
    /// @todo Remove.
    // print
    // cout << model_node->get_node()->get_id() <<": \t";
    // Containers::vector_print(model_node->get_FTable());
  }

  // Step 5.: Create model elements
  for (ElementIterator ei = pA->get_domain()->get_elements().begin();
              ei != pA->get_domain()->get_elements().end(); ei++) {
    Element* elem = ei->second;
    const std::vector<Node*> nodes = elem->get_nodes();
    // Copy ftable from the first node.
    std::vector<int> ftable = ftables[nodes[0]->get_id()];
    for (unsigned i = 1; i < nodes.size(); i++) {
      // Append ftables from the rest of the nodes
      const std::vector<int> next_ftable = ftables[nodes[i]->get_id()];
      ftable.insert(ftable.end(), next_ftable.begin(), next_ftable.end());
    }
    // Create model element
    ModelElement* model_elem;
    if (Containers::all_positive(ftable)) {
      // Create StandardModelElement
      model_elem = new StandardModelElement(ftable, elem);
    } else {
      // Create EliminationModelElement
      model_elem = new EliminationModelElement(ftable, elem);
    }
    // Add it to the model
    theModel->addModelElement(model_elem);
    /// @todo Remove.
    // print
    // cout << model_elem->get_element()->get_id() <<": \t";
    // Containers::vector_print(model_elem->get_FTable());
  }

  // Set the model as constrained
  theModel->set_constrained(true);

  return 0;
}
