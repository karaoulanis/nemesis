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

#include "imposer/elimination_imposer.h"
#include "analysis/analysis.h"
#include "constraints/constraint.h"
#include "domain/domain.h"
#include "model/elimination_model_element.h"
#include "model/elimination_model_node.h"
#include "model/model.h"
#include "model/standard_model_element.h"
#include "model/standard_model_node.h"
#include "node/node.h"

EliminationImposer::EliminationImposer()
  :Imposer() {
  myTag = TAG_IMPOSER_ELIMINATION;
}
EliminationImposer::~EliminationImposer() {
}
int EliminationImposer::impose() {
  if (theModel->isConstrained()) return 0;

  Node* pNode;
  Element* pElement;
  Constraint* pConstraint;
  StandardModelNode* pStdModelNode;
  StandardModelElement* pStdModelElement;
  EliminationModelNode* pElimModelNode;
  EliminationModelElement* pElimModelElement;

  // ---------------------------------------------------------------------------
  // Find nodal global numbering and store it
  // ---------------------------------------------------------------------------
  theModel->clear();
  int nDofs = this->createGlobalDofNumbering();

  // ---------------------------------------------------------------------------
  // Number constraints with global numbering
  // ---------------------------------------------------------------------------
  theNewDofs.resize(nDofs);
  for (int i = 0;i < nDofs;i++) theNewDofs[i]=1;
  for (ConstraintIterator cIter = theConstraints->begin();
              cIter != theConstraints->end(); cIter++) {
    pConstraint = cIter->second;
    if (pConstraint->get_num_cdofs() == 0) continue;
    if (pConstraint->get_num_cdofs() > 1)
      throw SException("[nemesis:%d] %s", 9999,
      "Elimination method cannot be used for MultiFreedom constraints.");
    int globalDof = this->get_global_dof(
      pConstraint->get_cdof(0).node->get_id(), pConstraint->get_cdof(0).dof);
    // Check if dof is activated
    if (globalDof < 0) continue;
    theNewDofs[globalDof]=-pConstraint->get_id();
    if (fabs(pConstraint->get_val())>1e-18)
      throw SException("[nemesis:%d] %s", 9999,
      "Elimination method cannot be used for non homogeneous constraints.");
  }
  int nEquations = 0;
  for (int i = 0;i < nDofs;i++) if (theNewDofs[i]>0) theNewDofs[i]=nEquations++;
  theModel->set_equations(nEquations);

  // ===========================================================================
  // Create ModelNodes
  // ===========================================================================
  for (NodeIterator nIter = pA->get_domain()->get_nodes().begin();
            nIter != pA->get_domain()->get_nodes().end(); nIter++) {
    static IDContainer oldFTable;
    static IDContainer newFTable;
    // Get a pointer to a node
    pNode = nIter->second;
    oldFTable = this->get_global_dofs(pNode->get_id());
    newFTable.resize(oldFTable.size());
    // Map old dofs to new ones
    for (unsigned i = 0;i < oldFTable.size();i++)
      newFTable[i]=theNewDofs[oldFTable[i]];
    // If constraints are found within the FTable
    if (Containers::all_positive(newFTable)) {
      // Create StandardModelNode
      pStdModelNode = new StandardModelNode(newFTable, pNode);
      theModel->addModelNode(pStdModelNode);
    } else {
      // Create EliminationModelNode
      pElimModelNode = new EliminationModelNode(newFTable, pNode);
      pElimModelNode->set_old_FTable(oldFTable);
      theModel->addModelNode(pElimModelNode);
    }
  }
  // ===========================================================================
  // Create ModelElements
  // ===========================================================================
  for (ElementIterator eIter = pA->get_domain()->get_elements().begin();
    eIter != pA->get_domain()->get_elements().end(); eIter++) {
    /// @todo Needs speed improvements
    // Get next (randomly chosen) element
    pElement = eIter->second;
    // Get the ids of the nodes
    IDContainer myNodalIDs = pElement->get_nodal_ids();
    // Get the local dofs of each node of the element
    IDContainer theNodalLocalDofs = pElement->get_local_nodal_dofs();
    // Create an element freedom table
    IDContainer elemFTable(myNodalIDs.size()*theNodalLocalDofs.size());

    for (unsigned j = 0; j < myNodalIDs.size(); j++) {
      int NodeID = myNodalIDs[j];
      for (unsigned k = 0; k < theNodalLocalDofs.size(); k++) {
        int localDof = theNodalLocalDofs[k];
        int globalDof = this->get_global_dof(NodeID, localDof);
        int newGlobalDof = theNewDofs[globalDof];
        elemFTable[j*theNodalLocalDofs.size()+k]=newGlobalDof;
      }
    }
    // Check if constraint found
    if (Containers::all_positive(elemFTable)) {
      // Create StandardModelElement
      pStdModelElement = new StandardModelElement(elemFTable, pElement);
      theModel->addModelElement(pStdModelElement);
    } else {
      // Create EliminationModelElement
      pElimModelElement = new EliminationModelElement(elemFTable, pElement);
      theModel->addModelElement(pElimModelElement);
    }
  }
  theModel->set_constrained(true);
  return 0;
}