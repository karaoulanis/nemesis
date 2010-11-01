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
* along with this program.  If not, see < http://www.gnu.org/licenses/>.        *
*******************************************************************************/

// *****************************************************************************
// $LastChangedDate$
// $LastChangedRevision$
// $LastChangedBy$
// $HeadURL$
// Author(s): F.E. Karaoulanis (fkar@nemesis-project.org)
// *****************************************************************************

#include "imposer/lagrange_imposer.h"

LagrangeImposer::LagrangeImposer()
  :Imposer() {
  myTag = TAG_IMPOSER_LAGRANGE;
}

int LagrangeImposer::impose() {
  if (theModel->isConstrained()) return 0;

  Node* pNode;                // Pointer to a domain node
  Element* pElement;              // Pointer to a domain element
  ModelNode* pStdModelNode;         // Pointer to a model node
  StandardModelElement* pStdModelElement;   // Pointer to a model element
  LagrangeModelElement* pLagModelElement;
  LagrangeModelNode* pLagModelNode;
  Constraint* pConstraint;

  //=========================================================================
  // Find nodal global numbering and store it
  //=========================================================================
  int nDofs = this->createGlobalDofNumbering();
  theModel->setEquations(nDofs);

  //=========================================================================
  // Create Standard ModelNodes
  //=========================================================================
  for (NodeIterator nIter = pA->getDomain()->getNodes().begin();
    nIter != pA->getDomain()->getNodes().end();nIter++) 
  {
    // Get a pointer to a node
    pNode = nIter->second;
    // Get the nodal id
    int nodeID = pNode->getID();
    IDContainer nodalFTable = this->getGlobalDofs(nodeID);
    pStdModelNode = new StandardModelNode(nodalFTable, pNode);
    theModel->addModelNode(pStdModelNode);
  }
  //=========================================================================
  // Create Standard ModelElements
  //=========================================================================
    for (ElementIterator eIter = pA->getDomain()->getElements().begin();
    eIter != pA->getDomain()->getElements().end();eIter++) 
  {
    // Get next (randomly chosen) element
    pElement = eIter->second;
    // Get the ids of the nodes
    IDContainer myNodalIDs = pElement->getNodalIDs();
    // Get the local dofs of each node of the element
    IDContainer theNodalLocalDofs = pElement->getLocalNodalDofs();
    // Create an element freedom table
    IDContainer elemFTable(myNodalIDs.size()*theNodalLocalDofs.size());

    int nElemNodes = myNodalIDs.size();
    for (int j = 0; j < nElemNodes; j++) {
      int NodeID = myNodalIDs[j];
      for (unsigned k = 0; k < theNodalLocalDofs.size(); k++) {
        int localDof = theNodalLocalDofs[k];
        int globalDof = this->getGlobalDof(NodeID, localDof);
        elemFTable[j*theNodalLocalDofs.size()+k]=globalDof;
      }
    }
    pStdModelElement = new StandardModelElement(elemFTable, pElement);
    theModel->addModelElement(pStdModelElement);
  }
  //=========================================================================
  // Constraints
  //=========================================================================
  int nEquations = theModel->getnEquations();
  for (ConstraintIterator cIter = theConstraints->begin();
    cIter != theConstraints->end();cIter++)
  {
    pConstraint = cIter->second;
    if (pConstraint->getncDofs()==0) continue;
    // Lagrange ModelNode
    static IDContainer nodalFTable(1);
    nodalFTable[0]=nEquations;
    pLagModelNode = new LagrangeModelNode(nodalFTable, 0, pConstraint);
    theModel->addModelNode(pLagModelNode);
    // Lagrange ModelElement
    int ncDofs = pConstraint->getncDofs();
    IDContainer cFTable(ncDofs+1);
    for (int j = 0;j < ncDofs;j++)
      cFTable[j]=this->getGlobalDof(
        pConstraint->getcDof(j).pNode->getID(), pConstraint->getcDof(j).dof);
    cFTable[ncDofs]=nEquations;
    pLagModelElement = new LagrangeModelElement(cFTable, pConstraint);
    theModel->addModelElement(pLagModelElement);
    nEquations++;
  }
  theModel->setEquations(nEquations);
  theModel->setConstrained(true);
  return 0;
}
