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

#include "model/model_node.h"
#include "node/node.h"

/**
 * Default constructor.
 */
ModelNode::ModelNode()
    : myNode(0),
      myConstraint(0) {
}
/**
 * Constructor.
 * Initializes the ModelObject, which in turn initializes the FEObject, passes
 * the FTable to the ModelObject and copies the address of it's node.
 */
ModelNode::ModelNode(const IDContainer& FTable, Node* node,
                     Constraint* constraint)
    : ModelObject(FTable),
      myNode(node),
      myConstraint(constraint) {
  // Set a pointer to the corresponding static Vector
  myVector = theStaticVectors[FTable.size()];
}
ModelNode::~ModelNode() {
}
void ModelNode::add_uTrial(double factor) {
  if (myNode != 0) myVector->add_cV(factor, myNode->get_disp_trial());
}
void ModelNode::add_vTrial(double factor) {
  if (myNode != 0) myVector->add_cV(factor, myNode->get_velc_trial());
}
void ModelNode::rollback() {
  if (myNode != 0) myNode->rollback();
}
