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

#include "model/standard_model_node.h"
#include "node/node.h"

/**
 * Default constructor.
 */
StandardModelNode::StandardModelNode() {
}
/**
 * Constructor.
 * Initializes the ModelObject, which in turn initializes the FEObject, passes
 * the FTable to the ModelObject and copies the address of it's node.
 */
StandardModelNode::StandardModelNode(const IDContainer& FTable, Node* pNode)
  :ModelNode(FTable, pNode) {
}
StandardModelNode::~StandardModelNode() {
}
void StandardModelNode::incTrialDisp(const Vector& du) {
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=du[theFTable[i]];
  myNode->incTrialDisp(*myVector);
}
void StandardModelNode::incTrialVecs(const Vector& du, const Vector& dv,
                                     const Vector& da) {
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=du[theFTable[i]];
  myNode->incTrialDisp(*myVector);
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=dv[theFTable[i]];
  myNode->addTrialVelc(*myVector);
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=da[theFTable[i]];
  myNode->addTrialAccl(*myVector);
}
void StandardModelNode::set_trial_disp(const Vector& u) {
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=u[theFTable[i]];
  myNode->set_trial_disp(*myVector);
}
void StandardModelNode::set_trial_vecs(const Vector& u, const Vector& v,
                                     const Vector& a) {
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=u[theFTable[i]];
  myNode->set_trial_disp(*myVector);
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=v[theFTable[i]];
  myNode->set_trial_velc(*myVector);
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=a[theFTable[i]];
  myNode->set_trial_accl(*myVector);
}
void StandardModelNode::Commit() {
  myNode->Commit();
}
void StandardModelNode::commitSens(const Vector& X, int param) {
  for (unsigned i = 0; i < theFTable.size(); i++)
    (*myVector)[i]=X[theFTable[i]];
  myNode->commitSens(*myVector, param);
}
void StandardModelNode::add_R(double factor) {
  if (myNode->existsLoad()) myVector->Add_cV(factor, myNode->get_R());
}
