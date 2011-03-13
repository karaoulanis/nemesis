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

#include "model/lagrange_model_node.h"
#include "constraints/constraint.h"
#include "node/node.h"

/**
 * Default constructor.
 */
LagrangeModelNode::LagrangeModelNode()
:ModelNode() {
}
/**
 * Constructor.
 */
LagrangeModelNode::LagrangeModelNode(const IDContainer& FTable, Node* pNode,
                                     Constraint* pConstraint)
  :ModelNode(FTable, pNode, pConstraint) {
}
void LagrangeModelNode::add_R(double /*factor*/) {
}
void LagrangeModelNode::incTrialDisp(const Vector& du) {
  myConstraint->inc_trial_force(du[theFTable[0]]);
}
void LagrangeModelNode::incTrialVecs(const Vector& /*du*/, const Vector& /*dv*/,
                                     const Vector& /*da*/) {
//  myConstraint->incTrialDisp(du[theFTable[0]]);
}
void LagrangeModelNode::set_trial_disp(const Vector& /*u*/) {
//  myConstraint->set_trial_disp(u[theFTable[0]]);
}
void LagrangeModelNode::set_trial_vecs(const Vector& /*u*/, const Vector& /*v*/,
                                     const Vector& /*a*/) {
//  myConstraint->set_trial_disp(u[theFTable[0]]);
}
void LagrangeModelNode::commit() {
  myConstraint->commit();
}
