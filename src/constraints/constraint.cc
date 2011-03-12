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

// Included files
#include "constraints/constraint.h"

// Initialize variables
int Constraint::num_constraints_ = 0;

/**
 * Constructor.
 */
Constraint::Constraint()
  :DomainObject(++num_constraints_), val_(0), f_trial_(0), f_convg_(0) {
}
/**
 * Create a constrained dof.
 * Check if \a node_id and \a dof is valid; otherwise throw an exception.
 * @param node_id The nodal id to be constrained.
 * @param dof The nodal dof to be constrained.
 * @param coeff The Constraint's coefficient.
 */
void Constraint::set_cdof(int node_id, int dof, double coeff) {
  static CDof cdof;
  // Get node; if does not exists re-throw the exception
  Node* pNode;
  try {
    pNode = pD->get<Node>(pD->get_nodes(), node_id);
  } catch(SException) {
    throw;
  }
  // Check if dof is active
  if (pNode->get_activated_dof(dof-1) < 0) {
    throw SException("[nemesis:%d] Node %d dof %d is not yet active.\n", 1110,
      pNode->get_id(), dof);
    // Containers::vector_print(pNode->get_connected_elements());cout << endl;
    // return;
  }
  // Now it is ok to continue
  cdof.pNode_ = pNode;
  cdof.dof_   = dof-1;
  // Add the coeffiecient
  cdof.coeff_ = coeff;
  // Keep constrained dof
  cdofs_.push_back(cdof);
}
const CDof& Constraint::get_cdof(int i) {
  return cdofs_[i];
}
void Constraint::set_val(double val) {
  val_ = val;
}
double Constraint::get_val(double /*time*/) {
  return val_;
}
int Constraint::get_num_cdofs() {
  return cdofs_.size();
}
double Constraint::get_disp(int i) {
  return cdofs_[i].pNode_->get_disp_trial_at_dof(cdofs_[i].dof_);
}
double Constraint::get_velc(int i) {
  return cdofs_[i].pNode_->get_velc_trial_at_dof(cdofs_[i].dof_);
}
double Constraint::get_disp_convg(int i) {
  return cdofs_[i].pNode_->get_disp_convg_at_dof(cdofs_[i].dof_);
}
double Constraint::get_disp_trial(int i) {
  return cdofs_[i].pNode_->get_disp_trial_at_dof(cdofs_[i].dof_);
}
double Constraint::get_velc_convg(int i) {
  return cdofs_[i].pNode_->get_velc_convg_at_dof(cdofs_[i].dof_);
}
double Constraint::get_accl_convg(int i) {
  return cdofs_[i].pNode_->get_accl_convg_at_dof(cdofs_[i].dof_);
}
void Constraint::inc_trial_force(double f) {
  f_trial_+=f;
}
void Constraint::update(double f)  {
  f_trial_+=f;
}
void Constraint::commit() {
  f_convg_ = f_trial_;
}
double Constraint::get_F() {
  return f_trial_;
}
