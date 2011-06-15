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

#ifndef SRC_MODEL_ELIMINATION_MODEL_NODE_H_
#define SRC_MODEL_ELIMINATION_MODEL_NODE_H_

#include "model/model_node.h"
#include "node/node.h"

class EliminationModelNode: public ModelNode {
 public:
  // Constructors
  EliminationModelNode();
  EliminationModelNode(const IDContainer& FTable, Node* pNode);
  ~EliminationModelNode();
  void add_R(double factor);
  void incTrialDisp(const Vector& du);
  void incTrialVecs(const Vector& du, const Vector& dv, const Vector& da);
  void set_trial_disp(const Vector& u);
  void set_trial_vecs(const Vector& u, const Vector& v, const Vector& a);
  void commit();
  void commitSens(const Vector& X, int param);
  /// @todo: Mark for removal
  void set_old_FTable(const IDContainer& ftab);
 private:
  /// @todo: Mark for removal
  IDContainer theOldFTable;
};

#endif  // SRC_MODEL_ELIMINATION_MODEL_NODE_H_
