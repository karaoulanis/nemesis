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

#include "model/elimination_model_node.h"

/**
 * Default constructor.
 */
EliminationModelNode::EliminationModelNode()
    : ModelNode(),
      theOldFTable(0) {
}

/**
 * Constructor.
 */
EliminationModelNode::EliminationModelNode(const IDContainer& FTable,
                                           Node* pNode)
    : ModelNode(FTable, pNode),
      theOldFTable(0) {
}

/**
 * Destructor.
 */
EliminationModelNode::~EliminationModelNode() {
}

void EliminationModelNode::set_old_FTable(const IDContainer& ftab) {
  theOldFTable = ftab;
}

void EliminationModelNode::add_R(double factor) {
  if (myNode->existsLoad())
    myVector->Add_cV(factor, myNode->get_R());
}

void EliminationModelNode::incTrialDisp(const Vector& du) {
  myVector->Clear();
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=du[theFTable[i]];
    }
  }
  myNode->incTrialDisp(*myVector);
}

void EliminationModelNode::incTrialVecs(const Vector& du, const Vector& dv,
                                        const Vector& da) {
  /// @todo Run these within a loop
  // Add displacements
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=du[theFTable[i]];
    } else {
      (*myVector)[i]=0.;
    }
  }
  myNode->incTrialDisp(*myVector);
  // Add velocities
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=dv[theFTable[i]];
    } else {
      (*myVector)[i]=0.;
    }
  }
  myNode->addTrialVelc(*myVector);
  // Add accelerations
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=da[theFTable[i]];
    } else {
      (*myVector)[i]=0.;
    }
  }
  myNode->addTrialAccl(*myVector);
}

void EliminationModelNode::set_trial_disp(const Vector& u) {
  myVector->Clear();
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=u[theFTable[i]];
    }
  }
  myNode->set_trial_disp(*myVector);
}

void EliminationModelNode::set_trial_vecs(const Vector& u, const Vector& v,
                                        const Vector& a) {
  /// @todo Run these within a loop
  // Set displacements
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=u[theFTable[i]];
    } else {
      (*myVector)[i]=0.;
    }
  }
  myNode->set_trial_disp(*myVector);
  // Set velocities
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=v[theFTable[i]];
    } else {
      (*myVector)[i]=0.;
    }
  }
  myNode->set_trial_velc(*myVector);
  // Set accelerations
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=a[theFTable[i]];
    } else {
      (*myVector)[i]=0.;
    }
  }
  myNode->set_trial_accl(*myVector);
}

void EliminationModelNode::Commit() {
  myNode->Commit();
}
void EliminationModelNode::commitSens(const Vector& X, int param) {
  myVector->Clear();
  for (unsigned i = 0; i < theFTable.size(); i++) {
    if (theFTable[i] >= 0) {
      (*myVector)[i]=X[theFTable[i]];
    }
  }
  myNode->commitSens(*myVector, param);
}
