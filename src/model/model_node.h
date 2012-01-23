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

#ifndef SRC_MODEL_MODEL_NODE_H_
#define SRC_MODEL_MODEL_NODE_H_

#include "model/model_object.h"

class Node;
class Constraint;

class ModelNode :public ModelObject {
 public:
  // Constructors
  ModelNode();
  ModelNode(const IDContainer& FTable, Node* pNode = 0,
            Constraint* pConstraint = 0);
  virtual ~ModelNode();

  // Access to data members
  Node* get_node()                {return myNode;}
  Constraint* get_constraint()    {return myConstraint;}

  virtual void add_R(double factor = 1.0)=0;
  virtual void add_uTrial(double factor = 1.0);
  virtual void add_vTrial(double factor = 1.0);

  virtual void incTrialDisp(const Vector& du)=0;
  virtual void incTrialVecs(const Vector& du, const Vector& dv,
                            const Vector& da)=0;
  virtual void set_trial_disp(const Vector& u)=0;
  virtual void set_trial_vecs(const Vector& u, const Vector& v,
                              const Vector& a)=0;

  virtual void Commit()=0;
  virtual void commitSens(const Vector& /*ds*/, int /*param*/) {}
  virtual void rollback();

 protected:
  Node* myNode;
  Constraint* myConstraint;

 private:
  // Dummy copy constructor and copy assignment as to explicitly disable them.
  // Only the declarations are provided and not the definitions.
  // When called a linking error will occur.
  ModelNode(const ModelNode&);
  void operator=(const ModelNode&);
};

#endif  // SRC_MODEL_MODEL_NODE_H_
