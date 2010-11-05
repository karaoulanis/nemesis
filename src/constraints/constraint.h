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

#ifndef SRC_CONSTRAINTS_CONSTRAINT_H_
#define SRC_CONSTRAINTS_CONSTRAINT_H_

// C++ system files
#include <vector>

// Project files
#include "domain/domain.h"
#include "domain/domain_object.h"
#include "node/node.h"

// Forward declarations
class Domain;
class Node;

struct CDof {
  Node* pNode_;
  int dof_;
  double coeff_;
};

/**
 * The Constraint class.
 */
class Constraint: public DomainObject {
 protected:
  static int num_constraints_;
  std::vector<CDof> cdofs_;
  double val_;
  double f_trial_;
  double f_convg_;
 public:
  // Constructors
  Constraint();

  void set_cdof(int NodeID, int dof, double coeff);
  void set_val(double val);

  int get_num_cdofs();
  const CDof& get_cdof(int i);
  virtual double get_val(double time = 0.);

  void inc_trial_force(double f);
  double get_disp(int i);
  double get_velc(int i);
  double get_disp_convg(int i);
  double get_disp_trial(int i);
  double get_velc_convg(int i);
  double get_accl_convg(int i);
  double get_F();

  void update(double f);
  void commit();
};
#endif  // SRC_CONSTRAINTS_CONSTRAINT_H_
