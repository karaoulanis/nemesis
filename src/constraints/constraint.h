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

#ifndef NEMESIS_CONSTRAINTS_CONSTRAINT_H_
#define NEMESIS_CONSTRAINTS_CONSTRAINT_H_

// C++ system files
#include < vector>

// Project files
#include "domain/domain.h"
#include "domain/domain_object.h"
#include "node/node.h"

// Forward declarations
class Domain;
class Node;

struct cDof {
  Node* pNode;
  int dof;
  double coeff;
};

/**
 * The Constraint class.
 */
class Constraint: public DomainObject {
 protected:
  static int nConstraints;

  std::vector < cDof > theCDofs;
  double cVal;
  double fTrial;
  double fConvg;
  public:
  // Constructors
  Constraint();

  void setcDof(int NodeID, int dof, double coeff);
  void setcVal(double val);

  int getncDofs();
  const cDof& getcDof(int i);
  virtual double getcVal(double time = 0.);

  void incTrialForce(double f);
  double getDisp(int i);
  double getVelc(int i);
  double getDispConvg(int i);
  double getDispTrial(int i);
  double getVelcConvg(int i);
  double getAcclConvg(int i);
  double getF();

  void update(double f);
  void commit();
};
#endif  // NEMESIS_CONSTRAINTS_CONSTRAINT_H_
