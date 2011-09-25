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

#ifndef SRC_CONTROL_SENSITIVITY_CONTROL_H_
#define SRC_CONTROL_SENSITIVITY_CONTROL_H_

#include "control/control.h"
#include "numeric/vector.h"

class SensitivityControl :public Control {
 protected:
  Vector ds;
  int currParameter;
  public:

  // Constructors and destructor
  SensitivityControl();
  virtual ~SensitivityControl();

  // Form tangent and residual element by element
  virtual void FormElementalTangent(ModelElement* pModelElement);
  virtual void FormElementalResidual(ModelElement* pModelElement,
                                     double time = 0.);

  // Form residual node by node
  void FormNodalResidual(ModelNode* /*pModelNode*/) {}

  // Methods that are used through analysis
  virtual void Init();
  virtual void Predict()              {}
  virtual void Correct()              {}
  virtual void Commit();
  virtual void Rollback()             {}

  virtual void FormResidual(double factor);
};

#endif  // SRC_CONTROL_SENSITIVITY_CONTROL_H_
