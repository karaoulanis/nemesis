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

#ifndef SRC_CONTROL_TRANSIENT_CONTROL_H_
#define SRC_CONTROL_TRANSIENT_CONTROL_H_

#include "control/control.h"
#include "numeric/vector.h"

class ModelElement;
class ModelNode;

class TransientControl :public Control {
 public:
  // Constructors and destructor
  TransientControl();
  virtual ~TransientControl();

  // Form tangent and residuals
  virtual void FormElementalTangent(ModelElement* pModelElement);
  virtual void FormElementalResidual(ModelElement* pModelElement,
                                     double time = 0.);
  virtual void FormNodalResidual(ModelNode* pModelNode);
  virtual void FormResidual(double factor);

  virtual void Init();
  virtual void Commit();
  virtual void Rollback();

 protected:
  Vector u, v, a, ut, vt, at;
  double c[3];
};
#endif  // SRC_CONTROL_TRANSIENT_CONTROL_H_
