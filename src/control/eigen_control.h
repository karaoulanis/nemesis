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

#ifndef SRC_CONTROL_EIGEN_CONTROL_H_
#define SRC_CONTROL_EIGEN_CONTROL_H_

#include "control/control.h"

class EigenControl :public Control {
  public:
  // Constructors and destructor
  EigenControl();
  ~EigenControl();

  void Init()     {}
  void Predict()    {}
  void Correct()    {}
  void Commit()   {}
  void FormResidual(double /*factor*/) {}
  void FormTangent();

  // Form tangent and residual element by element
  void FormElementalTangent(ModelElement* pModelElement);
  void FormElementalMassMatrix(ModelElement* pModelElement);
  void FormElementalResidual(ModelElement* pModelElement, double time = 0.);

  // Form residual node by node
  void FormNodalResidual(ModelNode* pModelNode);

  // Methods that are used through analysis
  void Rollback();
};
#endif  // SRC_CONTROL_EIGEN_CONTROL_H_
