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

#ifndef SRC_CONTROL_EIGEN_CONTROL_H_
#define SRC_CONTROL_EIGEN_CONTROL_H_

#include "control/control.h"

class EigenControl :public Control {
  public:

  // Constructors and destructor
  EigenControl();
  ~EigenControl();

  void init()     {}
  void predict()    {}
  void correct()    {}
  void commit()   {}
  void formResidual(double /*factor*/) {}
  void formTangent();

  // Form tangent and residual element by element
  void formElementalTangent(ModelElement* pModelElement);
  void formElementalMassMatrix(ModelElement* pModelElement);
  void formElementalResidual(ModelElement* pModelElement, double time = 0.);

  // Form residual node by node
  void formNodalResidual(ModelNode* pModelNode);

  // Methods that are used through analysis
  void rollback();
};
#endif  // SRC_CONTROL_EIGEN_CONTROL_H_
