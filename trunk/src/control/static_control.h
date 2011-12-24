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

#ifndef SRC_CONTROL_STATIC_CONTROL_H_
#define SRC_CONTROL_STATIC_CONTROL_H_

#include "control/control.h"

/**
 * The StaticControl Class.
 * Static Control is a class derived from the Control class. It should
 * control the model through a static analysis, i.e.
 * \li 1. It should form for each ModelElement and for each ModelNode its
 * tangent matrix and its residual vector and it should insert those into the
 * SOE. I.e. It must create the SOE for a static analysis.
 * \li 2. It should create a predictor step. (new step)
 * \li 3. It should be able to create corrective steps (update).
 * \li 4. It should update the model with the iterative displacements.
 * \li 5. It should commit the converged displacements to the model after
 * each successive increment.
 * \li 6. It should return to the previous converged step, if the current
 * step fails.
 */
class StaticControl :public Control {
  public:
  // Constructors and destructor
  StaticControl();
  StaticControl(double D0, double minD, double maxD,
                int IterDesired, double n, double DeltaTime);
  virtual ~StaticControl();

  void FormResidual(double factor);

  // Form tangent and residual element by element
  virtual void FormElementalTangent(ModelElement* pModelElement);
  virtual void FormElementalResidual(ModelElement* pModelElement,
                                     double time = 0.);

  // Form residual node by node
  void FormNodalResidual(ModelNode* pModelNode);

  // Methods that are used through analysis
  virtual void Init();
  virtual void Commit();
  virtual void Rollback();

 protected:
  double Delta0;        /// < Initial Delta
  double minDelta;      /// < Lower bound for Delta (absolute);
  double maxDelta;      /// < Upper bound for Delta (absolute);
  Vector du;            /// < Iterative displacements
  Vector duT;           /// < SOE solution for q
  Vector duBar;         /// < SOE solution for lambda*q
  Vector Du;            /// < Accumulative solution within each step
  int Io;               /// < Number of iterations in the last step
  int Id;               /// < Desired number of iterations in this step
  double nExp;          /// < Exponent for the auto incrementation
  double Dt;            /// < Timestep for viscoplastic solutions
};
#endif  // SRC_CONTROL_STATIC_CONTROL_H_
