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

#include "control/load_control.h"
#include < cmath>

/**
 * Constructor.
 * \param DL0     Initial Delta lambda.
 * \param minDL     Lower bound for Delta lambda.
 * \param maxDL     Upper bound for Delta lambda.
 * \param IterDesired Desired number of iterations.
 * \param n       Exponent parameter.
 * \param DeltaTime     Timestep for viscoplastic solutions.
*/
LoadControl::LoadControl(double DL0, double minDL, double maxDL,
             int IterDesired, double n, double DeltaTime)
:StaticControl(DL0, minDL, maxDL, IterDesired, n, DeltaTime) {
  myTag = TAG_CONTROL_LOAD;
}
/**
 * Destructor.
 */
LoadControl::~LoadControl() {
}
/**
 * Creates new step for a load control based static analysis.
 * It does two things:
 * \li Auto-incrementation \n
 * \f$\Delta\lambda\f$ at each new step is given as decribed in "Non-linear 
 * Finite Element Analysis of Solids and Strucures", Vol.1, p.287, by (9.40):
 * \f[\Delta\lambda_{n}=\Delta\lambda_{0}
 *            \left(
 *            \frac{I_d}{I_o}
 *                      \right)^n
 * \f]
 * \f$I_d\f$ is the desired number of iterations within each step (Crisfield 
 * suggests ~3), \f$I_o\f$ is the number of iterations in the last step and 
 * \f$n\f$ is an exponent, usually set to 0.5 as suggested by Ramm.\n
 * \f$\Delta\lambda\f$ is also limited within min\f$\Delta\lambda\f$ and
 * max\f$\Delta\lambda\f$. By setting those equal to \f$\Delta\lambda_{0}\f$,
 * then no auto-incrementation takes place and 
 * \f$\Delta\lambda=\Delta\lambda_{0}\f$. It should be noted that
 * min\f$\Delta\lambda\f$ and min\f$\Delta\lambda\f$ are given as absolute 
 * values.
 * \li Predictor Step\n
 * The predictor step is based on a forward Euler scheme. As soon this is
 * computed the domain is updated.
 */
void LoadControl::predict() {
  // Find lambdaTrial
  ///@todo Auto-incrementation involves abs() and this might be a problem...
  ///@todo Check defaults
  // DLambda*=pow(((double)Id/(double)Io), nExp);
  // if (DLambda < minDelta) DLambda = minDelta;
  // else if (DLambda>maxDelta) DLambda = maxDelta;

  lambdaTrial = lambdaConvg+DLambda;

  // Find duT (tangent du)
  pA->getSOE()->setB(qRef);
  pA->getSOE()->solve();
  duT = pA->getSOE()->getX();

  // Find du and update model (predictor)
  du = DLambda*duT;
  pA->getModel()->incTrialDisp(du);
  pA->getModel()->update();
  Du = du;

  // Set num of achieved iterations to one
  Io = 1;

  // Set initial displacement to the SOE so that the norm can find it there
  pA->getSOE()->setX(du);
}
/**
 * Updates that occur within each iterative step.
 * The nodal displacements are updated, while \f$\delta\lambda\f$ is assummed 
 * to be zero, so no updates for \f$\Delta\lambda\f$ take place in this step. 
 */
void LoadControl::correct() {
  // Find incremental and accumulative displacements
  du = pA->getSOE()->getX();
  Du+=du;

  // Set accumulative displacement to the SOE so that the norm can find it there
  pA->getSOE()->setX(Du);

  // Update displacements in the model
  pA->getModel()->incTrialDisp(du);
  pA->getModel()->update();

  // Increase number of iterations
  Io++;
}
