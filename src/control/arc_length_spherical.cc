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

#include "control/arc_length_spherical.h"
#include <cmath>
#include "analysis/analysis.h"
#include "domain/domain.h"
#include "model/model.h"
#include "model/model_element.h"
#include "model/model_node.h"
#include "soe/soe.h"

/**
 * Constructor.
 * \param DL0     Initial Dl.
 * \param minDL     Lower bound for Dl.
 * \param maxDL     Upper bound for Dl.
 * \param IterDesired Desired number of iterations.
 * \param n       Exponent parameter.
 * \param DeltaTime     Timestep for viscoplastic solutions.
 */
ArcLengthSpherical::ArcLengthSpherical(double DL0, double minDL, double maxDL,
                     int IterDesired, double n, double DeltaTime)
:StaticControl(DL0, minDL, maxDL, IterDesired, n, DeltaTime), DeltaL(DL0) {
  DeltaL = DL0;
  myTag = TAG_CONTROL_ARC_LENGTH_SPHERICAL;
}

/**
 * Destructor.
 */
ArcLengthSpherical::~ArcLengthSpherical() {
  // Does nothing
}

/**
 * Creates new step for a Arc Length Spherical control based static analysis.
 * It does two things:
 * \li Auto-incrementation \n
 * \f$\Delta l\f$ at each new step is given as decribed in "Non-linear
 * Finite Element Analysis of Solids and Strucures", Vol.1, p.287, by (9.40):
 * \f[\Delta l_{n}=\Delta l_{0}
 *            \left(
 *            \frac{I_d}{I_o}
 *                      \right)^n
 * \f]
 * \f$I_d\f$ is the desired number of iterations within each step (Crisfield
 * suggests ~3), \f$I_o\f$ is the number of iterations in the last step and
 * \f$n\f$ is an exponent, usually set to 0.5 as suggested by Ramm.\n
 * \f$\Delta l\f$ is also limited within min\f$\Delta l\f$ and
 * max\f$\Delta l\f$. By setting those equal to \f$\Delta l_{0}\f$,
 * then no auto-incrementation takes place and \f$\Delta l=\Delta l_{0}\f$.
 * It should be noted that min\f$\Delta l\f$ and min\f$\Delta l\f$ are given as
 * absolute values.
 * \li Predictor Step\n
 * The predictor step is based on a forward Euler scheme. It is based on
 * p.285-286 of Cridfield's book. As soon this is computed the domain is
 * updated.
 */
void ArcLengthSpherical::Predict() {
  // Find DeltaL increment
  /// @todo Auto-incrementation involves abs() and this might be a problem...
  DeltaL*=pow(static_cast<double>(Id)/static_cast<double>(Io), nExp);
  if (DeltaL < minDelta) DeltaL = minDelta;
  else if (DeltaL>maxDelta) DeltaL = maxDelta;

  // Find duT (tangent du)
  pA->get_soe()->set_B(qRef);
  pA->get_soe()->solve();
  duT = pA->get_soe()->get_X();

  // Now DLambda can be found, see Crisfield, (9.34).
  int sign = pA->get_soe()->get_eigen_sign();
  if (sign == 0) exit(-11111);
  DLambda = sign*DeltaL/sqrt(duT*duT);
  lambdaTrial+=DLambda;

  // Find displacement vectors and update the model
  du = DLambda*duT;
  Du = du;
  pA->get_model()->incTrialDisp(du);
  pA->get_model()->Update();

  // Set num of achieved iterations to one
  Io = 1;

  // Set du to the SOE so the norm can find it there
  pA->get_soe()->set_X(du);
}

/**
 * Updates that occur within each iterative step.
 * \f$\delta\lambda\f$ is the root of the equation
 * \f[
 *    a_1\delta\lambda^2+a_2\delta\lambda+a3 = 0.
 * \f]
 * This equation gives two roots. The minimum angle criterion is used to choose
 * which of these two roots will be used.
 */
void ArcLengthSpherical::Correct() {
  // Find du_bar
  duBar = pA->get_soe()->get_X();

  // Find duT
  pA->get_soe()->set_B(qRef);
  pA->get_soe()->solve();
  duT = pA->get_soe()->get_X();

  // Find a-coefficients
  double psi = 0;
  double psi2qTq = psi*psi*qRef*qRef;
  double a1 = duT*duT+psi2qTq;
  double a2 = 2*duT*(Du+duBar)+2*DLambda*psi2qTq;
  double a3=(Du+duBar)*(Du+duBar)-DeltaL*DeltaL;
  double a4 = Du*duBar+Du*Du;
  double a5 = Du*duT;

  // Solve a1*dLambda^2+a2*dLambda+a3 = 0
  double dLambda1, dLambda2;
  double det = a2*a2-4*a1*a3;
  if (det < 0) {
    exit(-8989);
  } else if (det>0) {
    dLambda1=(-a2+sqrt(det))/(2*a1);
    dLambda2=(-a2-sqrt(det))/(2*a1);
  } else {
    dLambda1 = dLambda2=-a3/a2;  /// @todo Check the linear case.
  }

  // Apply minimum angle criterion
  double DL2cosTheta1 = a4+a5*dLambda1;
  double DL2cosTheta2 = a4+a5*dLambda2;
  dLambda = dLambda1;
  if (DL2cosTheta2>DL2cosTheta1) dLambda = dLambda2;

  /// @todo Find more about the external work criterion.
  // Positive external work
  // dLambda = dLambda1;
  // double DW=(*q)*((*Du)+(*duBar)+dLambda*(*duT));
  // if (DW < 0) dLambda = dLambda2;

  // Update lambdaTrial and Du
  lambdaTrial+=dLambda;
  du = duBar+dLambda*duT;
  Du+=du;

  // Update displacements in the model
  pA->get_model()->incTrialDisp(du);
  pA->get_model()->Update();

  // Increase number of iterations
  Io++;
}
