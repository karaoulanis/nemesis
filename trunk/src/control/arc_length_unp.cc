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

#include "control/arc_length_unp.h"
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
 * \param DeltaTime     Timestep for viscoplastic solutions
 */
ArcLengthUNP::ArcLengthUNP(double DL0, double minDL, double maxDL,
                     int IterDesired, double n, double DeltaTime)
:StaticControl(DL0, minDL, maxDL, IterDesired, n, DeltaTime) {
  myTag = TAG_CONTROL_ARC_LENGTH_UNP;
}

/**
 * Destructor.
 */
ArcLengthUNP::~ArcLengthUNP() {
  // Does nothing
}

/**
 * Creates new step for a Arc Length UNP control based static analysis.
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
 * then no auto-incrementation takes place and \f$\Delta l=\Delta l_{0}\f$.
 * It should be noted that min\f$\Delta\lambda\f$ and min\f$\Delta\lambda\f$
 * are given as absolute values.
 * \li Predictor Step\n
 * The predictor step is based on a forward Euler scheme. It is based on
 * p.285-286 of Cridfield's book. As soon this is computed the domain is
 * updated.
 */
void ArcLengthUNP::Predict() {
  // Find DLambda
  /// @todo Auto-incrementation involves abs() and this might be a problem...
//  DLambda = DLambda*pow(((double)Id/(double)Io), nExp);
//  if (fabs(DLambda)<minDelta)    DLambda = minDelta;
//  else if (fabs(DLambda)>maxDelta) DLambda = maxDelta;

  // Find duT (tangent du)
  pA->get_soe()->set_B(qRef);
  pA->get_soe()->solve();
  duT = pA->get_soe()->get_X();

  // Now DLambda can be found
  /// @todo rewrite eigenSign
  int sign = pA->get_soe()->get_eigen_sign();
  // if (sign == 0) exit(-11111);
  DLambda*=sign;
  lambdaTrial+=DLambda;


  // Find displacement vectors and update the model
  du = DLambda*duT;
  Du = du;
  pA->get_model()->incTrialDisp(du);
  pA->get_model()->update();

  // Set num of achieved iterations to one
  Io = 1;

  // Set du to the SOE so the norm can find it there
  pA->get_soe()->set_X(du);
}

/**
 * Updates that occur within each iterative step.
 * In this case only one root is taken, so no criterion is needed.
 */
void ArcLengthUNP::Correct() {
  // Find du_bar
  duBar = pA->get_soe()->get_X();

  // Find duT
  pA->get_soe()->set_B(qRef);
  pA->get_soe()->solve();
  duT = pA->get_soe()->get_X();

  dLambda=-Du*duBar/(Du*duT+DLambda);
  DLambda+=dLambda;
  Du+=du;

  // Update lambdaTrial and Du
  lambdaTrial+=dLambda;
  du = duBar+dLambda*duT;

  // Update displacements in the model
  pA->get_model()->incTrialDisp(du);
  pA->get_model()->update();

  // Increase number of iterations
  Io++;
}
