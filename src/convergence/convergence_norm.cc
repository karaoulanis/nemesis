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

#include "convergence/convergence_norm.h"
#include <stdio.h>
#include "analysis/analysis.h"
#include "control/control.h"
#include "domain/domain.h"
#include "elements/element.h"
#include "soe/soe.h"

ConvergenceNorm::ConvergenceNorm()
    : tol(3),
      maxIter(0),
      iter(0),
      ro(0.),
      uo(0.),
      wo(0.),
      LC(0),
      nSteps(0),
      step(0) {
}

ConvergenceNorm::~ConvergenceNorm() {
}

void ConvergenceNorm::set_check(int maxIterations,
                 double tolRabs, double tolRrel, double tolWrel) {
  maxIter = maxIterations;
  tol[0] = tolRabs;
  tol[1] = tolRrel;
  tol[2] = tolWrel;
}

void ConvergenceNorm::Init(int LCid, int steps) {
  LC = LCid;
  nSteps = steps;
  step = 0;
  printf("_________________________________________________________________\n");
  printf("Step   Iter   Lambda     Time    R(abs)    R(rel)    E(rel)  P.P.\n");
}

void ConvergenceNorm::NewStep() {
  step++;
  iter = 0;
  ro = 0;
  uo = 0;
  wo = 0;
  const Vector& Ro = pA->get_soe()->get_B();
  const Vector& dUo = pA->get_soe()->get_X();
  for (int i = 0; i < Ro.get_size(); i++) {
    double c1 = Ro[i];
    double c2 = dUo[i];
    ro+=c1*c1;
    // uo+=c2*c2;
    wo+=c1*c2;
  }
  ro = sqrt(ro);
  // uo = sqrt(uo);
  wo = abs(wo);
}

int ConvergenceNorm::Update() {
  iter++;
  double ri = 0, wi = 0;
  double ridro, widwo;
  // Get lambda
  double lambda = pA->get_control()->get_lambda();
  double time = pA->get_domain()->get_time_curr();
  // Get number of plastic points
  int nPP = 0;
  for (ElementIterator eIter = pA->get_domain()->get_elements().begin();
            eIter != pA->get_domain()->get_elements().end(); eIter++)
    nPP+=eIter->second->get_num_plastic_points();


  const Vector& dR = pA->get_soe()->get_B();
  const Vector& dU = pA->get_soe()->get_X();
  // for (int i = 0;i < dR.size();i++) cout << dR[i] <<endl;
  for (int i = 0; i < dR.get_size(); i++) {
    double c1 = dR[i];
    double c2 = dU[i];
    ri+=c1*c1;
    // ui+=c2*c2;
    wi+=c1*c2;
  }
  ri = sqrt(ri);
  // ui = sqrt(ui);
  wi = abs(wi);

  /// @todo Check if by setting to zero has any problems.
  if (num::smaller(ro, 1e-9)) {
    ridro = 0.;
  } else {
    ridro = ri/ro;  // If ro, uo, wo very
  }
  //  if (num::smaller(uo, 1e-9)) {
  //    uiduo = 0.;
  //  } else {
  //    uiduo = ui/uo;  // small, a 'pseudo'
  //  }
  if (num::smaller(wo, 1e-9)) {
    widwo = 0.;
  } else {
    widwo = wi/wo;  // residual appears.
  }

  // cout<<' ';
  if (iter == 1) {
    num::print_i(step,  6);
  } else {
    printf("      ");
  }
  num::print_i(iter,  5);
  printf(" ");
  num::print_d(lambda, 8, 3);
  printf(" ");
  num::print_d(time,  8, 3);
  printf(" ");
  num::print_d(ri,    9, 3);
  printf(" ");
  num::print_d(ridro, 9, 6);
  printf(" ");
  // num::print_d(ui,    9, 3);
  // cout<<' ';
  // num::print_d(uiduo, 9, 6);
  // cout<<' ';
  //  num::print_d(wi,    9, 3);
  // cout<<' ';
  num::print_d(widwo, 9, 6);
  printf(" ");
  num::print_i(nPP, 5);
  printf("\n");

  if (iter <= maxIter &&
      ri    < tol[0]  &&
      ridro < tol[1]  &&
      widwo < tol[2]) {
        return  0;
  } else if (iter > maxIter-1) {
    return -2;
  } else {
    return  1;
  }
}
