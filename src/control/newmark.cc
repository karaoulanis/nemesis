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

#include "control/newmark.h"
#include <cmath>
#include "analysis/analysis.h"
#include "domain/domain.h"
#include "model/model.h"
#include "model/model_element.h"
#include "model/model_node.h"
#include "soe/soe.h"

/**
 * Constructor.
 */
Newmark::Newmark(double beta_, double gamma_, double dt_)
    : TransientControl(),
      beta(beta_),
      gamma(gamma_),
      dt(dt_) {
  myTag = TAG_CONTROL_LOAD;
  c[0]=1.0;
  c[1]=1.0/(beta*dt*dt);
  c[2]=gamma/(beta*dt);
}

/**
 * Destructor.
 */
Newmark::~Newmark() {
}

/**
 *
 */
void Newmark::Predict() {
  pA->get_domain()->incTime(dt);
  ut = u;
  vt = v;
  at = a;
  double c1 = 1.0-gamma/beta;
  double c2 = dt*(1.0-0.5*gamma/beta);
  double c3 = -1.0/(beta*dt);
  double c4 = 1.0-0.5/beta;
  v = c1*v+c2*at;
  a = c4*a+c3*vt;

  pA->get_model()->set_trial_vecs(u, v, a);
  pA->get_model()->Update();

  this->FormResidual(1.0);
  pA->get_soe()->solve();

  u = u+(pA->get_soe()->get_X());
  v = v+c[2]*(pA->get_soe()->get_X());
  a = a+c[1]*(pA->get_soe()->get_X());
  pA->get_model()->set_trial_vecs(u, v, a);
  pA->get_model()->Update();
}

/**
 *
 */
void Newmark::Correct() {
  u = u+(pA->get_soe()->get_X());
  v = v+c[2]*(pA->get_soe()->get_X());
  a = a+c[1]*(pA->get_soe()->get_X());

  pA->get_model()->set_trial_vecs(u, v, a);
  pA->get_model()->Update();
}
